"""Generalized RDF (Radial Distribution Function) analysis for MD trajectories."""

import argparse
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from ovito.io import import_file
from ovito.data import CutoffNeighborFinder, DataCollection, SimulationCell
from ovito.pipeline import ModifierInterface
from traits.api import Float, Int, List as TraitsList, Dict as TraitsDict
from tqdm import tqdm
from scipy.optimize import curve_fit
import scipy.special
import scipy.integrate

try:
    import tomllib
except ImportError:
    import tomli as tomllib


class PartialRDFModifier(ModifierInterface):
    """OVITO modifier to calculate partial Radial Distribution Functions (RDFs)."""
    
    num_bins = Int(400, label="Number of bins")
    r_max = Float(6.0, label="Cutoff radius")
    projection_requests = TraitsList(TraitsDict)
    cos_angle_threshold = Float(0.866)  # Default cos(30 deg)

    def modify(self, data: DataCollection, *, frame: int, **kwargs):
        """Calculate RDFs for the current frame."""
        # Define bin edges and centers
        bins = np.linspace(0, self.r_max, self.num_bins + 1)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        dr = bins[1] - bins[0]

        # Get particle data
        if not data.particles or not hasattr(data.particles, 'particle_types'):
            raise RuntimeError(f"Frame {frame}: Particle data missing")
        
        ptypes = data.particles.particle_types[:]
        type_table = data.particles.particle_types
        positions = data.particles.positions
        unique_type_ids = np.unique(ptypes)
        num_particles = data.particles.count

        # Cell information
        volume = data.cell.volume
        if volume <= 0:
            raise RuntimeError(f"Frame {frame}: Invalid cell volume: {volume}")

        finder = CutoffNeighborFinder(self.r_max, data)
        
        # Count particles by type
        central_counts = {t_id: np.sum(ptypes == t_id) for t_id in unique_type_ids}

        # Calculate standard RDFs
        histograms = []
        labels = []
        
        for central_type_id in sorted(unique_type_ids):
            if central_counts[central_type_id] == 0:
                continue
            central_type_name = type_table.type_by_id(central_type_id).name or f"Type_{central_type_id}"

            for neighbor_type_id in sorted(unique_type_ids):
                neighbor_type_name = type_table.type_by_id(neighbor_type_id).name or f"Type_{neighbor_type_id}"
                label = f"{central_type_name}-{neighbor_type_name}"
                labels.append(label)

                if central_counts[neighbor_type_id] == 0:
                    histograms.append(np.zeros(self.num_bins))
                    continue

                # Collect distances
                all_dists = []
                central_indices = np.where(ptypes == central_type_id)[0]
                for central_idx in central_indices:
                    for neigh in finder.find(central_idx):
                        if ptypes[neigh.index] == neighbor_type_id:
                            if central_idx != neigh.index:  # Avoid self-counting
                                all_dists.append(neigh.distance)
                
                # Calculate histogram and normalize
                hist, _ = np.histogram(all_dists, bins=bins)
                N_central = central_counts[central_type_id]
                N_neighbor = central_counts[neighbor_type_id]
                
                if N_central == 0 or N_neighbor == 0:
                    normalized_hist = np.zeros(self.num_bins)
                else:
                    shell_volumes = 4.0 * np.pi * bin_centers**2 * dr
                    shell_volumes[shell_volumes <= 1e-10] = 1e-10
                    rho_neighbor = N_neighbor / volume
                    norm_factor = N_central * rho_neighbor * shell_volumes
                    valid_norm = norm_factor > 1e-10
                    normalized_hist = np.zeros_like(hist, dtype=float)
                    normalized_hist[valid_norm] = hist[valid_norm] / norm_factor[valid_norm]
                
                histograms.append(normalized_hist)

        # Handle projected RDFs if requested
        projected_histograms = []
        projected_labels = []
        
        if self.projection_requests:
            projected_dists_collected = {req['label']: [] for req in self.projection_requests}
            
            for central_idx in range(num_particles):
                central_type_id = ptypes[central_idx]
                for neigh in finder.find(central_idx):
                    if central_idx == neigh.index:
                        continue
                    neighbor_type_id = ptypes[neigh.index]
                    bond_vector = neigh.delta
                    
                    for proj_req in self.projection_requests:
                        if (central_type_id == proj_req['c_type_id'] and 
                            neighbor_type_id == proj_req['n_type_id']):
                            current_distance = neigh.distance
                            if current_distance < 1e-9:
                                continue
                            
                            norm_bond_vector = bond_vector / current_distance
                            
                            # Check alignment with any symmetry-equivalent direction
                            bond_added = False
                            for direction in proj_req.get('directions', [proj_req.get('direction', np.array([1,0,0]))]):
                                if bond_added:
                                    break
                                dot_product = np.dot(norm_bond_vector, direction)
                                if abs(dot_product) >= self.cos_angle_threshold:
                                    projected_dists_collected[proj_req['label']].append(current_distance)
                                    bond_added = True
            
            # Process projected histograms
            for proj_req in self.projection_requests:
                proj_label = proj_req['label']
                projected_labels.append(proj_label)
                
                dists = projected_dists_collected[proj_label]
                hist_proj, _ = np.histogram(dists, bins=bins)
                
                N_central_proj = central_counts.get(proj_req['c_type_id'], 0)
                N_neighbor_proj = central_counts.get(proj_req['n_type_id'], 0)
                
                if N_central_proj == 0 or N_neighbor_proj == 0 or volume <= 1e-9:
                    normalized_hist_proj = np.zeros(self.num_bins)
                else:
                    shell_volumes = 4.0 * np.pi * bin_centers**2 * dr
                    shell_volumes[shell_volumes <= 1e-10] = 1e-10
                    rho_neighbor_proj = N_neighbor_proj / volume
                    norm_factor_proj = N_central_proj * rho_neighbor_proj * shell_volumes
                    valid_norm_proj = norm_factor_proj > 1e-10
                    normalized_hist_proj = np.zeros_like(hist_proj, dtype=float)
                    normalized_hist_proj[valid_norm_proj] = hist_proj[valid_norm_proj] / norm_factor_proj[valid_norm_proj]
                
                projected_histograms.append(normalized_hist_proj)

        # Combine all histograms
        all_histograms = histograms + projected_histograms
        all_labels = labels + projected_labels
        
        if not all_histograms:
            return

        all_histograms_np = np.array(all_histograms)

        # Store results in data table
        table_id = "rdf-data"
        if table_id in data.tables:
            data.tables.delete(table_id)

        table = data.tables.create(identifier=table_id, title="Radial Distribution Functions")
        table.x = table.create_property("distance", data=bin_centers)
        table.y = table.create_property("g(r)", data=all_histograms_np.T, components=all_labels)

        # Store metadata
        data.attributes['rdf_particle_counts'] = central_counts
        data.attributes['rdf_volume'] = volume
        data.attributes['rdf_type_map'] = {t.id: t.name for t in type_table.types}


def load_config(config_path: Path) -> Dict[str, Any]:
    """Load RDF configuration from TOML file."""
    with open(config_path, 'rb') as f:
        return tomllib.load(f)


def parse_vector(vector_str: str) -> np.ndarray:
    """Parse comma-separated vector string."""
    try:
        return np.array([float(x) for x in vector_str.split(',')])
    except ValueError:
        raise ValueError(f"Invalid vector format: '{vector_str}'")


def generate_cubic_symmetry_directions(direction: np.ndarray) -> List[np.ndarray]:
    """Generate all symmetry-equivalent directions for cubic system."""
    norm = np.linalg.norm(direction)
    if norm < 1e-9:
        return []
    
    normalized_dir = direction / norm
    x, y, z = normalized_dir
    
    # All permutations and sign changes
    symmetry_dirs = set()
    permutations = [
        [x, y, z], [x, z, y], [y, x, z], 
        [y, z, x], [z, x, y], [z, y, x]
    ]
    
    for perm in permutations:
        for sx in [1, -1]:
            for sy in [1, -1]:
                for sz in [1, -1]:
                    new_dir = np.array([sx * perm[0], sy * perm[1], sz * perm[2]])
                    new_dir = np.round(new_dir, decimals=6)
                    symmetry_dirs.add(tuple(new_dir))
    
    return [np.array(d) for d in symmetry_dirs]


def analyze_shell_direct(r: np.ndarray, g_r: np.ndarray, rho_neighbor: float, 
                        r_min: float, r_max: float) -> Dict[str, Optional[float]]:
    """Analyze RDF shell to calculate coordination number and moments."""
    mask = (r >= r_min) & (r <= r_max)
    r_shell = r[mask]
    g_r_shell = g_r[mask]

    if len(r_shell) == 0 or np.all(g_r_shell < 1e-9):
        return {'coordination': 0.0, 'mean': None, 'variance': None, 'third_cumulant': None}

    # Coordination number
    integrand_cn = 4.0 * np.pi * r_shell**2 * rho_neighbor * g_r_shell
    try:
        coordination = np.trapz(integrand_cn, r_shell)
    except Exception:
        coordination = 0.0

    # Calculate moments using 4πr²g(r) as weights
    weights = 4.0 * np.pi * r_shell**2 * g_r_shell
    sum_weights = np.sum(weights)

    if sum_weights < 1e-9:
        return {'coordination': coordination, 'mean': None, 'variance': None, 'third_cumulant': None}

    # Mean distance
    try:
        avg_r = np.average(r_shell, weights=weights)
    except Exception:
        avg_r = None

    # Variance (MSRD)
    msrd = None
    if avg_r is not None:
        try:
            msrd = np.average((r_shell - avg_r)**2, weights=weights)
        except Exception:
            pass

    # Third cumulant
    third_moment = None
    if avg_r is not None:
        try:
            third_moment = np.average((r_shell - avg_r)**3, weights=weights)
        except Exception:
            pass

    return {
        'coordination': coordination,
        'mean': avg_r,
        'variance': msrd,
        'third_cumulant': third_moment
    }


def gaussian(x: np.ndarray, amp: float, mu: float, sigma: float) -> np.ndarray:
    """Gaussian function for fitting."""
    sigma = max(abs(sigma), 1e-6)
    return amp * np.exp(-(x - mu)**2 / (2 * sigma**2))


def gaussian_constrained_mean(x: np.ndarray, amp: float, sigma: float, mu_fixed: float) -> np.ndarray:
    """Gaussian function with fixed mean for fitting."""
    sigma = max(abs(sigma), 1e-6)
    return amp * np.exp(-(x - mu_fixed)**2 / (2 * sigma**2))


def skewed_gaussian(x: np.ndarray, amp: float, mu: float, sigma: float, alpha: float) -> np.ndarray:
    """Skewed Gaussian function for fitting."""
    sigma = max(abs(sigma), 1e-6)
    t = (x - mu) / sigma
    skew_factor = np.exp(-t**2/2) * (1 + scipy.special.erf(alpha * t / np.sqrt(2)))
    return amp * skew_factor


def fit_peak(r: np.ndarray, g_r: np.ndarray, rho_neighbor: float, 
             r_min: float, r_max: float, shell_name: str, 
             fit_type: str = "both", constrain_gaussian_mean: Optional[float] = None) -> Dict[str, Any]:
    """Fit RDF peak with Gaussian and/or skewed Gaussian."""
    mask = (r >= r_min) & (r <= r_max)
    r_shell = r[mask]
    g_r_shell = g_r[mask]

    if len(r_shell) < 4 or np.all(g_r_shell < 1e-9):
        return {}

    # Initial parameter estimates
    max_idx = np.argmax(g_r_shell)
    max_g = g_r_shell[max_idx]
    max_r = r_shell[max_idx]
    sigma_guess = (r_max - r_min) / 4

    results = {'shell_name': shell_name, 'r_range': (r_min, r_max)}

    # Gaussian fit
    if fit_type in ["gaussian", "both"]:
        try:
            if constrain_gaussian_mean is not None:
                # Constrained fit: only fit amplitude and sigma
                gaussian_fixed = lambda x, amp, sigma: gaussian_constrained_mean(x, amp, sigma, constrain_gaussian_mean)
                p0_constrained = [max_g, sigma_guess]
                bounds_constrained = ([0, 1e-4], [max_g * 5, (r_max - r_min) * 2])
                popt_constrained, pcov_constrained = curve_fit(gaussian_fixed, r_shell, g_r_shell,
                                                               p0=p0_constrained, bounds=bounds_constrained)
                
                # Unpack parameters and reconstruct full parameter set
                fit_amp, fit_sigma = popt_constrained
                fit_mean = constrain_gaussian_mean
                fit_sigma = abs(fit_sigma)
                popt_gauss = [fit_amp, fit_mean, fit_sigma]
                pcov_gauss = pcov_constrained
                
                g_r_fit = gaussian(r_shell, *popt_gauss)
            else:
                # Standard fit: fit all three parameters
                p0_gauss = [max_g, max_r, sigma_guess]
                bounds_gauss = ([0, r_min, 1e-4], [max_g * 5, r_max, (r_max - r_min) * 2])
                popt_gauss, pcov_gauss = curve_fit(gaussian, r_shell, g_r_shell, 
                                                  p0=p0_gauss, bounds=bounds_gauss)
                
                g_r_fit = gaussian(r_shell, *popt_gauss)
                fit_amp, fit_mean, fit_sigma = popt_gauss
                fit_sigma = abs(fit_sigma)
            
            residuals = g_r_shell - g_r_fit
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((g_r_shell - np.mean(g_r_shell))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 1e-9 else 1.0
            
            # Coordination number
            integrand_cn = lambda r_val: 4.0 * np.pi * r_val**2 * rho_neighbor * gaussian(r_val, *popt_gauss)
            try:
                coord, _ = scipy.integrate.quad(integrand_cn, r_min, r_max)
            except Exception:
                coord = 0.0

            results['gaussian'] = {
                'parameters': {'amplitude': fit_amp, 'mean': fit_mean, 'sigma': fit_sigma},
                'covariance': pcov_gauss,
                'r_squared': r_squared,
                'fitted_y': g_r_fit,
                'derived': {
                    'coordination': coord,
                    'mean': fit_mean,
                    'variance': fit_sigma**2,
                    'third_cumulant': 0.0
                },
                'constrained_mean': constrain_gaussian_mean is not None
            }
        except Exception:
            results['gaussian'] = None

    # Skewed Gaussian fit
    if fit_type in ["skewed", "both"]:
        try:
            alpha_guess = 0.5
            p0_skew = [max_g / 2.0, max_r, sigma_guess, alpha_guess]
            bounds_skew = ([0, r_min, 1e-4, -15], [max_g * 10, r_max, (r_max - r_min) * 2, 15])
            popt_skew, pcov_skew = curve_fit(skewed_gaussian, r_shell, g_r_shell,
                                            p0=p0_skew, bounds=bounds_skew)
            
            g_r_fit = skewed_gaussian(r_shell, *popt_skew)
            residuals = g_r_shell - g_r_fit
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((g_r_shell - np.mean(g_r_shell))**2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 1e-9 else 1.0

            # Calculate derived properties
            fit_amp, fit_loc, fit_scale, fit_alpha = popt_skew
            fit_scale = abs(fit_scale)
            
            # Analytical moments for skewed normal
            delta = fit_alpha / np.sqrt(1 + fit_alpha**2)
            delta_term = delta * np.sqrt(2 / np.pi)
            analytical_mean = fit_loc + fit_scale * delta_term
            analytical_variance = fit_scale**2 * (1 - delta_term**2)
            
            # Third cumulant
            gamma1 = ((4 - np.pi) / 2) * (delta_term**3) / (1 - delta_term**2)**(3/2)
            third_cumulant = gamma1 * (analytical_variance**(3/2))

            # Coordination number
            integrand_cn = lambda r_val: 4.0 * np.pi * r_val**2 * rho_neighbor * skewed_gaussian(r_val, *popt_skew)
            try:
                coord, _ = scipy.integrate.quad(integrand_cn, r_min, r_max)
            except Exception:
                coord = 0.0

            results['skewed'] = {
                'parameters': {'amplitude': fit_amp, 'location': fit_loc, 
                              'scale': fit_scale, 'alpha': fit_alpha},
                'covariance': pcov_skew,
                'r_squared': r_squared,
                'fitted_y': g_r_fit,
                'derived': {
                    'coordination': coord,
                    'mean': analytical_mean,
                    'variance': analytical_variance,
                    'third_cumulant': third_cumulant
                }
            }
        except Exception:
            results['skewed'] = None

    return results


def save_results_to_file(filename: Path, analysis_results: Dict[str, Any], 
                        config: Dict[str, Any], projected_results: Dict[str, Any]) -> None:
    """Save analysis results to text file."""
    with open(filename, 'w', encoding='utf-8') as f:
        f.write("RDF Analysis Report\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Input file: {config['input']['trajectory_file']}\n")
        f.write(f"Frames analyzed: {config['analysis']['frame_start']} to "
                f"{config['analysis']['frame_end']} (step {config['analysis']['frame_step']})\n")
        f.write(f"RDF cutoff: {config['rdf']['cutoff_radius']} Å\n")
        f.write(f"Number of bins: {config['rdf']['num_bins']}\n")
        f.write("\n")

        # Standard RDF results
        f.write("STANDARD RDF ANALYSIS\n")
        f.write("-" * 80 + "\n\n")
        
        for peak_key, results in analysis_results.items():
            if not results:
                continue
            
            pair_name = results.get('shell_name', peak_key)
            f.write(f"Analysis for: {pair_name}\n")
            r_range = results.get('r_range', (None, None))
            if r_range[0] is not None:
                f.write(f"Range: {r_range[0]:.3f} - {r_range[1]:.3f} Å\n")
            f.write("-" * 40 + "\n")
            
            # Direct calculation results
            direct = results.get('direct', {})
            if direct:
                f.write("Direct calculation:\n")
                f.write(f"  Coordination number: {direct.get('coordination', 0):.4f}\n")
                if direct.get('mean') is not None:
                    f.write(f"  Mean distance: {direct['mean']:.4f} Å\n")
                if direct.get('variance') is not None:
                    f.write(f"  Variance (MSRD): {direct['variance']:.6f} Å²\n")
                if direct.get('third_cumulant') is not None:
                    f.write(f"  Third cumulant: {direct['third_cumulant']:.6f} Å³\n")
            
            # Gaussian fit results
            gaussian_fit = results.get('gaussian')
            if gaussian_fit and gaussian_fit != 'None':
                f.write("\nGaussian fit")
                if gaussian_fit.get('constrained_mean', False):
                    f.write(" (mean constrained):\n")
                else:
                    f.write(":\n")
                f.write(f"  R-squared: {gaussian_fit.get('r_squared', 0):.4f}\n")
                params = gaussian_fit.get('parameters', {})
                f.write(f"  Mean: {params.get('mean', 0):.4f} Å")
                if gaussian_fit.get('constrained_mean', False):
                    f.write(" (fixed)")
                f.write("\n")
                f.write(f"  Sigma: {params.get('sigma', 0):.4f} Å\n")
                f.write(f"  Variance: {params.get('sigma', 0)**2:.6f} Å²\n")
                derived = gaussian_fit.get('derived', {})
                f.write(f"  Coordination: {derived.get('coordination', 0):.4f}\n")
            
            # Skewed Gaussian fit results
            skewed_fit = results.get('skewed')
            if skewed_fit and skewed_fit != 'None':
                f.write("\nSkewed Gaussian fit:\n")
                f.write(f"  R-squared: {skewed_fit.get('r_squared', 0):.4f}\n")
                params = skewed_fit.get('parameters', {})
                f.write(f"  Location: {params.get('location', 0):.4f} Å\n")
                f.write(f"  Scale: {params.get('scale', 0):.4f} Å\n")
                f.write(f"  Alpha: {params.get('alpha', 0):.4f}\n")
                derived = skewed_fit.get('derived', {})
                f.write(f"  Coordination: {derived.get('coordination', 0):.4f}\n")
                f.write(f"  Mean: {derived.get('mean', 0):.4f} Å\n")
                f.write(f"  Variance: {derived.get('variance', 0):.6f} Å²\n")
                f.write(f"  Third cumulant: {derived.get('third_cumulant', 0):.6f} Å³\n")
            
            f.write("\n" + "=" * 80 + "\n\n")

        # Projected RDF results
        if projected_results:
            f.write("\nPROJECTED RDF ANALYSIS\n")
            f.write("-" * 80 + "\n\n")
            
            for proj_label, results in projected_results.items():
                if not results or not results.get('direct'):
                    continue
                
                f.write(f"Analysis for: {proj_label}\n")
                r_range = results.get('r_range', (None, None))
                if r_range[0] is not None:
                    f.write(f"Range: {r_range[0]:.3f} - {r_range[1]:.3f} Å\n")
                f.write("-" * 40 + "\n")
                
                direct = results.get('direct', {})
                if direct:
                    f.write(f"  Coordination number: {direct.get('coordination', 0):.4f}\n")
                    if direct.get('mean') is not None:
                        f.write(f"  Mean distance: {direct['mean']:.4f} Å\n")
                    if direct.get('variance') is not None:
                        f.write(f"  Variance (MSRD): {direct['variance']:.6f} Å²\n")
                
                f.write("\n" + "=" * 80 + "\n\n")


def plot_rdfs(r: np.ndarray, rdfs: Dict[str, np.ndarray], 
             analysis_results: Dict[str, Any], output_file: Path,
             config: Dict[str, Any]) -> None:
    """Create RDF plots with fits."""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot RDFs
    colors = plt.cm.tab10(np.linspace(0, 1, len(rdfs)))
    for (label, g_r), color in zip(rdfs.items(), colors):
        ax.plot(r, g_r, label=label, lw=2, color=color)
    
    # Plot fits if available
    for peak_key, results in analysis_results.items():
        if not results:
            continue
        
        # Extract the actual pair name from the results
        pair_name = results.get('shell_name', '')
        if pair_name not in rdfs:
            continue
        
        r_min, r_max = results.get('r_range', (None, None))
        if r_min is None:
            continue
        
        # Find color for this pair
        pair_idx = list(rdfs.keys()).index(pair_name)
        color = colors[pair_idx]
        
        # Plot fit region
        mask = (r >= r_min) & (r <= r_max)
        r_fit = r[mask]
        
        # Plot Gaussian fit
        gaussian_fit = results.get('gaussian')
        if gaussian_fit and 'fitted_y' in gaussian_fit:
            ax.plot(r_fit, gaussian_fit['fitted_y'], '--', color=color, 
                   label=f"{pair_name} Gaussian fit", lw=2)
        
        # Plot skewed Gaussian fit
        skewed_fit = results.get('skewed')
        if skewed_fit and 'fitted_y' in skewed_fit:
            ax.plot(r_fit, skewed_fit['fitted_y'], ':', color=color,
                   label=f"{pair_name} Skewed fit", lw=2.5)
        
        # Highlight analysis region
        ax.axvspan(r_min, r_max, alpha=0.1, color=color)
    
    ax.set_xlabel('Distance r (Å)', fontsize=12)
    ax.set_ylabel('g(r)', fontsize=12)
    ax.set_title('Radial Distribution Functions', fontsize=14)
    ax.axhline(1, color='grey', linestyle='--', linewidth=0.8)
    ax.legend(loc='best', fontsize='small')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0, right=config['rdf']['cutoff_radius'])
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


def main():
    """Main function for RDF analysis."""
    parser = argparse.ArgumentParser(
        description='Calculate and analyze RDFs from MD trajectories'
    )
    parser.add_argument('--rdf', type=str, required=True,
                       help='Path to RDF configuration TOML file')
    
    args = parser.parse_args()
    
    # Load configuration
    config_path = Path(args.rdf)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    config = load_config(config_path)
    
    # Set up paths
    traj_file = Path(config['input']['trajectory_file'])
    if not traj_file.exists():
        raise FileNotFoundError(f"Trajectory file not found: {traj_file}")
    
    output_dir = Path(config['output']['directory'])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load trajectory
    print(f"Loading trajectory: {traj_file}")
    pipeline = import_file(str(traj_file))
    num_frames = pipeline.source.num_frames
    
    # Set up simulation cell if provided
    if 'lattice' in config:
        from ovito.data import SimulationCell
        
        dynamic_lattice = config['lattice'].get('dynamic', False)
        
        if dynamic_lattice:
            # Dynamic lattice vectors - will be read from each frame
            print("Using dynamic lattice vectors from trajectory file")
            
            # Just set PBC if specified
            pbc = config['lattice'].get('pbc', [True, True, True])
            
            def set_pbc(frame, data):
                if data.cell is not None:
                    data.cell.pbc = pbc
                else:
                    # Create a dummy cell with PBC if none exists
                    # This should rarely happen with proper trajectory files
                    print(f"Warning: No cell in frame {frame}, cannot apply PBC")
            
            pipeline.modifiers.append(set_pbc)
            print(f"Set periodic boundary conditions: {pbc}")
        else:
            # Static lattice vectors - use provided values
            a = config['lattice']['a']
            b = config['lattice']['b']
            c = config['lattice']['c']
            pbc = config['lattice'].get('pbc', [True, True, True])
            
            cell_matrix = np.array([
                [a[0], a[1], a[2], 0.0],
                [b[0], b[1], b[2], 0.0],
                [c[0], c[1], c[2], 0.0]
            ])
            
            # Add modifier to set the cell
            def set_cell(frame, data):
                data.cell = SimulationCell(matrix=cell_matrix, pbc=pbc)
            
            pipeline.modifiers.append(set_cell)
            print(f"Set static simulation cell with PBC: {pbc}")
            print(f"  a = {a}")
            print(f"  b = {b}")
            print(f"  c = {c}")
    
    # Get frame range
    frame_start = config['analysis'].get('frame_start', 0)
    frame_end = min(config['analysis'].get('frame_end', num_frames), num_frames)
    frame_step = config['analysis'].get('frame_step', 1)
    
    print(f"Processing frames {frame_start} to {frame_end-1} (step {frame_step})")
    
    # Parse projection requests if any
    projection_requests = []
    projections_config = config.get('projections', {})
    cos_angle_threshold = np.cos(np.radians(projections_config.get('angle_threshold', 30.0)))
    
    if projections_config.get('bonds'):
        # Get type mapping from first frame
        data = pipeline.compute(0)
        type_map = {pt.name: pt.id for pt in data.particles.particle_types.types}
        
        for bond in projections_config['bonds']:
            center_name = bond['center']
            neighbor_name = bond['neighbor']
            direction = np.array(bond['direction'])
            
            c_id = type_map.get(center_name)
            n_id = type_map.get(neighbor_name)
            
            if c_id is None or n_id is None:
                print(f"Warning: Unknown atom type in projection {center_name}-{neighbor_name}")
                continue
            
            # Handle symmetry
            if projections_config.get('use_symmetry', False):
                sym_dirs = generate_cubic_symmetry_directions(direction)
                label = f"proj_{center_name}-{neighbor_name}_family"
                projection_requests.append({
                    'c_type_id': c_id,
                    'n_type_id': n_id,
                    'directions': sym_dirs,
                    'label': label,
                    'original_c_name': center_name,
                    'original_n_name': neighbor_name,
                    'original_dir_str': str(direction)
                })
            else:
                norm_dir = direction / np.linalg.norm(direction)
                label = f"proj_{center_name}-{neighbor_name}_{norm_dir[0]:.2f},{norm_dir[1]:.2f},{norm_dir[2]:.2f}"
                projection_requests.append({
                    'c_type_id': c_id,
                    'n_type_id': n_id,
                    'direction': norm_dir,
                    'label': label,
                    'original_c_name': center_name,
                    'original_n_name': neighbor_name,
                    'original_dir_str': str(direction)
                })
    
    # Add RDF modifier
    rdf_modifier = PartialRDFModifier(
        r_max=config['rdf']['cutoff_radius'],
        num_bins=config['rdf']['num_bins'],
        projection_requests=projection_requests,
        cos_angle_threshold=cos_angle_threshold
    )
    pipeline.modifiers.append(rdf_modifier)
    
    # Process frames
    frames_to_process = range(frame_start, frame_end, frame_step)
    all_rdf_data = []
    
    print(f"\nCalculating RDFs for {len(frames_to_process)} frames...")
    for frame_idx in tqdm(frames_to_process, desc="Processing frames"):
        try:
            data = pipeline.compute(frame_idx)
            table = data.tables.get('rdf-data')
            if table:
                all_rdf_data.append({
                    'y_data': table.y[:],
                    'x_data': table.x[:],
                    'counts': data.attributes['rdf_particle_counts'],
                    'volume': data.attributes['rdf_volume'],
                    'type_map': data.attributes['rdf_type_map'],
                    'labels': list(table.y.component_names)
                })
        except Exception as e:
            print(f"\nError processing frame {frame_idx}: {e}")
            continue
    
    if not all_rdf_data:
        raise RuntimeError("No valid RDF data computed")
    
    # Average RDFs
    print("\nAveraging RDFs...")
    avg_y_data = np.mean([item['y_data'] for item in all_rdf_data], axis=0)
    avg_volume = np.mean([item['volume'] for item in all_rdf_data])
    
    # Get distances and labels
    r_distances = all_rdf_data[-1]['x_data']
    labels = all_rdf_data[-1]['labels']
    type_map = all_rdf_data[-1]['type_map']
    
    # Create RDF dictionary
    rdfs = {label: avg_y_data[:, i] for i, label in enumerate(labels)}
    
    # Get particle counts
    all_type_ids = set().union(*(item['counts'].keys() for item in all_rdf_data))
    avg_counts = {tid: np.mean([item['counts'].get(tid, 0) for item in all_rdf_data]) 
                  for tid in all_type_ids}
    
    print(f"\nAverage volume: {avg_volume:.3f} Å³")
    print("Average particle counts:")
    for tid, count in avg_counts.items():
        print(f"  {type_map.get(tid, f'Type_{tid}')} (ID {tid}): {count:.2f}")
    
    # Analyze peaks
    print("\nAnalyzing RDF peaks...")
    analysis_results = {}
    
    for peak_config in config.get('peaks', []):
        pair_name = f"{peak_config['center']}-{peak_config['neighbor']}"
        if pair_name not in rdfs:
            print(f"Warning: RDF for {pair_name} not found")
            continue
        
        print(f"\nAnalyzing {pair_name} peak...")
        g_r = rdfs[pair_name]
        
        # Get neighbor density
        neighbor_type = peak_config['neighbor']
        neighbor_id = None
        for tid, tname in type_map.items():
            if tname == neighbor_type:
                neighbor_id = tid
                break
        
        if neighbor_id is None:
            print(f"Warning: Type {neighbor_type} not found")
            continue
        
        rho_neighbor = avg_counts[neighbor_id] / avg_volume
        
        # Create unique key for this peak (including range)
        peak_key = f"{pair_name}_{peak_config['r_min']:.1f}-{peak_config['r_max']:.1f}"
        
        # Analyze
        results = {
            'shell_name': pair_name,
            'r_range': (peak_config['r_min'], peak_config['r_max'])
        }
        
        # Direct calculation
        results['direct'] = analyze_shell_direct(
            r_distances, g_r, rho_neighbor,
            peak_config['r_min'], peak_config['r_max']
        )
        
        # Fitting
        fit_type = peak_config.get('fit_type', 'both')
        if fit_type != 'none':
            constrain_gaussian_mean = peak_config.get('constrain_gaussian_mean', None)
            fit_results = fit_peak(
                r_distances, g_r, rho_neighbor,
                peak_config['r_min'], peak_config['r_max'],
                pair_name, fit_type, constrain_gaussian_mean
            )
            results.update(fit_results)
        
        analysis_results[peak_key] = results
    
    # Analyze projected RDFs
    projected_results = {}
    if projection_requests:
        print("\nAnalyzing projected RDFs...")
        for proj_req in projection_requests:
            label = proj_req['label']
            if label in rdfs:
                g_r = rdfs[label]
                neighbor_id = proj_req['n_type_id']
                rho_neighbor = avg_counts[neighbor_id] / avg_volume
                
                # Find appropriate range from peaks config
                c_name = proj_req['original_c_name']
                n_name = proj_req['original_n_name']
                r_min, r_max = 0, config['rdf']['cutoff_radius']
                
                for peak in config.get('peaks', []):
                    if peak['center'] == c_name and peak['neighbor'] == n_name:
                        r_min, r_max = peak['r_min'], peak['r_max']
                        break
                
                results = {
                    'r_range': (r_min, r_max),
                    'direct': analyze_shell_direct(r_distances, g_r, rho_neighbor, r_min, r_max)
                }
                projected_results[label] = results
    
    # Save results
    report_file = output_dir / config['output']['report_file']
    print(f"\nSaving analysis report to: {report_file}")
    save_results_to_file(report_file, analysis_results, config, projected_results)
    
    # Create plots
    if config['output'].get('plot_file'):
        plot_file = output_dir / config['output']['plot_file']
        print(f"Creating RDF plot: {plot_file}")
        plot_rdfs(r_distances, rdfs, analysis_results, plot_file, config)
    
    print("\nRDF analysis complete!")


if __name__ == '__main__':
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    warnings.filterwarnings("ignore", category=scipy.integrate.IntegrationWarning)
    main()
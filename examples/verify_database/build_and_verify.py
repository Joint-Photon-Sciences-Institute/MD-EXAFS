#!/usr/bin/env python
"""
Build a database from multipath_folders and verify the chi(k) calculations
by comparing database results with FEFF chi.dat files.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from md_exafs.build_chi_database import build_database
from md_exafs.chi_database_query import ChiDatabaseQuery, PathQuery


def build_test_database():
    """Build database from reference multipath folders."""
    base_dir = Path("/home/nickj/claude-code/MD-EXAFS/reference_code_base/multipath_folders")
    db_path = Path("test_verification.db")
    
    if db_path.exists():
        print(f"Removing existing database: {db_path}")
        db_path.unlink()
    
    print(f"Building database from: {base_dir}")
    print("This may take a few minutes...")
    
    try:
        build_database(base_dir, db_path, num_workers=4, rebuild=True)
        print(f"Database built successfully: {db_path}")
        return db_path
    except Exception as e:
        print(f"Error building database: {e}")
        return None


def get_atom_folders(base_dir):
    """Find all atom folders containing chi.dat files."""
    atom_folders = []
    
    for chi_file in base_dir.rglob("chi.dat"):
        atom_dir = chi_file.parent
        # Extract frame and atom info
        frame = 0
        atom_id = 0
        
        # Try to extract atom_id
        if atom_dir.name.startswith("atom_"):
            try:
                atom_id = int(atom_dir.name.split('_')[1])
            except:
                pass
        
        # Try to find frame number
        for parent in atom_dir.parents:
            if parent.name.startswith("frame_"):
                try:
                    frame = int(parent.name.split('_')[1])
                    break
                except:
                    pass
        
        atom_folders.append((atom_dir, frame, atom_id))
    
    return sorted(atom_folders)


def sum_all_paths_from_database(db_path, frame, atom_id):
    """Sum all paths for a specific atom from the database."""
    with ChiDatabaseQuery(db_path) as db:
        # Query all paths for this atom
        query = PathQuery(frames=[frame], atom_ids=[atom_id])
        paths = db.query_paths(query)
        
        if not paths:
            return None, 0
        
        # Get all path IDs
        path_ids = [p['id'] for p in paths]
        
        # Sum chi data (includes degeneracy)
        summed_chi = db.sum_chi_data(path_ids)
        
        return summed_chi, len(paths)


def compare_chi_files(atom_folders, db_path):
    """Compare chi.dat files with database calculations."""
    results = []
    
    print(f"\nComparing {len(atom_folders)} atom folders...")
    
    for atom_dir, frame, atom_id in tqdm(atom_folders[:10], desc="Processing atoms"):  # Limit to first 10 for testing
        # Load FEFF chi.dat
        chi_dat_file = atom_dir / "chi.dat"
        if not chi_dat_file.exists():
            continue
        
        try:
            feff_data = np.loadtxt(chi_dat_file, comments='#')
            feff_k = feff_data[:, 0]
            feff_chi = feff_data[:, 1]
        except Exception as e:
            print(f"Error reading {chi_dat_file}: {e}")
            continue
        
        # Get chi from database
        db_chi_data, num_paths = sum_all_paths_from_database(db_path, frame, atom_id)
        
        if db_chi_data is None:
            print(f"No database data for frame {frame}, atom {atom_id}")
            continue
        
        db_k = db_chi_data[:, 0]
        db_chi = db_chi_data[:, 1]
        
        # Calculate differences
        # Ensure same length for comparison
        min_len = min(len(feff_chi), len(db_chi))
        diff = db_chi[:min_len] - feff_chi[:min_len]
        abs_diff = np.abs(diff)
        
        # Calculate metrics
        max_abs_diff = np.max(abs_diff)
        avg_abs_diff = np.mean(abs_diff)
        rms_diff = np.sqrt(np.mean(diff**2))
        
        results.append({
            'atom_dir': atom_dir,
            'frame': frame,
            'atom_id': atom_id,
            'num_paths': num_paths,
            'feff_k': feff_k,
            'feff_chi': feff_chi,
            'db_k': db_k,
            'db_chi': db_chi,
            'max_abs_diff': max_abs_diff,
            'avg_abs_diff': avg_abs_diff,
            'rms_diff': rms_diff
        })
    
    return results


def plot_comparisons(results, output_dir):
    """Create comparison plots."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Summary statistics
    max_diffs = [r['max_abs_diff'] for r in results]
    avg_diffs = [r['avg_abs_diff'] for r in results]
    
    print(f"\nSummary of {len(results)} comparisons:")
    print(f"  Average of max differences: {np.mean(max_diffs):.6f}")
    print(f"  Maximum of max differences: {np.max(max_diffs):.6f}")
    print(f"  Average of avg differences: {np.mean(avg_diffs):.6f}")
    
    # Plot first few comparisons
    n_plots = min(6, len(results))
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for i in range(n_plots):
        ax = axes[i]
        r = results[i]
        
        # Plot k^2 weighted chi
        k2_feff = r['feff_k']**2 * r['feff_chi']
        k2_db = r['db_k']**2 * r['db_chi']
        
        ax.plot(r['feff_k'], k2_feff, 'r-', label='FEFF chi.dat', linewidth=2)
        ax.plot(r['db_k'], k2_db, 'b--', label='Database', linewidth=2, alpha=0.8)
        
        ax.set_xlabel('k (Å⁻¹)')
        ax.set_ylabel('k² × χ(k)')
        ax.set_title(f"Frame {r['frame']}, Atom {r['atom_id']}\n" + 
                     f"Max diff: {r['max_abs_diff']:.2e}")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 20)
    
    plt.tight_layout()
    plot_file = output_dir / "chi_comparison_k2.png"
    plt.savefig(plot_file, dpi=150)
    print(f"\nComparison plot saved to: {plot_file}")
    plt.close()
    
    # Plot difference histogram
    plt.figure(figsize=(8, 6))
    plt.hist(max_diffs, bins=30, edgecolor='black', alpha=0.7)
    plt.xlabel('Maximum Absolute Difference')
    plt.ylabel('Count')
    plt.title('Distribution of Maximum Differences\nbetween FEFF and Database')
    plt.yscale('log')
    
    hist_file = output_dir / "difference_histogram.png"
    plt.savefig(hist_file, dpi=150)
    print(f"Histogram saved to: {hist_file}")
    plt.close()


def main():
    """Main verification workflow."""
    print("=" * 60)
    print("Database Verification Test")
    print("=" * 60)
    
    # Build database
    db_path = build_test_database()
    if db_path is None:
        return
    
    # Find atom folders
    base_dir = Path("/home/nickj/claude-code/MD-EXAFS/reference_code_base/multipath_folders")
    atom_folders = get_atom_folders(base_dir)
    print(f"\nFound {len(atom_folders)} atom folders with chi.dat files")
    
    # Compare results
    results = compare_chi_files(atom_folders, db_path)
    
    if not results:
        print("No comparison results obtained!")
        return
    
    # Create plots
    plot_comparisons(results, "verification_plots")
    
    # Check if differences are acceptable
    all_max_diffs = [r['max_abs_diff'] for r in results]
    if np.max(all_max_diffs) < 1e-3:
        print("\n✓ SUCCESS: Database chi calculations match FEFF within 0.001!")
    else:
        print(f"\n⚠ WARNING: Maximum difference is {np.max(all_max_diffs):.2e}")
    
    print("\nVerification complete!")


if __name__ == "__main__":
    main()
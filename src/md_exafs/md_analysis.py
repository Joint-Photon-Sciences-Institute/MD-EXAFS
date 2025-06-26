"""MD analysis tool for convergence checking and RDF calculations."""

import argparse
import sys
from pathlib import Path

# Import the specific analysis modules
from .convergence_check import parse_ener_file, plot_convergence
from .rdf_analysis import main as rdf_main


def convergence_analysis(args):
    """Handle convergence checking."""
    input_path = Path(args.input)
    output_path = Path(args.output)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    
    # Parse energy file
    print(f"Reading energy file: {input_path}")
    data = parse_ener_file(input_path)
    
    if len(data['step']) == 0:
        raise ValueError("No data found in energy file")
    
    # Create convergence plots
    print("Creating convergence plots...")
    plot_convergence(data, output_path)


def get_lattice_vectors(args):
    """Get lattice vectors from a specific frame."""
    from ovito.io import import_file
    import numpy as np
    
    traj_path = Path(args.trajectory)
    if not traj_path.exists():
        raise FileNotFoundError(f"Trajectory file not found: {traj_path}")
    
    frame = args.frame
    
    print(f"Loading trajectory: {traj_path}")
    pipeline = import_file(str(traj_path))
    
    if frame >= pipeline.source.num_frames:
        raise ValueError(f"Frame {frame} out of range. Trajectory has {pipeline.source.num_frames} frames.")
    
    # Compute the specified frame
    data = pipeline.compute(frame)
    
    if data.cell is None:
        print(f"\nWarning: No cell information found in frame {frame}")
        print("The trajectory file may not contain lattice vector information.")
        return
    
    # Get cell matrix
    matrix = data.cell.matrix
    a = matrix[0, :3]
    b = matrix[1, :3]
    c = matrix[2, :3]
    
    # Calculate cell parameters
    a_len = np.linalg.norm(a)
    b_len = np.linalg.norm(b)
    c_len = np.linalg.norm(c)
    
    # Calculate angles
    alpha = np.arccos(np.dot(b, c) / (b_len * c_len)) * 180 / np.pi
    beta = np.arccos(np.dot(a, c) / (a_len * c_len)) * 180 / np.pi
    gamma = np.arccos(np.dot(a, b) / (a_len * b_len)) * 180 / np.pi
    
    volume = data.cell.volume
    
    print(f"\nLattice vectors for frame {frame}:")
    print(f"  a = [{a[0]:.6f}, {a[1]:.6f}, {a[2]:.6f}]")
    print(f"  b = [{b[0]:.6f}, {b[1]:.6f}, {b[2]:.6f}]")
    print(f"  c = [{c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}]")
    print(f"\nCell parameters:")
    print(f"  |a| = {a_len:.6f} Å")
    print(f"  |b| = {b_len:.6f} Å")
    print(f"  |c| = {c_len:.6f} Å")
    print(f"  α = {alpha:.2f}°")
    print(f"  β = {beta:.2f}°")
    print(f"  γ = {gamma:.2f}°")
    print(f"  Volume = {volume:.3f} Å³")
    
    # Print PBC status
    pbc = data.cell.pbc
    print(f"\nPeriodic boundary conditions: [{pbc[0]}, {pbc[1]}, {pbc[2]}]")
    
    # Print TOML format for easy copying
    print("\n# For TOML configuration file:")
    print("[lattice]")
    print(f"a = [{a[0]:.6f}, {a[1]:.6f}, {a[2]:.6f}]")
    print(f"b = [{b[0]:.6f}, {b[1]:.6f}, {b[2]:.6f}]")
    print(f"c = [{c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}]")
    print(f"pbc = [{str(pbc[0]).lower()}, {str(pbc[1]).lower()}, {str(pbc[2]).lower()}]")


def main():
    """Main entry point for MD analysis tools."""
    parser = argparse.ArgumentParser(
        description='MD analysis tools for convergence checking and RDF calculations',
        epilog='Use one of: --check_convergence, --rdf, or --get_lattice_vectors'
    )
    
    # Convergence checking arguments
    parser.add_argument('--check_convergence', action='store_true',
                       help='Enable convergence checking mode')
    parser.add_argument('--input', type=str,
                       help='Path to .ener file (for convergence checking)')
    parser.add_argument('--output', type=str,
                       help='Path for output plot (PNG) (for convergence checking)')
    
    # RDF analysis arguments
    parser.add_argument('--rdf', type=str,
                       help='Path to RDF configuration TOML file')
    
    # Lattice vector extraction arguments
    parser.add_argument('--get_lattice_vectors', action='store_true',
                       help='Extract lattice vectors from a specific frame')
    parser.add_argument('--trajectory', type=str,
                       help='Path to trajectory file (for lattice vector extraction)')
    parser.add_argument('--frame', type=int, default=0,
                       help='Frame number to extract lattice vectors from (default: 0)')
    
    args = parser.parse_args()
    
    # Count how many modes are selected
    modes_selected = sum([args.check_convergence, bool(args.rdf), args.get_lattice_vectors])
    
    if modes_selected > 1:
        parser.error("Please specify only one mode: --check_convergence, --rdf, or --get_lattice_vectors")
    
    if modes_selected == 0:
        parser.error("Please specify one mode: --check_convergence, --rdf, or --get_lattice_vectors")
    
    # Handle convergence checking
    if args.check_convergence:
        if not args.input or not args.output:
            parser.error("--check_convergence requires --input and --output arguments")
        convergence_analysis(args)
    
    # Handle RDF analysis
    elif args.rdf:
        # Pass the argument to the RDF analysis module
        sys.argv = ['md-exafs-md', '--rdf', args.rdf]
        rdf_main()
    
    # Handle lattice vector extraction
    elif args.get_lattice_vectors:
        if not args.trajectory:
            parser.error("--get_lattice_vectors requires --trajectory argument")
        get_lattice_vectors(args)


if __name__ == '__main__':
    main()
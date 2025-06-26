#!/usr/bin/env python3
"""
Plot chi(k) spectra from the MD-EXAFS database.

This script allows plotting of:
- Individual path spectra
- Sum of paths for a specific atom
- Average of paths matching criteria
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import List, Optional, Tuple
import sqlite3

# Add the src directory to path to import md_exafs modules
import sys
sys.path.insert(0, str(Path(__file__).parent / "src"))

from md_exafs.chi_database_query import ChiDatabaseQuery, PathQuery


def plot_individual_paths(db: ChiDatabaseQuery, 
                         frame: int, 
                         atom_id: int, 
                         path_nums: Optional[List[int]] = None,
                         output_file: Optional[str] = None) -> None:
    """Plot individual path spectra for a specific atom."""
    # Query paths for this atom
    query = PathQuery(frames=[frame], atom_ids=[atom_id])
    paths = db.query_paths(query)
    
    if not paths:
        print(f"No paths found for frame {frame}, atom {atom_id}")
        return
    
    # Filter by path numbers if specified
    if path_nums:
        paths = [p for p in paths if p['path_num'] in path_nums]
    
    if not paths:
        print(f"No paths found matching path numbers {path_nums}")
        return
    
    # Get chi data
    path_ids = [p['id'] for p in paths]
    chi_data = db.get_chi_data(path_ids)
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    for path in paths:
        if path['id'] in chi_data:
            k_grid, chi_values = chi_data[path['id']]
            label = f"Path {path['path_num']}: {path['path_type']} (reff={path['reff']:.2f} Å)"
            plt.plot(k_grid, chi_values, label=label, alpha=0.7)
    
    plt.xlabel('k (Å⁻¹)')
    plt.ylabel('χ(k)')
    plt.title(f'Individual Path Spectra - Frame {frame}, Atom {atom_id}')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()


def plot_atom_sum(db: ChiDatabaseQuery,
                  frame: int,
                  atom_id: int,
                  path_types: Optional[List[str]] = None,
                  max_reff: Optional[float] = None,
                  output_file: Optional[str] = None) -> None:
    """Plot sum of paths for a specific atom."""
    # Query paths for this atom
    query = PathQuery(
        frames=[frame], 
        atom_ids=[atom_id],
        path_types=path_types,
        max_reff=max_reff
    )
    paths = db.query_paths(query)
    
    if not paths:
        print(f"No paths found for frame {frame}, atom {atom_id}")
        return
    
    # Get path IDs
    path_ids = [p['id'] for p in paths]
    
    # Sum chi data (paths are interpolated to uniform k-grid before summing)
    summed_data = db.sum_chi_data(path_ids)
    
    if summed_data is None:
        print("No chi data found for selected paths")
        return
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    k_grid = summed_data[:, 0]
    chi_sum = summed_data[:, 1]
    
    plt.plot(k_grid, chi_sum, 'b-', linewidth=2)
    
    # Add information about what was summed
    info_text = f"Frame {frame}, Atom {atom_id}\n"
    info_text += f"{len(paths)} paths summed"
    if path_types:
        info_text += f"\nPath types: {', '.join(path_types)}"
    if max_reff:
        info_text += f"\nMax reff: {max_reff} Å"
    
    plt.text(0.02, 0.98, info_text, transform=plt.gca().transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.xlabel('k (Å⁻¹)')
    plt.ylabel('χ(k)')
    plt.title(f'Sum of Paths - Frame {frame}, Atom {atom_id}')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()


def plot_averaged_atoms(db: ChiDatabaseQuery,
                       frame_range: Tuple[int, int],
                       path_types: Optional[List[str]] = None,
                       max_reff: Optional[float] = None,
                       output_file: Optional[str] = None) -> None:
    """Plot averaged chi(k) across multiple atoms."""
    # Get frames in range
    frame_counts = db.get_frame_atom_counts()
    frames = [f for f in frame_counts.keys() 
              if frame_range[0] <= f <= frame_range[1]]
    
    if not frames:
        print(f"No frames found in range {frame_range[0]}-{frame_range[1]}")
        return
    
    # Query paths
    query = PathQuery(
        frames=frames,
        path_types=path_types,
        max_reff=max_reff
    )
    paths = db.query_paths(query)
    
    if not paths:
        print("No paths found matching criteria")
        return
    
    # Get path IDs
    path_ids = [p['id'] for p in paths]
    
    # Average using correct method (interpolate to uniform k-grid, sum within atoms, then average)
    averaged_data, num_atoms = db.sum_chi_within_atoms_then_average(path_ids)
    
    if averaged_data is None:
        print("No chi data found for selected paths")
        return
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    k_grid = averaged_data[:, 0]
    chi_avg = averaged_data[:, 1]
    
    plt.plot(k_grid, chi_avg, 'b-', linewidth=2)
    
    # Add information
    info_text = f"Frames {frame_range[0]}-{frame_range[1]}\n"
    info_text += f"{num_atoms} atoms averaged\n"
    info_text += f"{len(paths)} total paths"
    if path_types:
        info_text += f"\nPath types: {', '.join(path_types)}"
    if max_reff:
        info_text += f"\nMax reff: {max_reff} Å"
    
    plt.text(0.02, 0.98, info_text, transform=plt.gca().transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.xlabel('k (Å⁻¹)')
    plt.ylabel('χ(k)')
    plt.title('Averaged χ(k) - Sum within atoms, then average')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()


def plot_path_comparison(db: ChiDatabaseQuery,
                        frame: int,
                        atom_id: int,
                        path_types: List[str],
                        output_file: Optional[str] = None) -> None:
    """Compare different path types for the same atom."""
    plt.figure(figsize=(10, 6))
    
    for path_type in path_types:
        query = PathQuery(
            frames=[frame],
            atom_ids=[atom_id],
            path_types=[path_type]
        )
        paths = db.query_paths(query)
        
        if paths:
            path_ids = [p['id'] for p in paths]
            summed_data = db.sum_chi_data(path_ids)
            
            if summed_data is not None:
                k_grid = summed_data[:, 0]
                chi_sum = summed_data[:, 1]
                plt.plot(k_grid, chi_sum, label=f"{path_type} ({len(paths)} paths)", 
                        linewidth=2, alpha=0.8)
    
    plt.xlabel('k (Å⁻¹)')
    plt.ylabel('χ(k)')
    plt.title(f'Path Type Comparison - Frame {frame}, Atom {atom_id}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Plot chi(k) spectra from MD-EXAFS database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot individual paths for a specific atom
  %(prog)s database.db --mode paths --frame 180000 --atom 0 --paths 1,2,3,4

  # Plot sum of all paths for an atom
  %(prog)s database.db --mode atom --frame 180000 --atom 0

  # Plot sum of U-O paths within 3.0 Å for an atom
  %(prog)s database.db --mode atom --frame 180000 --atom 0 --path-types U-O --max-reff 3.0

  # Plot averaged chi(k) for U-U paths across all atoms
  %(prog)s database.db --mode average --start-frame 180000 --end-frame 190000 --path-types U-U

  # Compare different path types for the same atom
  %(prog)s database.db --mode compare --frame 180000 --atom 0 --path-types U-O,U-U,U-O-O
        """
    )
    
    parser.add_argument('database', type=str, help='Path to chi database file')
    parser.add_argument('--mode', choices=['paths', 'atom', 'average', 'compare'], 
                       required=True, help='Plotting mode')
    parser.add_argument('--frame', type=int, help='Frame number (for paths/atom/compare modes)')
    parser.add_argument('--atom', type=int, help='Atom ID (for paths/atom/compare modes)')
    parser.add_argument('--paths', type=str, help='Comma-separated path numbers (for paths mode)')
    parser.add_argument('--path-types', type=str, 
                       help='Comma-separated path types (e.g., U-O,U-U)')
    parser.add_argument('--max-reff', type=float, help='Maximum reff distance in Angstroms')
    parser.add_argument('--start-frame', type=int, help='Start frame (for average mode)')
    parser.add_argument('--end-frame', type=int, help='End frame (for average mode)')
    parser.add_argument('--output', type=str, help='Output PNG file (optional, shows plot if not specified)')
    
    args = parser.parse_args()
    
    # Validate arguments based on mode
    if args.mode in ['paths', 'atom', 'compare']:
        if args.frame is None or args.atom is None:
            parser.error(f"--frame and --atom are required for {args.mode} mode")
    elif args.mode == 'average':
        if args.start_frame is None or args.end_frame is None:
            parser.error("--start-frame and --end-frame are required for average mode")
    
    # Parse path types
    path_types = None
    if args.path_types:
        path_types = [pt.strip() for pt in args.path_types.split(',')]
    
    # Open database
    db_path = Path(args.database)
    if not db_path.exists():
        print(f"Error: Database not found: {db_path}")
        return 1
    
    with ChiDatabaseQuery(db_path) as db:
        if args.mode == 'paths':
            # Parse path numbers
            path_nums = None
            if args.paths:
                path_nums = [int(p.strip()) for p in args.paths.split(',')]
            plot_individual_paths(db, args.frame, args.atom, path_nums, args.output)
            
        elif args.mode == 'atom':
            plot_atom_sum(db, args.frame, args.atom, path_types, args.max_reff, args.output)
            
        elif args.mode == 'average':
            plot_averaged_atoms(db, (args.start_frame, args.end_frame), 
                              path_types, args.max_reff, args.output)
            
        elif args.mode == 'compare':
            if not path_types:
                parser.error("--path-types is required for compare mode")
            plot_path_comparison(db, args.frame, args.atom, path_types, args.output)
    
    return 0


if __name__ == "__main__":
    exit(main())
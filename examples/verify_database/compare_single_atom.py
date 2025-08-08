#!/usr/bin/env python
"""
Compare chi calculations for a single atom folder using both methods:
1. Direct summation using xraylarch
2. From database
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from larch import Interpreter
from larch.xafs import feffpath, path2chi
from md_exafs.chi_database_query import ChiDatabaseQuery, PathQuery


def sum_paths_directly(atom_dir):
    """Sum all feff paths in an atom directory using xraylarch directly."""
    # Initialize xraylarch
    larch_interp = Interpreter()
    
    # Find all feff*.dat files
    feff_files = sorted(atom_dir.glob("feff*.dat"))
    if not feff_files:
        print(f"No feff*.dat files found in {atom_dir}")
        return None
    
    print(f"Found {len(feff_files)} feff*.dat files")
    
    # Use standard k-grid
    k_standard = np.arange(0, 20.05, 0.05)
    chi_sum = np.zeros_like(k_standard)
    
    # Process each feff file
    for feff_file in feff_files:
        try:
            # Create FeffPath object
            path = feffpath(str(feff_file), _larch=larch_interp)
            
            # Calculate chi using path2chi (includes degeneracy)
            path2chi(path, k=k_standard, _larch=larch_interp)
            
            if hasattr(path, 'chi') and path.chi is not None:
                chi_sum += path.chi
                if len(feff_files) <= 10:  # Show details for small sets
                    print(f"  {feff_file.name}: deg={path.degen}, reff={path.reff:.3f}")
        except Exception as e:
            print(f"  Error processing {feff_file.name}: {e}")
    
    return np.column_stack((k_standard, chi_sum))


def get_chi_from_database(db_path, atom_dir):
    """Get summed chi for an atom from the database."""
    # Extract frame and atom_id from directory
    frame = 0
    atom_id = 0
    
    if atom_dir.name.startswith("atom_"):
        try:
            atom_id = int(atom_dir.name.split('_')[1])
        except:
            pass
    
    for parent in atom_dir.parents:
        if parent.name.startswith("frame_"):
            try:
                frame = int(parent.name.split('_')[1])
                break
            except:
                pass
    
    print(f"Querying database for frame={frame}, atom_id={atom_id}")
    
    with ChiDatabaseQuery(db_path) as db:
        # Query all paths for this atom
        query = PathQuery(frames=[frame], atom_ids=[atom_id])
        paths = db.query_paths(query)
        
        if not paths:
            print("No paths found in database!")
            return None
        
        print(f"Found {len(paths)} paths in database")
        
        # Get all path IDs and sum
        path_ids = [p['id'] for p in paths]
        summed_chi = db.sum_chi_data(path_ids)
        
        return summed_chi


def main():
    if len(sys.argv) < 2:
        print("Usage: python compare_single_atom.py <atom_folder> [database.db]")
        print("Example: python compare_single_atom.py /path/to/frame_0/atom_0")
        sys.exit(1)
    
    atom_dir = Path(sys.argv[1])
    if not atom_dir.exists():
        print(f"Error: Directory not found: {atom_dir}")
        sys.exit(1)
    
    # Database path (optional)
    if len(sys.argv) > 2:
        db_path = Path(sys.argv[2])
    else:
        db_path = Path("test_verification.db")
    
    print(f"Comparing chi calculations for: {atom_dir}")
    print("=" * 60)
    
    # Method 1: Direct summation
    print("\n1. Direct summation using xraylarch:")
    direct_chi = sum_paths_directly(atom_dir)
    
    # Method 2: From database (if available)
    db_chi = None
    if db_path.exists():
        print(f"\n2. From database ({db_path.name}):")
        db_chi = get_chi_from_database(db_path, atom_dir)
    else:
        print(f"\n2. Database not found: {db_path}")
    
    # Method 3: FEFF chi.dat (if available)
    feff_chi = None
    chi_dat_file = atom_dir / "chi.dat"
    if chi_dat_file.exists():
        print(f"\n3. FEFF chi.dat file:")
        feff_data = np.loadtxt(chi_dat_file, comments='#')
        feff_chi = feff_data
        print(f"   Loaded FEFF chi.dat")
    
    # Create comparison plots
    plt.figure(figsize=(12, 8))
    
    # Plot 1: Chi(k)
    plt.subplot(2, 2, 1)
    if direct_chi is not None:
        plt.plot(direct_chi[:, 0], direct_chi[:, 1], 'b-', label='Direct sum', linewidth=2)
    if db_chi is not None:
        plt.plot(db_chi[:, 0], db_chi[:, 1], 'g--', label='Database', linewidth=2)
    if feff_chi is not None:
        plt.plot(feff_chi[:, 0], feff_chi[:, 1], 'r:', label='FEFF chi.dat', linewidth=2)
    
    plt.xlabel('k (Å⁻¹)')
    plt.ylabel('χ(k)')
    plt.title('Chi(k) Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 20)
    
    # Plot 2: k² weighted
    plt.subplot(2, 2, 2)
    if direct_chi is not None:
        k = direct_chi[:, 0]
        plt.plot(k, k**2 * direct_chi[:, 1], 'b-', label='Direct sum', linewidth=2)
    if db_chi is not None:
        k = db_chi[:, 0]
        plt.plot(k, k**2 * db_chi[:, 1], 'g--', label='Database', linewidth=2)
    if feff_chi is not None:
        k = feff_chi[:, 0]
        plt.plot(k, k**2 * feff_chi[:, 1], 'r:', label='FEFF chi.dat', linewidth=2)
    
    plt.xlabel('k (Å⁻¹)')
    plt.ylabel('k² × χ(k)')
    plt.title('k²-weighted Chi(k)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 20)
    
    # Plot 3: Differences
    plt.subplot(2, 2, 3)
    if direct_chi is not None and feff_chi is not None:
        min_len = min(len(direct_chi), len(feff_chi))
        diff_direct = direct_chi[:min_len, 1] - feff_chi[:min_len, 1]
        plt.plot(direct_chi[:min_len, 0], diff_direct, 'b-', label='Direct - FEFF', linewidth=2)
        print(f"\nDirect vs FEFF: max diff = {np.max(np.abs(diff_direct)):.2e}")
    
    if db_chi is not None and feff_chi is not None:
        min_len = min(len(db_chi), len(feff_chi))
        diff_db = db_chi[:min_len, 1] - feff_chi[:min_len, 1]
        plt.plot(db_chi[:min_len, 0], diff_db, 'g--', label='Database - FEFF', linewidth=2)
        print(f"Database vs FEFF: max diff = {np.max(np.abs(diff_db)):.2e}")
    
    plt.xlabel('k (Å⁻¹)')
    plt.ylabel('Δχ(k)')
    plt.title('Differences from FEFF')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 20)
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
    
    # Plot 4: Zoomed k² weighted (k=2-15)
    plt.subplot(2, 2, 4)
    if direct_chi is not None:
        k = direct_chi[:, 0]
        mask = (k >= 2) & (k <= 15)
        plt.plot(k[mask], k[mask]**2 * direct_chi[mask, 1], 'b-', label='Direct sum', linewidth=2)
    if db_chi is not None:
        k = db_chi[:, 0]
        mask = (k >= 2) & (k <= 15)
        plt.plot(k[mask], k[mask]**2 * db_chi[mask, 1], 'g--', label='Database', linewidth=2)
    if feff_chi is not None:
        k = feff_chi[:, 0]
        mask = (k >= 2) & (k <= 15)
        plt.plot(k[mask], k[mask]**2 * feff_chi[mask, 1], 'r:', label='FEFF chi.dat', linewidth=2)
    
    plt.xlabel('k (Å⁻¹)')
    plt.ylabel('k² × χ(k)')
    plt.title('k²-weighted Chi(k) (zoomed)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(2, 15)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = atom_dir / "chi_comparison.png"
    plt.savefig(plot_file, dpi=150)
    print(f"\nPlot saved to: {plot_file}")
    plt.show()


if __name__ == "__main__":
    main()
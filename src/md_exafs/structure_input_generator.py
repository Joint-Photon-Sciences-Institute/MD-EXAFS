#!/usr/bin/env python3
"""
Generate structure input files for MD simulations from CIF files.

This module provides functionality to read CIF files and generate
supercells in various MD simulation formats.
"""

import argparse
import sys
from pathlib import Path
from typing import Tuple, Optional

try:
    from ase.io import read, write
    from ase.build import make_supercell
    import numpy as np
except ImportError:
    print("Error: ASE (Atomic Simulation Environment) is required.")
    print("Please install it with: pip install ase")
    sys.exit(1)


def parse_size_argument(size_str: str) -> Tuple[int, int, int]:
    """
    Parse size argument from string format "nx,ny,nz" to tuple.
    
    Args:
        size_str: String in format "nx,ny,nz" where nx, ny, nz are integers
        
    Returns:
        Tuple of three integers representing supercell dimensions
        
    Raises:
        ValueError: If format is invalid
    """
    try:
        parts = size_str.strip().split(',')
        if len(parts) != 3:
            raise ValueError("Size must have exactly 3 components")
        nx, ny, nz = map(int, parts)
        if any(n <= 0 for n in (nx, ny, nz)):
            raise ValueError("All size components must be positive integers")
        return nx, ny, nz
    except ValueError as e:
        raise ValueError(f"Invalid size format '{size_str}': {e}")


def generate_supercell(cif_file: Path, size: Tuple[int, int, int]) -> 'ase.Atoms':
    """
    Read CIF file and generate supercell.
    
    Args:
        cif_file: Path to input CIF file
        size: Tuple of (nx, ny, nz) for supercell dimensions
        
    Returns:
        ASE Atoms object representing the supercell
        
    Raises:
        FileNotFoundError: If CIF file doesn't exist
        ValueError: If CIF file cannot be read
    """
    if not cif_file.exists():
        raise FileNotFoundError(f"CIF file not found: {cif_file}")
    
    try:
        # Read the CIF file
        structure = read(str(cif_file))
    except Exception as e:
        raise ValueError(f"Failed to read CIF file: {e}")
    
    # Create supercell matrix
    nx, ny, nz = size
    supercell_matrix = np.diag([nx, ny, nz])
    
    # Generate supercell
    supercell = make_supercell(structure, supercell_matrix)
    
    return supercell


def write_cp2k_format(atoms: 'ase.Atoms', output_file: Path) -> None:
    """
    Write structure in CP2K format (no header, element x y z).
    
    Args:
        atoms: ASE Atoms object
        output_file: Path to output file
    """
    # Ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        for atom in atoms:
            symbol = atom.symbol
            x, y, z = atom.position
            # CP2K format: element x y z (scientific notation)
            f.write(f"{symbol}  {x:24.16E}  {y:24.16E}  {z:24.16E}\n")


def main():
    """Main entry point for the structure input generator."""
    parser = argparse.ArgumentParser(
        description="Generate structure input files for MD simulations from CIF files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate a 3x3x3 supercell from a CIF file in CP2K format
  md-exafs-md-input-gen --cp2k --input structure.cif --size 3,3,3 --output supercell.xyz
  
  # Generate a 2x2x2 supercell
  md-exafs-md-input-gen --cp2k --input Au.cif --size 2,2,2 --output Au_2x2x2.xyz
        """
    )
    
    # Format options (currently only CP2K)
    format_group = parser.add_mutually_exclusive_group(required=True)
    format_group.add_argument(
        '--cp2k',
        action='store_true',
        help='Output in CP2K format (no header, element x y z)'
    )
    
    # Input/output options
    parser.add_argument(
        '--input',
        type=Path,
        required=True,
        help='Input CIF file path'
    )
    parser.add_argument(
        '--size',
        type=str,
        required=True,
        help='Supercell size as "nx,ny,nz" (e.g., "3,3,3" for 3x3x3)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output file path'
    )
    
    args = parser.parse_args()
    
    try:
        # Parse size argument
        size = parse_size_argument(args.size)
        
        # Generate supercell
        print(f"Reading CIF file: {args.input}")
        supercell = generate_supercell(args.input, size)
        print(f"Generated {size[0]}x{size[1]}x{size[2]} supercell with {len(supercell)} atoms")
        
        # Write output
        if args.cp2k:
            print(f"Writing CP2K format to: {args.output}")
            write_cp2k_format(supercell, args.output)
        
        print("Structure generation completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
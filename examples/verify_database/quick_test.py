#!/usr/bin/env python
"""Quick test to verify the refactored code works correctly."""

import sys
from pathlib import Path
import numpy as np

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from md_exafs.chi_averager import _convert_feffdat_to_chi
from larch import Interpreter


def test_single_path():
    """Test conversion of a single FEFF path."""
    # Use a known test file
    test_file = Path("/home/nickj/claude-code/MD-EXAFS/reference_code_base/UO_data/atom_0/feff0001.dat")
    
    if not test_file.exists():
        print(f"Test file not found: {test_file}")
        return
    
    print(f"Testing FEFF path conversion on: {test_file}")
    print("=" * 60)
    
    # Initialize larch
    larch_interp = Interpreter()
    
    # Convert using the refactored function
    chi_data = _convert_feffdat_to_chi(test_file, larch_interp)
    
    if chi_data is not None:
        k = chi_data[:, 0]
        chi = chi_data[:, 1]
        
        print(f"✓ Conversion successful!")
        print(f"  k-grid: {k[0]:.2f} to {k[-1]:.2f} Å⁻¹ ({len(k)} points)")
        print(f"  k-step: {k[1]-k[0]:.3f} Å⁻¹")
        print(f"  chi(k=5): {chi[100]:.6f}")
        print(f"  chi(k=10): {chi[200]:.6f}")
        
        # Check if it's on standard grid
        expected_k = np.arange(0, 20.05, 0.05)
        if len(k) == len(expected_k) and np.allclose(k, expected_k):
            print(f"✓ Using standard k-grid (0-20, step 0.05)")
        else:
            print(f"⚠ Not on standard k-grid!")
    else:
        print("✗ Conversion failed!")


if __name__ == "__main__":
    test_single_path()
#!/usr/bin/env python3
"""Test script for singular absorbing atom potential fix."""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from md_exafs.feff_generator import assign_potential_indices, generate_feff_input
import numpy as np


def test_singular_atom():
    """Test case where there's only one atom of the absorbing type (e.g., Fe impurity in Al matrix)."""
    
    # Configuration for Fe impurity in Al matrix
    config = {
        "feff": {
            "absorbing_element": "Fe",
            "header": """TITLE Fe impurity in Al matrix test
EDGE L3
S02 1.0
CONTROL 1 1 1 1 1 1
PRINT 0 0 0 0 0 0
RPATH 6.5
NLEG 8
FMS 10.0 0
EXCHANGE 0 0.0 0.0 2
SCF 7.0 0 100 0.20 3
XANES 4.0 0.07 0.0"""
        }
    }
    
    # Test 1: Fe atom with only Al neighbors (singular Fe case)
    print("Test 1: Singular Fe atom in Al matrix")
    print("-" * 50)
    
    # Central Fe atom at origin, surrounded by Al atoms
    element_symbols = ["Fe", "Al", "Al", "Al", "Al", "Al", "Al"]  # 1 Fe, 6 Al neighbors
    central_element = "Fe"
    
    # Get potential indices
    pot_indices = assign_potential_indices(element_symbols, central_element, config)
    
    print(f"Element symbols: {element_symbols}")
    print(f"Potential indices: {pot_indices}")
    
    # Expected: Fe gets 0 (central), all Al get 1 (only neighbor type)
    expected = [0, 1, 1, 1, 1, 1, 1]
    assert pot_indices == expected, f"Expected {expected}, got {pot_indices}"
    print("✓ Potential indices correct for singular Fe atom")
    
    # Generate FEFF input to verify POTENTIALS section
    atoms_data = [
        (np.array([0.0, 0.0, 0.0]), pot_indices[0], element_symbols[0]),
        (np.array([2.0, 0.0, 0.0]), pot_indices[1], element_symbols[1]),
        (np.array([-2.0, 0.0, 0.0]), pot_indices[2], element_symbols[2]),
        (np.array([0.0, 2.0, 0.0]), pot_indices[3], element_symbols[3]),
        (np.array([0.0, -2.0, 0.0]), pot_indices[4], element_symbols[4]),
        (np.array([0.0, 0.0, 2.0]), pot_indices[5], element_symbols[5]),
        (np.array([0.0, 0.0, -2.0]), pot_indices[6], element_symbols[6]),
    ]
    
    feff_content = generate_feff_input(atoms_data, config)
    print("\nGenerated FEFF POTENTIALS section:")
    for line in feff_content.split('\n'):
        if 'POTENTIALS' in line:
            for i in range(10):  # Print next few lines after POTENTIALS
                print(line)
                idx = feff_content.split('\n').index(line)
                if idx + i + 1 < len(feff_content.split('\n')):
                    line = feff_content.split('\n')[idx + i + 1]
                    if line.strip() == '' or 'ATOMS' in line:
                        break
    
    # Verify only potentials 0 and 1 exist (not potential 2 for Fe neighbors)
    assert "    0   26    Fe" in feff_content, "Potential 0 (Fe central) should exist"
    assert "    1   13    Al" in feff_content, "Potential 1 (Al neighbors) should exist"
    assert "    2   26    Fe" not in feff_content, "Potential 2 (Fe neighbors) should NOT exist"
    print("✓ POTENTIALS section correct - no spurious Fe neighbor potential")
    
    print("\n" + "=" * 50)
    
    # Test 2: Fe atom with mixed neighbors including other Fe atoms
    print("\nTest 2: Fe atom with both Al and Fe neighbors")
    print("-" * 50)
    
    element_symbols = ["Fe", "Al", "Al", "Fe", "Al", "Fe"]  # Central Fe with mixed neighbors
    
    pot_indices = assign_potential_indices(element_symbols, central_element, config)
    
    print(f"Element symbols: {element_symbols}")
    print(f"Potential indices: {pot_indices}")
    
    # Expected: Fe central gets 0, Al neighbors get 1, Fe neighbors get 2
    expected = [0, 1, 1, 2, 1, 2]
    assert pot_indices == expected, f"Expected {expected}, got {pot_indices}"
    print("✓ Potential indices correct for Fe with mixed neighbors")
    
    # Generate FEFF input for mixed case
    atoms_data = [
        (np.array([0.0, 0.0, 0.0]), pot_indices[0], element_symbols[0]),
        (np.array([2.0, 0.0, 0.0]), pot_indices[1], element_symbols[1]),
        (np.array([-2.0, 0.0, 0.0]), pot_indices[2], element_symbols[2]),
        (np.array([0.0, 2.0, 0.0]), pot_indices[3], element_symbols[3]),
        (np.array([0.0, -2.0, 0.0]), pot_indices[4], element_symbols[4]),
        (np.array([0.0, 0.0, 2.0]), pot_indices[5], element_symbols[5]),
    ]
    
    feff_content = generate_feff_input(atoms_data, config)
    
    # Verify all three potentials exist in this case
    assert "    0   26    Fe" in feff_content, "Potential 0 (Fe central) should exist"
    assert "    1   13    Al" in feff_content, "Potential 1 (Al neighbors) should exist"
    assert "    2   26    Fe" in feff_content, "Potential 2 (Fe neighbors) should exist"
    print("✓ POTENTIALS section correct - includes Fe neighbor potential when needed")
    
    print("\n" + "=" * 50)
    print("\nAll tests passed! ✓")
    print("The fix correctly handles singular absorbing atoms.")
    
    # Save FEFF input files for review
    print("\n" + "=" * 50)
    print("\nSaving FEFF input files for review...")
    
    # Test 1: Singular Fe atom case
    element_symbols = ["Fe", "Al", "Al", "Al", "Al", "Al", "Al"]
    pot_indices = assign_potential_indices(element_symbols, "Fe", config)
    atoms_data = [
        (np.array([0.0, 0.0, 0.0]), pot_indices[0], element_symbols[0]),
        (np.array([2.0, 0.0, 0.0]), pot_indices[1], element_symbols[1]),
        (np.array([-2.0, 0.0, 0.0]), pot_indices[2], element_symbols[2]),
        (np.array([0.0, 2.0, 0.0]), pot_indices[3], element_symbols[3]),
        (np.array([0.0, -2.0, 0.0]), pot_indices[4], element_symbols[4]),
        (np.array([0.0, 0.0, 2.0]), pot_indices[5], element_symbols[5]),
        (np.array([0.0, 0.0, -2.0]), pot_indices[6], element_symbols[6]),
    ]
    feff_content = generate_feff_input(atoms_data, config)
    
    with open("feff_singular_Fe.inp", "w") as f:
        f.write(feff_content)
    print("✓ Saved feff_singular_Fe.inp (Fe with only Al neighbors)")
    
    # Test 2: Fe with mixed neighbors
    element_symbols = ["Fe", "Al", "Al", "Fe", "Al", "Fe"]
    pot_indices = assign_potential_indices(element_symbols, "Fe", config)
    atoms_data = [
        (np.array([0.0, 0.0, 0.0]), pot_indices[0], element_symbols[0]),
        (np.array([2.0, 0.0, 0.0]), pot_indices[1], element_symbols[1]),
        (np.array([-2.0, 0.0, 0.0]), pot_indices[2], element_symbols[2]),
        (np.array([0.0, 2.0, 0.0]), pot_indices[3], element_symbols[3]),
        (np.array([0.0, -2.0, 0.0]), pot_indices[4], element_symbols[4]),
        (np.array([0.0, 0.0, 2.0]), pot_indices[5], element_symbols[5]),
    ]
    feff_content = generate_feff_input(atoms_data, config)
    
    with open("feff_mixed_neighbors.inp", "w") as f:
        f.write(feff_content)
    print("✓ Saved feff_mixed_neighbors.inp (Fe with Al and Fe neighbors)")
    
    print("\nYou can review the generated FEFF input files:")
    print("  - feff_singular_Fe.inp: Singular Fe case (should have only potentials 0 and 1)")
    print("  - feff_mixed_neighbors.inp: Mixed case (should have potentials 0, 1, and 2)")


if __name__ == "__main__":
    test_singular_atom()
#!/usr/bin/env python3
"""Test script for user-defined POTENTIALS handling."""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from md_exafs.feff_generator import assign_potential_indices, generate_feff_input, _parse_potentials_from_header
import numpy as np


def test_user_defined_potentials():
    """Test that user-defined POTENTIALS are correctly parsed and used."""
    
    print("=" * 60)
    print("Testing User-Defined POTENTIALS Handling")
    print("=" * 60)
    
    # Test 1: Parse user-defined POTENTIALS from header
    print("\nTest 1: Parse POTENTIALS from header")
    print("-" * 50)
    
    header_with_potentials = """TITLE   EvAX_FEFF_inp
  
EDGE      K 
S02       1.000000000000

POTENTIALS
*   ipot   z [ label   l_scmt  l_fms  stoichiometry ]
    0   26    Fe     -1      -1       1
    1   13    Al     -1      -1       6
    2   26    Fe     -1      -1       2

ATOMS"""
    
    mapping = _parse_potentials_from_header(header_with_potentials)
    print(f"Parsed mapping: {mapping}")
    
    # Should have Al -> 1, Fe -> 2 (excluding potential 0)
    assert mapping == {"Al": 1, "Fe": 2}, f"Expected {{'Al': 1, 'Fe': 2}}, got {mapping}"
    print("✓ Correctly parsed user-defined POTENTIALS")
    
    # Test 2: Use user-defined POTENTIALS for singular Fe case
    print("\nTest 2: User-defined POTENTIALS with singular Fe atom")
    print("-" * 50)
    
    config_with_user_potentials = {
        "feff": {
            "absorbing_element": "Fe",
            "header": """TITLE Fe in Al with user-defined potentials
EDGE L3

POTENTIALS
*   ipot   z [ label   l_scmt  l_fms  stoichiometry ]
    0   26    Fe     -1      -1       1
    1   13    Al     -1      -1       6"""
        }
    }
    
    # Fe atom with only Al neighbors
    element_symbols = ["Fe", "Al", "Al", "Al", "Al", "Al", "Al"]
    central_element = "Fe"
    
    pot_indices = assign_potential_indices(element_symbols, central_element, config_with_user_potentials)
    print(f"Element symbols: {element_symbols}")
    print(f"Potential indices: {pot_indices}")
    
    # Should use the user-defined mapping: Fe -> 0 (central), Al -> 1
    expected = [0, 1, 1, 1, 1, 1, 1]
    assert pot_indices == expected, f"Expected {expected}, got {pot_indices}"
    print("✓ Correctly used user-defined POTENTIALS for singular Fe")
    
    # Test 3: User-defined POTENTIALS with mixed neighbors
    print("\nTest 3: User-defined POTENTIALS with Fe neighbors")
    print("-" * 50)
    
    config_with_fe_neighbors = {
        "feff": {
            "absorbing_element": "Fe",
            "header": """TITLE Fe in Al with Fe neighbors
EDGE L3

POTENTIALS
*   ipot   z [ label   l_scmt  l_fms  stoichiometry ]
    0   26    Fe     -1      -1       1
    1   13    Al     -1      -1       4
    2   26    Fe     -1      -1       2"""
        }
    }
    
    # Fe atom with both Al and Fe neighbors
    element_symbols = ["Fe", "Al", "Al", "Fe", "Al", "Fe"]
    
    pot_indices = assign_potential_indices(element_symbols, central_element, config_with_fe_neighbors)
    print(f"Element symbols: {element_symbols}")
    print(f"Potential indices: {pot_indices}")
    
    # Should use: Fe central -> 0, Al -> 1, Fe neighbors -> 2
    expected = [0, 1, 1, 2, 1, 2]
    assert pot_indices == expected, f"Expected {expected}, got {pot_indices}"
    print("✓ Correctly used user-defined POTENTIALS with Fe neighbors")
    
    # Test 4: Auto-generated POTENTIALS (no user definition)
    print("\nTest 4: Auto-generated POTENTIALS")
    print("-" * 50)
    
    config_auto_potentials = {
        "feff": {
            "absorbing_element": "Fe",
            "header": """TITLE Fe in Al auto-generated
EDGE L3
S02 1.0"""  # No POTENTIALS section
        }
    }
    
    # Fe with only Al neighbors - should auto-generate
    element_symbols = ["Fe", "Al", "Al", "Al", "Al"]
    
    pot_indices = assign_potential_indices(element_symbols, central_element, config_auto_potentials)
    print(f"Element symbols: {element_symbols}")
    print(f"Potential indices: {pot_indices}")
    
    # Auto-generated should be: Fe -> 0, Al -> 1
    expected = [0, 1, 1, 1, 1]
    assert pot_indices == expected, f"Expected {expected}, got {pot_indices}"
    print("✓ Correctly auto-generated POTENTIALS")
    
    # Test 5: Error handling - missing element in user-defined POTENTIALS
    print("\nTest 5: Error handling for missing elements")
    print("-" * 50)
    
    config_missing_element = {
        "feff": {
            "absorbing_element": "Fe",
            "header": """TITLE Missing O potential
POTENTIALS
    0   26    Fe     -1      -1       1
    1   13    Al     -1      -1       6"""
        }
    }
    
    # Fe with Al and O neighbors, but O not defined in POTENTIALS
    element_symbols = ["Fe", "Al", "Al", "O", "Al"]
    
    try:
        pot_indices = assign_potential_indices(element_symbols, central_element, config_missing_element)
        assert False, "Should have raised ValueError for missing element"
    except ValueError as e:
        print(f"✓ Correctly raised error: {e}")
        assert "Element 'O' found in neighbors but not defined in POTENTIALS" in str(e)
    
    print("\n" + "=" * 60)
    print("All tests passed! ✓")
    print("User-defined POTENTIALS handling works correctly.")


if __name__ == "__main__":
    test_user_defined_potentials()
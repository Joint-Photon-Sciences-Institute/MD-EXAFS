"""FEFF input file generation utilities."""

from typing import List, Tuple, Dict, Any
import numpy as np
from .config import get_atomic_number


def generate_feff_input(
    atoms_data: List[Tuple[np.ndarray, int, str]],
    config: Dict[str, Any]
) -> str:
    """
    Generate FEFF input file content.
    
    Args:
        atoms_data: List of (position, potential_index, element_symbol) tuples
        config: Configuration dictionary
        
    Returns:
        Complete FEFF input file content as string
    """
    # Get header from config
    header = config["feff"]["header"]
    
    # Check if POTENTIALS section is already in the header
    if "POTENTIALS" in header:
        # Header already contains POTENTIALS section, just add ATOMS section
        atoms_section = _build_atoms_section(atoms_data)
        feff_content = f"{header}\n\n{atoms_section}\nEND\n"
    else:
        # Build potentials section dynamically
        potentials_section = _build_potentials_section(atoms_data, config)
        atoms_section = _build_atoms_section(atoms_data)
        feff_content = f"{header}\n\n{potentials_section}\n\n{atoms_section}\nEND\n"
    
    return feff_content


def _build_potentials_section(
    atoms_data: List[Tuple[np.ndarray, int, str]],
    config: Dict[str, Any]
) -> str:
    """Build the POTENTIALS section of FEFF input."""
    # Get unique elements and their potential indices
    unique_potentials = {}
    for _, pot_idx, element in atoms_data:
        if pot_idx not in unique_potentials:
            unique_potentials[pot_idx] = element
    
    # Sort by potential index
    sorted_potentials = sorted(unique_potentials.items())
    
    # Build potentials lines
    lines = ["POTENTIALS", "*   ipot   z [ label   l_scmt  l_fms  stoichiometry ]"]
    
    for pot_idx, element in sorted_potentials:
        atomic_num = get_atomic_number(element)
        # Count atoms with this potential (for stoichiometry)
        count = sum(1 for _, p, _ in atoms_data if p == pot_idx)
        lines.append(f"    {pot_idx}   {atomic_num}    {element}     -1      -1       {count}")
    
    return "\n".join(lines)


def _build_atoms_section(atoms_data: List[Tuple[np.ndarray, int, str]]) -> str:
    """Build the ATOMS section of FEFF input."""
    lines = ["ATOMS"]
    
    for pos, pot_idx, element in atoms_data:
        lines.append(
            f"{pos[0]:.6f}\t{pos[1]:.6f}\t{pos[2]:.6f}\t{pot_idx}\t{element}"
        )
    
    return "\n".join(lines)


def assign_potential_indices(
    element_symbols: List[str],
    central_element: str,
    config: Dict[str, Any]
) -> List[int]:
    """
    Assign potential indices to atoms based on element type.
    
    FEFF Convention:
    - Potential 0: Always the central (absorbing) atom
    - Potentials 1, 2, 3, ...: Different atomic species
    
    The absorbing element always gets potential 1 when it appears as a neighbor.
    Other elements are assigned sequential potentials starting from 2.
    
    Args:
        element_symbols: List of element symbols for all atoms
        central_element: Element symbol of the central atom
        config: Configuration dictionary
        
    Returns:
        List of potential indices corresponding to each atom
    """
    absorbing_element = config["feff"]["absorbing_element"]
    
    # Build the species mapping
    # First, get all unique elements except the absorbing element
    unique_elements = sorted(set(element_symbols) - {absorbing_element})
    
    # Create potential mapping
    # Absorbing element always gets index 1 when it's a neighbor
    species_to_potential = {absorbing_element: 1}
    
    # Assign sequential indices to other elements
    next_potential = 2
    for element in unique_elements:
        if element not in species_to_potential:
            species_to_potential[element] = next_potential
            next_potential += 1
    
    # Now assign potential indices
    potential_indices = []
    for i, element in enumerate(element_symbols):
        if i == 0:  # Central atom
            potential_indices.append(0)
        else:  # Neighbor atom
            potential_indices.append(species_to_potential[element])
    
    return potential_indices
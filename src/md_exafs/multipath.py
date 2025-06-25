"""Multipath analysis for FEFF calculations."""

import re
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

try:
    from larch import Interpreter
    from larch.xafs import feffpath
    LARCH_AVAILABLE = True
except ImportError:
    LARCH_AVAILABLE = False


def parse_feff_header(feff_file: Path) -> Dict[str, Any]:
    """
    Parse FEFF path file header to extract path characteristics.
    
    Returns:
        Dictionary containing:
        - nleg: Number of legs in path
        - deg: Degeneracy
        - reff: Effective path length (Angstroms)
        - atoms: List of atom types in path (e.g., ['U', 'O', 'U'])
        - path_type: Path type string (e.g., 'uo', 'uu', 'uoo')
    """
    with open(feff_file, 'r') as f:
        lines = f.readlines()
    
    # Parse header information
    info = {}
    atom_types = []
    in_atom_list = False
    
    for i, line in enumerate(lines):
        # Look for the line with nleg, deg, reff info
        if 'nleg, deg, reff,' in line:
            # The values are at the beginning of the same line
            parts = line.strip().split()
            if len(parts) >= 4:
                info['nleg'] = int(parts[0])
                info['deg'] = float(parts[1])
                info['reff'] = float(parts[2])
            in_atom_list = True
        # Parse atom positions
        elif in_atom_list and not line.strip().startswith('k'):
            # Skip header line
            if 'pot at#' in line:
                continue
            parts = line.strip().split()
            # Atom lines have format: x y z pot at# label [optional comments]
            if len(parts) >= 6:
                try:
                    # Check if first element is a number (x coordinate)
                    float(parts[0])
                    # Check if the 4th and 5th elements are integers (pot and at#)
                    int(parts[3])
                    int(parts[4])
                    # The 6th element is the atom label
                    atom_label = parts[5].upper()
                    atom_types.append(atom_label)
                except (ValueError, IndexError):
                    continue
        elif line.strip().startswith('k') and 'real' in line:
            # End of atom list
            break
    
    # Build path type string with full notation (e.g., 'U-O', 'U-U', 'U-O-O')
    # Include all atoms in the path to make it clear
    info['atoms'] = atom_types
    
    if len(atom_types) >= 2:
        # For FEFF paths, we need to handle the fact that paths are stored differently
        # A 2-leg path (single scattering) has 2 atoms: absorber and scatterer
        # A 3-leg path (double scattering) has 3 atoms: absorber, scatterer1, scatterer2
        # The path doesn't explicitly list the return to absorber
        
        # Build the path type showing the full scattering sequence
        info['path_type'] = '-'.join(atom_types)
    else:
        info['path_type'] = 'unknown'
    
    return info


def characterize_paths_in_atom_folder(atom_dir: Path) -> Dict[int, Dict[str, Any]]:
    """
    Characterize all FEFF paths in an atom folder.
    
    Returns:
        Dictionary mapping path number to path characteristics
    """
    path_info = {}
    
    # Find all feff####.dat files
    feff_files = sorted(atom_dir.glob("feff*.dat"))
    
    for feff_file in feff_files:
        # Extract path number
        match = re.match(r'feff(\d{4})\.dat', feff_file.name)
        if match:
            path_num = int(match.group(1))
            try:
                info = parse_feff_header(feff_file)
                path_info[path_num] = info
            except Exception as e:
                print(f"Error parsing {feff_file}: {e}")
    
    return path_info


def filter_paths_by_criteria(
    path_info: Dict[int, Dict[str, Any]],
    path_types: Optional[List[str]] = None,
    max_distance: Optional[float] = None
) -> List[int]:
    """
    Filter paths based on type and distance criteria.
    
    Args:
        path_info: Dictionary of path characteristics
        path_types: List of path types to include (e.g., ['uo', 'uoo'])
        max_distance: Maximum effective path distance (reff)
    
    Returns:
        List of path numbers that match criteria
    """
    matching_paths = []
    
    for path_num, info in path_info.items():
        # Check path type
        if path_types and info.get('path_type') not in path_types:
            continue
        
        # Check distance
        if max_distance and info.get('reff', float('inf')) > max_distance:
            continue
        
        matching_paths.append(path_num)
    
    return sorted(matching_paths)




def process_atom_folder_multipath(
    atom_dir: Path,
    path_types: Optional[List[str]],
    max_distance: Optional[float],
    larch_interp: Optional[Any] = None
) -> Optional[np.ndarray]:
    """
    Process an atom folder using multipath criteria.
    
    Returns:
        Summed chi(k) data for paths matching criteria
    """
    from .chi_averager import _convert_feffdat_to_chi, _sum_chi_data
    
    # Characterize all paths in this atom folder
    path_info = characterize_paths_in_atom_folder(atom_dir)
    
    # Filter paths based on criteria
    matching_paths = filter_paths_by_criteria(path_info, path_types, max_distance)
    
    if not matching_paths:
        return None
    
    # Initialize larch if not provided
    if larch_interp is None and LARCH_AVAILABLE:
        larch_interp = Interpreter()
    
    # Convert and sum matching paths
    chi_data_list = []
    
    for path_num in matching_paths:
        feff_file = atom_dir / f"feff{path_num:04d}.dat"
        if feff_file.exists():
            chi_data = _convert_feffdat_to_chi(feff_file, larch_interp)
            if chi_data is not None:
                chi_data_list.append(chi_data)
    
    if chi_data_list:
        return _sum_chi_data(chi_data_list)
    
    return None


def process_atom_folder_wrapper(args: Tuple) -> Tuple[Path, Optional[np.ndarray]]:
    """Wrapper function for multiprocessing."""
    atom_dir, path_types, max_distance = args
    try:
        result = process_atom_folder_multipath(atom_dir, path_types, max_distance)
        return atom_dir, result
    except Exception as e:
        print(f"Error processing {atom_dir}: {e}")
        return atom_dir, None
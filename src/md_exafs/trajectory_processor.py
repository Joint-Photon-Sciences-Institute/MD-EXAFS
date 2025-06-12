"""Process MD trajectories to generate FEFF input files."""

import os
import random
import argparse
from pathlib import Path
from typing import Dict, Any, List, Tuple

import numpy as np
from ovito.io import import_file
from ovito.data import CutoffNeighborFinder, SimulationCell
from tqdm import tqdm

from .config import load_config, ConfigError
from .feff_generator import generate_feff_input, assign_potential_indices


def process_trajectory(config: Dict[str, Any]) -> None:
    """
    Process trajectory and generate FEFF input files based on configuration.
    
    Args:
        config: Configuration dictionary from TOML file
    """
    # Extract configuration values
    system = config["system"]
    processing = config["processing"]
    
    xyz_file = system["trajectory_file"]
    base_dir = Path(system["output_directory"])
    
    frame_range = range(
        processing["start_frame"], 
        processing["end_frame"], 
        processing["frame_step"]
    )
    num_atoms = processing["atoms_per_frame"]
    
    print(f"Current working directory: {os.getcwd()}")
    print(f"Will save files in: {base_dir.absolute()}")
    print(f"Reading trajectory from: {xyz_file}")
    
    # Create base directory
    base_dir.mkdir(parents=True, exist_ok=True)
    
    # Import trajectory and check number of frames
    pipeline = import_file(xyz_file)
    num_frames = pipeline.source.num_frames
    
    # Convert frame_range to list and validate
    frame_list = list(frame_range)
    max_frame = max(frame_list)
    min_frame = min(frame_list)
    
    if max_frame >= num_frames:
        raise ValueError(f"End frame {max_frame} exceeds the number of frames in trajectory ({num_frames})")
    if min_frame < 0:
        raise ValueError(f"Start frame {min_frame} cannot be negative")
    
    print(f"Trajectory contains {num_frames} frames")
    print(f"Processing {len(frame_list)} frames from {min_frame} to {max_frame}")
    
    # Get central atom type from config
    central_atom_type = _get_central_atom_type(config)
    
    total_files = 0
    total_iterations = len(frame_list) * num_atoms
    
    # Process frames
    with tqdm(total=total_iterations, desc="Processing FEFF inputs") as pbar:
        for frame in frame_list:
            data = pipeline.compute(frame)
            
            # Process central atoms
            central_indices = _get_central_atom_indices(
                data, central_atom_type, num_atoms, config
            )
            
            for atom_idx, central_index in enumerate(central_indices):
                dir_path = _create_directory_structure(
                    base_dir, frame, atom_idx
                )
                feff_path = dir_path / "feff.inp"
                
                # Generate FEFF input
                feff_content = _generate_single_feff_input(
                    data, central_index, config
                )
                
                with open(feff_path, 'w') as f:
                    f.write(feff_content)
                
                total_files += 1
                pbar.update(1)
    
    print(f"\nFinished! Created {total_files} FEFF input files in {base_dir}")
    print(f"Total frames processed: {len(frame_list)}")
    print(f"FEFF inputs per frame: {num_atoms}")



def _get_central_atom_type(config: Dict[str, Any]) -> int:
    """Determine the central atom type from configuration."""
    # For now, use the first atom type in the atoms section
    # In future, this could be made configurable
    atoms = config["atoms"]
    # Get the atom with the lowest type ID (typically the central atom)
    return min(atoms.values())


def _get_central_atom_indices(
    data: Any, 
    central_atom_type: int, 
    num_atoms: int,
    config: Dict[str, Any]
) -> List[int]:
    """Get indices of central atoms to process."""
    # Find all atoms of the central type
    atom_types = data.particles.particle_types[:]
    central_indices = np.where(atom_types == central_atom_type)[0]
    
    if len(central_indices) == 0:
        raise ValueError(f"No atoms of type {central_atom_type} found in trajectory")
    
    # Randomly select the requested number of atoms
    if num_atoms > len(central_indices):
        raise ValueError(
            f"Requested {num_atoms} central atoms but only "
            f"{len(central_indices)} atoms of type {central_atom_type} found"
        )
    
    return random.sample(list(central_indices), num_atoms)


def _create_directory_structure(
    base_dir: Path, 
    frame: int, 
    atom_idx: int
) -> Path:
    """Create nested directory structure for FEFF input."""
    dir_path = base_dir / f"frame_{frame}" / f"atom_{atom_idx}"
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def _set_simulation_cell(data: Any, config: Dict[str, Any]) -> SimulationCell:
    """Set up the simulation cell with lattice vectors and PBC from config."""
    lattice = config["lattice"]
    
    # Create the cell matrix
    cell_matrix = np.array([
        lattice["a"] + [0.0],
        lattice["b"] + [0.0],
        lattice["c"] + [0.0]
    ])
    
    # Set the simulation cell
    data.cell = SimulationCell(matrix=cell_matrix)
    data.cell.pbc = tuple(lattice["pbc"])
    
    return data.cell


def _wrap_position(pos: np.ndarray, cell: SimulationCell) -> np.ndarray:
    """Wrap atomic position into primary cell."""
    # Convert to fractional coordinates
    inv_cell = np.linalg.inv(cell.matrix[:3,:3])
    frac = np.dot(inv_cell, pos)
    
    # Wrap to [0,1)
    frac = frac % 1.0
    
    # Convert back to Cartesian
    return np.dot(cell.matrix[:3,:3], frac)


def _get_minimum_image_vector(
    pos1: np.ndarray, 
    pos2: np.ndarray, 
    cell: SimulationCell
) -> np.ndarray:
    """Get minimum image vector between two positions."""
    # Convert to fractional coordinates
    inv_cell = np.linalg.inv(cell.matrix[:3,:3])
    frac1 = np.dot(inv_cell, pos1)
    frac2 = np.dot(inv_cell, pos2)
    
    # Get difference and apply minimum image convention
    diff = frac1 - frac2
    diff = diff - np.round(diff)
    
    # Convert back to Cartesian
    return np.dot(cell.matrix[:3,:3], diff)


def _generate_single_feff_input(
    frame_data: Any,
    central_atom_index: int,
    config: Dict[str, Any]
) -> str:
    """Generate FEFF input content for a specific central atom."""
    # Set up simulation cell with PBC
    cell = _set_simulation_cell(frame_data, config)
    
    # Get configuration values
    cutoff_radius = config["processing"]["cutoff_radius"]
    atoms_config = config["atoms"]
    
    # Create neighbor finder with PBC enabled
    finder = CutoffNeighborFinder(cutoff_radius, frame_data)
    finder.pbc = True
    
    positions = frame_data.particles.positions[:]
    types = frame_data.particles.particle_types[:]
    
    # Get central atom position and wrap it into primary cell
    central_pos = _wrap_position(positions[central_atom_index], cell)
    central_type = types[central_atom_index]
    
    # Get element symbol for central atom
    central_element = _get_element_from_type(central_type, atoms_config)
    
    # Start with central atom (always potential 0)
    atoms_data = [(central_pos, 0, central_element)]
    element_symbols = [central_element]
    
    # Find neighbors and wrap their positions relative to central atom
    for neigh in finder.find(central_atom_index):
        # Get neighbor position relative to central atom using PBC-aware delta vector
        rel_pos = central_pos + neigh.delta
        neighbor_type = types[neigh.index]
        neighbor_element = _get_element_from_type(neighbor_type, atoms_config)
        
        atoms_data.append((rel_pos, -1, neighbor_element))  # -1 as placeholder
        element_symbols.append(neighbor_element)
    
    # Assign potential indices
    potential_indices = assign_potential_indices(element_symbols, central_element, config)
    
    # Update atoms_data with correct potential indices
    atoms_data_final = []
    for i, (pos, _, element) in enumerate(atoms_data):
        atoms_data_final.append((pos, potential_indices[i], element))
    
    # Generate FEFF input
    return generate_feff_input(atoms_data_final, config)


def _get_element_from_type(atom_type: int, atoms_config: Dict[str, int]) -> str:
    """Get element symbol from atom type ID."""
    # Reverse the mapping
    type_to_element = {v: k for k, v in atoms_config.items()}
    
    if atom_type not in type_to_element:
        raise ValueError(f"Unknown atom type: {atom_type}")
    
    return type_to_element[atom_type]


def main():
    """Command-line interface for trajectory processing."""
    parser = argparse.ArgumentParser(
        description='Generate FEFF input files from MD trajectory.'
    )
    parser.add_argument(
        'config', 
        type=str, 
        help='Path to TOML configuration file'
    )
    parser.add_argument(
        '--create-tar', 
        action='store_true',
        help='Create a tar.gz archive of the output directory'
    )
    
    args = parser.parse_args()
    
    # Load configuration
    try:
        config = load_config(args.config)
    except ConfigError as e:
        print(f"Configuration error: {e}")
        return 1
    except Exception as e:
        print(f"Failed to load configuration: {e}")
        return 1
    
    # Process trajectory
    try:
        process_trajectory(config)
        
        # Create tar archive if requested
        if args.create_tar:
            output_dir = config["system"]["output_directory"]
            _create_tarfile(output_dir)
            
    except Exception as e:
        print(f"Error processing trajectory: {e}")
        return 1
    
    return 0


def _create_tarfile(source_dir: str) -> None:
    """Create a tar.gz archive from the source directory."""
    import tarfile
    
    output_filename = f"{source_dir}.tar.gz"
    print(f"\nCreating tar.gz archive: {output_filename}")
    
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
    
    print(f"Archive created successfully at: {output_filename}")


if __name__ == "__main__":
    exit(main())
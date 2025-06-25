"""Average chi(k) data from FEFF calculations."""

import os
import argparse
from pathlib import Path
from typing import List, Tuple, Optional, Any, Dict

import numpy as np
from tqdm import tqdm

try:
    from larch import Interpreter
    from larch.xafs import feffpath
    LARCH_AVAILABLE = True
except ImportError:
    LARCH_AVAILABLE = False

from .config import load_config, ConfigError
from .multipath import (
    process_atom_folder_multipath,
    process_atom_folder_wrapper
)


def average_chi_files(
    input_dir: str,
    frame_range: Tuple[int, int],
    output_file: str,
    paths: Optional[List[int]] = None,
    multipath_config: Optional[Dict[str, Any]] = None
) -> None:
    """
    Average chi.dat files from FEFF calculations within specified frame range.
    
    Args:
        input_dir: Directory containing FEFF calculation results
        frame_range: Tuple of (start_frame, end_frame) to include
        output_file: Path to save averaged chi data
        paths: Optional list of path numbers to process (e.g., [1, 2, 3, 4, 5, 6, 7, 8])
               If provided, will:
               - Convert feffxxxx.dat files to chi(k)
               - Sum paths within each atom folder (saved as chi_partial_0.dat)
               - Average the sums across all atoms (saved as output_file)
        multipath_config: Optional multipath configuration dictionary with keys:
                         - paths: List of path types (e.g., ["U-O", "U-O-O"])
                         - max_distance: Maximum path distance in Angstroms
                         - num_processes: Number of parallel processes
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        raise ValueError(f"Input directory not found: {input_dir}")
    
    start_frame, end_frame = frame_range
    
    # Check if we're processing paths or existing chi.dat files
    if multipath_config:
        if not LARCH_AVAILABLE:
            raise ImportError("xraylarch is required for processing FEFF paths. "
                            "Please install it with: conda install -c conda-forge xraylarch")
        
        # Process using multipath configuration
        _average_multipath(input_path, frame_range, output_file, multipath_config)
    elif paths:
        if not LARCH_AVAILABLE:
            raise ImportError("xraylarch is required for processing FEFF paths. "
                            "Please install it with: conda install -c conda-forge xraylarch")
        
        # Process FEFF paths
        _average_feff_paths(input_path, frame_range, output_file, paths)
    else:
        # Original behavior: average existing chi.dat files
        chi_files = _find_chi_files(input_path, start_frame, end_frame)
        
        if not chi_files:
            print(f"No chi.dat files found in frame range {start_frame}-{end_frame}")
            return
        
        print(f"Found {len(chi_files)} chi.dat files in the frame range")
        
        # Load and average the data
        all_data = []
        
        for filepath in tqdm(chi_files, desc="Processing chi.dat files"):
            try:
                # Load chi.dat file (k and chi(k) columns)
                data = np.loadtxt(filepath, comments="#", usecols=(0, 1))
                all_data.append(data)
            except Exception as e:
                print(f"Error reading {filepath}: {e}")
                continue
        
        if not all_data:
            print("No valid chi.dat files could be read")
            return
        
        # Average the data
        averaged_data = _average_chi_data(all_data)
        
        # Save the averaged data
        try:
            np.savetxt(output_file, averaged_data, fmt='%.6f')
            print(f"Averaged data saved to {output_file}")
            print(f"Averaged over {len(all_data)} files")
        except Exception as e:
            print(f"Error saving averaged data: {e}")
            raise


def _find_chi_files(
    root_dir: Path,
    start_frame: int,
    end_frame: int
) -> List[Path]:
    """Find all chi.dat files within the specified frame range."""
    chi_files = []
    
    # Walk through directory structure
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename == "chi.dat":
                # Extract frame number from path
                path_parts = Path(dirpath).parts
                frame_num = _extract_frame_number(path_parts)
                
                if frame_num is not None and start_frame <= frame_num < end_frame:
                    chi_files.append(Path(dirpath) / filename)
    
    # Sort by frame number for consistency
    chi_files.sort(key=lambda p: _extract_frame_number(p.parts) or 0)
    
    return chi_files


def _extract_frame_number(path_parts: Tuple[str, ...]) -> Optional[int]:
    """Extract frame number from path components."""
    for part in path_parts:
        if part.startswith("frame_"):
            try:
                return int(part.replace("frame_", ""))
            except ValueError:
                continue
    return None


def _average_chi_data(data_list: List[np.ndarray]) -> np.ndarray:
    """Average multiple chi(k) datasets."""
    # Determine the length of the shortest data array
    min_length = min(len(data) for data in data_list)
    
    # Truncate all arrays to the minimum length
    truncated_data = [data[:min_length] for data in data_list]
    
    # Stack and average
    stacked_data = np.stack(truncated_data)
    averaged_data = np.mean(stacked_data, axis=0)
    
    return averaged_data


def _sum_chi_data(data_list: List[np.ndarray]) -> np.ndarray:
    """Sum multiple chi(k) datasets."""
    # Determine the length of the shortest data array
    min_length = min(len(data) for data in data_list)
    
    # Truncate all arrays to the minimum length
    truncated_data = [data[:min_length] for data in data_list]
    
    # Use the k-grid from the first dataset
    k_grid = truncated_data[0][:, 0]
    
    # Sum only the chi values (column 1), not the k values
    chi_sum = np.zeros(min_length)
    for data in truncated_data:
        chi_sum += data[:, 1]
    
    # Return k-grid with summed chi
    return np.column_stack((k_grid, chi_sum))


def _convert_feffdat_to_chi(feff_file: Path, larch_interp: Any) -> Optional[np.ndarray]:
    """Convert a single feff####.dat file to chi(k) using xraylarch."""
    try:
        # Create a FeffPath object
        path = feffpath(str(feff_file), _larch=larch_interp)
        
        # The path needs to be properly initialized first
        path.create_path_params()
        
        # Get k grid from the path data itself
        if hasattr(path, 'k') and path.k is not None:
            k = np.array(path.k)
        else:
            # If k is not available, use the _feffdat data
            if hasattr(path, '_feffdat') and hasattr(path._feffdat, 'k'):
                k = np.array(path._feffdat.k)
            else:
                print(f"Warning: No k grid found in path data for {feff_file}")
                return None
        
        # Calculate chi for this k grid
        path._calc_chi(k)
        
        # Get the chi array
        if hasattr(path, 'chi') and path.chi is not None:
            chi = np.array(path.chi)
        else:
            print(f"Warning: Could not calculate chi for {feff_file}")
            return None
        
        # Ensure k and chi have the same length
        if len(k) != len(chi):
            min_len = min(len(k), len(chi))
            k = k[:min_len]
            chi = chi[:min_len]
        
        return np.column_stack((k, chi))
    except Exception as e:
        print(f"Error converting {feff_file}: {e}")
        return None


def _process_paths_in_atom_folder(
    atom_dir: Path,
    path_numbers: List[int],
    larch_interp: Any
) -> Optional[np.ndarray]:
    """Process specific feff paths in an atom folder and return summed chi."""
    chi_data_list = []
    
    for path_num in path_numbers:
        feff_file = atom_dir / f"feff{path_num:04d}.dat"
        if feff_file.exists():
            chi_data = _convert_feffdat_to_chi(feff_file, larch_interp)
            if chi_data is not None:
                chi_data_list.append(chi_data)
                # Save individual chi file
                chi_file = atom_dir / f"chi{path_num:04d}.dat"
                np.savetxt(chi_file, chi_data, fmt='%.6f')
        else:
            print(f"Warning: {feff_file} not found")
    
    if chi_data_list:
        # Sum the chi data for this atom
        summed_data = _sum_chi_data(chi_data_list)
        return summed_data
    return None


def _find_atom_directories(
    root_dir: Path,
    start_frame: int,
    end_frame: int
) -> List[Path]:
    """Find all atom directories within the specified frame range."""
    atom_dirs = []
    
    for dirpath, _, _ in os.walk(root_dir):
        path_obj = Path(dirpath)
        path_parts = path_obj.parts
        
        # Check if this is an atom directory within a frame
        for i, part in enumerate(path_parts):
            if part.startswith("frame_"):
                try:
                    frame_num = int(part.replace("frame_", ""))
                    if start_frame <= frame_num < end_frame:
                        # Check if this is an atom directory
                        if i + 1 < len(path_parts) and path_parts[i + 1].startswith("atom_"):
                            atom_dirs.append(path_obj)
                            break
                except ValueError:
                    continue
    
    return sorted(atom_dirs)


def _average_feff_paths(
    input_dir: Path,
    frame_range: Tuple[int, int],
    output_file: str,
    paths: List[int]
) -> None:
    """Process FEFF path files, sum within atoms, then average across atoms."""
    # Initialize xraylarch
    larch_interp = Interpreter()
    
    start_frame, end_frame = frame_range
    
    # Find all atom directories in the frame range
    atom_dirs = _find_atom_directories(input_dir, start_frame, end_frame)
    
    if not atom_dirs:
        print(f"No atom directories found in frame range {start_frame}-{end_frame}")
        return
    
    print(f"Found {len(atom_dirs)} atom directories in frame range {start_frame}-{end_frame}")
    print(f"Processing paths: {paths}")
    
    all_partial_chi = []
    
    # Process each atom directory
    for atom_dir in tqdm(atom_dirs, desc="Processing atom directories"):
        # This returns the SUM of paths for this atom
        summed_chi = _process_paths_in_atom_folder(atom_dir, paths, larch_interp)
        if summed_chi is not None:
            # Save summed chi for this atom (chi_partial_0.dat is the sum of paths)
            partial_file = atom_dir / "chi_partial_0.dat"
            np.savetxt(partial_file, summed_chi, fmt='%.6f', 
                      header="k(A^-1)  chi(k) - Sum of selected paths for this atom")
            all_partial_chi.append(summed_chi)
    
    if all_partial_chi:
        # Average all partial chi files (average of sums across atoms)
        final_averaged = _average_chi_data(all_partial_chi)
        
        # Interpolate to standard k-grid (0 to 20 with step 0.05)
        from scipy.interpolate import interp1d
        
        # Create standard k-grid
        k_standard = np.arange(0, 20.05, 0.05)
        
        # Extract k and chi from averaged data
        k_orig = final_averaged[:, 0]
        chi_orig = final_averaged[:, 1]
        
        # Create interpolation function
        f_interp = interp1d(k_orig, chi_orig, kind='cubic', 
                           bounds_error=False, fill_value='extrapolate')
        
        # Interpolate to standard grid
        chi_interp = f_interp(k_standard)
        
        # Create final data array
        final_data = np.column_stack((k_standard, chi_interp))
        
        # Save the final averaged and interpolated data
        np.savetxt(output_file, final_data, fmt='%.6f',
                  header="k(A^-1)  chi(k) - Average of path sums across all atoms (interpolated to k=0:20:0.05)")
        print(f"Final averaged data saved to {output_file}")
        print(f"Averaged over {len(all_partial_chi)} atom folders")
        print(f"Interpolated to standard k-grid: 0 to 20 Å⁻¹ with step 0.05")
        print(f"Note: Each atom's chi_partial_0.dat contains the SUM of paths {paths}")
        print(f"      The final output is the AVERAGE of these sums across atoms")
    else:
        print("No data to average")


def _average_multipath(
    input_dir: Path,
    frame_range: Tuple[int, int],
    output_file: str,
    multipath_config: Dict[str, Any]
) -> None:
    """Process FEFF paths using multipath configuration with optional multiprocessing."""
    from concurrent.futures import ProcessPoolExecutor, as_completed
    
    # Extract configuration
    path_types = multipath_config.get('paths', [])
    max_distance = multipath_config.get('max_distance', None)
    num_processes = multipath_config.get('num_processes', 1)
    
    start_frame, end_frame = frame_range
    
    # Find all atom directories in the frame range
    atom_dirs = _find_atom_directories(input_dir, start_frame, end_frame)
    
    if not atom_dirs:
        print(f"No atom directories found in frame range {start_frame}-{end_frame}")
        return
    
    print(f"Found {len(atom_dirs)} atom directories in frame range {start_frame}-{end_frame}")
    print(f"Processing with multipath criteria:")
    print(f"  Path types: {path_types if path_types else 'all'}")
    print(f"  Max distance: {max_distance if max_distance else 'unlimited'} Å")
    print(f"  Using {num_processes} process(es)")
    
    all_partial_chi = []
    
    if num_processes > 1:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            # Submit all tasks
            future_to_dir = {
                executor.submit(process_atom_folder_wrapper, (atom_dir, path_types, max_distance)): atom_dir 
                for atom_dir in atom_dirs
            }
            
            # Process results as they complete
            for future in tqdm(as_completed(future_to_dir), total=len(atom_dirs), 
                              desc="Processing atom directories"):
                atom_dir = future_to_dir[future]
                try:
                    _, summed_chi = future.result()
                    if summed_chi is not None:
                        # Save summed chi for this atom
                        partial_file = atom_dir / "chi_multipath.dat"
                        np.savetxt(partial_file, summed_chi, fmt='%.6f',
                                  header="k(A^-1)  chi(k) - Sum of paths matching multipath criteria")
                        all_partial_chi.append(summed_chi)
                except Exception as e:
                    print(f"Error processing {atom_dir}: {e}")
    else:
        # Sequential processing
        larch_interp = Interpreter()
        
        for atom_dir in tqdm(atom_dirs, desc="Processing atom directories"):
            summed_chi = process_atom_folder_multipath(atom_dir, path_types, max_distance, larch_interp)
            if summed_chi is not None:
                # Save summed chi for this atom
                partial_file = atom_dir / "chi_multipath.dat"
                np.savetxt(partial_file, summed_chi, fmt='%.6f',
                          header="k(A^-1)  chi(k) - Sum of paths matching multipath criteria")
                all_partial_chi.append(summed_chi)
    
    if all_partial_chi:
        # Average all partial chi files
        final_averaged = _average_chi_data(all_partial_chi)
        
        # Interpolate to standard k-grid
        from scipy.interpolate import interp1d
        
        # Create standard k-grid
        k_standard = np.arange(0, 20.05, 0.05)
        
        # Extract k and chi from averaged data
        k_orig = final_averaged[:, 0]
        chi_orig = final_averaged[:, 1]
        
        # Create interpolation function
        f_interp = interp1d(k_orig, chi_orig, kind='cubic',
                           bounds_error=False, fill_value='extrapolate')
        
        # Interpolate to standard grid
        chi_interp = f_interp(k_standard)
        
        # Create final data array
        final_data = np.column_stack((k_standard, chi_interp))
        
        # Save the final averaged and interpolated data
        np.savetxt(output_file, final_data, fmt='%.6f',
                  header="k(A^-1)  chi(k) - Average of multipath sums across all atoms (interpolated to k=0:20:0.05)")
        print(f"Final averaged data saved to {output_file}")
        print(f"Averaged over {len(all_partial_chi)} atom folders")
        print(f"Interpolated to standard k-grid: 0 to 20 Å⁻¹ with step 0.05")
    else:
        print("No data to average")


def _parse_paths_arg(paths_str: str) -> Optional[List[int]]:
    """Parse paths argument like '[1,8]' to list of integers."""
    if not paths_str:
        return None
    try:
        # Remove brackets and whitespace
        paths_str = paths_str.strip('[]')
        path_nums = [int(x.strip()) for x in paths_str.split(',')]
        if len(path_nums) == 2 and path_nums[1] > path_nums[0]:
            # Range format [start, end]
            return list(range(path_nums[0], path_nums[1] + 1))
        else:
            # Individual paths
            return path_nums
    except ValueError:
        raise ValueError(f"Invalid paths format: {paths_str}. "
                        "Use format like '[1,8]' for range or '[1,3,5,7]' for specific paths")


def main():
    """Command-line interface for chi averaging."""
    parser = argparse.ArgumentParser(
        description='Average chi(k) data from FEFF calculations.'
    )
    
    # Support both TOML config and CLI arguments
    parser.add_argument(
        '--config',
        type=str,
        help='Path to TOML configuration file for averaging'
    )
    
    # CLI argument alternatives
    parser.add_argument(
        '--input-dir',
        type=str,
        help='Directory containing FEFF calculation results'
    )
    parser.add_argument(
        '--start',
        type=int,
        help='Start frame for averaging'
    )
    parser.add_argument(
        '--end',
        type=int,
        help='End frame for averaging'
    )
    parser.add_argument(
        '--output',
        type=str,
        help='Output file for averaged chi data'
    )
    parser.add_argument(
        '--paths',
        type=str,
        help='Path numbers to process, e.g., "[1,8]" for paths 1-8 or "[1,3,5,7]" for specific paths'
    )
    
    args = parser.parse_args()
    
    # Determine configuration source
    if args.config:
        # Load from TOML
        try:
            config = load_config(args.config)
            if "averaging" not in config:
                print("Error: Configuration file must contain [averaging] section")
                return 1
                
            averaging = config["averaging"]
            input_dir = averaging["input_directory"]
            frame_range = tuple(averaging["frame_range"])
            output_file = averaging["output_file"]
            
            # Check for paths in config
            paths = None
            multipath_config = None
            
            if "paths" in averaging:
                paths_config = averaging["paths"]
                if isinstance(paths_config, list):
                    paths = paths_config
                elif isinstance(paths_config, str):
                    paths = _parse_paths_arg(paths_config)
                else:
                    print(f"Warning: Invalid paths format in config: {paths_config}")
            
            # Check for multipath configuration
            if "multipath" in averaging:
                multipath_config = averaging["multipath"]
            
        except ConfigError as e:
            print(f"Configuration error: {e}")
            return 1
        except Exception as e:
            print(f"Failed to load configuration: {e}")
            return 1
    
    elif all([args.input_dir, args.start is not None, 
              args.end is not None, args.output]):
        # Use CLI arguments
        input_dir = args.input_dir
        frame_range = (args.start, args.end)
        output_file = args.output
        paths = _parse_paths_arg(args.paths) if args.paths else None
        multipath_config = None  # Multipath only available via TOML
        
        # Validate frame range
        if args.start >= args.end:
            print("Error: Start frame must be less than end frame")
            return 1
    
    else:
        # Missing arguments
        print("Error: Must provide either --config or all of: "
              "--input-dir, --start, --end, --output")
        parser.print_help()
        return 1
    
    # Run averaging
    try:
        average_chi_files(input_dir, frame_range, output_file, paths, multipath_config)
    except Exception as e:
        print(f"Error during averaging: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
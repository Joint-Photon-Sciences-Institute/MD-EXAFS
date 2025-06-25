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
from .build_chi_database import build_database
from .chi_database_query import ChiDatabaseQuery, PathQuery


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


def average_chi_from_database(
    db_path: Path,
    output_file: str,
    multipath_config: Dict[str, Any],
    frame_range: Optional[Tuple[int, int]] = None
) -> None:
    """
    Average chi(k) data using the pre-built database.
    
    Args:
        db_path: Path to the chi database
        output_file: Path to save averaged chi data
        multipath_config: Multipath configuration dictionary with keys:
                         - paths: List of path types (e.g., ["U-O", "U-O-O"])
                         - max_distance: List of maximum distances for each path type
        frame_range: Optional tuple of (start_frame, end_frame) to include
    """
    print(f"Using database: {db_path}")
    
    # Parse multipath configuration
    if "paths" not in multipath_config:
        raise ValueError("multipath_config must contain 'paths' key")
    
    path_types = multipath_config["paths"]
    
    # Handle max_distance configuration
    if "max_distance" in multipath_config:
        max_distances = multipath_config["max_distance"]
        if isinstance(max_distances, (int, float)):
            # Single value applies to all path types
            max_distances = [max_distances] * len(path_types)
        elif not isinstance(max_distances, list):
            max_distances = [max_distances]
    else:
        # No distance limit
        max_distances = [None] * len(path_types)
    
    # Ensure we have a max_distance for each path type
    if len(max_distances) < len(path_types):
        # Pad with the last value or None
        last_val = max_distances[-1] if max_distances else None
        max_distances.extend([last_val] * (len(path_types) - len(max_distances)))
    
    # Prepare frames list if frame_range is specified
    frames = None
    if frame_range:
        start_frame, end_frame = frame_range
        # Query database to find which frames actually exist
        with ChiDatabaseQuery(db_path) as db:
            frame_counts = db.get_frame_atom_counts()
            available_frames = sorted(frame_counts.keys())
            frames = [f for f in available_frames if start_frame <= f <= end_frame]
            print(f"Found {len(frames)} frames in range {start_frame}-{end_frame}")
    
    # Query and average data
    with ChiDatabaseQuery(db_path) as db:
        all_path_ids = []
        
        # Query each path type
        for path_type, max_dist in zip(path_types, max_distances):
            query = PathQuery(
                path_types=[path_type],
                max_reff=max_dist,
                frames=frames
            )
            
            paths = db.query_paths(query)
            path_ids = [p['id'] for p in paths]
            all_path_ids.extend(path_ids)
            
            print(f"Found {len(path_ids)} paths of type '{path_type}'" + 
                  (f" with reff <= {max_dist}" if max_dist else ""))
        
        if not all_path_ids:
            print("No paths found matching criteria")
            return
        
        # Average chi data using correct method (sum within atoms, then average)
        print(f"Processing {len(all_path_ids)} total paths...")
        averaged_data = db.sum_chi_within_atoms_then_average(all_path_ids)
        
        if averaged_data is None:
            print("No chi data found for selected paths")
            return
        
        # Save the averaged data
        np.savetxt(output_file, averaged_data, fmt='%.6f',
                  header="k(A^-1)  chi(k) - Sum within atoms then average across atoms (k=0:20:0.05)")
        print(f"Averaged data saved to {output_file}")
        
        # Print some statistics
        unique_path_types = db.get_unique_path_types()
        print(f"Database contains {len(unique_path_types)} unique path types")


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
    
    # Database-related arguments
    parser.add_argument(
        '--build-database',
        action='store_true',
        help='Build chi(k) database from FEFF calculations'
    )
    parser.add_argument(
        '--database',
        type=str,
        help='Path to chi(k) database file'
    )
    parser.add_argument(
        '--use-database',
        action='store_true',
        help='Use database for averaging instead of file-based processing'
    )
    parser.add_argument(
        '--rebuild',
        action='store_true',
        help='Rebuild database even if it exists'
    )
    parser.add_argument(
        '--num-workers',
        type=int,
        default=4,
        help='Number of parallel workers for database building (default: 4)'
    )
    
    args = parser.parse_args()
    
    # Handle database building first (doesn't need averaging parameters)
    if args.build_database:
        if not args.database:
            print("Error: --database required when using --build-database")
            return 1
        
        # For database building, we only need input directory
        if args.input_dir:
            base_dir = Path(args.input_dir)
        elif args.config:
            # Try to get input_dir from config
            try:
                config = load_config(args.config)
                if "averaging" in config and "input_directory" in config["averaging"]:
                    base_dir = Path(config["averaging"]["input_directory"])
                else:
                    print("Error: Input directory not found in config file")
                    return 1
            except Exception as e:
                print(f"Error loading config: {e}")
                return 1
        else:
            print("Error: --input-dir required for database building")
            return 1
        
        try:
            db_path = Path(args.database)
            build_database(base_dir, db_path, 
                         num_workers=args.num_workers,
                         rebuild=args.rebuild)
            print(f"Database successfully built at {db_path}")
            return 0
        except Exception as e:
            print(f"Error building database: {e}")
            return 1
    
    # For non-database operations, determine configuration source
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
    
    # Handle database-based averaging
    use_database = args.use_database
    if args.config and 'config' in locals():
        use_database = use_database or ("database" in config.get("averaging", {}))
    
    if use_database:
        # Determine database path
        if args.database:
            db_path = Path(args.database)
        elif args.config and 'config' in locals() and "database" in config.get("averaging", {}):
            db_config = config["averaging"]["database"]
            if isinstance(db_config, dict) and "path" in db_config:
                db_path = Path(db_config["path"])
            elif isinstance(db_config, str):
                db_path = Path(db_config)
            else:
                print("Error: Invalid database configuration")
                return 1
        else:
            print("Error: --database required when using --use-database")
            return 1
        
        if not db_path.exists():
            print(f"Error: Database not found: {db_path}")
            print("Run with --build-database first to create the database")
            return 1
        
        if not multipath_config:
            print("Error: Multipath configuration required for database-based averaging")
            print("Use --config with a TOML file containing [averaging.multipath] section")
            return 1
        
        try:
            average_chi_from_database(db_path, output_file, multipath_config, frame_range)
            return 0
        except Exception as e:
            print(f"Error during database averaging: {e}")
            return 1
    
    # Run file-based averaging
    try:
        average_chi_files(input_dir, frame_range, output_file, paths, multipath_config)
    except Exception as e:
        print(f"Error during averaging: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
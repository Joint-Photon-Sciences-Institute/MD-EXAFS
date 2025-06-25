"""
Build chi(k) database from FEFF calculation results.

This module scans the feff_calculations directory, processes all FEFF path files,
converts them to chi(k) using xraylarch, and stores the results in a SQLite database
for fast subsequent access.
"""

import argparse
import sqlite3
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from tqdm import tqdm
import json

from .db_schema import (
    create_database, validate_database, insert_path_metadata,
    insert_chi_data, insert_atom_info, get_database_stats,
    clear_database, optimize_database
)
from .multipath import parse_feff_header

try:
    from larch import Interpreter
    from larch.xafs import feffpath
    LARCH_AVAILABLE = True
except ImportError:
    LARCH_AVAILABLE = False
    print("Warning: xraylarch not available. Database building requires xraylarch.")

logger = logging.getLogger(__name__)


def _convert_feffdat_to_chi(feff_file: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Convert a single feff####.dat file to chi(k) using xraylarch.
    
    Args:
        feff_file: Path to the feff####.dat file
        
    Returns:
        Tuple of (k_grid, chi_values) as numpy arrays, or (None, None) on error
    """
    if not LARCH_AVAILABLE:
        return None, None
        
    try:
        # Create a new larch interpreter for this process
        larch_interp = Interpreter()
        
        # Create a FeffPath object
        path = feffpath(str(feff_file), _larch=larch_interp)
        
        # Initialize path parameters
        path.create_path_params()
        
        # Get k grid from the path data
        if hasattr(path, 'k') and path.k is not None:
            k = np.array(path.k)
        elif hasattr(path, '_feffdat') and hasattr(path._feffdat, 'k'):
            k = np.array(path._feffdat.k)
        else:
            logger.warning(f"No k grid found in path data for {feff_file}")
            return None, None
        
        # Calculate chi for this k grid
        path._calc_chi(k)
        
        # Get the chi array
        if hasattr(path, 'chi') and path.chi is not None:
            chi = np.array(path.chi)
        else:
            logger.warning(f"Could not calculate chi for {feff_file}")
            return None, None
        
        # Ensure k and chi have the same length
        if len(k) != len(chi):
            min_len = min(len(k), len(chi))
            k = k[:min_len]
            chi = chi[:min_len]
        
        return k, chi
        
    except Exception as e:
        logger.error(f"Error converting {feff_file}: {e}")
        return None, None


def process_atom_folder(atom_dir: Path, frame: int, atom_id: int) -> List[Dict[str, Any]]:
    """
    Process all FEFF paths in a single atom folder.
    
    Args:
        atom_dir: Path to the atom directory
        frame: Frame number
        atom_id: Atom ID
        
    Returns:
        List of dictionaries containing path data and chi(k) results
    """
    results = []
    
    # Find all feff####.dat files
    feff_files = sorted(atom_dir.glob("feff*.dat"))
    
    for feff_file in feff_files:
        # Extract path number from filename
        match = feff_file.name.replace('feff', '').replace('.dat', '')
        if match.isdigit():
            path_num = int(match)
        else:
            continue
            
        # Parse header to get path characteristics
        try:
            path_info = parse_feff_header(feff_file)
        except Exception as e:
            logger.warning(f"Failed to parse {feff_file}: {e}")
            continue
            
        # Convert to chi(k)
        k_grid, chi_values = _convert_feffdat_to_chi(feff_file)
        
        if k_grid is not None and chi_values is not None:
            results.append({
                'frame': frame,
                'atom_id': atom_id,
                'path_num': path_num,
                'path_info': path_info,
                'k_grid': k_grid,
                'chi_values': chi_values
            })
    
    return results


def scan_feff_calculations(base_dir: Path) -> List[Tuple[Path, int, int]]:
    """
    Scan the feff_calculations directory structure to find all atom folders.
    
    Args:
        base_dir: Base directory containing working_X folders
        
    Returns:
        List of (atom_dir, frame, atom_id) tuples
    """
    atom_folders = []
    
    # Find all working_X directories
    working_dirs = sorted([d for d in base_dir.glob("working_*") if d.is_dir()])
    
    for working_dir in working_dirs:
        # Find all frame_X directories
        frame_dirs = sorted([d for d in working_dir.glob("frame_*") if d.is_dir()])
        
        for frame_dir in frame_dirs:
            # Extract frame number
            try:
                frame = int(frame_dir.name.split('_')[1])
            except (IndexError, ValueError):
                continue
                
            # Find all atom_X directories
            atom_dirs = sorted([d for d in frame_dir.glob("atom_*") if d.is_dir()])
            
            for atom_dir in atom_dirs:
                # Extract atom ID
                try:
                    atom_id = int(atom_dir.name.split('_')[1])
                except (IndexError, ValueError):
                    continue
                    
                atom_folders.append((atom_dir, frame, atom_id))
    
    return atom_folders


def build_database(base_dir: Path, db_path: Path, num_workers: int = 4, 
                  rebuild: bool = False) -> None:
    """
    Build the chi(k) database from FEFF calculations.
    
    Args:
        base_dir: Base directory containing feff_calculations
        db_path: Path to the database file to create/update
        num_workers: Number of parallel workers
        rebuild: If True, clear existing database before building
    """
    if not LARCH_AVAILABLE:
        raise RuntimeError("xraylarch is required for database building")
    
    # Create or validate database
    if db_path.exists():
        if rebuild:
            logger.info(f"Clearing existing database at {db_path}")
            clear_database(db_path)
        elif not validate_database(db_path):
            raise ValueError(f"Invalid database at {db_path}")
    else:
        create_database(db_path)
    
    # Scan for all atom folders
    logger.info("Scanning for atom folders...")
    atom_folders = scan_feff_calculations(base_dir)
    logger.info(f"Found {len(atom_folders)} atom folders to process")
    
    if not atom_folders:
        logger.warning("No atom folders found to process")
        return
    
    # Process atom folders in parallel
    logger.info(f"Processing with {num_workers} workers...")
    
    # Open database connection
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode = WAL")  # Enable write-ahead logging
    conn.execute("PRAGMA synchronous = NORMAL")  # Faster writes
    
    try:
        # Process in batches to avoid memory issues
        batch_size = 100
        total_paths_processed = 0
        
        with tqdm(total=len(atom_folders), desc="Processing atoms") as pbar:
            for i in range(0, len(atom_folders), batch_size):
                batch = atom_folders[i:i + batch_size]
                
                # Process batch in parallel
                with ProcessPoolExecutor(max_workers=num_workers) as executor:
                    futures = {
                        executor.submit(process_atom_folder, atom_dir, frame, atom_id): (atom_dir, frame, atom_id)
                        for atom_dir, frame, atom_id in batch
                    }
                    
                    for future in as_completed(futures):
                        atom_dir, frame, atom_id = futures[future]
                        
                        try:
                            results = future.result()
                            
                            # Insert results into database
                            for result in results:
                                # Insert atom info (if not already present)
                                insert_atom_info(conn, result['frame'], result['atom_id'], 
                                               'U')  # Assuming uranium for now
                                
                                # Insert path metadata
                                path_id = insert_path_metadata(
                                    conn, result['frame'], result['atom_id'],
                                    result['path_num'], result['path_info']
                                )
                                
                                # Insert chi data
                                insert_chi_data(conn, path_id, result['k_grid'], 
                                              result['chi_values'])
                                
                                total_paths_processed += 1
                            
                        except Exception as e:
                            logger.error(f"Failed to process {atom_dir}: {e}")
                        
                        pbar.update(1)
                
                # Commit after each batch
                conn.commit()
        
        logger.info(f"Processed {total_paths_processed} paths total")
        
        # Optimize database
        logger.info("Optimizing database...")
        optimize_database(db_path)
        
        # Print statistics
        stats = get_database_stats(db_path)
        logger.info(f"Database statistics:")
        logger.info(f"  Total paths: {stats['total_paths']}")
        logger.info(f"  Total frames: {stats['total_frames']}")
        logger.info(f"  Total atoms: {stats['total_atoms']}")
        logger.info(f"  Database size: {stats['size_bytes'] / 1024 / 1024:.1f} MB")
        
    finally:
        conn.close()


def main():
    """Command-line interface for database building."""
    parser = argparse.ArgumentParser(
        description="Build chi(k) database from FEFF calculations"
    )
    parser.add_argument(
        "base_dir",
        type=Path,
        help="Base directory containing feff_calculations"
    )
    parser.add_argument(
        "database",
        type=Path,
        help="Path to output database file"
    )
    parser.add_argument(
        "-n", "--num-workers",
        type=int,
        default=4,
        help="Number of parallel workers (default: 4)"
    )
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Clear existing database before building"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Build database
    try:
        build_database(
            args.base_dir,
            args.database,
            num_workers=args.num_workers,
            rebuild=args.rebuild
        )
    except Exception as e:
        logger.error(f"Database building failed: {e}")
        raise


if __name__ == "__main__":
    main()
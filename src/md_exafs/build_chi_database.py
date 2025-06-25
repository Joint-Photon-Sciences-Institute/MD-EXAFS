"""
Build chi(k) database from FEFF calculation results.

This module scans the feff_calculations directory, processes all FEFF path files,
converts them to chi(k) using xraylarch, and stores the results in a SQLite database
for fast subsequent access.
"""

import argparse
import sqlite3
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from concurrent.futures import ProcessPoolExecutor, as_completed, TimeoutError
import numpy as np
from tqdm import tqdm
import json
import time

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

# Global variable to store Larch interpreter per process
_LARCH_INTERPRETER = None


def _get_larch_interpreter():
    """Get or create a Larch interpreter for this process."""
    global _LARCH_INTERPRETER
    if _LARCH_INTERPRETER is None:
        logger.debug(f"Creating Larch interpreter for process {os.getpid()}")
        _LARCH_INTERPRETER = Interpreter()
    return _LARCH_INTERPRETER


def _init_worker():
    """Initialize worker process with Larch interpreter."""
    # Pre-create the interpreter when worker starts
    _get_larch_interpreter()
    logger.debug(f"Worker {os.getpid()} initialized with Larch")


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
        # Get cached interpreter for this process
        larch_interp = _get_larch_interpreter()
        
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
    Scan directory to find all folders containing FEFF calculation files.
    
    Looks for directories containing feff*.dat files and attempts to extract
    frame and atom information from the directory structure.
    
    Args:
        base_dir: Base directory to search for FEFF calculations
        
    Returns:
        List of (atom_dir, frame, atom_id) tuples
    """
    atom_folders = []
    
    # Find all directories containing feff*.dat files
    logger.info(f"Scanning {base_dir} for FEFF calculations...")
    logger.info("This may take a while for large datasets...")
    
    # Use rglob to find all feff*.dat files recursively
    start_time = time.time()
    feff_files = []
    file_count = 0
    
    # Show progress during scanning
    for feff_file in base_dir.rglob("feff*.dat"):
        feff_files.append(feff_file)
        file_count += 1
        
        # Log progress every 10000 files
        if file_count % 10000 == 0:
            elapsed = time.time() - start_time
            rate = file_count / elapsed
            logger.info(f"  Scanned {file_count} files so far ({rate:.0f} files/sec)...")
    
    scan_time = time.time() - start_time
    logger.info(f"Scan complete: found {len(feff_files)} files in {scan_time:.1f} seconds")
    
    if not feff_files:
        logger.warning(f"No feff*.dat files found in {base_dir}")
        return atom_folders
    
    # Group files by their parent directory
    from collections import defaultdict
    dirs_with_feff = defaultdict(list)
    
    for feff_file in feff_files:
        dirs_with_feff[feff_file.parent].append(feff_file)
    
    logger.info(f"Found {len(dirs_with_feff)} directories containing FEFF files")
    
    # Process each directory containing FEFF files
    for atom_dir, files in sorted(dirs_with_feff.items()):
        # Try to extract frame and atom information from path
        frame = 0  # Default frame
        atom_id = 0  # Default atom ID
        
        # Try to extract atom_id from directory name (e.g., "atom_5")
        if atom_dir.name.startswith("atom_"):
            try:
                atom_id = int(atom_dir.name.split('_')[1])
            except (IndexError, ValueError):
                pass
        
        # Try to find frame number from parent directories
        for parent in atom_dir.parents:
            if parent == base_dir:
                break
            if parent.name.startswith("frame_"):
                try:
                    frame = int(parent.name.split('_')[1])
                    break
                except (IndexError, ValueError):
                    pass
        
        # If we still don't have IDs, generate them based on directory position
        if atom_id == 0:
            # Use directory index as atom_id
            atom_id = len(atom_folders)
        
        atom_folders.append((atom_dir, frame, atom_id))
        
        # Log first few directories found
        if len(atom_folders) <= 3:
            logger.info(f"  Found: {atom_dir.relative_to(base_dir)} (frame={frame}, atom={atom_id})")
    
    logger.info(f"Total atom folders found: {len(atom_folders)}")
    
    return atom_folders


def run_diagnostic(base_dir: Path, num_workers: int = 4) -> bool:
    """
    Run diagnostic tests before database building.
    
    Args:
        base_dir: Base directory containing FEFF calculations
        num_workers: Number of workers to test
        
    Returns:
        True if all tests pass
    """
    logger.info("Running diagnostic tests...")
    logger.info("=" * 60)
    
    # Test 1: Check xraylarch availability
    logger.info("Test 1: Checking xraylarch availability...")
    if not LARCH_AVAILABLE:
        logger.error("✗ xraylarch is not available. Install with: conda install -c conda-forge xraylarch")
        return False
    logger.info("✓ xraylarch is available")
    
    # Test 2: Check directory exists and is readable
    logger.info("\nTest 2: Checking directory access...")
    if not base_dir.exists():
        logger.error(f"✗ Directory not found: {base_dir}")
        return False
    if not os.access(base_dir, os.R_OK):
        logger.error(f"✗ Directory not readable: {base_dir}")
        return False
    logger.info(f"✓ Directory is accessible: {base_dir}")
    
    # Test 3: Quick file scan
    logger.info("\nTest 3: Quick file scan...")
    feff_count = 0
    for feff_file in base_dir.rglob("feff*.dat"):
        feff_count += 1
        if feff_count >= 10:
            break
    
    if feff_count == 0:
        logger.error("✗ No feff*.dat files found")
        return False
    logger.info(f"✓ Found FEFF files (sampled {feff_count})")
    
    # Test 4: Test multiprocessing with Larch
    logger.info(f"\nTest 4: Testing multiprocessing with {num_workers} workers...")
    try:
        def test_worker(dummy):
            _init_worker()
            return os.getpid()
        
        with ProcessPoolExecutor(max_workers=min(2, num_workers)) as executor:
            futures = [executor.submit(test_worker, i) for i in range(2)]
            pids = [f.result(timeout=30) for f in futures]
        
        logger.info(f"✓ Multiprocessing works with PIDs: {pids}")
    except Exception as e:
        logger.error(f"✗ Multiprocessing test failed: {e}")
        return False
    
    # Test 5: Test database creation
    logger.info("\nTest 5: Testing database creation...")
    test_db = base_dir / "test_diagnostic.db"
    try:
        create_database(test_db)
        if test_db.exists():
            logger.info("✓ Database creation works")
            test_db.unlink()
        else:
            logger.error("✗ Database file not created")
            return False
    except Exception as e:
        logger.error(f"✗ Database creation failed: {e}")
        return False
    
    logger.info("\n" + "=" * 60)
    logger.info("✓ All diagnostic tests passed!")
    return True


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
        
        # Create a persistent pool with worker initialization
        with ProcessPoolExecutor(
            max_workers=num_workers,
            initializer=_init_worker
        ) as executor:
            
            with tqdm(total=len(atom_folders), desc="Processing atoms") as pbar:
                for i in range(0, len(atom_folders), batch_size):
                    batch = atom_folders[i:i + batch_size]
                    
                    # Submit batch for processing
                    futures = {}
                    for atom_dir, frame, atom_id in batch:
                        # Submit with timeout
                        future = executor.submit(process_atom_folder, atom_dir, frame, atom_id)
                        futures[future] = (atom_dir, frame, atom_id)
                    
                    # Process results as they complete
                    for future in as_completed(futures, timeout=300):  # 5 minute timeout per batch
                        atom_dir, _, _ = futures[future]
                        
                        try:
                            results = future.result(timeout=60)  # 1 minute timeout per atom
                            
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
                            
                        except TimeoutError:
                            logger.error(f"Timeout processing {atom_dir}")
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
    parser.add_argument(
        "--diagnostic",
        action="store_true",
        help="Run diagnostic tests before building database"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Run diagnostic if requested
    if args.diagnostic:
        if not run_diagnostic(args.base_dir, args.num_workers):
            logger.error("Diagnostic tests failed. Fix issues before building database.")
            return 1
        
        # Ask if user wants to continue
        response = input("\nContinue with database building? (y/n): ")
        if response.lower() != 'y':
            logger.info("Database building cancelled.")
            return 0
    
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
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
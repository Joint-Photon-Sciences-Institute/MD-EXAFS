#!/usr/bin/env python
"""
Run FEFF calculations locally with parallel processing.

This script provides local parallel execution of FEFF calculations
as an alternative to HPC/SLURM submission.
"""

import os
import sys
import argparse
import subprocess
import multiprocessing
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm


def find_feff_executable():
    """Find the FEFF executable, preferring the bundled version."""
    # First check if md_exafs is installed
    try:
        import md_exafs
        bundled_feff = Path(md_exafs.__file__).parent / "feff8_lite" / "feff85L.exe"
        if bundled_feff.exists():
            return str(bundled_feff)
    except ImportError:
        pass
    
    # Check common locations
    possible_paths = [
        "feff85L.exe",
        "feff85L",
        "../src/md_exafs/feff8_lite/feff85L.exe",
        os.path.expanduser("~/bin/feff85L.exe"),
        "/usr/local/bin/feff85L.exe",
    ]
    
    for path in possible_paths:
        if Path(path).exists():
            return path
    
    raise FileNotFoundError(
        "Could not find FEFF executable. Please specify with --feff-path"
    )


def run_single_feff(args):
    """Run FEFF calculation in a single directory."""
    directory, feff_exe = args
    
    # Change to the directory
    original_dir = os.getcwd()
    os.chdir(directory)
    
    try:
        # Run FEFF
        with open("feff.out", "w") as outfile:
            result = subprocess.run(
                [feff_exe],
                stdout=outfile,
                stderr=subprocess.STDOUT,
                timeout=300  # 5 minute timeout
            )
        
        # Check if chi.dat was created
        if Path("chi.dat").exists():
            # Clean up intermediate files (keep only essential ones)
            for file in Path(".").glob("*"):
                if file.name not in ["feff.inp", "feff.out", "chi.dat"]:
                    try:
                        file.unlink()
                    except:
                        pass
            return True, directory
        else:
            return False, directory
            
    except subprocess.TimeoutExpired:
        return False, f"{directory} (timeout)"
    except Exception as e:
        return False, f"{directory} ({str(e)})"
    finally:
        os.chdir(original_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Run FEFF calculations locally with parallel processing"
    )
    parser.add_argument(
        "--base-dir",
        type=str,
        default="feff_calculations",
        help="Base directory containing FEFF input files (default: feff_calculations)"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: number of CPU cores)"
    )
    parser.add_argument(
        "--feff-path",
        type=str,
        default=None,
        help="Path to FEFF executable (default: auto-detect)"
    )
    
    args = parser.parse_args()
    
    # Find FEFF executable
    if args.feff_path:
        feff_exe = args.feff_path
        if not Path(feff_exe).exists():
            print(f"Error: FEFF executable not found at {feff_exe}")
            sys.exit(1)
    else:
        try:
            feff_exe = find_feff_executable()
        except FileNotFoundError as e:
            print(f"Error: {e}")
            sys.exit(1)
    
    print(f"Using FEFF executable: {feff_exe}")
    
    # Find all directories with feff.inp files
    base_path = Path(args.base_dir)
    if not base_path.exists():
        print(f"Error: Base directory {base_path} does not exist")
        sys.exit(1)
    
    feff_dirs = []
    for feff_inp in base_path.rglob("feff.inp"):
        feff_dirs.append(str(feff_inp.parent))
    
    if not feff_dirs:
        print(f"No feff.inp files found in {base_path}")
        sys.exit(1)
    
    print(f"Found {len(feff_dirs)} directories to process")
    
    # Count working directories
    working_dirs = set()
    for feff_dir in feff_dirs:
        # Extract working_X from path
        parts = Path(feff_dir).parts
        for part in parts:
            if part.startswith("working_"):
                working_dirs.add(part)
                break
    
    num_working_dirs = len(working_dirs)
    print(f"Found {num_working_dirs} working directories")
    
    # Determine number of workers
    if args.workers is None:
        cpu_count = multiprocessing.cpu_count()
        print(f"\nError: --workers must be specified")
        print(f"Your system has {cpu_count} CPU cores")
        print(f"The feff_calculations directory has {num_working_dirs} working directories")
        print(f"For optimal performance, set --workers equal to the number of working directories")
        print(f"Example: python run_feff_local.py --workers {num_working_dirs}")
        sys.exit(1)
    else:
        workers = args.workers
    
    print(f"Using {workers} parallel workers")
    
    # Process in parallel
    successful = 0
    failed = []
    
    # Prepare arguments for parallel execution
    work_items = [(dir_path, feff_exe) for dir_path in feff_dirs]
    
    with ProcessPoolExecutor(max_workers=workers) as executor:
        # Submit all tasks
        futures = {executor.submit(run_single_feff, item): item[0] 
                   for item in work_items}
        
        # Process results with progress bar
        with tqdm(total=len(feff_dirs), desc="Running FEFF calculations") as pbar:
            for future in as_completed(futures):
                success, info = future.result()
                if success:
                    successful += 1
                else:
                    failed.append(info)
                pbar.update(1)
    
    # Report results
    print(f"\nCompleted: {successful}/{len(feff_dirs)} successful")
    
    if failed:
        print(f"\nFailed calculations ({len(failed)}):")
        for fail_info in failed[:10]:  # Show first 10 failures
            print(f"  - {fail_info}")
        if len(failed) > 10:
            print(f"  ... and {len(failed) - 10} more")
    
    # Verify chi.dat files
    chi_files = list(base_path.rglob("chi.dat"))
    print(f"\nTotal chi.dat files created: {len(chi_files)}")
    
    return 0 if len(failed) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
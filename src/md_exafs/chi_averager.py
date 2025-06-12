"""Average chi(k) data from FEFF calculations."""

import os
import argparse
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any

import numpy as np
from tqdm import tqdm

from .config import load_config, ConfigError


def average_chi_files(
    input_dir: str,
    frame_range: Tuple[int, int],
    output_file: str
) -> None:
    """
    Average chi.dat files from FEFF calculations within specified frame range.
    
    Args:
        input_dir: Directory containing FEFF calculation results
        frame_range: Tuple of (start_frame, end_frame) to include
        output_file: Path to save averaged chi data
    """
    input_path = Path(input_dir)
    if not input_path.exists():
        raise ValueError(f"Input directory not found: {input_dir}")
    
    start_frame, end_frame = frame_range
    
    # Find all chi.dat files in the specified frame range
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
    for dirpath, dirnames, filenames in os.walk(root_dir):
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
        average_chi_files(input_dir, frame_range, output_file)
    except Exception as e:
        print(f"Error during averaging: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
#!/usr/bin/env python
"""
Example of using MD-EXAFS as a Python library.
"""

from md_exafs import load_config, process_trajectory, average_chi_files

def main():
    print("MD-EXAFS Python API Example")
    print("==========================")
    
    # Load configuration
    print("\n1. Loading configuration...")
    config = load_config("uo2_config.toml")
    print(f"   System: {config['system']['name']}")
    print(f"   Frames: {config['processing']['start_frame']}-{config['processing']['end_frame']}")
    print(f"   Step: {config['processing']['frame_step']}")
    
    # Process trajectory
    print("\n2. Processing trajectory...")
    try:
        process_trajectory(config)
        print("   ✓ Processing complete")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return
    
    # After FEFF calculations are done, average chi data
    print("\n3. Averaging chi(k) data...")
    print("   Note: This requires completed FEFF calculations")
    
    # Example of programmatic averaging
    try:
        average_chi_files(
            input_dir="feff_calculations",
            frame_range=(90000, 100000),
            output_file="chi_avg_api.dat"
        )
        print("   ✓ Averaging complete")
    except Exception as e:
        print(f"   ℹ Info: {e}")
        print("   (This is expected if FEFF calculations haven't been run)")

if __name__ == "__main__":
    main()
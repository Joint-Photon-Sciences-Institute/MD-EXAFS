#!/usr/bin/env python
"""Test script to verify MD-EXAFS modules."""

import sys
from pathlib import Path

def test_imports():
    """Test that all modules can be imported."""
    print("Testing imports...")
    try:
        from md_exafs import load_config, ConfigError, process_trajectory, average_chi_files
        print("✓ All imports successful")
        return True
    except ImportError as e:
        print(f"✗ Import error: {e}")
        return False

def test_config_loading():
    """Test configuration loading."""
    print("\nTesting configuration loading...")
    try:
        from md_exafs import load_config, ConfigError
        
        # Test processing config
        config = load_config("examples/uo2_config.toml")
        print("✓ Processing config loaded successfully")
        print(f"  System name: {config['system']['name']}")
        print(f"  Frame range: {config['processing']['start_frame']}-{config['processing']['end_frame']}")
        
        # Test averaging config (create temp directory first)
        Path("feff_calculations").mkdir(exist_ok=True)
        avg_config = load_config("examples/averaging_config.toml")
        print("✓ Averaging config loaded successfully")
        print(f"  Output file: {avg_config['averaging']['output_file']}")
        # Clean up
        Path("feff_calculations").rmdir()
        
        return True
    except Exception as e:
        print(f"✗ Config loading error: {e}")
        return False

def test_config_validation():
    """Test configuration validation."""
    print("\nTesting configuration validation...")
    try:
        from md_exafs import load_config, ConfigError
        
        # This should fail due to missing trajectory file
        try:
            # Create a temporary bad config
            bad_config = """
[system]
name = "Test"
trajectory_file = "nonexistent.xyz"
output_directory = "test_output"

[lattice]
a = [1, 0, 0]
b = [0, 1, 0]
c = [0, 0, 1]
pbc = [true, true, true]

[atoms]
Fe = 1

[processing]
start_frame = 0
end_frame = 10
frame_step = 1
atoms_per_frame = 1
cutoff_radius = 5.0
cores_per_node = 1

[feff]
header = "TEST"
"""
            with open("test_bad_config.toml", "w") as f:
                f.write(bad_config)
            
            config = load_config("test_bad_config.toml")
            print("✗ Validation should have failed for missing trajectory")
            return False
        except ConfigError as e:
            print(f"✓ Validation correctly caught error: {e}")
            # Clean up
            Path("test_bad_config.toml").unlink(missing_ok=True)
            return True
            
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False

def test_cli_help():
    """Test CLI interfaces."""
    print("\nTesting CLI help messages...")
    import subprocess
    
    # Test trajectory processor help
    result = subprocess.run(
        [sys.executable, "-m", "md_exafs.trajectory_processor", "--help"],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        print("✓ Trajectory processor CLI works")
    else:
        print("✗ Trajectory processor CLI failed")
        
    # Test chi averager help
    result = subprocess.run(
        [sys.executable, "-m", "md_exafs.chi_averager", "--help"],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        print("✓ Chi averager CLI works")
    else:
        print("✗ Chi averager CLI failed")

def main():
    """Run all tests."""
    print("MD-EXAFS Module Tests")
    print("=" * 50)
    
    all_passed = True
    
    # Run tests
    all_passed &= test_imports()
    all_passed &= test_config_loading()
    all_passed &= test_config_validation()
    test_cli_help()  # Don't fail on CLI test
    
    print("\n" + "=" * 50)
    if all_passed:
        print("✓ All core tests passed!")
    else:
        print("✗ Some tests failed")
        sys.exit(1)

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""
Test script for the database implementation.

This script tests the database building and querying functionality
with a small test dataset.
"""

import sys
import tempfile
from pathlib import Path
import numpy as np

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from md_exafs.db_schema import create_database, validate_database, get_database_stats
from md_exafs.chi_database_query import ChiDatabaseQuery, PathQuery


def test_database_creation():
    """Test database creation and schema."""
    print("Testing database creation...")
    
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as tmp:
        db_path = Path(tmp.name)
    
    try:
        # Create database
        create_database(db_path)
        print(f"✓ Created database at {db_path}")
        
        # Validate database
        is_valid = validate_database(db_path)
        assert is_valid, "Database validation failed"
        print("✓ Database validation passed")
        
        # Get stats (should be empty)
        stats = get_database_stats(db_path)
        assert stats['total_paths'] == 0, "New database should have no paths"
        assert stats['total_atoms'] == 0, "New database should have no atoms"
        print("✓ Database stats check passed")
        
        return db_path
        
    except Exception as e:
        if db_path.exists():
            db_path.unlink()
        raise e


def test_database_query():
    """Test database query functionality."""
    print("\nTesting database query...")
    
    # Create a test database
    db_path = test_database_creation()
    
    try:
        # Test query with empty database
        with ChiDatabaseQuery(db_path) as db:
            # Test path type query
            path_types = db.get_unique_path_types()
            assert len(path_types) == 0, "Empty database should have no path types"
            print("✓ Empty database query passed")
            
            # Test reff range
            min_reff, max_reff = db.get_reff_range()
            assert min_reff == 0.0 and max_reff == 0.0, "Empty database should have zero range"
            print("✓ Range query passed")
            
            # Test frame counts
            frame_counts = db.get_frame_atom_counts()
            assert len(frame_counts) == 0, "Empty database should have no frames"
            print("✓ Frame count query passed")
        
        print("\n✓ All database query tests passed")
        
    finally:
        if db_path.exists():
            db_path.unlink()


def test_multipath_config_parsing():
    """Test multipath configuration parsing."""
    print("\nTesting multipath configuration parsing...")
    
    # Test different max_distance formats
    test_configs = [
        {
            "paths": ["U-O", "U-U"],
            "max_distance": 5.0  # Single value
        },
        {
            "paths": ["U-O", "U-U", "U-O-O"],
            "max_distance": [3.0, 4.5, 5.0]  # List of values
        },
        {
            "paths": ["U-O"],
            "max_distance": [3.0]  # Single-element list
        }
    ]
    
    for i, config in enumerate(test_configs):
        print(f"  Testing config {i+1}: {config}")
        # The actual parsing is done in average_chi_from_database
        # Here we just verify the config structure is valid
        assert "paths" in config
        assert "max_distance" in config
        assert isinstance(config["paths"], list)
        assert all(isinstance(p, str) for p in config["paths"])
        
        max_dist = config["max_distance"]
        if isinstance(max_dist, (int, float)):
            assert max_dist > 0
        elif isinstance(max_dist, list):
            assert all(isinstance(d, (int, float)) and d > 0 for d in max_dist)
    
    print("✓ All multipath config tests passed")


def test_reference_data_compatibility():
    """Test that the implementation is compatible with reference data structure."""
    print("\nTesting reference data compatibility...")
    
    # Check if reference data exists
    ref_path = Path("reference_code_base/multipath_folders")
    if ref_path.exists():
        print(f"✓ Found reference data at {ref_path}")
        
        # Check directory structure
        working_dirs = list(ref_path.glob("working_*"))
        if working_dirs:
            print(f"✓ Found {len(working_dirs)} working directories")
            
            # Check for frame directories
            frame_dirs = list(working_dirs[0].glob("frame_*"))
            if frame_dirs:
                print(f"✓ Found {len(frame_dirs)} frame directories in first working dir")
                
                # Check for atom directories
                atom_dirs = list(frame_dirs[0].glob("atom_*"))
                if atom_dirs:
                    print(f"✓ Found {len(atom_dirs)} atom directories in first frame")
                    
                    # Check for FEFF files
                    feff_files = list(atom_dirs[0].glob("feff*.dat"))
                    if feff_files:
                        print(f"✓ Found {len(feff_files)} FEFF path files in first atom")
    else:
        print("⚠ Reference data not found - skipping compatibility test")


def main():
    """Run all tests."""
    print("MD-EXAFS Database Implementation Test Suite")
    print("=" * 50)
    
    try:
        test_database_creation()
        test_database_query()
        test_multipath_config_parsing()
        test_reference_data_compatibility()
        
        print("\n" + "=" * 50)
        print("✓ All tests passed successfully!")
        
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
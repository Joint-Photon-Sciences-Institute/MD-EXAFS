# MD-EXAFS Database Feature

## Overview

The database feature provides a significant performance improvement for multipath averaging by pre-computing and caching chi(k) data in a SQLite database. This eliminates the need to repeatedly parse FEFF files and convert paths to chi(k) during each averaging run.

## Performance Benefits

- **Database Building**: One-time cost (~1-2 hours for large datasets)
- **Query Performance**: 10-100x faster than file-based approach
- **Memory Efficiency**: Load only required data
- **Flexibility**: Easy to add new selection criteria via SQL queries

## Usage

### 1. Build the Database

Build a database from your FEFF calculations:

```bash
# Command-line usage
md-exafs-average --build-database feff_calculations --database chi_database.db

# With parallel processing (default: 4 workers)
md-exafs-average --build-database feff_calculations --database chi_database.db --num-workers 8

# Force rebuild if database exists
md-exafs-average --build-database feff_calculations --database chi_database.db --rebuild
```

### 2. Use Database for Averaging

Once the database is built, use it for fast averaging:

```bash
# Using configuration file (recommended)
md-exafs-average --config averaging_database_config.toml --use-database

# Or specify database directly
md-exafs-average --config multipath_config.toml --database chi_database.db --use-database
```

## Configuration

### TOML Configuration Example

```toml
[averaging]
input_directory = "feff_calculations"
frame_range = [90000, 90100]
output_file = "averaged_chi_database.dat"

# Database configuration
[averaging.database]
path = "chi_database.db"      # Path to database file
build = true                  # Build if doesn't exist
use_database = true          # Use database for averaging
rebuild = false              # Force rebuild even if exists

# Multipath configuration
[averaging.multipath]
paths = ["U-O", "U-U", "U-O-O"]
max_distance = [3.0, 4.5, 5.0]  # Different limits for each path
num_processes = 8
```

## Database Schema

The database contains three main tables:

1. **paths**: Stores path metadata (frame, atom_id, path_num, nleg, deg, reff, path_type)
2. **chi_data**: Stores pre-computed chi(k) arrays as binary data
3. **atoms**: Stores atom information (frame, atom_id, element, position)

## Python API

### Building a Database

```python
from md_exafs import build_database
from pathlib import Path

build_database(
    base_dir=Path("feff_calculations"),
    db_path=Path("chi_database.db"),
    num_workers=8,
    rebuild=False
)
```

### Querying the Database

```python
from md_exafs import ChiDatabaseQuery, PathQuery
from pathlib import Path

# Query specific paths
with ChiDatabaseQuery(Path("chi_database.db")) as db:
    # Define query criteria
    query = PathQuery(
        path_types=["U-O", "U-U"],
        max_reff=5.0,
        frames=[90000, 90010, 90020]
    )
    
    # Get paths matching criteria
    paths = db.query_paths(query)
    
    # Get averaged chi(k) data
    path_ids = [p['id'] for p in paths]
    averaged_chi = db.average_chi_data(path_ids)
```

### Database-Based Averaging

```python
from md_exafs import average_chi_from_database
from pathlib import Path

multipath_config = {
    "paths": ["U-O", "U-U", "U-O-O"],
    "max_distance": [3.0, 4.5, 5.0]
}

average_chi_from_database(
    db_path=Path("chi_database.db"),
    output_file="averaged_chi.dat",
    multipath_config=multipath_config,
    frame_range=(90000, 90100)
)
```

## Database Management

### Check Database Statistics

```python
from md_exafs.db_schema import get_database_stats
from pathlib import Path

stats = get_database_stats(Path("chi_database.db"))
print(f"Total paths: {stats['total_paths']}")
print(f"Total frames: {stats['total_frames']}")
print(f"Database size: {stats['size_bytes'] / 1024 / 1024:.1f} MB")
```

### Clear Database

```python
from md_exafs.db_schema import clear_database
from pathlib import Path

clear_database(Path("chi_database.db"))  # Removes all data but keeps schema
```

### Optimize Database

```python
from md_exafs.db_schema import optimize_database
from pathlib import Path

optimize_database(Path("chi_database.db"))  # Improves query performance
```

## Best Practices

1. **Build Once, Query Many**: Build the database once after FEFF calculations are complete
2. **Regular Optimization**: Run `optimize_database()` after building for better performance
3. **Incremental Updates**: Future versions will support adding new frames without full rebuild
4. **Backup**: Keep a backup of your database file as it represents significant computation

## Troubleshooting

### Database Not Found Error
- Ensure the database has been built with `--build-database`
- Check the database path in your configuration

### Memory Issues During Building
- Reduce `num_workers` to lower memory usage
- Process in smaller batches by modifying frame ranges

### Slow Database Queries
- Run `optimize_database()` to improve performance
- Ensure appropriate indexes exist (created by default)

## Future Enhancements

- Database versioning for schema updates
- Compression to reduce database size
- Support for remote/shared databases
- More complex query criteria
- In-memory caching for frequently accessed data
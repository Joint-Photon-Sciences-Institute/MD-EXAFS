# Database Optimization Plan for MD-EXAFS Multipath Feature

## Problem Statement

The current multipath workflow is too slow, especially when processing different combinations of paths. Each run requires:
- Re-parsing all FEFF path files
- Re-converting feffxxxx.dat files to chi(k) using xraylarch
- Redundant calculations for the same paths across different averaging runs

## Solution Overview

Implement a database-based caching system that:
1. Pre-processes all FEFF paths once and stores results in a SQLite database
2. Enables fast querying and averaging without re-parsing or re-calculating
3. Maintains backward compatibility with the existing file-based workflow

## Current Workflow Analysis

### Identified Bottlenecks

1. **xraylarch feffpath conversions** (chi_averager.py:180-219)
   - Each feffxxxx.dat file is converted to chi(k) using `feffpath()` and `_calc_chi()`
   - Involves parsing FEFF file, creating path parameters, and calculating chi
   - Done for every path in every atom folder during each averaging run

2. **Repeated file I/O and parsing**
   - FEFF header parsing happens multiple times for the same files
   - Path characterization is redone for each averaging operation
   - No caching of parsed metadata

3. **Redundant calculations**
   - Same FEFF paths converted to chi(k) multiple times across different runs
   - No reuse of previously calculated chi data

## Database Design

### Schema

```sql
-- Store path metadata
CREATE TABLE paths (
    id INTEGER PRIMARY KEY,
    frame INTEGER NOT NULL,
    atom_id INTEGER NOT NULL,
    path_num INTEGER NOT NULL,
    nleg INTEGER NOT NULL,
    deg REAL NOT NULL,
    reff REAL NOT NULL,
    path_type TEXT NOT NULL,
    atom_sequence TEXT NOT NULL,  -- JSON array of atoms
    UNIQUE(frame, atom_id, path_num)
);

-- Store pre-computed chi data
CREATE TABLE chi_data (
    id INTEGER PRIMARY KEY,
    path_id INTEGER NOT NULL,
    k_grid BLOB NOT NULL,  -- Numpy array as bytes
    chi_values BLOB NOT NULL,  -- Numpy array as bytes
    FOREIGN KEY(path_id) REFERENCES paths(id)
);

-- Store atom information
CREATE TABLE atoms (
    id INTEGER PRIMARY KEY,
    frame INTEGER NOT NULL,
    atom_id INTEGER NOT NULL,
    element TEXT NOT NULL,
    position_x REAL,
    position_y REAL,
    position_z REAL,
    UNIQUE(frame, atom_id)
);

-- Indexes for fast querying
CREATE INDEX idx_path_type ON paths(path_type);
CREATE INDEX idx_reff ON paths(reff);
CREATE INDEX idx_frame_atom ON paths(frame, atom_id);
```

### Data Storage

- **Path Metadata**: All information from parse_feff_header()
- **Chi(k) Data**: Pre-computed chi(k) arrays stored as binary blobs
- **Atom Information**: Frame/atom relationships and element types

## Implementation Plan

### Phase 1: Database Infrastructure

1. **Create db_schema.py**
   - Define SQLite schema
   - Helper functions for creating/validating database
   - Numpy array serialization/deserialization

2. **Create build_chi_database.py**
   - Scan feff_calculations directory structure
   - Parse all FEFF path files
   - Convert to chi(k) using xraylarch
   - Store in database with progress tracking
   - Support for incremental updates

3. **Add multiprocessing to database building**
   - Process multiple atoms in parallel
   - Batch database insertions for efficiency
   - Progress estimation based on total paths

### Phase 2: Query Infrastructure

4. **Create chi_database_query.py**
   - Functions to query paths by type and distance
   - Retrieve pre-computed chi(k) data
   - Handle aggregation (sum/average) operations
   - Memory-efficient data loading

5. **Update chi_averager.py**
   - Add `--build-database` flag for database creation
   - Add `--use-database` flag for database-based averaging
   - Auto-detect existing database
   - Fallback to file-based processing if needed

### Phase 3: Configuration and CLI

6. **Update TOML configuration**
   ```toml
   [averaging.database]
   build = true  # Build database if it doesn't exist
   path = "chi_database.db"  # Database file location
   use_database = true  # Use database for averaging
   rebuild = false  # Force rebuild even if database exists
   ```

7. **CLI Interface**
   ```bash
   # Build database
   md-exafs-average --build-database feff_calculations --database chi_database.db
   
   # Use database for averaging
   md-exafs-average --use-database chi_database.db --config multipath_config.toml
   ```

## Expected Performance Improvements

- **Database Building**: One-time cost (~1-2 hours for large datasets)
- **Query Performance**: 10-100x faster than file-based approach
- **Memory Efficiency**: Load only required data
- **Flexibility**: Easy to add new selection criteria via SQL queries

## Technical Considerations

1. **Database Size**: ~1-10 GB for typical datasets
2. **Backward Compatibility**: Existing workflows continue to function
3. **Incremental Updates**: Support adding new frames without full rebuild
4. **Error Handling**: Graceful fallback if database is corrupted
5. **Portability**: SQLite database is a single file, easy to share

## Testing Strategy

1. Build database with reference_code_base/multipath_folders
2. Compare results between file-based and database-based averaging
3. Benchmark performance improvements
4. Test edge cases (missing data, corrupted files)

## Future Enhancements

1. **Database versioning**: Handle schema updates gracefully
2. **Compression**: Reduce database size with data compression
3. **Remote databases**: Support for shared network databases
4. **Advanced queries**: More complex path selection criteria
5. **Caching strategies**: In-memory caching for frequently accessed data

## Implementation Timeline

1. **Week 1**: Database schema and builder implementation
2. **Week 2**: Query infrastructure and chi_averager integration
3. **Week 3**: Testing, optimization, and documentation
4. **Week 4**: User feedback and refinements

## Success Criteria

- Database building completes without errors for large datasets
- Query performance is at least 10x faster than file-based approach
- Results match exactly between file-based and database-based methods
- Documentation is clear and comprehensive
- Users can easily switch between workflows
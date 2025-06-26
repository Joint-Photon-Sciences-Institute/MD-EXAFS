# MD-EXAFS Database Implementation Continuation

## Current Status

### Database Building
- **Issue**: Database building was hanging at 0% on HPC with 10,000 atoms (3.6M files)
- **Root Causes Identified**:
  1. Multiprocessing pickle errors with nested functions
  2. Larch interpreter initialization overhead (12s per file)
  3. No progress feedback during 167-second file scan
  4. ProcessPoolExecutor being recreated for each batch

- **Fixes Implemented** (commit 2a3d65b):
  1. Moved worker functions to module level to fix pickle errors
  2. Cache Larch interpreter per worker process (created once, reused for all files)
  3. Added progress logging during file scanning
  4. Use persistent ProcessPoolExecutor with initializer
  5. Added timeouts (5 min/batch, 1 min/atom) to prevent hanging
  6. Added `--diagnostic` mode to test setup

### Database Averaging
- **Critical Bug Found**: Database was averaging all paths directly instead of sum-within-atoms then average
- **Fix Implemented** (commit dbb8130):
  1. Added `sum_chi_within_atoms_then_average()` method to `chi_database_query.py`
  2. Updated `average_chi_from_database()` to use corrected method
  3. Now correctly: Groups by atom → Sums paths within atom → Averages across atoms

## Current Database Building Command
```bash
md-exafs-average --build-database --input-dir feff_calculations_keep --database 298_database.db --num-workers 112
```

## Key Implementation Details

### File Structure Found
- 10,000 directories with FEFF files
- ~360 feff*.dat files per directory
- Total: 3.6 million files
- Example: `frame_180000/atom_0/` contains 368 files

### Database Schema
- `paths` table: frame, atom_id, path_num, nleg, deg, reff, path_type, atom_sequence
- `chi_data` table: k_grid and chi_values as binary blobs
- `atoms` table: frame, atom_id, element, position

### Critical Functions
1. `_init_worker()` - Initializes Larch once per worker
2. `_get_larch_interpreter()` - Returns cached interpreter
3. `scan_feff_calculations()` - Now works with any directory structure
4. `sum_chi_within_atoms_then_average()` - Correct physical averaging

## Testing Status
- SQLite works on NFS filesystem
- Multiprocessing with Larch works when functions at module level
- File scanning takes ~167 seconds for 3.6M files
- Database creation successful with WAL mode

## Next Steps if Database Building Still Has Issues
1. Monitor progress during file scan (should show rate)
2. Check if workers are initializing (should see progress after scan)
3. If still hanging, try with fewer workers first (`--num-workers 4`)
4. Use `--diagnostic` flag to test setup before full run

## Configuration Files
- Using `average_multi.toml` for multipath averaging
- Database path: `298_database.db`
- Frame range: 180000-200000 (1000 frames)
- Path types being averaged: U-U with reff <= 4.0

## Important Notes
- The database builder now processes ALL files found, not filtered by frame range
- Frame filtering happens during averaging, not building
- Each atom must be treated as an independent measurement
- The averaging must be: sum within atoms, then average across atoms

## Git Status
- On branch: main
- Latest commits include database building and averaging fixes
- All changes committed and ready for testing

## Environment
- Running on HPC with NFS filesystem
- Python 3.10.18, SQLite 3.50.1
- Using conda environment: md-exafs
- 112 workers available for parallel processing
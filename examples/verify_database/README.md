# Database Verification Example

This example demonstrates how to verify that the chi(k) calculations from the database match the original FEFF calculations.

## Scripts

### 1. build_and_verify.py
Builds a database from the reference multipath folders and compares the results with FEFF chi.dat files.

```bash
python build_and_verify.py
```

This script will:
- Build a SQLite database from `/home/nickj/claude-code/MD-EXAFS/reference_code_base/multipath_folders`
- Compare chi(k) from the database with FEFF chi.dat files
- Create comparison plots showing k²-weighted EXAFS
- Generate a histogram of differences

### 2. compare_single_atom.py
Compares chi(k) calculations for a single atom folder using three methods:
1. Direct summation using xraylarch (current implementation)
2. From the database
3. FEFF chi.dat (ground truth)

```bash
# Compare a specific atom folder
python compare_single_atom.py /home/nickj/claude-code/MD-EXAFS/reference_code_base/multipath_folders/working_0/frame_0/atom_0

# Use a specific database
python compare_single_atom.py /path/to/atom_folder /path/to/database.db
```

## Expected Results

With the corrected implementation using `path2chi()`, you should see:
- Maximum differences < 0.001 between all methods
- Excellent visual agreement in k²-weighted plots
- The database results should match both direct summation and FEFF chi.dat

## Troubleshooting

If you see large differences (> 0.01):
1. Check that `path2chi()` is being used (not `_calc_chi()`)
2. Verify that the k-grid is consistent (0-20 Å⁻¹, step 0.05)
3. Ensure xraylarch is properly installed via conda

## Files Generated

- `test_verification.db` - SQLite database with chi(k) data
- `verification_plots/chi_comparison_k2.png` - Comparison plots
- `verification_plots/difference_histogram.png` - Distribution of differences
- Individual `chi_comparison.png` files in atom folders when using compare_single_atom.py
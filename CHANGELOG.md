# Changelog

All notable changes to MD-EXAFS will be documented in this file.

## [Unreleased]

### Added
- `plot_chi_from_database.py` - Flexible plotting script for visualizing chi(k) spectra from database
  - Multiple plotting modes: individual paths, atom sums, averages, and path comparisons
  - PNG export capability for all plot types

### Changed
- Database query methods now use cubic spline interpolation instead of linear interpolation
  - Fixes jagged peaks in chi(k) spectra when interpolating from 0.1 to 0.05 Å⁻¹ grid
  - Produces smooth curves consistent with file-based averaging
  - Affects `average_chi_data()`, `sum_chi_data()`, and `sum_chi_within_atoms_then_average()`

### Fixed
- Database averaging now displays the number of unique atoms being processed
- Output explicitly states the averaging method (sum within atoms, then average)
- File headers clarify that data is interpolated to uniform k-grid

## Previous Updates

### Database Feature
- Pre-compute chi(k) data for 10-100x faster multipath averaging
- SQLite database with optimized schema for path queries
- Support for building and querying large datasets (millions of paths)

### Multipath Analysis
- Selective path processing by scattering type and distance
- Parallel processing with configurable workers
- TOML-based configuration for complex path selections
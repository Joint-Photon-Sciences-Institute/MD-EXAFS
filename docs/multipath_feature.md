# Multipath Feature Documentation

## Overview

The multipath feature in MD-EXAFS allows selective processing and averaging of FEFF scattering paths based on:
- **Path type**: The atoms involved in the scattering path (e.g., U-O, U-U, U-O-O)
- **Distance**: The effective path length (reff) in Angstroms

This feature enables detailed analysis of specific scattering contributions to the total EXAFS signal.

## Configuration

The multipath feature is configured exclusively through TOML files. It cannot be accessed via command-line arguments.

### TOML Structure

```toml
[averaging]
input_directory = "feff_calculations"
frame_range = [90000, 100000]
output_file = "chi_multipath.dat"

[averaging.multipath]
paths = ["U-O", "U-O-O"]     # Path types to include
max_distance = 3.0            # Maximum path distance in Angstroms (optional)
num_processes = 4             # Number of parallel processes (optional)
```

### Configuration Parameters

- **paths**: List of path types using full atom notation
  - Single scattering: `"U-O"`, `"U-U"`
  - Double scattering: `"U-O-O"`, `"U-U-O"`, `"U-O-U"`
  - Triple scattering: `"U-O-U-O"`, etc.
- **max_distance**: Maximum effective path length in Angstroms (optional)
  - If not specified, no distance limit is applied
- **num_processes**: Number of parallel processes for faster computation (optional)
  - Default: 1 (sequential processing)

## Path Type Notation

Path types show the complete scattering sequence:
- The first atom is always the absorber (e.g., U)
- Subsequent atoms show the scattering path
- Examples:
  - `U-O`: U → O → U (single scattering through oxygen)
  - `U-U`: U → U → U (single scattering through uranium)
  - `U-O-O`: U → O → O → U (double scattering through two oxygens)

## Usage Example

1. Create a TOML configuration file (e.g., `multipath_config.toml`):

```toml
[averaging]
input_directory = "reference_code_base/multipath_folders"
frame_range = [180360, 180400]
output_file = "chi_uo_paths.dat"

[averaging.multipath]
paths = ["U-O", "U-O-O"]
max_distance = 4.0
num_processes = 4
```

2. Run the averaging:

```bash
md-exafs-average --config multipath_config.toml
```

## Output

The multipath averaging produces:
1. **Per-atom files**: `chi_multipath.dat` in each atom folder containing the sum of selected paths
2. **Final averaged file**: Specified by `output_file`, containing the average across all atoms
3. **Standard k-grid**: Output is interpolated to k = 0 to 20 Å⁻¹ with step 0.05

## Performance Considerations

- Multiprocessing scales linearly with available cores
- Processing time depends on:
  - Number of atoms
  - Number of paths selected
  - Frame range
- Typical performance: ~1 atom/second per core

## Requirements

- xraylarch must be installed: `conda install -c conda-forge xraylarch`
- FEFF calculations must include individual path files (feffxxxx.dat)
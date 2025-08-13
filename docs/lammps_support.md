# LAMMPS Trajectory Support

## Overview

The md-exafs package now includes native support for LAMMPS trajectory files (`.lammpstrj` and other LAMMPS formats). LAMMPS uses numeric atom type IDs instead of element symbols, and the package automatically detects and handles this format difference.

## Configuration Changes

### Atom Type Mapping

For LAMMPS trajectories, the `[atoms]` section in your TOML configuration file should map numeric type IDs (as strings) to element symbols:

```toml
[atoms]
# For LAMMPS: map type ID -> element symbol
# Type IDs must be quoted strings
"1" = "Zn"
"2" = "Fe"
"3" = "O"
```

This differs from Cp2k/XYZ format where element symbols are mapped to type IDs:

```toml
[atoms]
# For Cp2k/XYZ: map element symbol -> type ID
Zn = 1
Fe = 2
O = 3
```

## Example Configuration

Here's a complete example configuration for analyzing RDF from a LAMMPS trajectory:

```toml
[input]
trajectory_file = "zfo_unitcell.lammpstrj"

[atoms]
# LAMMPS type ID mapping
1 = "Zn"
2 = "Fe"
3 = "O"

[lattice]
a = [8.435505, 0.000000, 0.000000]
b = [0.000000, 8.435505, 0.000000]
c = [0.000000, 0.000000, 8.435505]
pbc = [true, true, true]

[output]
directory = "."
report_file = "rdf_analysis.txt"
plot_file = "rdf_plot.png"

[analysis]
frame_start = 9000
frame_end = 10000
frame_step = 20

[rdf]
cutoff_radius = 6.0
num_bins = 400

[[peaks]]
center = "Fe"
neighbor = "O"
r_min = 1
r_max = 3.2
fit_type = "both"
```

## How It Works

1. **Automatic Detection**: The package automatically detects LAMMPS format by checking if particle type names in the trajectory are empty or numeric.

2. **Type Mapping**: For LAMMPS trajectories, the package creates an internal mapping from type IDs to element names and uses this throughout the analysis pipeline.

3. **Compatibility**: Full backward compatibility with Cp2k/XYZ format is maintained - the package automatically detects the format and handles it appropriately.

## Supported Features

- RDF analysis with LAMMPS trajectories
- FEFF input generation from LAMMPS trajectories
- All existing features of md-exafs package

## Usage

### RDF Analysis

To analyze RDF from a LAMMPS trajectory:

```bash
python -m md_exafs.rdf_analysis --rdf config.toml
```

### FEFF Input Generation

To generate FEFF input files from a LAMMPS trajectory:

```bash
md-exafs-process config.toml
```

## Important Notes

- **Type IDs must be quoted**: In the TOML configuration, type IDs must be strings (e.g., `"1"` not `1`)
- **Type IDs must match**: Ensure type IDs in the configuration match those in your LAMMPS dump file
- **Use element symbols**: When specifying `absorbing_element` or peak centers/neighbors, use element symbols (e.g., "Fe"), not type IDs
- **Lattice parameters**: Should match your LAMMPS simulation box dimensions
- **File formats**: Supports standard LAMMPS dump formats readable by OVITO
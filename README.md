# MD-EXAFS

A Python package for calculating Extended X-ray Absorption Fine Structure (EXAFS) spectra from Molecular Dynamics (MD) trajectories.

## Overview

MD-EXAFS processes MD simulation trajectories to generate FEFF input files and compute averaged EXAFS χ(k) spectra. The package is designed for high-performance computing environments and supports parallel processing of large trajectory datasets.

## Key Features

- **Flexible Configuration**: TOML-based configuration for all parameters
- **Multi-Element Support**: Works with any atomic system
- **Parallel Processing**: Optimized for HPC environments
- **Database Caching**: Pre-compute chi(k) data for 10-100x faster multipath averaging
- **Multipath Analysis**: Advanced path selection by scattering type and distance
- **PBC Support**: Full periodic boundary condition handling
- **Multiple Formats**: Compatible with OVITO-supported trajectories (XYZ, LAMMPS, VASP, etc.)

## Installation

### Using Conda (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/md-exafs.git
cd md-exafs

# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate md-exafs

# Install package in development mode
pip install -e .
```

### Manual Installation

```bash
# IMPORTANT: First install OVITO via conda (required)
conda install --strict-channel-priority -c https://conda.ovito.org -c conda-forge ovito=3.12.4

# Install xraylarch (required for path processing)
conda install -c conda-forge xraylarch

# Clone and install
git clone https://github.com/yourusername/md-exafs.git
cd md-exafs
pip install -e .
```

**Requirements**: Python ≥ 3.8, OVITO (conda-only), xraylarch, numpy, tqdm

## Quick Start

### 1. Configure Your System

Create a TOML configuration file (see `examples/uo2_config.toml`):

```toml
[system]
name = "UO2"
trajectory_file = "path/to/trajectory.xyz"
output_directory = "feff_calculations"

[lattice]
a = [21.875, 0.328, 0.0]
b = [0.328, 21.875, 0.0]
c = [0.0, 0.0, 21.880]
pbc = [true, true, true]

[atoms]
U = 1  # Atom type mapping from trajectory
O = 2

[processing]
start_frame = 90000
end_frame = 100000
frame_step = 10
atoms_per_frame = 10
cutoff_radius = 10.0

[feff]
absorbing_element = "U"
```

### 2. Standard Workflow

```bash
# Step 1: Process trajectory to generate FEFF inputs
md-exafs-process config.toml

# Step 2: Run FEFF calculations (choose one)
# Option A: Local execution
md-exafs-feff-local --base-dir feff_calculations --workers 10

# Option B: HPC execution
sbatch examples/run_slurm.sh

# Step 3: Average chi(k) data
md-exafs-average --config averaging_config.toml
```

### 3. Database-Optimized Workflow (Recommended for Multipath)

For faster multipath averaging, use the database feature:

```bash
# Build database once (processes ALL frames/atoms found)
md-exafs-average --build-database --input-dir feff_calculations --database chi_database.db

# Then use database for fast averaging with different configurations
md-exafs-average --config multipath_config.toml --use-database --database chi_database.db
```

The database provides 10-100x speedup for multipath averaging by eliminating redundant calculations.

## Detailed Usage

### Processing MD Trajectories

The `md-exafs-process` command reads your trajectory and generates FEFF input files:

```bash
md-exafs-process config.toml
```

This creates a directory structure optimized for parallel FEFF execution:
```
feff_calculations/
├── working_0/
│   ├── frame_90000/
│   │   ├── atom_0/feff.inp
│   │   ├── atom_1/feff.inp
│   │   └── ...
│   └── frame_90010/
└── working_1/
```

### Running FEFF Calculations

#### Local Execution
For testing or small datasets:

```bash
md-exafs-feff-local --base-dir feff_calculations --workers 10
```

#### HPC Execution
For production runs, modify the SLURM script for your system:

```bash
sbatch examples/run_slurm.sh
```

### Averaging Chi(k) Data

#### Basic Averaging
Average all chi.dat files in a frame range:

```bash
# Using configuration file
md-exafs-average --config averaging_config.toml

# Using command-line arguments
md-exafs-average --input-dir feff_calculations --start 90000 --end 100000 --output chi_avg.dat
```

#### Path-Specific Averaging
Process specific FEFF paths:

```bash
# Paths 1-8
md-exafs-average --input-dir feff_calculations --start 90000 --end 100000 --paths "[1,8]" --output chi_partial.dat

# Specific paths only
md-exafs-average --input-dir feff_calculations --start 90000 --end 100000 --paths "[1,3,5,7]" --output chi_selected.dat
```

#### Multipath Analysis
For advanced path selection by type and distance:

```toml
[averaging]
input_directory = "feff_calculations"
frame_range = [90000, 100000]
output_file = "chi_multipath.dat"

[averaging.multipath]
paths = ["U-O", "U-U", "U-O-O"]    # Path types to include
max_distance = [3.0, 4.5, 5.0]     # Max distance for each type
num_processes = 8                   # Parallel workers
```

### Database Feature

For repeated averaging with different multipath configurations:

#### Build Database
```bash
# Build from all FEFF calculations (one-time operation)
md-exafs-average --build-database --input-dir feff_calculations --database chi_database.db --num-workers 8

# Force rebuild
md-exafs-average --build-database --input-dir feff_calculations --database chi_database.db --rebuild
```

#### Use Database
```toml
[averaging]
input_directory = "feff_calculations"
frame_range = [90000, 100000]  # Only these frames will be averaged
output_file = "chi_database_avg.dat"

[averaging.database]
path = "chi_database.db"
use_database = true

[averaging.multipath]
paths = ["U-O", "U-U"]
max_distance = [3.0, 5.0]
```

See `docs/database_feature_README.md` for detailed documentation.

## Python API

### Basic Usage

```python
from md_exafs import load_config, process_trajectory, average_chi_files

# Process trajectory
config = load_config("config.toml")
process_trajectory(config)

# Average chi data
average_chi_files(
    input_dir="feff_calculations",
    frame_range=(90000, 100000),
    output_file="chi_avg.dat"
)
```

### Database Usage

```python
from md_exafs import build_database, average_chi_from_database
from pathlib import Path

# Build database
build_database(
    base_dir=Path("feff_calculations"),
    db_path=Path("chi_database.db"),
    num_workers=8
)

# Use database for averaging
multipath_config = {
    "paths": ["U-O", "U-U"],
    "max_distance": [3.0, 5.0]
}

average_chi_from_database(
    db_path=Path("chi_database.db"),
    output_file="chi_avg.dat",
    multipath_config=multipath_config,
    frame_range=(90000, 100000)
)
```

## Utility Tools

### Structure Generator

Generate MD input files from CIF structures:

```bash
# Create 3×3×3 supercell in CP2K format
md-exafs-md-input-gen --cp2k --input structure.cif --size 3,3,3 --output supercell.xyz
```

### MD Convergence Check

Check MD simulation convergence by analyzing energy evolution from CP2K .ener files:

```bash
# Generate convergence plots
md-exafs-md --check_convergence --input simulation.ener --output convergence.png
```

This creates a comprehensive convergence analysis including:
- Temperature evolution with running average
- Total energy conservation
- Kinetic and potential energy trends
- Energy drift analysis (in kcal/mol)

The tool also prints convergence statistics including average temperature, energy drift, and relative drift percentage.

### Plotting Chi(k) Spectra

Plot chi(k) spectra directly from the database:

```bash
# Plot individual paths for a specific atom
./plot_chi_from_database.py database.db --mode paths --frame 180000 --atom 0 --paths 1,2,3,4

# Plot sum of paths for an atom
./plot_chi_from_database.py database.db --mode atom --frame 180000 --atom 0 --path-types U-O --max-reff 3.0

# Plot averaged chi(k) across multiple atoms
./plot_chi_from_database.py database.db --mode average --start-frame 180000 --end-frame 190000 --path-types U-U

# Save plot as PNG
./plot_chi_from_database.py database.db --mode atom --frame 180000 --atom 0 --output chi_plot.png
```

## Examples

The `examples/` directory contains:

- `uo2_single_frame.xyz` - Test trajectory for quick validation
- `uo2_config.toml` - Complete processing configuration
- `averaging_config.toml` - Basic averaging configuration
- `multipath_averaging.toml` - Multipath feature example
- `averaging_database_config.toml` - Database-based averaging
- `workflow_example.sh` - Complete workflow script
- `run_slurm.sh` - HPC batch script template

## Documentation

- `docs/multipath_feature.md` - Detailed multipath analysis guide
- `docs/database_feature_README.md` - Database optimization guide
- `examples/README.md` - Example walkthrough

## HPC Considerations

1. **Directory Distribution**: Files are distributed across `working_*` directories based on `cores_per_node`
2. **Independent Calculations**: Each FEFF run is independent, enabling perfect parallel scaling
3. **Memory Management**: Database building processes atoms in batches to manage memory
4. **Load Balancing**: Workers dynamically pull from calculation queue

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.
# MD-EXAFS

A Python package for calculating Extended X-ray Absorption Fine Structure (EXAFS) spectra from Molecular Dynamics (MD) trajectories.

## Overview

MD-EXAFS processes MD simulation trajectories to generate FEFF input files and compute averaged EXAFS χ(k) spectra. The package is designed for high-performance computing environments and supports parallel processing of large trajectory datasets.

The package has been tested with CP2K trajectory files but should work with any trajectory format compatible with OVITO (XYZ, LAMMPS, VASP, etc.).

## Features

- **Flexible Configuration**: TOML-based configuration for all parameters
- **Multi-Element Support**: Works with any atomic system (not limited to specific elements)
- **Parallel Processing**: Optimized for HPC environments with configurable core distribution
- **PBC Support**: Full periodic boundary condition handling with minimum image convention
- **Independent Averaging**: Chi averaging can be run separately with different frame ranges
- **Multipath Analysis**: Advanced path selection by scattering type and distance
- **CLI and API**: Both command-line tools and Python API access

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

### Using pip

```bash
# IMPORTANT: First install OVITO via conda (required)
conda install --strict-channel-priority -c https://conda.ovito.org -c conda-forge ovito=3.12.4

# Install xraylarch (required for path processing)
pip install xraylarch

# Clone and install
git clone https://github.com/yourusername/md-exafs.git
cd md-exafs
pip install -e .
```

**Note**: Requires Python ≥ 3.8. OVITO must be installed via conda (not pip) for trajectory processing. Xraylarch is required for FEFF path processing features.

## Quick Start

### Testing with Example Data

The package includes a single-frame UO₂ trajectory (`examples/uo2_single_frame.xyz`) for quick testing without needing large trajectory files. The example configuration is pre-configured for this test data.

### 1. Prepare Configuration

Create a TOML configuration file for your system (see `examples/uo2_config.toml` for the test example):

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
absorbing_element = "U"  # Which element is the X-ray absorber
```

### 2. Process Trajectory

```bash
md-exafs-process config.toml
```

This generates FEFF input files organized in directories optimized for parallel processing.

### 3. Run FEFF Calculations

#### Option A: Local Execution (for testing/small datasets)

Run FEFF calculations locally with parallel processing:

```bash
# Check how many FEFF calculations need to be run
find feff_calculations -name "feff.inp" | wc -l

# Run with appropriate number of workers
md-exafs-feff-local --base-dir feff_calculations --workers 10
```

**Note:** The `--workers` parameter is required. Each FEFF calculation (one per atom) runs independently, so you can use multiple workers up to the number of calculations.

#### Option B: HPC Execution (for production)

Submit the SLURM script to your HPC queue:

```bash
sbatch examples/run_slurm.sh
```

### 4. Average Chi(k) Data

After FEFF calculations complete:

```bash
# Using TOML config
md-exafs-average --config averaging_config.toml

# Or using CLI arguments
md-exafs-average --input-dir feff_calculations --start 90000 --end 100000 --output chi_avg.dat
```

## Command-Line Tools

MD-EXAFS provides several command-line tools for different stages of the workflow:

### md-exafs-process

Process MD trajectories to generate FEFF input files.

```bash
md-exafs-process config.toml
```

This tool reads the trajectory file and creates FEFF input files for each selected atom in each frame, organized into directories for parallel processing.

### md-exafs-feff-local

Run FEFF calculations locally with parallel processing.

```bash
md-exafs-feff-local --base-dir feff_calculations --workers 10
```

Options:
- `--base-dir`: Directory containing FEFF input files (required)
- `--workers`: Number of parallel workers (required)

### md-exafs-average

Average chi(k) data from completed FEFF calculations.

```bash
# Using TOML config
md-exafs-average --config averaging_config.toml

# Using CLI arguments
md-exafs-average --input-dir feff_calculations --start 90000 --end 100000 --output chi_avg.dat

# Process specific FEFF paths
md-exafs-average --input-dir feff_calculations --start 90000 --end 100000 --paths "[1,8]" --output chi_partial.dat
```

Options:
- `--config`: Path to TOML configuration file
- `--input-dir`: Directory containing FEFF calculations
- `--start`: Starting frame number
- `--end`: Ending frame number
- `--output`: Output file for averaged chi(k) data
- `--paths`: Process specific FEFF paths (e.g., "[1,8]" for paths 1-8, or "[1,3,5,7]" for specific paths)

When using `--paths`:
- Converts feffxxxx.dat files to chi(k) data
- Sums paths within each atom folder (saved as chi_partial_0.dat)
- Averages the sums across all atoms (saved as output file)
- Automatically interpolates output to standard k-grid (0 to 20 Å⁻¹ with step 0.05)

#### Multipath Feature

For advanced path selection based on scattering type and distance, use the multipath feature via TOML configuration:

```toml
[averaging.multipath]
paths = ["U-O", "U-O-O"]     # Select specific scattering paths
max_distance = 3.0            # Maximum path distance in Angstroms
num_processes = 4             # Parallel processing
```

This feature allows filtering paths by:
- **Path type**: Full atom notation (e.g., "U-O" for single scattering, "U-O-O" for double scattering)
- **Distance**: Maximum effective path length (reff)

See `examples/multipath_averaging.toml` and `docs/multipath_feature.md` for detailed documentation.

### md-exafs-md-input-gen

Generate structure input files for MD simulations from CIF files.

```bash
md-exafs-md-input-gen --cp2k --input structure.cif --size 3,3,3 --output supercell.xyz
```

This tool reads CIF (Crystallographic Information File) format structures and generates supercells in formats suitable for MD simulations.

Options:
- `--cp2k`: Output in CP2K format (no header, element x y z)
- `--input`: Input CIF file path (required)
- `--size`: Supercell size as "nx,ny,nz" (e.g., "3,3,3" for 3×3×3) (required)
- `--output`: Output file path (required)

Example:
```bash
# Generate a 3×3×3 supercell of gold in CP2K format
md-exafs-md-input-gen --cp2k --input Au.cif --size 3,3,3 --output Au_3x3x3.xyz

# Generate a 2×2×2 supercell
md-exafs-md-input-gen --cp2k --input UO2.cif --size 2,2,2 --output UO2_2x2x2.xyz
```

The CP2K format output consists of:
- No header lines
- One atom per line
- Format: `Element  X_coordinate  Y_coordinate  Z_coordinate`
- Coordinates in scientific notation (e.g., 1.2345678901234567E+00)

## Configuration

### Processing Configuration

The main configuration file controls trajectory processing:

- **[system]**: Project name, trajectory file, output directory
- **[lattice]**: Lattice vectors and periodic boundary conditions
- **[atoms]**: Element to atom type mapping
- **[processing]**: Frame range, sampling, and parallelization settings
- **[feff]**: FEFF executable path and input header template

### Averaging Configuration

Separate configuration for chi averaging:

```toml
[averaging]
input_directory = "feff_calculations"
frame_range = [90000, 100000]
output_file = "chi_avg_100k.dat"
```

## Python API

```python
from md_exafs import load_config, process_trajectory, average_chi_files

# Load and process
config = load_config("config.toml")
process_trajectory(config)

# Average chi data
average_chi_files(
    input_dir="feff_calculations",
    frame_range=(90000, 100000),
    output_file="chi_avg.dat"
)
```

## HPC Usage

The package is designed for HPC environments:

1. **SLURM Integration**: Example script provided (modify for your HPC system)
2. **Parallel FEFF**: Each FEFF calculation runs independently
3. **Dynamic Load Balancing**: Workers process calculations from a shared queue
4. **Scalable**: Works efficiently from small local runs to large HPC jobs

## Examples

See the `examples/` directory for:

- `README.md` - Instructions for running the examples
- `uo2_single_frame.xyz` - Single-frame UO₂ trajectory for testing
- `uo2_config.toml` - UO₂ system configuration (configured for test data)
- `averaging_config.toml` - Configuration for chi(k) averaging
- `workflow_example.sh` - Complete workflow demonstration with local execution
- `python_example.py` - Python API usage
- `run_slurm.sh` - HPC batch script

## Requirements

- Python ≥ 3.8
- NumPy ≥ 1.20.0
- OVITO ≥ 3.0.0
- tqdm ≥ 4.60.0
- tomli ≥ 2.0.0 (Python < 3.11)
- ASE ≥ 3.22.0
- xraylarch (for FEFF path processing)

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.
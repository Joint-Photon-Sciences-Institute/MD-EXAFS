# MD-EXAFS Usage Guide

## Configuration Schema

### Processing Configuration

```toml
[system]
name = "string"              # Project name
trajectory_file = "path"     # Path to MD trajectory file
output_directory = "path"    # Where to save FEFF inputs

[lattice]
a = [float, float, float]    # Lattice vector a (Angstroms)
b = [float, float, float]    # Lattice vector b
c = [float, float, float]    # Lattice vector c
pbc = [bool, bool, bool]     # Periodic boundaries [x, y, z]

[atoms]
Element1 = integer           # Map element symbol to type ID
Element2 = integer           # Must match trajectory file

[processing]
start_frame = integer        # First frame to process
end_frame = integer          # Last frame (exclusive)
frame_step = integer         # Process every Nth frame
atoms_per_frame = integer    # Central atoms per frame
cutoff_radius = float        # Neighbor search radius (Angstroms)

[feff]
absorbing_element = "symbol" # Element that absorbs X-rays (required)
feff_executable = "path"     # Optional: custom FEFF path
header = """                 # FEFF input header template
multiline string
"""
```

### Averaging Configuration

```toml
[averaging]
input_directory = "path"     # Directory with chi.dat files
frame_range = [start, end]   # Frames to include in average
output_file = "path"         # Output file for averaged data
paths = [1, 2, 3, 4, 5, 6, 7, 8]  # Optional: specific FEFF paths to process
# Alternative: paths = "[1,8]"    # Range notation as string
```

## Command Line Interface

### Processing Trajectories

```bash
md-exafs-process [--create-tar] config.toml
```

Options:
- `--create-tar`: Create compressed archive of output

### Averaging Chi Data

```bash
# Using TOML configuration
md-exafs-average --config averaging_config.toml

# Using command line arguments
md-exafs-average --input-dir DIR --start FRAME --end FRAME --output FILE

# Process specific FEFF paths
md-exafs-average --input-dir DIR --start FRAME --end FRAME --paths "[1,8]" --output FILE
```

Options:
- `--paths`: Process specific FEFF paths instead of existing chi.dat files
  - Format: "[start,end]" for range (e.g., "[1,8]") or "[1,3,5,7]" for specific paths
  - Requires xraylarch package
  - Sums paths within each atom, then averages across atoms

## Workflow Example

1. **Prepare your trajectory file** - The package supports any format compatible with OVITO (tested with CP2K XYZ files)
2. **Create configuration** based on your system
3. **Process trajectory** to generate FEFF inputs
4. **Run FEFF calculations** locally or on HPC
5. **Average results** to get final χ(k) spectrum

### Running FEFF Calculations

#### Local Execution

For testing or small datasets, you can run FEFF calculations locally:

```bash
md-exafs-feff-local --base-dir feff_calculations --workers 4
```

Options:
- `--base-dir`: Directory containing FEFF inputs (default: feff_calculations)
- `--workers`: Number of parallel processes (REQUIRED)
- `--feff-path`: Path to FEFF executable (default: auto-detect bundled version)

**Important:** Each FEFF calculation (one per atom) runs independently. You can use multiple workers to parallelize these calculations, up to the total number of calculations.

#### HPC Execution

For large production runs, submit to your HPC queue:

```bash
sbatch examples/run_slurm.sh
```

### Parallelization Strategy

Both local and HPC execution use the same parallelization approach:

1. **All FEFF calculations are independent** - Each atom in each frame can be processed in parallel
2. **Dynamic work queue** - Workers take the next available calculation when ready
3. **No frame boundaries** - A worker might process atoms from different frames
4. **Optimal efficiency** - No worker sits idle if there's work available

The directory structure (`frame_X/atom_Y/`) is purely for organization, not for work distribution.

## Tips

- Start with a small frame range for testing
- Verify atom type mappings match your trajectory
- Check that lattice vectors are correct for your system
- Use appropriate cutoff radius for your element
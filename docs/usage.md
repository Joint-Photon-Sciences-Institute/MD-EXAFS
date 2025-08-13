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
# For Cp2k/XYZ format: Element symbol -> Type ID
Element1 = integer           # Map element symbol to type ID
Element2 = integer           # Must match trajectory file

# For LAMMPS format: Type ID -> Element symbol (note the quotes)
# "1" = "Zn"
# "2" = "Fe"
# "3" = "O"

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
  - Sums paths within each atom, then averages across atoms
  - Automatically interpolates output to standard k-grid (0 to 20 Å⁻¹ with step 0.05)

## Workflow Example

1. **Prepare your trajectory file** - The package supports any format compatible with OVITO (tested with CP2K XYZ files and LAMMPS dump files)
2. **Create configuration** based on your system (see LAMMPS example below)
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

## LAMMPS Example

When using LAMMPS dump files, the atom type mapping uses numeric IDs as keys:

```toml
[system]
name = "ZnFeO"
trajectory_file = "trajectory.lammpstrj"
output_directory = "feff_calculations"

[lattice]
a = [8.435505, 0.0, 0.0]
b = [0.0, 8.435505, 0.0]
c = [0.0, 0.0, 8.435505]
pbc = [true, true, true]

[atoms]
# LAMMPS format: Type ID -> Element symbol
# Type IDs must be quoted strings
"1" = "Zn"
"2" = "Fe"
"3" = "O"

[processing]
start_frame = 0
end_frame = 100
frame_step = 10
atoms_per_frame = 5
cutoff_radius = 6.0

[feff]
absorbing_element = "Fe"  # Use element symbol, not type ID
# ... rest of FEFF configuration
```

The package automatically detects LAMMPS format by checking if particle types have empty or numeric names, and applies the mapping accordingly.

## RDF Analysis Options

The RDF analysis module (`python -m md_exafs.rdf_analysis`) supports the following options:

```bash
# Basic RDF analysis
python -m md_exafs.rdf_analysis --rdf config.toml

# Save RDF data to columnar text file
python -m md_exafs.rdf_analysis --rdf config.toml --rdf-save rdf_data.txt
```

The `--rdf-save` option creates a tab-delimited text file with:
- First column: r (distance in Angstroms)
- Individual columns for each RDF pair (e.g., Fe-O, Zn-O, etc.)
- **Average_RDF**: Simple arithmetic average of all partial RDFs
- **Weighted_RDF_Xray**: Properly weighted RDF for X-ray scattering comparison
  - Weighted by atomic concentrations (c_i × c_j)
  - Weighted by X-ray form factors (f_i × f_j, approximated as atomic numbers)
  - Suitable for comparison with experimental X-ray PDF data

This format is convenient for further analysis in tools like Excel, Origin, or Python/NumPy.

**Note**: For accurate comparison with experimental X-ray or neutron PDF data, the Weighted_RDF_Xray column should be used rather than the simple average.

## Tips

- Start with a small frame range for testing
- Verify atom type mappings match your trajectory
- For LAMMPS: Ensure type IDs in config match those in your dump file
- Check that lattice vectors are correct for your system
- Use appropriate cutoff radius for your element
- Use `--rdf-save` to export RDF data for custom analysis
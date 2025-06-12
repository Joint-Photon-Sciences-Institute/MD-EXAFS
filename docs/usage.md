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
cores_per_node = integer     # For parallel distribution

[feff]
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
```

## Workflow Example

1. **Prepare your trajectory file** in XYZ format
2. **Create configuration** based on your system
3. **Process trajectory** to generate FEFF inputs
4. **Run FEFF calculations** on HPC
5. **Average results** to get final χ(k) spectrum

## Tips

- Start with a small frame range for testing
- Verify atom type mappings match your trajectory
- Check that lattice vectors are correct for your system
- Use appropriate cutoff radius for your element
- Adjust `cores_per_node` to match your HPC system
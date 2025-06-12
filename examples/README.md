# MD-EXAFS Examples

This directory contains example files and scripts to help you get started with MD-EXAFS.

## Quick Test Example

We provide a single-frame UO₂ trajectory for quick testing of the complete workflow.

### Files Included

- `uo2_single_frame.xyz` - Single frame from a UO₂ MD simulation (96 atoms)
- `uo2_config.toml` - Configuration file pre-configured for the test trajectory
- `averaging_config.toml` - Configuration for averaging the chi(k) results
- `workflow_example.sh` - Complete workflow script demonstrating all steps
- `python_example.py` - Example of using MD-EXAFS as a Python library
- `run_slurm.sh` - SLURM script for HPC execution (modify for your system)

### Running the Test Example

1. **Activate your conda environment:**
   ```bash
   conda activate md-exafs
   ```

2. **Process the trajectory to generate FEFF inputs:**
   ```bash
   md-exafs-process uo2_config.toml
   ```
   
   This will create a `feff_calculations` directory with 10 FEFF input files (10 uranium atoms from 1 frame).

3. **Check the number of working directories created:**
   ```bash
   ls feff_calculations/
   ```
   
   For this example with 10 atoms and 4 cores configured, you should see 4 working directories (working_0 through working_3).

4. **Run FEFF calculations locally:**
   ```bash
   md-exafs-feff-local --base-dir feff_calculations --workers 4
   ```
   
   **Important:** Set `--workers` equal to the number of working directories for optimal performance.

5. **Average the chi(k) results:**
   ```bash
   md-exafs-average --config averaging_config.toml
   ```
   
   This creates `chi_avg_example.dat` with the averaged EXAFS spectrum.

### Automated Workflow

You can run all steps automatically:

```bash
./workflow_example.sh
```

This script will:
- Process the single-frame trajectory
- Run FEFF calculations locally
- Average the results
- Report the status of each step

### Adapting for Your System

To use MD-EXAFS with your own data:

1. **Prepare your trajectory** in XYZ format
2. **Copy and modify `uo2_config.toml`**:
   - Update `trajectory_file` path
   - Set appropriate lattice vectors
   - Map your atom types
   - Adjust frame range and processing parameters
   - Modify the FEFF header for your element and edge

3. **Run the workflow** as shown above

### Notes

- The example uses only 1 frame for quick testing. Real calculations typically use thousands of frames.
- Local execution is suitable for testing and small datasets. For production runs with many frames, use HPC.
- The number of workers should match the number of working directories for optimal load balancing.
- FEFF calculations can take several minutes even for this small example.
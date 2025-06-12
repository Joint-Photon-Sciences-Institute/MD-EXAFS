# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MD-EXAFS is a computational framework for calculating Extended X-ray Absorption Fine Structure (EXAFS) spectra from Molecular Dynamics (MD) trajectories, specifically designed for UO₂ (Uranium Dioxide) systems. The project processes MD simulation data to generate FEFF input files and compute averaged EXAFS chi(k) spectra.

## Architecture

The workflow consists of three main stages:

1. **Input Generation** (`make_inputs.py`): Processes MD trajectory files to create FEFF input files for multiple uranium atoms across selected frames, implementing periodic boundary conditions (PBC) for accurate neighbor calculations.

2. **FEFF Calculations** (`run.sh`): SLURM script that executes FEFF8L calculations in parallel across 112 cores, processing all generated input files.

3. **Data Averaging** (`build_chi.py`): Collects and averages chi(k) data from completed FEFF calculations.

## Key Commands

### Generate FEFF Input Files
```bash
python reference_code_base/tools/make_inputs.py <trajectory.xyz> [options]

# Common options:
# --start-frame 90000  # Starting frame (default: 90000)
# --end-frame 100000   # Ending frame (default: 100000)
# --everyother 10      # Process every nth frame (default: 10)
# --num-u-atoms 10     # U atoms per frame (default: 10)
# --cores-per-node 112 # For directory distribution (default: 40)
# --create-tar         # Create tar.gz archive of results

# Example:
python reference_code_base/tools/make_inputs.py reference_code_base/md-trajectories/UO2-pos-1.xyz --cores-per-node 112
```

### Run FEFF Calculations (HPC cluster)
```bash
sbatch reference_code_base/slurm_scripts/run.sh
```

### Average Chi(k) Data
```bash
python reference_code_base/tools/build_chi.py
```

## Critical Implementation Details

### Periodic Boundary Conditions
- The project uses specific lattice vectors for UO₂ (defined in `make_inputs.py:69-78`)
- PBC is crucial for correct neighbor identification in FEFF calculations
- All atomic positions are wrapped using minimum image convention

### FEFF Configuration
- Uses FEFF8 Lite with L3 edge for uranium
- Cutoff radius: 10.0 Å for neighbor search
- RPATH: 6.5 Å for scattering paths
- Three potentials: U (central), U (neighbor), O

### Directory Structure
- Input files are distributed across `working_X` directories for parallel processing
- Each directory contains: `frame_Y/U_Z/feff.inp`
- After calculation: only `feff.inp`, `feff.out`, and `chi.dat` are retained

## Dependencies

- **Python**: ovito, numpy, tqdm (for trajectory processing)
- **FEFF8L**: Located at `reference_code_base/feff8_lite/feff85L.exe`
- **HPC**: SLURM scheduler for parallel execution

## Development Notes

- No formal test suite exists
- No configuration files; parameters are hardcoded or passed via CLI
- Frame range 90000-100000 is commonly used for production runs
- The project lacks package management files (requirements.txt, setup.py)
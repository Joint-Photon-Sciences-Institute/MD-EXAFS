#!/bin/bash
# Example workflow for MD-EXAFS processing with single frame

echo "MD-EXAFS Example Workflow (Single Frame)"
echo "========================================"
echo "Note: This example uses a single-frame trajectory for quick testing"
echo ""

# Step 1: Generate FEFF input files from trajectory
echo "Step 1: Processing trajectory..."
md-exafs-process uo2_config.toml

# Check if processing was successful
if [ $? -eq 0 ]; then
    echo "✓ Trajectory processing complete"
else
    echo "✗ Trajectory processing failed"
    exit 1
fi

# Step 2: Run FEFF calculations locally
echo ""
echo "Step 2: Running FEFF calculations locally..."
num_inputs=$(find feff_calculations -name "feff.inp" | wc -l)
echo "Found $num_inputs FEFF input files to process"

# Count working directories
num_working_dirs=$(ls -d feff_calculations/working_* 2>/dev/null | wc -l)
echo "Found $num_working_dirs working directories"

# Run locally with the CLI tool
echo "Running parallel FEFF calculations with $num_working_dirs workers..."
md-exafs-feff-local --base-dir feff_calculations --workers $num_working_dirs

if [ $? -eq 0 ]; then
    echo "✓ FEFF calculations complete"
else
    echo "✗ Some FEFF calculations failed"
fi

# Note about HPC option
echo ""
echo "Note: For large production runs, use HPC with: sbatch run_slurm.sh"

# Step 3: Average chi(k) data
echo ""
echo "Step 3: Averaging chi(k) data..."

# Check if chi.dat files exist
num_chi=$(find feff_calculations -name "chi.dat" | wc -l)
if [ $num_chi -gt 0 ]; then
    echo "Found $num_chi chi.dat files"
    
    # Run averaging
    md-exafs-average --config averaging_config.toml
    
    if [ -f "chi_avg_example.dat" ]; then
        echo "✓ Averaging complete: chi_avg_example.dat created"
    fi
else
    echo "No chi.dat files found - FEFF calculations may have failed"
fi

echo ""
echo "Workflow example complete!"
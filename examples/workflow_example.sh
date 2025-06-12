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

# Run locally with Python script
echo "Running parallel FEFF calculations..."
python run_feff_local.py --base-dir feff_calculations

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
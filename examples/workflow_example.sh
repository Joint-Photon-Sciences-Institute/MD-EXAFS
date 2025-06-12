#!/bin/bash
# Example workflow for MD-EXAFS processing

echo "MD-EXAFS Example Workflow"
echo "========================="

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

# Step 2: Run FEFF calculations (locally for testing)
echo ""
echo "Step 2: Running FEFF calculations..."
echo "Note: In production, submit run_slurm.sh to your HPC queue instead"

# For demonstration, just show what would be done
num_inputs=$(find feff_calculations -name "feff.inp" | wc -l)
echo "Found $num_inputs FEFF input files ready for processing"
echo "To run on HPC: sbatch run_slurm.sh"

# Step 3: Average chi(k) data (after FEFF calculations complete)
echo ""
echo "Step 3: Averaging chi(k) data..."
echo "This step should be run after all FEFF calculations are complete"

# Show both TOML and CLI options
echo ""
echo "Option A - Using TOML config:"
echo "  md-exafs-average --config averaging_config.toml"

echo ""
echo "Option B - Using CLI arguments:"
echo "  md-exafs-average --input-dir feff_calculations --start 90000 --end 100000 --output chi_avg.dat"

echo ""
echo "Workflow example complete!"
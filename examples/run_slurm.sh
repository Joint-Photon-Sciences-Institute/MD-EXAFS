#!/bin/bash
#SBATCH --job-name=MD-EXAFS
#SBATCH --nodes=1       
#SBATCH --ntasks-per-node=112
#SBATCH --mem=0         # Requested Memory
#SBATCH --time=5:00:00
#SBATCH -o job_%j.out   # %j = job ID
#SBATCH -e job_%j.err   # %j = job ID

# NOTE: This is an example SLURM script. You may need to modify:
# - Number of cores (--ntasks-per-node) to match your HPC system
# - Memory requirements
# - Time limits
# - Add partition/account information as needed
# - Module loads for your environment

NCORES=112

# Load conda environment (adjust path as needed)
# module load anaconda3
# conda activate md-exafs

# Directory containing all working directories
BASE_DIR="feff_calculations"

# Print configuration
echo "Running with ${NCORES} cores"

# Function to run FEFF calculation and cleanup
run_feff() {
    local dir=$1
    cd "$dir" || exit

    # Run FEFF using the configured executable
    # The path to FEFF can be configured in your TOML file
    # Default uses the bundled feff85L.exe
    feff85L.exe > feff.out

    # Important: Wait for the process to finish
    wait $!

    # Check if chi.dat was created
    if [ -f "chi.dat" ]; then
        # Keep only the essential files
        find . -type f ! -name 'feff.inp' ! -name 'feff.out' ! -name 'chi.dat' -delete
        echo "Completed FEFF calculation in $dir"
    else
        echo "Warning: FEFF calculation failed in $dir"
    fi

    cd - > /dev/null
}

# Find all directories containing feff.inp
mapfile -t FEFF_DIRS < <(find "$BASE_DIR" -type f -name "feff.inp" -exec dirname {} \;)
total_dirs=${#FEFF_DIRS[@]}
echo "Found $total_dirs directories to process"

# Initialize arrays for job tracking
declare -A active_jobs
current_jobs=0
completed_jobs=0

# Process all directories
for dir in "${FEFF_DIRS[@]}"; do
    # Wait if we're at max jobs
    while [ $current_jobs -ge $NCORES ]; do
        # Check each running job
        for pid in "${!active_jobs[@]}"; do
            if ! kill -0 $pid 2>/dev/null; then
                wait $pid
                ((completed_jobs++))
                ((current_jobs--))
                unset active_jobs[$pid]
                echo "Completed job $completed_jobs of $total_dirs (Running: $current_jobs)"
            fi
        done
        sleep 2
    done

    # Start new job
    (run_feff "$dir") &
    pid=$!
    active_jobs[$pid]=$dir
    ((current_jobs++))
    echo "Started FEFF calculation in $dir (Running: $current_jobs, Total completed: $completed_jobs)"
done

# Wait for all remaining jobs
echo "Waiting for remaining jobs to complete..."
for pid in "${!active_jobs[@]}"; do
    wait $pid
    dir=${active_jobs[$pid]}
    echo "Completed calculation in $dir"
done

echo "All FEFF calculations completed!"

# Verify results
total_inp=$(find "$BASE_DIR" -type f -name "feff.inp" | wc -l)
total_chi=$(find "$BASE_DIR" -type f -name "chi.dat" | wc -l)
echo "Total feff.inp files: $total_inp"
echo "Total chi.dat files: $total_chi"

if [ $total_inp -ne $total_chi ]; then
    echo "Warning: Number of chi.dat files ($total_chi) does not match number of input files ($total_inp)"
    echo "Failed calculations in directories:"
    find "$BASE_DIR" -type f -name "feff.inp" -exec sh -c '
        dir=$(dirname "{}")
        if [ ! -f "$dir/chi.dat" ]; then
            echo "$dir"
        fi
    ' \;
fi
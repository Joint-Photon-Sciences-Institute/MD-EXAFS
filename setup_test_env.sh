#!/bin/bash

echo "Setting up MD-EXAFS test environment..."

# Create conda environment
echo "Creating conda environment from environment.yml..."
conda env create -f environment.yml

echo ""
echo "To activate the environment, run:"
echo "  conda activate md-exafs"
echo ""
echo "Then install the package in development mode:"
echo "  pip install -e ."
echo ""
echo "To test the package:"
echo "  python -m md_exafs.trajectory_processor examples/uo2_config.toml"
echo "  python -m md_exafs.chi_averager --help"
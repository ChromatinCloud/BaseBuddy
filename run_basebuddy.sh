#!/bin/bash
# Script to run BaseBuddy GUI with the correct environment

echo "Starting BaseBuddy GUI..."

# Initialize conda from miniforge3
eval "$(~/miniforge3/bin/conda shell.bash hook)"

# Activate the basebuddy environment from miniconda3
conda activate /Users/lauferva/miniconda3/envs/basebuddy

# Verify environment
echo "Using Python: $(which python)"
echo "Python version: $(python --version)"

# Run the GUI
exec python run_gui.py
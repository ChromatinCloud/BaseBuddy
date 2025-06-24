#!/bin/bash
# Setup script for BaseBuddy GUI environment

echo "Setting up BaseBuddy GUI environment..."

# Create conda environment with Python 3.12
echo "Creating conda environment with Python 3.12..."
conda create -n bb_gui python=3.12 tk -y

# Activate the environment
echo "Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bb_gui

# Install Python dependencies
echo "Installing Python dependencies..."
pip install customtkinter typer pysam

echo ""
echo "Setup complete! To use the GUI:"
echo "  conda activate bb_gui"
echo "  python run_gui.py"
echo ""
echo "Or run directly:"
echo "  conda run -n bb_gui python run_gui.py"
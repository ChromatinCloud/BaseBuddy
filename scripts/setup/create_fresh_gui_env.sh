#!/bin/bash
# Create a fresh GUI environment for BaseBuddy

echo "Creating fresh BaseBuddy GUI environment..."
echo "========================================="

# Remove the broken environment
echo "Removing broken bb_gui environment..."
conda env remove -n bb_gui -y 2>/dev/null || true

# Create new environment with Python 3.12
echo "Creating new bb_gui environment..."
conda create -n bb_gui python=3.12 -y

# Activate and install packages
echo "Installing packages..."
conda activate bb_gui

# Install GUI dependencies
pip install customtkinter typer pysam

# Install bioinformatics tools separately
conda install -n bb_gui -c conda-forge -c bioconda samtools art gsl -y

# Install aligners last to avoid conflicts
conda install -n bb_gui -c conda-forge -c bioconda bwa -y
conda install -n bb_gui -c conda-forge -c bioconda minimap2 -y

echo ""
echo "Environment setup complete!"
echo ""
echo "To use the GUI:"
echo "  conda activate bb_gui"
echo "  python run_gui.py"
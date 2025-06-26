#!/bin/bash
# Script to install BAMSurgeon for variant spiking functionality

echo "Installing BAMSurgeon for BaseBuddy variant spiking..."
echo "======================================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: Conda is not installed or not in PATH"
    exit 1
fi

# Check if basebuddy-gui environment exists
if ! conda env list | grep -q "basebuddy-gui"; then
    echo "Error: basebuddy-gui conda environment not found"
    echo "Please create it first or activate the correct environment"
    exit 1
fi

# Install BAMSurgeon
echo ""
echo "Installing BAMSurgeon in basebuddy-gui environment..."
echo "This will install addsnv.py and related tools."
echo ""

# Run the installation
conda install -n basebuddy-gui -c bioconda bamsurgeon -y

# Check if installation was successful
echo ""
echo "Checking installation..."

# Activate environment and check for addsnv.py
source $(conda info --base)/etc/profile.d/conda.sh
conda activate basebuddy-gui

if command -v addsnv.py &> /dev/null; then
    echo "✓ Success! addsnv.py is now available"
    echo "  Location: $(which addsnv.py)"
    echo ""
    echo "BAMSurgeon installation complete!"
    echo "You can now use the Variant Spiking feature in BaseBuddy GUI."
else
    echo "✗ Installation may have failed. addsnv.py not found in PATH"
    echo "  Try running: conda activate basebuddy-gui && which addsnv.py"
fi
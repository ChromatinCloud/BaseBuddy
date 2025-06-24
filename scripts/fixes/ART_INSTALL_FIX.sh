#!/bin/bash
# Fix ART installation for BaseBuddy

echo "Installing ART with all dependencies..."

# Make sure we're in the bb_gui environment
conda activate bb_gui

# Install ART and GSL together from bioconda
conda install -c bioconda -c conda-forge art gsl -y

# Test the installation
echo "Testing ART installation..."
if command -v art_illumina &> /dev/null; then
    echo "✅ ART installed successfully!"
    art_illumina -h | head -5
else
    echo "❌ ART installation failed"
fi

echo ""
echo "To use the GUI:"
echo "  conda activate bb_gui"
echo "  python run_gui.py"
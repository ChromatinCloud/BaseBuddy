#!/bin/bash
# Install aligners for BaseBuddy auto-align feature

echo "Installing aligners for BaseBuddy..."
echo "=================================="

# Make sure we're in the bb_gui environment
echo "Activating bb_gui environment..."
conda activate bb_gui

# Install BWA
echo ""
echo "Installing BWA..."
conda install -c bioconda bwa -y

# Verify installation
echo ""
echo "Verifying installations..."
echo -n "BWA: "
if command -v bwa &> /dev/null; then
    bwa 2>&1 | head -1
else
    echo "NOT FOUND"
fi

echo -n "samtools: "
if command -v samtools &> /dev/null; then
    samtools --version | head -1
else
    echo "NOT FOUND - Installing..."
    conda install -c bioconda samtools -y
fi

echo ""
echo "Installation complete!"
echo ""
echo "To use auto-align in BaseBuddy:"
echo "1. conda activate bb_gui"
echo "2. python run_gui.py"
echo "3. Check 'Auto-align to BAM' in Short Read Sim tab"
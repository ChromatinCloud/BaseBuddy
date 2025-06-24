#!/bin/bash
# Fix aligner installation with proper channels

echo "Installing aligners with all required channels..."

# Add conda-forge channel if not already added
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults

# Install with all channels specified
conda install -c conda-forge -c bioconda bwa minimap2 -y

# If that fails, try installing libzlib first
if [ $? -ne 0 ]; then
    echo "Trying alternative approach..."
    conda install -c conda-forge libzlib -y
    conda install -c conda-forge -c bioconda bwa minimap2 -y
fi

echo "Installation complete!"
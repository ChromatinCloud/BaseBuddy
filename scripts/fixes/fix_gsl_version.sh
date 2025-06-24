#!/bin/bash
# Fix GSL version mismatch for ART

echo "Fixing GSL library version mismatch..."

# Create symlink for GSL version compatibility
cd $CONDA_PREFIX/lib
ln -sf libgsl.27.dylib libgsl.25.dylib

echo "GSL symlink created."
echo ""
echo "Testing ART..."
art_illumina -h > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✅ ART is now working!"
else
    echo "❌ Still having issues. Try reinstalling ART:"
    echo "  conda remove art -y"
    echo "  conda install -c bioconda art=2016.06.05 -y"
fi
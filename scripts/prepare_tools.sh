#!/bin/bash
# Prepare external tools for PyInstaller packaging

set -e

echo "Preparing external tools for BaseBuddy packaging"
echo "=============================================="
echo ""

# Create tools directory
mkdir -p tools

# Function to copy tool if it exists
copy_tool() {
    local tool=$1
    local tool_path=$(which $tool 2>/dev/null || echo "")
    
    if [ -n "$tool_path" ]; then
        echo "Copying $tool from $tool_path"
        cp "$tool_path" tools/
        # Also copy any adjacent files (for tools with multiple executables)
        local tool_dir=$(dirname "$tool_path")
        local tool_base=$(basename "$tool" | cut -d. -f1)
        for related in "$tool_dir/$tool_base"*; do
            if [ -f "$related" ] && [ "$related" != "$tool_path" ]; then
                echo "  Also copying $(basename $related)"
                cp "$related" tools/
            fi
        done
    else
        echo "Warning: $tool not found in current environment"
    fi
}

# Ensure we're in the packaging environment
if [[ "$CONDA_DEFAULT_ENV" != "bb-pack" ]]; then
    echo "Error: Please activate the bb-pack environment first"
    echo "Run: conda activate bb-pack"
    exit 1
fi

echo "Copying tools from conda environment..."
echo ""

# BamSurgeon tools
copy_tool addsnv.py
copy_tool addindel.py
copy_tool addsv.py

# Core tools
copy_tool bwa
copy_tool samtools
copy_tool bcftools
copy_tool bedtools

# Read simulators
copy_tool art_illumina
copy_tool wgsim
copy_tool nanosim
copy_tool nanosim-h

# Variant manipulation dependencies
copy_tool exonerate
copy_tool velvet
copy_tool velveth
copy_tool velvetg

# Other tools
copy_tool pigz
copy_tool minimap2

# Copy SimuG if available
if [ -f "external/simuG/simuG.pl" ]; then
    echo "Copying SimuG..."
    cp external/simuG/simuG.pl tools/
    # Also copy SimuG modules if they exist
    if [ -d "external/simuG/modules" ]; then
        cp -r external/simuG/modules tools/
    fi
else
    echo "Warning: SimuG not found at external/simuG/simuG.pl"
fi

# Make all tools executable
chmod +x tools/* 2>/dev/null || true

# Copy any required libraries
echo ""
echo "Checking for required dynamic libraries..."

# Function to copy dylibs for a binary
copy_dylibs() {
    local binary=$1
    if [[ -f "tools/$binary" ]]; then
        echo "Checking dependencies for $binary..."
        # Get non-system dylibs
        otool -L "tools/$binary" 2>/dev/null | grep -v "/usr/lib" | grep -v "/System" | \
            awk '{print $1}' | while read lib; do
            if [[ -f "$lib" ]] && [[ ! -f "tools/$(basename $lib)" ]]; then
                echo "  Copying library: $(basename $lib)"
                cp "$lib" tools/
            fi
        done
    fi
}

# Check key binaries for dependencies
for tool in bwa samtools exonerate wgsim; do
    copy_dylibs $tool
done

echo ""
echo "Tools prepared in ./tools directory"
echo ""
echo "Contents:"
ls -la tools/ | head -20
echo ""
echo "Total size: $(du -sh tools | cut -f1)"
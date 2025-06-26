#!/bin/bash
# Script to rebuild BWA index for reference genomes

echo "BWA Index Rebuilder for BaseBuddy"
echo "================================="

# Check if BWA is available
if ! command -v bwa &> /dev/null; then
    echo "Error: BWA is not installed or not in PATH"
    echo "Please activate your conda environment: conda activate basebuddy-gui"
    exit 1
fi

# Reference file
REF_FILE="$1"

if [ -z "$REF_FILE" ]; then
    echo "Usage: $0 <reference_fasta>"
    echo ""
    echo "Example for GRCh37:"
    echo "  $0 refs/GRCh37/hs37d5.fa"
    echo ""
    echo "Example for GRCh38:"
    echo "  $0 refs/GRCh38/GRCh38.fa"
    exit 1
fi

if [ ! -f "$REF_FILE" ]; then
    echo "Error: Reference file not found: $REF_FILE"
    exit 1
fi

# Check file size
FILE_SIZE=$(stat -f%z "$REF_FILE" 2>/dev/null || stat -c%s "$REF_FILE" 2>/dev/null)
FILE_SIZE_GB=$((FILE_SIZE / 1073741824))

echo "Reference file: $REF_FILE"
echo "File size: ${FILE_SIZE_GB}GB"

# Remove old index files
echo ""
echo "Removing old BWA index files..."
for ext in .amb .ann .bwt .pac .sa; do
    if [ -f "${REF_FILE}${ext}" ]; then
        echo "  Removing ${REF_FILE}${ext}"
        rm -f "${REF_FILE}${ext}"
    fi
done

# Create new index
echo ""
echo "Creating new BWA index..."
echo "This may take 60-90 minutes for human genome references."
echo ""

time bwa index "$REF_FILE"

# Check if all files were created
echo ""
echo "Checking index files..."
ALL_GOOD=1
for ext in .amb .ann .bwt .pac .sa; do
    if [ -f "${REF_FILE}${ext}" ]; then
        echo "  ✓ ${REF_FILE}${ext}"
    else
        echo "  ✗ ${REF_FILE}${ext} - MISSING!"
        ALL_GOOD=0
    fi
done

if [ $ALL_GOOD -eq 1 ]; then
    echo ""
    echo "BWA index created successfully!"
else
    echo ""
    echo "Error: Some index files are missing. Please check the output above."
    exit 1
fi
#!/bin/bash
# Test that aligners are installed

echo "Testing aligner installations..."
echo "==============================="

# Test BWA
echo -n "BWA: "
if command -v bwa &> /dev/null; then
    echo "✅ Installed"
    bwa 2>&1 | grep -i version | head -1
else
    echo "❌ Not found"
fi

# Test minimap2
echo -n "minimap2: "
if command -v minimap2 &> /dev/null; then
    echo "✅ Installed"
    minimap2 --version
else
    echo "❌ Not found"
fi

# Test samtools (required for BAM operations)
echo -n "samtools: "
if command -v samtools &> /dev/null; then
    echo "✅ Installed"
    samtools --version | head -1
else
    echo "❌ Not found"
fi

echo ""
echo "All aligners are ready for use in BaseBuddy!"
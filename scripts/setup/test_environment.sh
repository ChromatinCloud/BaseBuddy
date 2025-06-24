#!/bin/bash
# Test that the environment can be created successfully

echo "Testing BaseBuddy environment creation..."
echo "========================================"

# Remove any existing test environment
conda env remove -n basebuddy-test -y 2>/dev/null || true

# Try to create from environment-gui.yml
echo "Creating environment from environment-gui.yml..."
if conda env create -f environment-gui.yml -n basebuddy-test; then
    echo "✅ Environment created successfully!"
    
    # Test imports
    echo ""
    echo "Testing Python imports..."
    conda run -n basebuddy-test python -c "
import sys
print(f'Python: {sys.version}')
try:
    import customtkinter
    print('✅ customtkinter')
except ImportError as e:
    print(f'❌ customtkinter: {e}')
try:
    import typer
    print('✅ typer')
except ImportError as e:
    print(f'❌ typer: {e}')
try:
    import pysam
    print('✅ pysam')
except ImportError as e:
    print(f'❌ pysam: {e}')
"
    
    # Test tools
    echo ""
    echo "Testing tools..."
    conda run -n basebuddy-test which bwa && echo "✅ bwa" || echo "❌ bwa"
    conda run -n basebuddy-test which samtools && echo "✅ samtools" || echo "❌ samtools"
    conda run -n basebuddy-test which art_illumina && echo "✅ art" || echo "❌ art"
    
    # Clean up
    conda env remove -n basebuddy-test -y
else
    echo "❌ Failed to create environment"
    exit 1
fi

echo ""
echo "Environment file validated successfully!"
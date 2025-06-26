#!/bin/bash
# Quick setup script for development environment

echo "BaseBuddy Development Environment Setup"
echo "======================================"
echo ""

# Check what environments exist
echo "Current conda environments:"
conda env list
echo ""

# Ask user which approach they want
echo "Which setup would you like?"
echo "1) Unified environment (all tools, some via Rosetta) - RECOMMENDED"
echo "2) Hybrid ARM/x86 environments (native performance where possible)"
echo "3) Simple environment (from original environment.yml)"
echo ""
read -p "Enter choice (1-3): " choice

case $choice in
    1)
        echo ""
        echo "Setting up unified environment..."
        if [ -f "scripts/setup_unified_env.sh" ]; then
            ./scripts/setup_unified_env.sh
        else
            echo "Creating unified environment directly..."
            mamba env create -f environment-unified.yml || conda env create -f environment-unified.yml
        fi
        echo ""
        echo "✓ Setup complete!"
        echo "To use: conda activate basebuddy"
        echo "Then: python run_gui.py"
        ;;
    
    2)
        echo ""
        echo "Setting up hybrid environments..."
        if [ -f "scripts/setup_hybrid_env.sh" ]; then
            ./scripts/setup_hybrid_env.sh
        else
            echo "Error: setup_hybrid_env.sh not found"
            exit 1
        fi
        echo ""
        echo "✓ Setup complete!"
        echo "To use: source activate_basebuddy.sh"
        ;;
    
    3)
        echo ""
        echo "Setting up from original environment.yml..."
        if [ -f "environment.yml" ]; then
            mamba env create -f environment.yml || conda env create -f environment.yml
        else
            echo "Error: environment.yml not found"
            exit 1
        fi
        echo ""
        echo "✓ Setup complete!"
        echo "To use: conda activate basebuddy"
        ;;
    
    *)
        echo "Invalid choice"
        exit 1
        ;;
esac

echo ""
echo "Next steps:"
echo "1. Activate the environment (see above)"
echo "2. Run: python run_gui.py"
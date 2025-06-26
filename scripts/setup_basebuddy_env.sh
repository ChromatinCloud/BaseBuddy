#!/bin/bash
# Complete setup script for BaseBuddy environment on macOS (Intel & Apple Silicon)
# Based on comprehensive guide for reproducible Conda environments

set -e  # Exit on error

echo "BaseBuddy Environment Setup"
echo "==========================="
echo ""

# Color codes for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}✓${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

print_error() {
    echo -e "${RED}✗${NC} $1"
}

# Step 1: Check for Conda/Miniforge
echo "Step 1: Checking for Conda installation..."
if ! command -v conda &> /dev/null; then
    print_error "Conda not found! Please install Miniforge first."
    echo ""
    echo "To install Miniforge on Apple Silicon:"
    echo "1. Download from: https://github.com/conda-forge/miniforge#download"
    echo "2. Run: bash Miniforge3-MacOSX-arm64.sh"
    echo "3. Restart your terminal"
    exit 1
else
    CONDA_VERSION=$(conda --version)
    print_status "Found $CONDA_VERSION"
fi

# Step 2: Configure Conda channels
echo ""
echo "Step 2: Configuring Conda channels..."

# Check current channels
CURRENT_CHANNELS=$(conda config --show channels 2>/dev/null | grep -E "conda-forge|bioconda" | wc -l)

if [ "$CURRENT_CHANNELS" -lt 2 ]; then
    print_warning "Setting up conda-forge and bioconda channels..."
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --set channel_priority strict
    print_status "Channels configured with strict priority"
else
    print_status "Channels already configured"
fi

# Step 3: Install mamba for faster resolution (optional)
echo ""
echo "Step 3: Checking for mamba..."
if ! command -v mamba &> /dev/null; then
    print_warning "Mamba not found. Installing for faster dependency resolution..."
    conda install mamba -n base -c conda-forge -y
    print_status "Mamba installed"
else
    print_status "Mamba already installed"
fi

# Step 4: Create environment from YAML
echo ""
echo "Step 4: Creating BaseBuddy environment..."

ENV_FILE="environment-complete.yml"
if [ ! -f "$ENV_FILE" ]; then
    ENV_FILE="environment.yml"
fi

if [ ! -f "$ENV_FILE" ]; then
    print_error "No environment.yml file found!"
    exit 1
fi

# Check if environment already exists
if conda env list | grep -q "^basebuddy "; then
    print_warning "Environment 'basebuddy' already exists."
    read -p "Do you want to remove and recreate it? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n basebuddy -y
        print_status "Removed existing environment"
    else
        echo "Exiting without changes."
        exit 0
    fi
fi

# Create environment using mamba if available, otherwise conda
print_status "Creating environment from $ENV_FILE..."
if command -v mamba &> /dev/null; then
    mamba env create -f "$ENV_FILE"
else
    conda env create -f "$ENV_FILE"
fi

# Step 5: Post-installation setup
echo ""
echo "Step 5: Post-installation setup..."

# Activate the environment for setup
source $(conda info --base)/etc/profile.d/conda.sh
conda activate basebuddy

# Install SimuG if not present
if ! command -v simuG.pl &> /dev/null; then
    print_warning "Installing SimuG..."
    if [ ! -d "external/simuG" ]; then
        mkdir -p external
        cd external
        git clone https://github.com/yjx1217/simuG.git
        cd ..
    fi
    # Add to PATH or create symlink
    if [ -d "external/simuG" ]; then
        chmod +x external/simuG/simuG.pl
        print_status "SimuG cloned to external/simuG/"
        echo "Note: Add $(pwd)/external/simuG to your PATH or use full path to simuG.pl"
    fi
fi

# Step 6: Verification
echo ""
echo "Step 6: Verifying installation..."

# Check critical tools
TOOLS_TO_CHECK="python bwa samtools art_illumina"
ALL_GOOD=true

for tool in $TOOLS_TO_CHECK; do
    if command -v $tool &> /dev/null; then
        # Check if it's ARM64 native (on Apple Silicon)
        if [[ $(uname -m) == "arm64" ]]; then
            FILE_INFO=$(file $(which $tool) 2>/dev/null | grep -o "arm64\|x86_64" | head -1)
            if [[ "$FILE_INFO" == "arm64" ]]; then
                print_status "$tool ✓ (ARM64 native)"
            elif [[ "$FILE_INFO" == "x86_64" ]]; then
                print_warning "$tool ✓ (x86_64 - will use Rosetta 2)"
            else
                print_status "$tool ✓"
            fi
        else
            print_status "$tool ✓"
        fi
    else
        print_error "$tool ✗ (not found)"
        ALL_GOOD=false
    fi
done

# Check Python packages
echo ""
echo "Checking Python packages..."
python -c "import pysam; print(f'  ✓ pysam {pysam.__version__}')" 2>/dev/null || print_error "pysam not found"
python -c "import customtkinter; print('  ✓ customtkinter')" 2>/dev/null || print_error "customtkinter not found"

# Step 7: Create activation script
echo ""
echo "Step 7: Creating convenience scripts..."

cat > activate_basebuddy.sh << 'EOF'
#!/bin/bash
# Convenience script to activate BaseBuddy environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate basebuddy
echo "BaseBuddy environment activated!"
echo "To run the GUI: python run_gui.py"
echo "To run the CLI: basebuddy --help"
EOF

chmod +x activate_basebuddy.sh
print_status "Created activate_basebuddy.sh"

# Step 8: Final instructions
echo ""
echo "================================"
echo "Setup Complete!"
echo "================================"
echo ""
echo "To activate the environment:"
echo "  conda activate basebuddy"
echo ""
echo "Or use the convenience script:"
echo "  source ./activate_basebuddy.sh"
echo ""

if [ "$ALL_GOOD" = true ]; then
    print_status "All critical tools verified successfully!"
else
    print_warning "Some tools are missing. Check the output above."
fi

echo ""
echo "Optional: To integrate with Jupyter:"
echo "  python -m ipykernel install --user --name basebuddy --display-name 'BaseBuddy'"
echo ""
echo "For troubleshooting, see: docs/environment_setup.md"
#!/bin/bash
# Unified environment setup for seamless GUI experience
# Creates a single environment with all tools (some via Rosetta 2)

set -e

echo "BaseBuddy Unified Environment Setup"
echo "==================================="
echo ""
echo "This setup creates a single environment where all GUI features work"
echo "seamlessly. Some tools will run via Rosetta 2 on Apple Silicon."
echo ""

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

print_status() { echo -e "${GREEN}✓${NC} $1"; }
print_warning() { echo -e "${YELLOW}⚠${NC} $1"; }
print_error() { echo -e "${RED}✗${NC} $1"; }

# Check prerequisites
if ! command -v conda &> /dev/null; then
    print_error "Conda not found! Please install Miniforge first."
    echo "Run: ./scripts/install_miniforge.sh"
    exit 1
fi

# Configure channels
echo "Configuring conda channels..."
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict

# Install mamba if needed
if ! command -v mamba &> /dev/null; then
    print_warning "Installing mamba for faster dependency resolution..."
    conda install -n base -c conda-forge mamba -y
fi

# Remove existing environment if present
if conda env list | grep -q "^basebuddy "; then
    print_warning "Removing existing basebuddy environment..."
    conda env remove -n basebuddy -y
fi

# Create unified environment
echo ""
echo "Creating unified environment..."
echo "This may take 10-20 minutes..."

# For Apple Silicon, we need to allow x86_64 packages
if [[ $(uname -m) == "arm64" ]]; then
    print_warning "Apple Silicon detected. Some tools will use Rosetta 2."
    export CONDA_SUBDIR=osx-64
fi

mamba env create -f environment-unified.yml

# Reset subdir for future operations
unset CONDA_SUBDIR

# Install SimuG
echo ""
echo "Setting up SimuG..."
if [ ! -d "external/simuG" ]; then
    mkdir -p external
    cd external
    git clone https://github.com/yjx1217/simuG.git
    cd ..
    chmod +x external/simuG/simuG.pl
fi

# Create launcher script
cat > run_basebuddy_gui.sh << 'EOF'
#!/bin/bash
# Launch BaseBuddy GUI with all features enabled

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate basebuddy

# Add SimuG to PATH
if [ -d "external/simuG" ]; then
    export PATH="$PWD/external/simuG:$PATH"
fi

# Launch GUI
echo "Starting BaseBuddy GUI..."
echo "All features should work seamlessly."
python run_gui.py
EOF

chmod +x run_basebuddy_gui.sh

# Create desktop app (macOS)
if [[ $(uname) == "Darwin" ]]; then
    mkdir -p BaseBuddy.app/Contents/MacOS
    cat > BaseBuddy.app/Contents/MacOS/BaseBuddy << 'EOF'
#!/bin/bash
cd "$(dirname "$0")/../../.."
./run_basebuddy_gui.sh
EOF
    chmod +x BaseBuddy.app/Contents/MacOS/BaseBuddy
    
    cat > BaseBuddy.app/Contents/Info.plist << 'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleExecutable</key>
    <string>BaseBuddy</string>
    <key>CFBundleIdentifier</key>
    <string>com.basebuddy.gui</string>
    <key>CFBundleName</key>
    <string>BaseBuddy</string>
    <key>CFBundleVersion</key>
    <string>1.0</string>
</dict>
</plist>
EOF
    print_status "Created BaseBuddy.app - you can double-click to launch!"
fi

echo ""
echo "========================================"
echo "Setup Complete!"
echo "========================================"
echo ""
echo "To run the GUI with all features:"
echo "  ./run_basebuddy_gui.sh"
echo ""
if [[ $(uname) == "Darwin" ]]; then
    echo "Or double-click BaseBuddy.app"
    echo ""
fi
echo "Note: First launch may be slow as Rosetta 2 initializes."
echo "Subsequent launches will be faster."
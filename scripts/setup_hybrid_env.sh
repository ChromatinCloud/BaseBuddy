#!/bin/bash
# Hybrid environment setup for BaseBuddy on Apple Silicon
# Creates both ARM64 and x86_64 environments

set -e

echo "BaseBuddy Hybrid Environment Setup for Apple Silicon"
echo "===================================================="
echo ""

# Color codes
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

print_status() { echo -e "${GREEN}✓${NC} $1"; }
print_warning() { echo -e "${YELLOW}⚠${NC} $1"; }
print_error() { echo -e "${RED}✗${NC} $1"; }
print_info() { echo -e "${BLUE}ℹ${NC} $1"; }

# Check for Apple Silicon
if [[ $(uname -m) != "arm64" ]]; then
    print_warning "This script is optimized for Apple Silicon. Detected: $(uname -m)"
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    [[ ! $REPLY =~ ^[Yy]$ ]] && exit 1
fi

# Step 1: Prerequisites check
echo "Step 1: Checking prerequisites..."

# Check Conda
if ! command -v conda &> /dev/null; then
    print_error "Conda not found! Please install Miniforge first."
    echo "Download from: https://github.com/conda-forge/miniforge#download"
    exit 1
fi
print_status "Conda found: $(conda --version)"

# Check Homebrew
if ! command -v brew &> /dev/null; then
    print_warning "Homebrew not found. Installing..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
fi
print_status "Homebrew found: $(brew --version | head -1)"

# Step 2: Configure Conda
echo ""
echo "Step 2: Configuring Conda channels..."
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
print_status "Channels configured"

# Install mamba if not present
if ! command -v mamba &> /dev/null; then
    print_info "Installing mamba for faster dependency resolution..."
    conda install -n base -c conda-forge mamba -y
fi

# Step 3: Create ARM64 environment
echo ""
echo "Step 3: Creating ARM64-native environment..."

if conda env list | grep -q "^basebuddy-arm "; then
    print_warning "basebuddy-arm environment exists. Removing..."
    conda env remove -n basebuddy-arm -y
fi

mamba env create -f environment-arm.yml
print_status "ARM64 environment created"

# Step 4: Create x86_64 environment
echo ""
echo "Step 4: Creating x86_64 environment for legacy tools..."

if conda env list | grep -q "^basebuddy-x86 "; then
    print_warning "basebuddy-x86 environment exists. Removing..."
    conda env remove -n basebuddy-x86 -y
fi

print_info "Creating x86_64 environment (will use Rosetta 2)..."
CONDA_SUBDIR=osx-64 mamba create -n basebuddy-x86 python=3.10 -y
conda activate basebuddy-x86
conda config --env --set subdir osx-64
mamba env update -f environment-x86.yml
conda deactivate
print_status "x86_64 environment created"

# Step 5: Install Homebrew packages
echo ""
echo "Step 5: Installing native ARM64 tools via Homebrew..."

BREW_PACKAGES=(
    "exonerate"     # Some ARM64 support
    "blast"         # Optional, for sequence analysis
)

for pkg in "${BREW_PACKAGES[@]}"; do
    if brew list --formula | grep -q "^${pkg}$"; then
        print_status "$pkg already installed"
    else
        print_info "Installing $pkg..."
        brew install "$pkg" || print_warning "Failed to install $pkg"
    fi
done

# Step 6: Install SimuG
echo ""
echo "Step 6: Setting up SimuG..."

if [ ! -d "external/simuG" ]; then
    mkdir -p external
    cd external
    git clone https://github.com/yjx1217/simuG.git
    cd ..
    chmod +x external/simuG/simuG.pl
    print_status "SimuG installed to external/simuG/"
else
    print_status "SimuG already installed"
fi

# Step 7: Create helper scripts
echo ""
echo "Step 7: Creating helper scripts..."

# Main activation script
cat > activate_basebuddy.sh << 'EOF'
#!/bin/bash
# Activate BaseBuddy environment with architecture detection

print_env_info() {
    echo "Environment: $1"
    echo "Architecture: $2"
    echo "Tools available: $3"
    echo ""
}

if [[ "$1" == "x86" ]] || [[ "$1" == "legacy" ]]; then
    conda activate basebuddy-x86
    print_env_info "basebuddy-x86" "x86_64 (Rosetta 2)" "BAMSurgeon, exonerate, velvet, wgsim, NanoSim"
elif [[ "$1" == "arm" ]] || [[ "$1" == "native" ]] || [[ -z "$1" ]]; then
    conda activate basebuddy-arm
    print_env_info "basebuddy-arm" "ARM64 (Native)" "BWA, samtools, ART, most Python tools"
else
    echo "Usage: source activate_basebuddy.sh [arm|x86]"
    echo "  arm (default) - Native ARM64 environment"
    echo "  x86          - Legacy tools via Rosetta 2"
    return 1
fi

# Add SimuG to PATH if available
if [ -d "$(pwd)/external/simuG" ]; then
    export PATH="$(pwd)/external/simuG:$PATH"
fi
EOF

# Tool wrapper for automatic environment switching
cat > basebuddy-tool << 'EOF'
#!/bin/bash
# Smart wrapper that switches environments based on tool

TOOL="$1"
shift

# Tools that need x86_64 environment
X86_TOOLS="addsnv.py addindel.py addsv.py exonerate velvet wgsim nanosim"

# Check if tool needs x86 environment
for x86_tool in $X86_TOOLS; do
    if [[ "$TOOL" == "$x86_tool" ]]; then
        echo "Switching to x86_64 environment for $TOOL..."
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate basebuddy-x86
        exec "$TOOL" "$@"
    fi
done

# Otherwise use ARM environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate basebuddy-arm
exec "$TOOL" "$@"
EOF

chmod +x activate_basebuddy.sh basebuddy-tool
print_status "Helper scripts created"

# Step 8: Verification
echo ""
echo "Step 8: Verifying installation..."

# Test ARM environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate basebuddy-arm

echo "ARM64 Environment:"
for tool in python bwa samtools art_illumina; do
    if command -v $tool &> /dev/null; then
        arch=$(file $(which $tool) 2>/dev/null | grep -o "arm64\|x86_64" | head -1)
        print_status "$tool (${arch:-unknown})"
    else
        print_error "$tool not found"
    fi
done

# Test x86 environment
conda activate basebuddy-x86
echo ""
echo "x86_64 Environment:"
for tool in addsnv.py exonerate wgsim; do
    if command -v $tool &> /dev/null; then
        print_status "$tool (x86_64/Rosetta)"
    else
        print_error "$tool not found"
    fi
done

conda deactivate

# Final instructions
echo ""
echo "========================================"
echo "Setup Complete!"
echo "========================================"
echo ""
echo "Usage:"
echo "  Native ARM64 environment (default):"
echo "    source activate_basebuddy.sh"
echo ""
echo "  Legacy x86_64 environment:"
echo "    source activate_basebuddy.sh x86"
echo ""
echo "  Smart tool wrapper:"
echo "    ./basebuddy-tool <command> [args]"
echo ""
echo "Examples:"
echo "  # Use native BWA"
echo "  source activate_basebuddy.sh"
echo "  bwa mem ref.fa reads.fq > output.sam"
echo ""
echo "  # Use BAMSurgeon (automatic env switch)"
echo "  ./basebuddy-tool addsnv.py -v chr1,100,A,T -f input.bam -o output.bam"
echo ""
print_info "SimuG available at: external/simuG/simuG.pl"
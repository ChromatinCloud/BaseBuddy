#!/bin/bash
# Verify packaging setup before building

echo "BaseBuddy Packaging Setup Verification"
echo "====================================="
echo ""

ERRORS=0
WARNINGS=0

# Color codes
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

check_pass() { echo -e "${GREEN}✓${NC} $1"; }
check_warn() { echo -e "${YELLOW}⚠${NC} $1"; ((WARNINGS++)); }
check_fail() { echo -e "${RED}✗${NC} $1"; ((ERRORS++)); }

echo "1. Checking prerequisites..."
echo "----------------------------"

# Check OS
if [[ $(uname) == "Darwin" ]]; then
    check_pass "Running on macOS"
else
    check_fail "Not running on macOS"
fi

# Check Python
if command -v python3 &> /dev/null; then
    check_pass "Python3 found: $(python3 --version)"
else
    check_fail "Python3 not found"
fi

# Check conda
if command -v conda &> /dev/null; then
    check_pass "Conda found: $(conda --version)"
else
    check_fail "Conda not found"
fi

# Check for GUI file
if [ -f "run_gui.py" ]; then
    check_pass "run_gui.py exists"
else
    check_fail "run_gui.py not found"
fi

# Check for spec file
if [ -f "basebuddy.spec" ]; then
    check_pass "basebuddy.spec exists"
else
    check_fail "basebuddy.spec not found"
fi

echo ""
echo "2. Checking Python imports..."
echo "-----------------------------"

# Test basic imports
python3 -c "import sys; print(f'Python path: {sys.executable}')" 2>/dev/null || check_fail "Python broken"

# Check if we can import tkinter
python3 -c "import tkinter" 2>/dev/null && check_pass "tkinter available" || check_fail "tkinter not available"

# Check source structure
if [ -d "src/basebuddy" ]; then
    check_pass "Source directory structure correct"
    
    # Check for key files
    [ -f "src/basebuddy/__init__.py" ] && check_pass "__init__.py exists" || check_warn "Missing __init__.py"
    [ -f "src/basebuddy/gui/main_app.py" ] && check_pass "main_app.py exists" || check_fail "Missing main_app.py"
    [ -d "src/basebuddy/data/signatures" ] && check_pass "Signature data exists" || check_warn "Missing signature data"
else
    check_fail "Source directory 'src/basebuddy' not found"
fi

echo ""
echo "3. Checking build environment..."
echo "--------------------------------"

# Check if bb-pack environment exists
if conda env list | grep -q "^bb-pack "; then
    check_pass "bb-pack environment exists"
    
    # Check if it's x86_64
    conda run -n bb-pack python -c "import platform; print(f'Architecture: {platform.machine()}')" 2>/dev/null
else
    check_warn "bb-pack environment not found (will be created during build)"
fi

echo ""
echo "4. Checking tools availability..."
echo "---------------------------------"

# Create tools directory if needed
mkdir -p tools

# Check for critical tools in current environment
for tool in bwa samtools; do
    if command -v $tool &> /dev/null; then
        check_pass "$tool available"
    else
        check_warn "$tool not found in current environment"
    fi
done

echo ""
echo "5. Checking potential issues..."
echo "-------------------------------"

# Check for spaces in path
if [[ "$PWD" =~ " " ]]; then
    check_fail "Current path contains spaces - this may cause issues"
else
    check_pass "No spaces in current path"
fi

# Check disk space
AVAILABLE=$(df -H . | awk 'NR==2 {print $4}' | sed 's/G//')
if (( $(echo "$AVAILABLE > 5" | bc -l) )); then
    check_pass "Sufficient disk space: ${AVAILABLE}G available"
else
    check_warn "Low disk space: ${AVAILABLE}G available (need ~5GB)"
fi

# Check for Rosetta on Apple Silicon
if [[ $(uname -m) == "arm64" ]]; then
    if /usr/bin/pgrep oahd >/dev/null 2>&1; then
        check_pass "Rosetta 2 installed (required for x86_64 build)"
    else
        check_fail "Rosetta 2 not installed - run: softwareupdate --install-rosetta"
    fi
fi

echo ""
echo "6. Testing GUI launch..."
echo "------------------------"

# Try to import the GUI
python3 -c "
import sys
sys.path.insert(0, 'src')
try:
    from basebuddy.gui.main_app import BaseBuddyGUI
    print('GUI imports successfully')
except ImportError as e:
    print(f'GUI import failed: {e}')
    sys.exit(1)
" && check_pass "GUI imports work" || check_fail "GUI import failed"

echo ""
echo "======================================="
echo "Verification Summary"
echo "======================================="
echo "Errors: $ERRORS"
echo "Warnings: $WARNINGS"
echo ""

if [ $ERRORS -eq 0 ]; then
    echo -e "${GREEN}Ready to build!${NC}"
    echo ""
    echo "Next steps:"
    echo "1. Run: ./scripts/build_app.sh"
    echo "2. Wait for build to complete (10-20 minutes)"
    echo "3. Test the app in dist/BaseBuddy.app"
else
    echo -e "${RED}Please fix errors before building${NC}"
    exit 1
fi

if [ $WARNINGS -gt 0 ]; then
    echo ""
    echo -e "${YELLOW}Note: Some warnings were found but build may still succeed${NC}"
fi
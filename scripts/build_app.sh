#!/bin/bash
# Build BaseBuddy.app for distribution

set -e

echo "BaseBuddy App Builder"
echo "===================="
echo ""

# Check if we're on macOS
if [[ $(uname) != "Darwin" ]]; then
    echo "Error: This script is for macOS only"
    exit 1
fi

# Step 1: Ensure Rosetta is installed
echo "Step 1: Checking Rosetta 2..."
if /usr/bin/pgrep oahd >/dev/null 2>&1; then
    echo "✓ Rosetta 2 is installed"
else
    echo "Installing Rosetta 2..."
    softwareupdate --install-rosetta --agree-to-license
fi

# Step 2: Create x86_64 packaging environment
echo ""
echo "Step 2: Setting up packaging environment..."

# Check if bb-pack environment exists
if conda env list | grep -q "^bb-pack "; then
    echo "✓ bb-pack environment exists"
else
    echo "Creating bb-pack environment..."
    export CONDA_SUBDIR=osx-64
    conda create -n bb-pack python=3.10 \
        pysam numpy scipy pandas customtkinter tk \
        bwa samtools minimap2 art \
        nanosim wgsim bamsurgeon exonerate velvet picard \
        biopython matplotlib pybedtools \
        bcftools bedtools pigz \
        -c conda-forge -c bioconda -y
fi

# Step 3: Prepare build in x86_64 mode
echo ""
echo "Step 3: Preparing build environment..."
echo "Switching to x86_64 mode and activating environment..."

# Create a temporary script to run in x86_64 mode
cat > build_in_x86.sh << 'EOF'
#!/bin/bash
set -e

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bb-pack

# Install PyInstaller if needed
if ! pip show pyinstaller >/dev/null 2>&1; then
    echo "Installing PyInstaller..."
    pip install pyinstaller
fi

# Prepare tools
echo ""
echo "Preparing external tools..."
./scripts/prepare_tools.sh

# Create resources directory and icon
echo ""
echo "Creating resources..."
mkdir -p resources

# Create a simple icon if it doesn't exist
if [ ! -f resources/basebuddy.icns ]; then
    echo "Creating placeholder icon..."
    # Create a simple PNG icon
    cat > create_icon.py << 'PYEOF'
from PIL import Image, ImageDraw, ImageFont
import os

# Create a simple icon
img = Image.new('RGBA', (512, 512), (52, 152, 219, 255))  # Blue background
draw = ImageDraw.Draw(img)

# Draw "BB" text
try:
    font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 200)
except:
    font = None

draw.text((128, 156), "BB", fill=(255, 255, 255, 255), font=font)

# Save as PNG
img.save("resources/basebuddy.png")

# Convert to ICNS (macOS icon format)
os.system("sips -s format icns resources/basebuddy.png --out resources/basebuddy.icns")
PYEOF
    
    python create_icon.py || echo "Warning: Could not create icon"
    rm -f create_icon.py resources/basebuddy.png
fi

# Update paths in main application
echo ""
echo "Updating application paths for bundled tools..."

# Create a launcher wrapper that sets up paths correctly
cat > run_gui_bundled.py << 'PYEOF'
#!/usr/bin/env python3
"""
Bundled launcher for BaseBuddy GUI
Sets up paths for bundled tools
"""

import os
import sys
from pathlib import Path

# Determine if we're running as a bundled app
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    # Running as bundled app
    bundle_dir = Path(sys._MEIPASS)
    tools_dir = bundle_dir / 'tools'
    
    # Add tools directory to PATH
    os.environ['PATH'] = f"{tools_dir}{os.pathsep}{os.environ.get('PATH', '')}"
    
    # Set environment variables for tools
    os.environ['BASEBUDDY_TOOLS_DIR'] = str(tools_dir)
    os.environ['BASEBUDDY_BUNDLED'] = '1'
else:
    # Running in development
    os.environ['BASEBUDDY_BUNDLED'] = '0'

# Import and run the main GUI
from run_gui import main

if __name__ == '__main__':
    main()
PYEOF

# Build the app
echo ""
echo "Building BaseBuddy.app..."
pyinstaller --clean --noconfirm basebuddy.spec

echo ""
echo "Build complete!"
EOF

chmod +x build_in_x86.sh

# Run the build in x86_64 mode
echo ""
echo "Running build in x86_64 mode..."
arch -x86_64 /bin/bash build_in_x86.sh

# Cleanup
rm -f build_in_x86.sh

# Step 4: Create DMG for distribution
echo ""
echo "Step 4: Creating DMG installer..."

if [ -d "dist/BaseBuddy.app" ]; then
    # Create a DMG
    echo "Creating DMG..."
    
    # Create temporary directory for DMG
    mkdir -p dist/dmg
    cp -r dist/BaseBuddy.app dist/dmg/
    ln -s /Applications dist/dmg/Applications
    
    # Create DMG
    hdiutil create -volname "BaseBuddy" \
        -srcfolder dist/dmg \
        -ov -format UDZO \
        dist/BaseBuddy.dmg
    
    # Cleanup
    rm -rf dist/dmg
    
    echo ""
    echo "✓ Created dist/BaseBuddy.dmg"
fi

# Step 5: Code signing (optional)
echo ""
echo "Step 5: Code signing..."

if security find-identity -p codesigning -v | grep -q "Developer ID"; then
    echo "Developer ID certificate found. Signing app..."
    
    # Get the first Developer ID
    IDENTITY=$(security find-identity -p codesigning -v | grep "Developer ID Application" | head -1 | awk '{print $2}')
    
    codesign --deep --force --verify --verbose \
        --sign "$IDENTITY" \
        --options runtime \
        dist/BaseBuddy.app
    
    echo "✓ App signed successfully"
else
    echo "No Developer ID certificate found. Skipping code signing."
    echo "Note: Users may see security warnings when opening the app."
fi

echo ""
echo "========================================"
echo "Build Complete!"
echo "========================================"
echo ""
echo "Output files:"
echo "  - dist/BaseBuddy.app (Application bundle)"
echo "  - dist/BaseBuddy.dmg (Installer image)"
echo ""
echo "To distribute:"
echo "1. Share the DMG file with users"
echo "2. Users drag BaseBuddy.app to Applications"
echo "3. Double-click to run - no installation needed!"
echo ""
echo "Note: First launch may trigger macOS security dialog."
echo "Users should right-click and select 'Open' to bypass."
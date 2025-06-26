#!/bin/bash
# Install Miniforge on Apple Silicon Mac

echo "Miniforge Installation for Apple Silicon"
echo "========================================"
echo ""

# Detect architecture
ARCH=$(uname -m)
if [[ "$ARCH" == "arm64" ]]; then
    echo "Detected Apple Silicon (ARM64)"
    INSTALLER="Miniforge3-MacOSX-arm64.sh"
    URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
elif [[ "$ARCH" == "x86_64" ]]; then
    echo "Detected Intel Mac (x86_64)"
    INSTALLER="Miniforge3-MacOSX-x86_64.sh"
    URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh"
else
    echo "Unknown architecture: $ARCH"
    exit 1
fi

# Download installer
echo ""
echo "Downloading Miniforge installer..."
curl -L -o "$INSTALLER" "$URL"

if [ ! -f "$INSTALLER" ]; then
    echo "Error: Failed to download installer"
    exit 1
fi

# Make executable and run
echo ""
echo "Running Miniforge installer..."
echo "Note: Follow the prompts to complete installation"
echo "  - Press ENTER to continue"
echo "  - Type 'yes' to accept the license"
echo "  - Press ENTER to accept default location or specify custom path"
echo "  - Type 'yes' when asked to initialize Miniforge3"
echo ""
bash "$INSTALLER"

# Clean up
rm -f "$INSTALLER"

echo ""
echo "Installation complete!"
echo ""
echo "IMPORTANT: You need to restart your terminal or run:"
echo "  source ~/.bashrc  (if using bash)"
echo "  source ~/.zshrc   (if using zsh)"
echo ""
echo "Then verify installation with:"
echo "  conda --version"
echo ""
echo "After that, you can run:"
echo "  ./scripts/setup_hybrid_env.sh"
# BaseBuddy Setup from Scratch on Apple Silicon

This guide will walk you through the complete setup process, starting with no prerequisites installed.

## Step 1: Install Miniforge

Miniforge is a minimal conda installer that uses conda-forge as the default channel.

```bash
# Run the installer script
./scripts/install_miniforge.sh
```

### Manual Installation (if script fails)

```bash
# Download Miniforge for Apple Silicon
curl -L -o Miniforge3-MacOSX-arm64.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh

# Run the installer
bash Miniforge3-MacOSX-arm64.sh

# Follow the prompts:
# - Press ENTER to continue
# - Type 'yes' to accept license
# - Press ENTER for default install location (or specify custom)
# - Type 'yes' to initialize Miniforge3

# Clean up
rm Miniforge3-MacOSX-arm64.sh
```

## Step 2: Restart Your Terminal

**This is crucial!** The conda command won't be available until you start a new terminal session.

```bash
# Close your terminal and open a new one, or run:
source ~/.zshrc  # for zsh (default on macOS)
# or
source ~/.bashrc  # for bash
```

## Step 3: Verify Conda Installation

```bash
# This should now work
conda --version
# Expected output: conda 24.x.x or similar

# Check that conda-forge is the default channel
conda config --show channels
```

## Step 4: Install Homebrew (if not installed)

```bash
# Check if Homebrew is installed
brew --version

# If not installed, install it:
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Follow the instructions at the end to add Homebrew to your PATH
```

## Step 5: Run the BaseBuddy Setup

Now you're ready to set up BaseBuddy:

```bash
# Navigate to BaseBuddy directory
cd /Users/lauferva/Desktop/Professional/Projects/2025/GOAL/BaseBuddy

# Run the hybrid setup script
./scripts/setup_hybrid_env.sh
```

This will:
1. Configure conda channels
2. Install mamba for faster package resolution
3. Create two environments (ARM64 and x86_64)
4. Install necessary Homebrew packages
5. Set up SimuG
6. Create helper scripts

## Step 6: Using BaseBuddy

After setup completes, you'll have two helper scripts in your BaseBuddy directory:

### For native ARM64 environment (most tools):
```bash
source activate_basebuddy.sh
# or
source activate_basebuddy.sh arm
```

### For x86_64 environment (BAMSurgeon, etc.):
```bash
source activate_basebuddy.sh x86
```

## Step 7: Test the Installation

### Test ARM64 environment:
```bash
source activate_basebuddy.sh
python --version  # Should show Python 3.10+
bwa  # Should show BWA help
samtools --version  # Should show samtools version
```

### Test x86_64 environment:
```bash
source activate_basebuddy.sh x86
which addsnv.py  # Should find BAMSurgeon's addsnv.py
```

### Test the GUI:
```bash
source activate_basebuddy.sh
python run_gui.py
```

## Troubleshooting

### "conda: command not found" after installing Miniforge
- Make sure you restarted your terminal or sourced your shell config
- Check if conda was installed to a custom location
- Try: `~/miniforge3/bin/conda --version`

### "source: no such file or directory: activate_basebuddy.sh"
- Make sure you're in the BaseBuddy directory when running the command
- The script is created by setup_hybrid_env.sh, so run that first
- Check if the file exists: `ls -la activate_basebuddy.sh`

### Homebrew installation issues
- On Apple Silicon, Homebrew installs to `/opt/homebrew`
- Make sure `/opt/homebrew/bin` is in your PATH
- Add to ~/.zshrc: `eval "$(/opt/homebrew/bin/brew shellenv)"`

### Permission denied errors
- Make sure scripts are executable: `chmod +x scripts/*.sh`
- Some Homebrew installations may require sudo

## Quick Reference

Once everything is set up:

```bash
# Daily workflow
cd /path/to/BaseBuddy
source activate_basebuddy.sh      # Use native ARM64 env
python run_gui.py                  # Run the GUI

# When you need BAMSurgeon
source activate_basebuddy.sh x86   # Switch to x86 env
python run_gui.py                  # Variant spiking will work

# Use specific tools
./basebuddy-tool bwa mem ...       # Auto-selects correct env
./basebuddy-tool addsnv.py ...     # Auto-switches to x86
```
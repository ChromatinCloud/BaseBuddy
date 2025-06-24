# Scripts Directory

This directory contains various utility scripts for BaseBuddy.

## Directory Structure

- `setup/` - Installation and environment setup scripts
  - `install_aligners.sh` - Install BWA and other aligners
  - `setup_gui_env.sh` - Set up conda environment for GUI
  - `create_fresh_gui_env.sh` - Create a fresh GUI environment
  - `test_environment.sh` - Test that the environment works correctly
  - `test_aligners.sh` - Test aligner installations

- `fixes/` - Scripts to fix common issues
  - `ART_INSTALL_FIX.sh` - Fix ART installation issues
  - `fix_aligner_install.sh` - Fix aligner installation problems
  - `fix_gsl_version.sh` - Fix GSL library version issues

- `run_gui.py` - Main GUI launcher script
- `launch_gui_macos.py` - macOS-specific GUI launcher with environment activation

## Usage

Most scripts can be run directly from the BaseBuddy root directory:
```bash
# From BaseBuddy root
./scripts/setup/test_environment.sh
python scripts/run_gui.py
```
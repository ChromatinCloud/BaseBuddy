# Fix ART GSL Library Error

The error shows ART can't find the GSL library (libgsl.25.dylib). Here's how to fix it:

## Solution 1: Install GSL in your conda environment

```bash
# Activate your environment
conda activate bb_gui

# Install GSL
conda install -c conda-forge gsl -y

# Verify the library is installed
ls $CONDA_PREFIX/lib/libgsl*

# Restart the GUI
python run_gui.py
```

## Solution 2: Reinstall ART with dependencies

```bash
# Activate environment
conda activate bb_gui

# Remove and reinstall ART with all dependencies
conda remove art -y
conda install -c bioconda -c conda-forge art gsl -y

# Test ART
art_illumina --help

# Restart GUI
python run_gui.py
```

## Solution 3: Use system-wide installation

If conda continues to have issues, use Homebrew:

```bash
# Exit conda environment
conda deactivate

# Install with Homebrew
brew install gsl
brew tap brewsci/bio
brew install art

# Use system Python for GUI
/usr/bin/python3 run_gui.py
```

## Quick Test

After fixing, test ART directly:

```bash
# Should show help without errors
art_illumina -h
```

## Alternative: Use Docker

If library issues persist, the Docker image has all dependencies:

```bash
# In the BaseBuddy directory
docker build -t basebuddy .
docker run -it -v $(pwd):/work basebuddy bash
# Then run commands inside container
```

The GSL library issue is common on macOS with conda. Solution 1 or 2 should resolve it quickly.
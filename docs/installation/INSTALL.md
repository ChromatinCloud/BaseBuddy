# BaseBuddy Installation Guide

## Quick Start (GUI)

For GUI usage with core features:

```bash
# Clone the repository
git clone https://github.com/yourusername/BaseBuddy.git
cd BaseBuddy

# Create environment from tested configuration
conda env create -f environment-gui.yml

# Activate environment
conda activate basebuddy-gui

# Run GUI
python run_gui.py
```

## Full Installation (All Features)

For complete functionality including all simulation tools:

```bash
# Create environment
conda env create -f environment.yml

# Activate environment
conda activate basebuddy

# Run CLI
basebuddy --help
```

## Features by Installation Type

### GUI Environment (environment-gui.yml)
✅ GUI interface  
✅ Short read simulation (ART)  
✅ Auto-align with BWA  
✅ Variant spiking  
✅ Quality control  
❌ Long read simulation  
❌ minimap2 aligner  

### Full Environment (environment.yml)
✅ All GUI features  
✅ Long read simulation (NanoSim)  
✅ Mutational signatures  
✅ BAMSurgeon  
✅ FastQC  
❌ minimap2 (conflicts with Python 3.12)  

## Known Issues

1. **minimap2 conflicts**: Due to conda dependency conflicts with Python 3.12, minimap2 is not included. Use BWA for alignment.

2. **macOS tkinter**: If you encounter tkinter issues on macOS, ensure you're using the conda-provided Python, not system Python.

## Manual Installation

If you prefer to install components separately:

```bash
# Create base environment
conda create -n basebuddy python=3.12 tk pip -y
conda activate basebuddy

# Install Python packages
pip install customtkinter typer pysam

# Install bioinformatics tools
conda install -c conda-forge -c bioconda samtools art gsl bwa -y

# Install BaseBuddy
pip install -e .
```

## Docker Installation

For a consistent environment across platforms:

```bash
# Build Docker image
docker build -t basebuddy .

# Run with Docker
docker run -v $(pwd):/data basebuddy --help
```

## Testing Your Installation

```bash
# Test CLI
basebuddy --version

# Test short read simulation
basebuddy short --reference test.fa --depth 30

# Test GUI (requires display)
python run_gui.py
```

## Troubleshooting

### "Module not found" errors
- Ensure you've activated the conda environment
- Check with: `conda info --envs`

### GUI won't start
- Verify tk is installed: `python -c "import tkinter"`
- Check customtkinter: `python -c "import customtkinter"`

### Alignment fails
- Ensure BWA is installed: `which bwa`
- Check samtools: `which samtools`
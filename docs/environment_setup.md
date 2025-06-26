# BaseBuddy Environment Setup Guide

This guide provides comprehensive instructions for setting up a reproducible Conda environment for BaseBuddy on macOS (both Intel and Apple Silicon).

## Quick Start

```bash
# Run the automated setup script
./scripts/setup_basebuddy_env.sh
```

## Manual Setup Steps

### 1. Install Miniforge (if not already installed)

For Apple Silicon (M1/M2/M3):
```bash
curl -L -o Miniforge3-MacOSX-arm64.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh
```

For Intel Macs:
```bash
curl -L -o Miniforge3-MacOSX-x86_64.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh
bash Miniforge3-MacOSX-x86_64.sh
```

### 2. Configure Conda Channels

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

### 3. Install Mamba (Optional but Recommended)

```bash
conda install mamba -n base -c conda-forge
```

### 4. Create the Environment

Using mamba (faster):
```bash
mamba env create -f environment-complete.yml
```

Or using conda:
```bash
conda env create -f environment-complete.yml
```

### 5. Activate the Environment

```bash
conda activate basebuddy
```

## Environment Contents

The environment includes:

### Core Tools
- **Python 3.10-3.11**: Compatible with all required packages
- **pysam**: For BAM/SAM file manipulation
- **numpy, scipy, pandas**: Scientific computing
- **customtkinter**: Modern GUI framework

### Alignment Tools
- **BWA**: Burrows-Wheeler Aligner
- **samtools**: SAM/BAM utilities
- **minimap2**: Versatile sequence aligner

### Read Simulation
- **ART**: Illumina read simulator
- **NanoSim**: Nanopore read simulator
- **wgsim**: Simple read simulator

### Variant Manipulation
- **BAMSurgeon**: For spiking variants into BAM files
- **exonerate**: Sequence alignment (BAMSurgeon dependency)
- **velvet**: De novo assembler (BAMSurgeon dependency)
- **Picard**: Java-based BAM utilities

### Additional Tools
- **bcftools**: VCF file manipulation
- **bedtools**: Genomic interval operations
- **SimuG**: Germline variant simulation (requires manual installation)

## Platform-Specific Notes

### Apple Silicon (ARM64)

Many tools will run natively on ARM64, providing better performance. Some tools (marked with ⚠️) only have x86_64 builds and will run under Rosetta 2:

- ⚠️ **exonerate**: No ARM64 build available
- ⚠️ **ART**: May require Rosetta 2
- ⚠️ **velvet**: May require Rosetta 2

To check if a tool is running natively:
```bash
file $(which <tool_name>)
```

### Intel Macs

All tools run natively without emulation.

## Post-Installation Setup

### SimuG Installation

SimuG needs to be cloned separately:
```bash
mkdir -p external
cd external
git clone https://github.com/yjx1217/simuG.git
cd ..
# Add to PATH or use full path: external/simuG/simuG.pl
```

### Jupyter Integration

To use BaseBuddy in Jupyter notebooks:
```bash
python -m ipykernel install --user --name basebuddy --display-name "BaseBuddy"
```

## Troubleshooting

### Common Issues

1. **"Package not found" errors**
   - Ensure channels are configured correctly
   - Try using mamba instead of conda

2. **"Solving environment" hangs**
   - Use mamba for faster resolution
   - Consider creating environment with fewer packages first

3. **Tools not found after installation**
   - Ensure environment is activated: `conda activate basebuddy`
   - Check PATH: `echo $PATH`

4. **x86_64 tools on Apple Silicon**
   - These will run automatically under Rosetta 2
   - Performance may be slightly reduced but functionality is preserved

### Verifying Installation

Run the verification script:
```bash
conda activate basebuddy
python -c "import pysam; print(f'pysam version: {pysam.__version__}')"
bwa || echo "BWA installed"
samtools --version
```

## Maintenance

### Updating the Environment

To update all packages:
```bash
mamba env update -f environment-complete.yml
```

### Exporting for Reproducibility

To create a locked environment file:
```bash
conda env export --no-builds > environment-locked.yml
```

### Removing the Environment

If you need to start fresh:
```bash
conda env remove -n basebuddy
```

## Docker Alternative

For maximum portability, use the Docker setup:
```bash
docker build -t basebuddy .
docker run -it basebuddy
```

## Additional Resources

- [Miniforge Documentation](https://github.com/conda-forge/miniforge)
- [Bioconda Documentation](https://bioconda.github.io/)
- [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/)
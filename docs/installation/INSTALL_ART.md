# Installing ART for BaseBuddy

ART (Artificial Read Transcript) simulator is required for short read simulation.

## Option 1: Install via Conda (Recommended)

```bash
# In your base environment
conda install -c bioconda art

# Or in the bb_gui environment
conda activate bb_gui
conda install -c bioconda art
```

## Option 2: Install via Homebrew (macOS)

```bash
brew tap brewsci/bio
brew install art
```

## Option 3: Download Pre-compiled Binary

```bash
# Download for macOS
cd ~/Downloads
wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05macos.tgz
tar -xzf artbinmountrainier2016.06.05macos.tgz

# Add to PATH (add to ~/.zshrc for permanent)
export PATH="$HOME/Downloads/art_bin_MountRainier:$PATH"
```

## Option 4: Use Docker (All tools included)

```bash
# Build the BaseBuddy Docker image
docker build -t basebuddy .

# Run with Docker
docker run -v $(pwd):/data basebuddy short --reference /data/refs/GRCh37/hs37d5.fa --depth 100
```

## Quick Test After Installation

```bash
# Check if ART is installed
which art_illumina

# Test with CLI
basebuddy short --reference refs/test/chr12_kras_region.fa --depth 30
```

## For GUI Usage

After installing ART, restart the GUI:
```bash
conda activate bb_gui  # if using conda environment
python run_gui.py
```

The error should be resolved and you can simulate reads around BRAF.
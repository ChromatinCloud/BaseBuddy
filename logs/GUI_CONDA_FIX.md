# Fixing GUI Launch with Conda Python 3.13

## The Issue

Python 3.13 in conda currently has compatibility issues with tkinter on macOS. This prevents the GUI from launching.

## Solutions

### Option 1: Create a new conda environment with Python 3.12

```bash
# Create new environment with Python 3.12
conda create -n basebuddy python=3.12 tk -y

# Activate the environment
conda activate basebuddy

# Install dependencies
pip install customtkinter typer pysam

# Launch the GUI
python run_gui.py
```

### Option 2: Use the CLI instead

The CLI works perfectly with your current setup:

```bash
# View all available commands
basebuddy --help

# Example: Simulate short reads
basebuddy short --reference genome.fa --depth 30 --profile HS25

# Example: Apply signature
basebuddy apply-signature genome.fa --bundled-type sbs --signature-name SBS1 --num-mutations 5000
```

### Option 3: Install Python with tkinter support

Install Python 3.12 from python.org or Homebrew:

```bash
# Via Homebrew
brew install python@3.12
brew install python-tk@3.12

# Then install dependencies
/opt/homebrew/opt/python@3.12/bin/python3.12 -m pip install customtkinter typer pysam

# Launch GUI
/opt/homebrew/opt/python@3.12/bin/python3.12 run_gui.py
```

## Quick Fix for Now

Since you already have all the CLI functionality working, I recommend using Option 1 (new conda environment) only if you specifically need the GUI. The CLI provides all the same functionality and is already fully operational.

To quickly test if the GUI would work in a new environment:

```bash
conda create -n bb_gui python=3.12 tk customtkinter typer pysam -c conda-forge -y
conda activate bb_gui
python run_gui.py
```
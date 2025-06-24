# Running BaseBuddy GUI

## Prerequisites

The GUI requires the `customtkinter` library for the modern interface:

```bash
pip install customtkinter
```

## Launch Methods

### Method 1: Using the launch script (Recommended)
```bash
python run_gui.py
```

### Method 2: Direct module execution
```bash
PYTHONPATH=src python -m basebuddy.gui.main_app
```

### Method 3: From Python
```python
import sys
sys.path.insert(0, 'src')
from basebuddy.gui.main_app import run_gui
run_gui()
```

## GUI Features

The BaseBuddy GUI provides tabs for:

1. **Short Read Sim** - Simulate Illumina reads with ART
   - Reference file selection
   - Depth and read length configuration
   - Profile selection (HS25, HSXt, etc.)
   - Paired/single-end options

2. **Variant Spiking** - Add mutations to BAM files
   - Multiple BAM input support
   - SNP and Indel VCF input
   - VAF configuration
   - IGV session generation

3. **Long Read Sim** - Simulate Nanopore/PacBio reads
   - Model selection
   - Depth or exact read count
   - Reference selection

4. **Apply Signature** - Apply mutational signatures to FASTA
   - Bundled signature selection (SBS1, SBS5)
   - Custom signature file support
   - Number of mutations control

5. **Germline Simulation** - Germline variant workflow
   - Combined germline + read simulation
   - Short or long read options

## Troubleshooting

### Import Errors
If you see import errors, ensure you're running from the BaseBuddy root directory:
```bash
cd /path/to/BaseBuddy
python run_gui.py
```

### Missing Dependencies
The GUI may require additional dependencies for full functionality:
```bash
pip install customtkinter typer pysam
```

### Display Issues
On some systems, you may need to set the display:
```bash
export DISPLAY=:0  # Linux/WSL
python run_gui.py
```

### macOS Security
On macOS, you may need to allow the terminal to control system events for file dialogs to work properly.

## GUI Workflow Example

1. Launch the GUI
2. Select "Short Read Sim" tab
3. Browse and select a reference FASTA
4. Set depth to 50x
5. Choose HS25 profile
6. Click "Run Short Read Simulation"
7. Monitor progress in the status area
8. Click "Open Output Directory" when complete

## Notes

- The GUI runs tasks in separate threads to keep the interface responsive
- Output is displayed in the status text area
- All runs create timestamped directories with manifest files
- The GUI uses the same backend runners as the CLI

## Known Limitations

- Some features may require external tools (ART, BAMSurgeon, etc.)
- Custom signature paths in Apply Signature tab need full paths
- The GUI doesn't currently validate all tool dependencies upfront
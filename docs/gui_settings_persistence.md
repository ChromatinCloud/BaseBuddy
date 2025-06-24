# GUI Settings Persistence

**Date**: 2025-06-24  
**Feature**: Per-tab settings persistence for BaseBuddy GUI

## Overview

The BaseBuddy GUI now automatically saves and restores your settings for each tab. This means you won't have to re-type file paths, parameters, and options every time you launch the GUI, even after updating to a new version.

## How It Works

### Settings Storage

- Settings are saved in individual JSON files in the `./refs/gui_settings/` directory
- Each tab has its own settings file:
  - `short_read_sim.settings.json` - Short Read Simulation tab
  - `variant_spiking.settings.json` - Variant Spiking tab
  - `long_read_sim.settings.json` - Long Read Simulation tab
  - `apply_signature.settings.json` - Apply Signature tab
  - `germline_simulation.settings.json` - Germline Simulation tab

### Automatic Save

Settings are automatically saved when:
- You close the GUI window
- The window close event is triggered

### Automatic Load

Settings are automatically loaded when:
- The GUI starts up
- All tabs are created and initialized

## What Gets Saved

### Short Read Simulation Tab
- Reference FASTA path
- Output root directory
- Run name
- ART platform and profile
- Genome build selection
- Read length, depth, fragment sizes
- Paired-end option
- Overwrite and auto-index options
- Auto-align settings
- Genomic ranges
- Aligner selection
- Output BAM path

### Variant Spiking Tab
- Reference FASTA path
- Output root directory
- Run name
- VCF file path
- Spike percentage
- Overwrite and auto-index options

### Long Read Simulation Tab
- Reference FASTA path
- Output root directory
- Run name
- NanoSim model
- Depth and number of reads
- Overwrite and auto-index options

### Apply Signature Tab
- Reference FASTA path
- Output FASTA path
- Signature ID selection
- Custom VCF path
- Number of mutations
- Random seed
- Genomic ranges

### Germline Simulation Tab
- Reference FASTA path
- Germline VCF path
- Output root directory
- Run name
- Read simulation type (short/long)
- Read-specific parameters
- Overwrite, auto-index, and keep intermediates options

## Benefits

1. **Time Saving**: No need to re-enter frequently used paths and parameters
2. **Version Persistence**: Settings survive GUI updates and new releases
3. **Project Continuity**: Easy to resume work on the same project
4. **Multiple Projects**: Can manually save/swap settings files for different projects

## Technical Details

- Settings are stored in JSON format for easy manual editing if needed
- Each tab's settings are independent - updating one tab doesn't affect others
- Missing or corrupted settings files are handled gracefully
- New fields added in updates will use default values if not in saved settings

## Manual Management

If you want to:
- **Reset settings**: Delete the specific `.settings.json` file in `refs/gui_settings/`
- **Share settings**: Copy the settings files to share project configurations
- **Switch projects**: Rename settings files to save multiple project configurations

## Example Settings File

Here's what a typical `short_read_sim.settings.json` might look like:

```json
{
  "reference_fasta": "/path/to/refs/GRCh38/GRCh38.fa",
  "output_root": "/path/to/output",
  "run_name": "BRAF_simulation",
  "art_platform": "illumina",
  "art_profile": "HS25",
  "genome_build": "hg38/GRCh38",
  "read_length": "150",
  "depth": "30",
  "mean_fragment": "400",
  "std_dev_fragment": "50",
  "is_paired_end": true,
  "overwrite": false,
  "auto_index": true,
  "auto_align": true,
  "genomic_ranges": "chr7:140753336-140753436",
  "aligner": "bwa",
  "output_bam": ""
}
```
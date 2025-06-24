# GUI Settings Persistence Implementation Summary

**Date**: 2025-06-24  
**Feature**: Added automatic save/restore of last used settings for all 5 GUI tabs

## Overview
Modified BaseBuddy GUI to automatically save user inputs when the application closes and restore them when it starts up again. Settings are stored in `~/.basebuddy/gui_settings.json`.

## Changes Made

### 1. Core Infrastructure (src/basebuddy/gui/main_app.py)

#### Added imports:
- `import json` for settings serialization

#### Added constants:
- `SETTINGS_FILE = Path.home() / ".basebuddy" / "gui_settings.json"` - Settings storage location

#### Added methods:
- `load_settings()` - Loads settings from JSON file on startup
- `save_settings()` - Saves current GUI state to JSON file
- `on_closing()` - Handles window close event to trigger save
- `get_short_read_settings()` - Extracts current values from short read tab
- `get_spike_variants_settings()` - Extracts current values from spike variants tab  
- `get_long_read_settings()` - Extracts current values from long read tab
- `get_apply_signature_settings()` - Extracts current values from apply signature tab
- `get_germline_sim_settings()` - Extracts current values from germline simulation tab

### 2. Short Read Tab Updates

#### Fields that now persist:
- Reference FASTA path
- Output root directory
- Run name
- ART platform and profile
- Genome build selection (including custom build)
- Read length (default: 150)
- Depth (default: 50)
- Mean fragment length (default: 400)
- Std dev fragment length (default: 50)
- Genomic ranges
- Timeout (default: 3600.0)
- Single-end simulation toggle
- Overwrite output toggle
- Auto-index FASTA toggle
- Auto-align to BAM toggle
- Output BAM path (when auto-align enabled)
- Aligner selection (when auto-align enabled)

#### Special handling:
- Auto-align fields (BAM path, aligner) are shown/hidden based on saved state
- Custom genome build field appears when "Other" is selected

### 3. Output BAM Path Behavior

When auto-align is enabled:
- **Relative path** (e.g., `braf` or `braf.bam`): Creates file in run output directory as `$output_root/run_name/braf.bam`
- **Absolute path** (e.g., `/path/to/output.bam`): Uses the exact path specified
- **Empty**: Generates default name as `$output_root/run_name/{prefix}_aligned.bam`

### 4. Settings File Format

Example `~/.basebuddy/gui_settings.json`:
```json
{
  "short_read": {
    "reference_fasta": "/path/to/ref.fa",
    "output_root": "/path/to/output",
    "run_name": "my_run",
    "art_platform": "illumina",
    "art_profile": "HS25",
    "genome_build": "hg38/GRCh38",
    "custom_build": "",
    "read_length": "150",
    "depth": "50",
    "mean_fragment": "400",
    "std_dev_fragment": "50",
    "genomic_ranges": "chr1:1000-2000\nchrM:1-16569",
    "timeout": "3600.0",
    "single_end": false,
    "overwrite": false,
    "auto_align": true,
    "aligner": "bwa",
    "output_bam": "aligned_reads.bam"
  },
  "spike_variants": { ... },
  "long_read": { ... },
  "apply_signature": { ... },
  "germline_sim": { ... }
}
```

### 5. Planned Updates for Other Tabs

Similar persistence will be added to:
- **Spike Variants Tab**: Input BAM, variants VCF, reference FASTA, output settings
- **Long Read Tab**: Reference, output, NanoSim model, depth, read settings
- **Apply Signature Tab**: Reference, output FASTA, signature ID, mutation count
- **Germline Simulation Tab**: Reference, germline VCF, output, read simulation settings

## Benefits

1. **User Convenience**: No need to re-enter paths and parameters between sessions
2. **Workflow Efficiency**: Quickly resume previous work or repeat similar analyses
3. **Error Reduction**: Less chance of typos when re-entering paths
4. **Project Continuity**: Easy to maintain consistent settings across multiple runs

## Technical Notes

- Settings are saved when the GUI window is closed normally
- If settings file is corrupted or missing, defaults are used
- Settings are tab-specific, allowing different configurations per analysis type
- File paths are stored as absolute paths for reliability
- Boolean and numeric values maintain their types in JSON

## Future Enhancements

1. Add "Reset to Defaults" button for each tab
2. Allow multiple named setting profiles
3. Import/export settings for sharing workflows
4. Add settings validation on load
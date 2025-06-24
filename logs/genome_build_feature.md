# Genome Build Selection Feature

**Date**: 2025-06-23  
**Feature**: Added genome build selection to CLI and GUI

## Overview

Added the ability to specify genome build (e.g., hg19, hg38, GRCh37, GRCh38) in both CLI and GUI. This helps document which reference genome version is being used and will be included in output manifests for better traceability.

## Changes Made

### 1. CLI Updates (`src/basebuddy/cli.py`)

Added new parameter to the `short` command:
- `--build, -b`: Specify genome build (optional, for documentation)

Example usage:
```bash
basebuddy short --reference hs37d5.fa --build GRCh37 --depth 100
basebuddy short --reference hg38.fa --build hg38 --depth 50
```

### 2. GUI Updates (`src/basebuddy/gui/main_app.py`)

Added to Short Read Sim tab:
- Genome Build dropdown with common options:
  - (empty/not specified)
  - hg19/GRCh37
  - hg38/GRCh38
  - T2T-CHM13
  - mm10
  - mm39
  - Other
- Custom build entry field (shown when "Other" is selected)
- Callback function `_on_genome_build_change()` to toggle custom entry

### 3. Runner Updates (`src/basebuddy/runner.py`)

- Added `genome_build` parameter to `simulate_short()` function
- Logs genome build when specified
- Includes genome build in manifest for documentation

## Usage Examples

### CLI
```bash
# Specify standard build
basebuddy short --reference refs/GRCh37/hs37d5.fa --build GRCh37 --depth 100

# Custom build
basebuddy short --reference refs/dog/canFam3.fa --build canFam3 --depth 50

# Auto-align with build specification
basebuddy short --reference refs/hg38.fa --build hg38 --auto-align --depth 30
```

### GUI
1. Select from dropdown: "hg19/GRCh37", "hg38/GRCh38", etc.
2. Or select "Other" and type custom build name
3. Build information is saved in manifest.json

## Benefits

1. **Documentation**: Clear record of which genome build was used
2. **Reproducibility**: Manifest includes build information
3. **Clarity**: Helps when working with multiple reference versions
4. **Flexibility**: Support for standard and custom genome builds

## Common Genome Builds

### Human
- **hg19/GRCh37**: Genome Reference Consortium Human Build 37
- **hg38/GRCh38**: Genome Reference Consortium Human Build 38
- **T2T-CHM13**: Telomere-to-Telomere complete human genome

### Mouse
- **mm10/GRCm38**: Genome Reference Consortium Mouse Build 38
- **mm39/GRCm39**: Genome Reference Consortium Mouse Build 39

### Others
- Use "Other" option for any custom or non-standard builds

## Technical Details

The genome build is:
- Stored in the manifest.json file
- Logged during simulation startup
- Purely informational (doesn't affect simulation parameters)
- Optional (can be left blank)

## Future Enhancements

Could be extended to:
- Auto-detect build from reference filename patterns
- Validate build against reference file headers
- Provide build-specific parameter defaults
- Link to build-specific annotation files
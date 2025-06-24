# Genome Build to Reference Linking Feature

**Date**: 2025-06-24  
**Feature**: Automatically link genome build selection to appropriate reference FASTA files

## Problem Solved

Users were experiencing issues when:
- Using GRCh38 coordinates with GRCh37 reference (or vice versa)
- BWA trying to index the wrong reference genome
- Confusion about which reference file matches which genome build

## Solution Implemented

### 1. Default Reference Paths

Added mapping of genome builds to expected reference file locations:

```python
DEFAULT_REFERENCES = {
    "hg19/GRCh37": [
        "refs/GRCh37/hs37d5.fa",
        "refs/GRCh37/human_g1k_v37.fasta",
        "refs/GRCh37/GRCh37.fa"
    ],
    "hg38/GRCh38": [
        "refs/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "refs/GRCh38/GRCh38.fa",
        "refs/GRCh38/hg38.fa"
    ],
    # ... other builds
}
```

### 2. Build-Specific Coordinate Examples

Added appropriate genomic range examples for each build:

```python
EXAMPLE_RANGES = {
    "hg19/GRCh37": "chr7:140453136-140453236\nchrM:1-16569",  # BRAF V600E in GRCh37
    "hg38/GRCh38": "chr7:140753336-140753436\nchrM:1-16569",  # BRAF V600E in GRCh38
    # ... other builds
}
```

### 3. Automatic Reference Suggestion

When user selects a genome build:
1. Searches for matching reference files in expected locations
2. If found and current reference is empty or mismatched, updates the reference path
3. Shows appropriate coordinate examples in the genomic ranges textbox
4. Warns if no matching reference is found

### 4. Mismatch Detection

The system detects when reference doesn't match build:
- GRCh37 selected but reference contains "grch38" or "hg38"
- GRCh38 selected but reference contains "grch37", "hg19", or "hs37d5"
- Mouse genome mismatches (mm10 vs mm39)

## User Experience Improvements

1. **Prevents coordinate/reference mismatches** - Major source of errors eliminated
2. **Auto-fills correct reference** - Saves time and reduces errors
3. **Shows appropriate examples** - Build-specific coordinate examples
4. **Clear warnings** - Notifies when references don't match builds

## Example Workflow

1. User selects "hg38/GRCh38" from genome build dropdown
2. System automatically:
   - Finds and fills in `/path/to/refs/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa`
   - Updates genomic ranges example to `chr7:140753336-140753436` (BRAF V600E in GRCh38)
   - Shows status message confirming reference update
3. User can now run simulation with correct coordinate/reference pairing

## Technical Notes

- Reference paths are resolved relative to BaseBuddy installation directory
- Only updates reference if field is empty or contains mismatched reference
- Preserves user's custom reference paths if they don't appear mismatched
- Works with settings persistence - saves and restores genome build selection
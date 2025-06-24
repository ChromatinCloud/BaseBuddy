# Auto-Align Feature for Short Read Simulation

**Date**: 2025-06-23  
**Feature**: Automatic alignment of simulated reads to BAM format

## Overview

Added automatic alignment capability to the short read simulation module. After ART generates FASTQ files, BaseBuddy can now automatically align them to the reference genome and produce sorted, indexed BAM files.

## Changes Made

### 1. CLI Updates (`src/basebuddy/cli.py`)

Added two new parameters to the `short` command:
- `--auto-align`: Boolean flag to enable automatic alignment
- `--aligner`: Choice of aligner (bwa or minimap2, default: bwa)

Example usage:
```bash
basebuddy short --reference genome.fa --depth 100 --auto-align --aligner bwa
```

### 2. GUI Updates (`src/basebuddy/gui/main_app.py`)

Added to Short Read Sim tab:
- "Auto-align to BAM" checkbox
- Output BAM file path entry (shown when auto-align is checked)
- Aligner selection dropdown (bwa/minimap2)
- Toggle function `_on_auto_align_toggle()` to show/hide options

### 3. Runner Updates (`src/basebuddy/runner.py`)

Added to `simulate_short()` function:
- New parameters: `auto_align`, `aligner`, `output_bam`
- Auto-alignment logic after ART completion:
  - Checks for aligner availability (bwa/minimap2)
  - Creates BWA index if needed
  - Aligns reads to reference
  - Sorts SAM to BAM using samtools
  - Indexes final BAM file
  - Cleans up intermediate SAM files
- Updates manifest with BAM/BAI files

## Technical Details

### Alignment Workflow

1. **BWA alignment**:
   ```
   bwa index reference.fa (if needed)
   bwa mem -t 4 reference.fa reads_1.fq reads_2.fq -o aligned.sam
   samtools sort -o aligned.bam aligned.sam
   samtools index aligned.bam
   ```

2. **Minimap2 alignment**:
   ```
   minimap2 -ax sr -t 4 reference.fa reads_1.fq reads_2.fq -o aligned.sam
   samtools sort -o aligned.bam aligned.sam
   samtools index aligned.bam
   ```

### Error Handling

- Checks if aligners are installed before attempting alignment
- Provides helpful error messages with installation instructions
- Validates aligner choice (only bwa/minimap2 supported)
- Cleans up temporary files even on failure

## Dependencies

New tool requirements when auto-align is enabled:
- BWA or minimap2 (based on user choice)
- samtools (for sorting and indexing)

Installation:
```bash
conda install -c bioconda bwa samtools
# or
conda install -c bioconda minimap2 samtools
```

## Usage Examples

### CLI
```bash
# Basic usage with auto-align
basebuddy short --reference refs/GRCh37/hs37d5.fa --depth 100 --auto-align

# Specify output BAM location
basebuddy short --reference genome.fa --depth 50 --auto-align --output-bam results/aligned.bam

# Use minimap2 instead of BWA
basebuddy short --reference genome.fa --depth 30 --auto-align --aligner minimap2
```

### GUI
1. Check "Auto-align to BAM" checkbox
2. Optionally specify output BAM file path
3. Select aligner (BWA or minimap2)
4. Run simulation

## Benefits

1. **Convenience**: One-step process from reference to aligned BAM
2. **Integration**: BAM files can be directly used in Variant Spiking tab
3. **Visualization**: BAM files can be loaded in IGV
4. **Efficiency**: Avoids manual alignment steps

## Testing

To test with BRAF example:
```bash
basebuddy short \
  --reference refs/GRCh37/hs37d5.fa \
  --depth 100 \
  --readlen 150 \
  --auto-align \
  --outdir test_output/braf_aligned
```

This produces:
- `gui_sim_reads1.fq` and `gui_sim_reads2.fq` (FASTQ files)
- `gui_sim_reads_aligned.bam` (sorted BAM)
- `gui_sim_reads_aligned.bam.bai` (BAM index)
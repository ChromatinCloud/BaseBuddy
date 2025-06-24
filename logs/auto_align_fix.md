# Auto-Align Reference Fix

**Date**: 2025-06-24  
**Issue**: BWA was trying to index the full genome instead of the subset reference, causing timeouts

## Problem

When using genomic ranges (e.g., chr7:140753336-140753436), BaseBuddy creates a small subset FASTA file containing only that region. However, the auto-align feature was incorrectly using the original full genome reference for alignment instead of the subset, causing:

1. BWA to attempt indexing the entire genome (hs37d5.fa - ~3GB)
2. Timeout after 600 seconds (10 minutes)
3. No BAM file generation

## Solution

Modified `src/basebuddy/runner.py` to use `current_reference_for_art` (the actual reference used for simulation) instead of `reference_fasta` (the original full genome) for alignment.

### Changes Made

1. **BWA alignment** - Lines 305-325:
   - Now uses `align_reference = str(current_reference_for_art)`
   - Indexes only the subset reference (much faster)
   - Aligns reads to the correct reference

2. **Minimap2 alignment** - Lines 363-368:
   - Same fix applied for consistency

## Benefits

1. **Speed**: Indexing a small subset (e.g., 100bp) takes <1 second vs 10+ minutes for full genome
2. **Correctness**: Reads are aligned to the same reference they were simulated from
3. **Resource efficiency**: Less memory and disk space required

## Example

For a BRAF region (chr7:140753336-140753436):
- Before: Tried to index 3GB hs37d5.fa → timeout
- After: Indexes ~100bp subset file → completes in seconds

## Testing

To verify the fix:
1. Select a genomic range in the GUI
2. Enable "Auto-align to BAM"
3. Run simulation
4. Should see BAM file created with name like `chr7_140753336_140753436_aligned.bam`
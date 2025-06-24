# IGV Viewing Guide for BaseBuddy BAM Files

**Date**: 2025-06-24  
**Topic**: How to view BaseBuddy-generated BAM files in IGV

## The Challenge

When simulating reads from a specific genomic region (e.g., BRAF), BaseBuddy:
1. Extracts just that region from the reference genome
2. Simulates reads from this small subset
3. Aligns reads back to the subset (fast!)
4. Creates a BAM where the reference name is the region (e.g., `7:140453136-140453236`)

IGV expects standard chromosome names (e.g., `chr7` or `7`), not region names.

## Solutions

### Option 1: Load the Subset Reference into IGV (Recommended)

1. In IGV, go to **Genomes → Load Genome from File**
2. Navigate to your BaseBuddy output directory
3. Select the subset reference file: `temp_subset_ref_[run_name].fa`
4. Load your BAM file - it will now display correctly!

**Advantages:**
- Fast loading (small reference)
- Perfect alignment visualization
- No conversion needed

### Option 2: Use Full Genome Alignment

When running BaseBuddy:
1. Ensure the full reference genome is BWA-indexed (one-time, ~60-90 min)
2. The alignment will be slower but produce IGV-ready BAMs

### Option 3: Convert Coordinates Post-Alignment

For existing BAM files, you can create an IGV session that remaps coordinates:
1. The subset reference and BAM are listed in your output files
2. Create a custom IGV session using the subset reference

## Quick Reference

### Output Files for IGV
After a successful run with genomic ranges, you'll see:
```
Generated Files:
  FASTQ Files:
    - gui_sim_reads1.fq
    - gui_sim_reads2.fq
  Alignment Files:
    - 7_140453136_140453236_aligned.bam (BAM)
    - 7_140453136_140453236_aligned.bam.bai (BAI)
  Other Files:
    - temp_subset_ref_braf.fa (Subset Reference for IGV)
    - temp_subset_ref_braf.fa.fai (FAI)
```

### Loading in IGV
1. **First**: Load `temp_subset_ref_braf.fa` as genome
2. **Then**: Load the BAM file
3. **Navigate**: The entire reference is your region of interest

## Why This Approach?

- **Performance**: Aligning to a 100bp region is instant vs. hours for full genome
- **Accuracy**: Reads align perfectly to their source sequence
- **Flexibility**: Can view in IGV with subset reference or convert later
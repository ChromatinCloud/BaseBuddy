# BAM File Naming Update

**Date**: 2025-06-24  
**Feature**: Enhanced auto-generated BAM filename to use genomic ranges

## Change Summary

Modified the auto-align feature in `src/basebuddy/runner.py` to generate more descriptive BAM filenames when no output name is specified.

### Previous Behavior
- Empty output BAM name → `{prefix}_aligned.bam`

### New Behavior
- Empty output BAM name + genomic ranges → `{first_range}_aligned.bam`
- Empty output BAM name + no ranges → `{prefix}_aligned.bam` (fallback)

### Examples

1. **With genomic range `chr7:10000-20000`**:
   - Output: `chr7_10000_20000_aligned.bam`

2. **With genomic range `chrX:1234567-2345678`**:
   - Output: `chrX_1234567_2345678_aligned.bam`

3. **With multiple ranges** (uses first range):
   - Ranges: `["chr1:1000-2000", "chr2:3000-4000"]`
   - Output: `chr1_1000_2000_aligned.bam`

4. **No ranges provided** (fallback):
   - Output: `gui_sim_reads_aligned.bam` (or similar prefix-based name)

### Implementation Details

The update sanitizes the genomic range string for filesystem compatibility:
- Colons (`:`) → underscores (`_`)
- Hyphens (`-`) → underscores (`_`)
- Commas (`,`) → removed (for ranges with comma separators like `1,000`)

This makes the output files more descriptive and easier to identify when working with specific genomic regions.
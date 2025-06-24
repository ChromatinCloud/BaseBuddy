# IGV Compatibility Fix for BAM Files

**Date**: 2025-06-24  
**Issue**: BAM files created from subset regions were not compatible with IGV

## Problem

When simulating reads from a genomic range (e.g., `7:140753300-140753900`):
1. BaseBuddy extracts just that region from the reference
2. Simulates reads from the subset
3. Aligns reads back to the subset
4. Creates BAM with non-standard reference names like `7:140753300-140753900`
5. IGV expects standard chromosome names like `chr7` or `7`

Error message:
```
File does not contain any sequence names which match the current genome.
Sequence names in 'filename': 7:140753300-140753900
Sequence names in genome reference sequence: chr1, chr2, chr3, ...
```

## Solution

Changed alignment strategy to use the full reference genome instead of the subset:

```python
# Before: Aligned to subset reference
align_reference = str(current_reference_for_art)  # subset with just the range

# After: Align to full reference
align_reference = str(original_reference_path_obj)  # full genome
```

## Benefits

1. **IGV Compatible**: BAM files now have standard chromosome names
2. **Correct Coordinates**: Reads are mapped to their correct genomic positions
3. **Better Visualization**: Can view reads in context of full chromosome
4. **No Post-Processing**: No need to fix headers after alignment

## Performance Considerations

- **BWA Indexing**: Full genome indexing takes 60+ minutes (one-time)
- **Alignment Speed**: Slightly slower but negligible for small read sets
- **Index Reuse**: Once indexed, the BWA index is reused for all future runs

## Workflow

1. Simulate reads from genomic subset (fast, small reference)
2. Align reads to full genome (correct coordinates)
3. Output BAM with standard chromosome names
4. Direct loading into IGV works

## Note on Reference Indexing

The fix includes intelligent timeout adjustment:
- Small references (<1GB): Standard timeout
- Large references (>1GB): Extended timeout (2 hours) for indexing

Users will see: "This may take a while for large genomes..." when indexing starts.
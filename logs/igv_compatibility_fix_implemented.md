# IGV Compatibility Fix Implementation

**Date**: 2025-06-24  
**Issue**: BAM files created from subset regions were not compatible with IGV

## Problems Fixed

### 1. IGV File Path Construction
- **Issue**: IGV loading was using relative paths (just filenames) instead of full paths
- **Fix**: Modified `on_igv_button_short_click()` and `on_igv_button_long_click()` to construct full paths by joining output directory with filename

### 2. BAM Chromosome Naming
- **Issue**: BAM files had subset reference names like `7:140753300-140753900` instead of standard chromosome names
- **Fix**: Changed alignment to use full original reference instead of subset reference

## Code Changes

### GUI Path Construction (main_app.py)
```python
# Before: 
bam_files = [f["path"] for f in output_files if f["type"] == "BAM"]

# After:
bam_files = []
for f in output_files:
    if f["type"] == "BAM":
        file_path = f["path"]
        if not Path(file_path).is_absolute():
            file_path = str(Path(output_dir) / file_path)
        bam_files.append(file_path)
```

### Alignment Reference (runner.py)
```python
# Before:
align_reference = str(current_reference_for_art)  # subset reference

# After:
align_reference = str(original_reference_path_obj)  # full reference
```

### BWA Index Check Fix
```python
# Before:
bwt_file = Path(align_reference).with_suffix(".fa.bwt")  # Wrong!

# After:
bwt_file = Path(str(align_reference) + ".bwt")  # Correct
```

### Large Genome Indexing
Added intelligent timeout adjustment for BWA indexing:
- Detects reference size
- Uses extended timeout (2 hours) for references >1GB
- Provides user feedback about large genome indexing

## Benefits

1. **IGV Compatible**: BAM files now have standard chromosome names
2. **Correct File Paths**: IGV can find and load BAM files
3. **Proper Coordinates**: Reads mapped to correct genomic positions
4. **Better Performance**: Smart timeout handling for large genomes

## Workflow Now

1. Simulate reads from genomic subset (fast)
2. Align reads to full genome (standard chromosome names)
3. IGV loads BAM files directly without issues
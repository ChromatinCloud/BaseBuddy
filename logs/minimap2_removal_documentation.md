# Minimap2 Removal Documentation

**Date**: 2025-06-24  
**Reason**: Minimap2 was causing conda environment conflicts with Python 3.12, removing Python from the environment when installed.

## Code Removed

### 1. From `src/basebuddy/runner.py` (lines 357-394)

```python
elif aligner.lower() == "minimap2":
    # Check if minimap2 is available
    try:
        mm2_cmd = bb_utils.Command("minimap2")
    except bb_utils.BaseBuddyConfigError:
        raise bb_utils.BaseBuddyConfigError(
            f"minimap2 not found in PATH but required for auto-align",
            details="Install minimap2 with: conda install -c bioconda minimap2"
        )
    
    # Run minimap2 alignment
    if is_paired_end:
        r1_path = art_output_prefix_path.with_name(art_output_prefix_path.name + "1.fq")
        r2_path = art_output_prefix_path.with_name(art_output_prefix_path.name + "2.fq")
        mm2 = bb_utils.Command("minimap2")
        mm2.parts.extend(["-ax", "sr", "-t", "4", align_reference, str(r1_path), str(r2_path)])
        mm2_cmd_parts = mm2.get_command_parts()
    else:
        r_path = art_output_prefix_path.with_suffix(".fq")
        mm2 = bb_utils.Command("minimap2")
        mm2.parts.extend(["-ax", "sr", "-t", "4", align_reference, str(r_path)])
        mm2_cmd_parts = mm2.get_command_parts()
    
    # Pipe to samtools sort
    try:
        samtools_cmd = bb_utils.Command("samtools")
    except bb_utils.BaseBuddyConfigError:
        raise bb_utils.BaseBuddyConfigError(
            f"samtools not found in PATH but required for auto-align",
            details="Install samtools with: conda install -c bioconda samtools"
        )
    
    # Create unsorted SAM first, then sort
    sam_path = bam_path.with_suffix(".sam")
    bb_utils.run_external_cmd(mm2_cmd_parts + ["-o", str(sam_path)], timeout_seconds=timeout, stream_output=True, cwd=run_output_dir)
    
    # Sort SAM to BAM
    samtools_sort = bb_utils.Command("samtools")
    samtools_sort.parts.extend(["sort", "-o", str(bam_path), str(sam_path)])
    bb_utils.run_external_cmd(samtools_sort.get_command_parts(), timeout_seconds=timeout, stream_output=True, cwd=run_output_dir)
    
    # Clean up SAM file
    sam_path.unlink(missing_ok=True)
```

### 2. From `src/basebuddy/cli.py`

The `--aligner` option allowed choosing between "bwa" and "minimap2". This will be modified to only accept "bwa".

### 3. From Installation Scripts

- `scripts/setup/install_aligners.sh` - Had minimap2 installation commands
- `scripts/setup/test_aligners.sh` - Had minimap2 testing
- `scripts/fixes/fix_aligner_install.sh` - Had minimap2 fixes
- `scripts/setup/create_fresh_gui_env.sh` - Mentioned minimap2 exclusion

### 4. From Documentation

- `docs/installation/INSTALL.md` - Mentioned minimap2 as an option
- `logs/auto_align_feature.md` - Documented minimap2 support
- `logs/auto_align_fix.md` - Had minimap2 fixes

## Future Considerations

If minimap2 support is needed in the future:
1. Wait for conda package compatibility with Python 3.12+
2. Consider using pip installation instead of conda
3. Use containerized environment (Docker) to avoid conflicts
4. Create separate environment specifically for minimap2

## User Feedback

User explicitly stated: "no. we need the align feature. but we dont have to have minimap2" and "minimap2 was ultimately not used because of issues with env."
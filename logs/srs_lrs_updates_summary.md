# BaseBuddy SRS and LRS Updates Summary

**Date**: 2025-06-23  
**Developer**: Claude (with human guidance)  
**Scope**: Short Read Simulation (SRS) and Long Read Simulation (LRS) improvements

## Overview

This document summarizes all improvements made to BaseBuddy's Short Read Simulation (ART) and Long Read Simulation (NanoSim) functionality, including enhanced error handling, parameter validation, and output detection.

---

## Short Read Simulation (SRS) Updates

### 1. Enhanced Error Messages

**Location**: `src/basebuddy/runner.py` - `simulate_short()` function

**Improvements**:
- Added detailed error messages with typical value ranges for all parameters
- Implemented comprehensive parameter validation with helpful context

**Specific Changes**:
```python
# Before: Basic error
if depth <= 0: raise BaseBuddyInputError(f"Sequencing depth must be a positive integer, got {depth}.")

# After: Detailed error with guidance
if depth <= 0:
    raise BaseBuddyInputError(
        f"Invalid sequencing depth: {depth}",
        details=f"Sequencing depth must be a positive integer (typically 1-100 for WGS, 100-1000 for targeted sequencing). Got {depth}."
    )
```

### 2. Parameter Validation Enhancements

**New Validations Added**:
- Read length validation with typical ranges (50-300bp for Illumina)
- Fragment length validation for paired-end reads
- Mean fragment length must be greater than read length for paired-end
- Platform validation (illumina, 454, solid)
- Profile validation with helpful list of valid profiles

**Example**:
```python
# Fragment length validation for paired-end
if mean_fragment_length <= read_length:
    raise BaseBuddyInputError(
        f"Fragment length too short for paired-end reads",
        details=f"Mean fragment length ({mean_fragment_length}bp) must be greater than read length ({read_length}bp) for paired-end sequencing."
    )
```

### 3. Tool Installation Guidance

**Improvement**: Better error messages when ART simulator is not found

```python
try:
    art_cmd_builder = bb_utils.Command(art_exe_name)
except bb_utils.BaseBuddyConfigError as e:
    raise bb_utils.BaseBuddyConfigError(
        f"ART simulator not found for platform '{art_platform}'",
        details=f"The executable '{art_exe_name}' is required but not found in PATH. "
               f"Please install ART (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) "
               f"or use the Docker image which includes all dependencies."
    )
```

### 4. Reference FASTA Error Handling

**Improvements**:
- Clear error when reference file is not accessible
- Helpful hints about using --reference flag or cached GRCh38
- Better error messages for indexing failures

---

## Long Read Simulation (LRS) Updates

### 1. num_reads Parameter Implementation

**Location**: `src/basebuddy/runner.py` - `simulate_long()` function  
**CLI**: `src/basebuddy/cli.py` - Added `--num-reads` option

**Implementation Details**:
- Added proper handling of `num_reads` parameter from CLI
- When specified, `num_reads` overrides `depth` parameter
- Proper validation to ensure either depth or num_reads is positive
- Clear error messages when both are invalid

**Code Logic**:
```python
if num_reads_param is not None:
    try:
        num_reads_val = int(num_reads_param)
        if num_reads_val > 0:
            nanosim_cmd_builder.add_option("-N", str(num_reads_val))
            # num_reads takes precedence over depth
        elif depth > 0:
            nanosim_cmd_builder.add_option("-c", str(depth))
            logger.warning(f"num_reads parameter '{num_reads_param}' is not positive. Using depth.")
```

### 2. Enhanced Output Detection

**Problem**: NanoSim generates different output file patterns depending on version and settings

**Solution**: Comprehensive output detection checking multiple patterns

**Patterns Checked**:
- `{prefix}_aligned_reads.fastq`
- `{prefix}_aligned_reads.fasta`
- `{prefix}_unaligned_reads.fastq`
- `{prefix}_unaligned_reads.fasta`
- `{prefix}.fastq`
- `{prefix}.fasta`

**Error Handling**: If no output found, lists all files in directory for debugging

### 3. Model Validation

**Added Validation**: List of valid NanoSim models with helpful error messages

```python
valid_models = [
    "nanopore_R9.4.1", "nanopore_R10.3", "nanopore_RNA",
    "pacbio_CLR", "pacbio_CCS", "pacbio_sequel",
    "dna_r9.4.1_e8.1", "dna_r10.3_e8.2", "rna_r9.4.1_e8.1"
]
```

### 4. Parameter Error Messages

**Improved Error**:
```python
if depth <= 0 and not (manifest_params.get("num_reads") and manifest_params["num_reads"] and int(manifest_params["num_reads"]) > 0):
    raise bb_utils.BaseBuddyInputError(
        f"Invalid sequencing parameters",
        details=f"Either depth (got {depth}) or num_reads must be positive. Use --depth for coverage-based simulation or --num-reads for exact read count."
    )
```

---

## Testing

### Unit Tests Created

**File**: `tests/test_new_features.py`

**Test Coverage**:
1. `TestLongReadSimulation`:
   - `test_num_reads_parameter_validation()` - Validates parameter checking
   - `test_num_reads_overrides_depth()` - Ensures num_reads takes precedence
   - `test_output_detection_multiple_patterns()` - Tests flexible output detection

2. `TestImprovedErrorMessages`:
   - `test_short_read_error_messages()` - Fragment length validation
   - `test_invalid_platform_error()` - Platform validation

---

## Code Changes Summary

### Files Modified:

1. **src/basebuddy/runner.py**:
   - Lines 86-91: Enhanced parameter validation for SRS
   - Lines 93-118: Improved tool path and platform validation
   - Lines 103-124: Better reference FASTA handling
   - Lines 162-183: Enhanced ART error handling
   - Lines 625-705: num_reads parameter handling for LRS
   - Lines 711-748: Improved output detection for LRS

2. **src/basebuddy/cli.py**:
   - Line 91: Added `num_reads` parameter to long read command
   - Already had proper typing imports fixed

### Key Improvements:

1. **User Experience**:
   - Clear, actionable error messages
   - Typical value ranges provided in errors
   - Installation hints for missing tools
   - Better guidance for parameter selection

2. **Robustness**:
   - Flexible output detection for tool version differences
   - Comprehensive parameter validation
   - Proper error propagation with context

3. **Functionality**:
   - num_reads parameter fully functional for exact read count simulation
   - Better handling of edge cases
   - Improved logging for debugging

---

## Usage Examples

### Short Read Simulation
```bash
# Will show helpful error about fragment length
basebuddy short --reference genome.fa --depth 30 --read-len 150 --mean-frag-len 100 --paired

# Will show platform validation error
basebuddy short --reference genome.fa --depth 30 --platform nanopore
```

### Long Read Simulation
```bash
# Use exact read count (overrides depth)
basebuddy long --reference genome.fa --num-reads 10000 --model nanopore_R10.3

# Coverage-based simulation
basebuddy long --reference genome.fa --depth 20 --model pacbio_CCS
```

---

## Impact

These updates significantly improve the usability and reliability of BaseBuddy's read simulation features:

1. **Reduced user frustration** through clear, helpful error messages
2. **Better parameter guidance** with typical values and constraints
3. **More flexible tool integration** with various NanoSim versions
4. **Enhanced functionality** with num_reads parameter for precise control
5. **Improved debugging** with comprehensive output detection and logging

The changes maintain backward compatibility while adding new capabilities and improving the overall user experience.
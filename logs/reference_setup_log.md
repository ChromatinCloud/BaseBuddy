# Reference Setup Log
**Date**: 2025-06-23  
**Task**: Configure reference genome handling for BaseBuddy

## Actions Taken

### 1. Updated .gitignore
- Added `/refs/` to ignore reference genome files (line 225)
- Prevents large reference files from being tracked in git

### 2. Analyzed Broad b37 Bundle
- User downloaded Broad Institute b37 reference bundle to `refs/GRCh37/`
- Bundle contains VCF files but main FASTA (hs37d5.fa) still downloading
- Identified essential files needed:
  - `hs37d5.fa` (main reference)
  - `hs37d5.fa.fai` (index)
  - `hs37d5.dict` (sequence dictionary)

### 3. Created Documentation
- Created `refs/REFERENCE_SETUP.md` with:
  - File descriptions and usage in BaseBuddy
  - Download commands for missing files
  - Setup instructions for CLI and GUI
  - Storage requirements (~50GB full bundle, ~3GB minimum)

### 4. Created Test Reference
- Generated `refs/test/chr12_kras_region.fa` (1.5kb test file)
- Indexed with samtools
- Immediately usable for testing BaseBuddy features

## Summary
User's reference download is ongoing. Provided interim test reference and complete documentation for reference setup. BaseBuddy can now run with test data while full GRCh37 reference completes downloading.
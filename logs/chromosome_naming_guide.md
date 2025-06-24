# Chromosome Naming Guide for BaseBuddy

**Date**: 2025-06-24  
**Issue**: Different reference genomes use different chromosome naming conventions

## Quick Reference

### GRCh37 (hg19) - NO "chr" prefix
- Autosomes: `1`, `2`, `3`, ... `22`
- Sex chromosomes: `X`, `Y`
- Mitochondrial: `MT`
- Example: `7:140453136-140453236` (BRAF V600E)

### GRCh38 (hg38) - WITH "chr" prefix
- Autosomes: `chr1`, `chr2`, `chr3`, ... `chr22`
- Sex chromosomes: `chrX`, `chrY`
- Mitochondrial: `chrM`
- Example: `chr7:140753336-140753436` (BRAF V600E)

## Common Files and Their Naming

| Reference File | Build | Chr Naming |
|---------------|-------|------------|
| hs37d5.fa | GRCh37 | No prefix (1, 2, X, MT) |
| human_g1k_v37.fasta | GRCh37 | No prefix |
| GRCh38.fa | GRCh38 | With prefix (chr1, chr2, chrX, chrM) |
| hg38.fa | GRCh38 | With prefix |

## Important Coordinates

### BRAF V600E
- **GRCh37**: `7:140453136`
- **GRCh38**: `chr7:140753336`

### KRAS G12C
- **GRCh37**: `12:25398284`
- **GRCh38**: `chr12:25398284`

## Troubleshooting

### Error: "Reference chr7:... not found in FASTA file"
This means you're using "chr" prefix with a reference that doesn't use it (like GRCh37).
- Solution: Remove "chr" prefix (use `7:...` instead of `chr7:...`)

### Error: "Reference 7:... not found in FASTA file"
This means you're NOT using "chr" prefix with a reference that requires it (like GRCh38).
- Solution: Add "chr" prefix (use `chr7:...` instead of `7:...`)

## How to Check Your Reference

To see what chromosome names your reference uses:
```bash
grep "^>" your_reference.fa | head
```

If you see:
- `>1` → No prefix (GRCh37 style)
- `>chr1` → With prefix (GRCh38 style)
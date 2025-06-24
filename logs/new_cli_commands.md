# New CLI Commands Added to BaseBuddy

**Date**: 2025-06-23  
**Developer**: Claude (with human guidance)

## Overview

Two new CLI commands have been added to BaseBuddy to expose existing functionality that was previously only available through the API:

1. **download-reference** - Download reference genomes and verify checksums
2. **apply-signature** - Apply mutational signatures to FASTA files

---

## 1. download-reference Command

### Purpose
Downloads reference files (genomes, annotations, etc.) from URLs with optional checksum verification and automatic indexing for FASTA files.

### Implementation Details
- **Location**: `src/basebuddy/cli.py` lines 460-532
- **Runner**: Uses existing `download_reference_runner` from `bb_runners.py`
- **Features**:
  - Automatic filename detection from URL
  - Checksum verification (SHA256, MD5, SHA1)
  - Automatic FASTA indexing with samtools
  - Progress streaming during download
  - Timeout configuration

### Usage Examples

```bash
# Download human reference genome
basebuddy download-reference \
  https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  --output-dir references \
  --checksum abc123...

# Download with custom filename
basebuddy download-reference \
  https://example.com/genome.fa \
  --filename my_genome.fa \
  --output-dir my_refs

# Download without checksum verification
basebuddy download-reference \
  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz \
  --output-dir ncbi_refs
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| url | str | Required | URL of file to download |
| --output-dir, -o | Path | references | Output directory |
| --filename, -f | str | Auto from URL | Save as filename |
| --checksum | str | None | Expected checksum |
| --checksum-type | str | sha256 | Algorithm (sha256, md5, sha1) |
| --run-name | str | Auto-generated | Custom run name |
| --timeout | int | 10800 | Timeout in seconds |
| --overwrite | bool | False | Overwrite existing |

---

## 2. apply-signature Command

### Purpose
Applies mutational signatures to reference FASTA files, creating mutated genomes for simulation studies.

### Implementation Details
- **Location**: `src/basebuddy/cli.py` lines 534-654
- **Runner**: Uses existing `apply_signature_to_fasta` from `runner.py`
- **Features**:
  - Support for bundled COSMIC signatures (SBS, DBS, ID)
  - Custom signature file support
  - Reproducible with seed
  - Automatic output indexing
  - Signature validation

### Bundled Signature Files
Located in `src/basebuddy/data/signatures/`:
- `sbs_grch37_cosmic_v3.3.tsv` - Single base substitutions
- `dbs_grch37_cosmic_v3.3.tsv` - Doublet base substitutions  
- `id_grch37_cosmic_v3.3.tsv` - Insertions/deletions

### Usage Examples

```bash
# Apply SBS1 signature (clock-like)
basebuddy apply-signature \
  genome.fa \
  --bundled-type sbs \
  --signature-name SBS1 \
  --num-mutations 5000 \
  --output mutated_sbs1.fa

# Apply UV signature (SBS7a)
basebuddy apply-signature \
  genome.fa \
  --bundled-type sbs \
  --signature-name SBS7a \
  --num-mutations 10000 \
  --seed 42

# Use custom signature file
basebuddy apply-signature \
  genome.fa \
  --signature-file my_signatures.tsv \
  --signature-name CustomSig1 \
  --num-mutations 3000 \
  --output custom_mutated.fa

# Apply doublet substitutions
basebuddy apply-signature \
  genome.fa \
  --bundled-type dbs \
  --signature-name DBS1 \
  --num-mutations 500
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| reference | Path | Required | Input FASTA file |
| --signature-file, -s | Path | None | Custom signature TSV |
| --signature-name, -n | str | Required | Signature to apply |
| --output, -o | Path | mutated.fa | Output FASTA |
| --num-mutations, -m | int | 1000 | Number of mutations |
| --seed | int | 42 | Random seed |
| --bundled-type | str | None | Use bundled (sbs, dbs, id) |

### Signature File Format
TSV files with mutation types as rows and signatures as columns:
```
MutationType	SBS1	SBS2	SBS3
A[C>A]A	0.011	0.00018	0.022
A[C>A]C	0.009	0.00021	0.017
...
```

---

## Integration with Existing Features

### Workflow Example: Signature-based Simulation

```bash
# 1. Download reference genome
basebuddy download-reference \
  https://example.com/chr12.fa.gz \
  --output-dir refs

# 2. Apply smoking signature
basebuddy apply-signature \
  refs/chr12.fa \
  --bundled-type sbs \
  --signature-name SBS4 \
  --num-mutations 10000 \
  --output refs/chr12_smoking.fa

# 3. Simulate reads from mutated genome
basebuddy short \
  --reference refs/chr12_smoking.fa \
  --depth 100 \
  --profile HS25

# 4. Run QC
basebuddy qc *.fastq.gz --output-dir qc_results
```

---

## Error Handling

Both commands include comprehensive error handling:

1. **File validation** - Check existence and permissions
2. **Network errors** - Timeout and connection handling  
3. **Checksum failures** - Clear error messages
4. **Signature validation** - List available signatures on error
5. **Tool dependencies** - Helpful messages if curl/samtools missing

---

## Testing

### Test download-reference:
```bash
# Test with small file
basebuddy download-reference \
  https://raw.githubusercontent.com/samtools/samtools/develop/test/ce.fa \
  --output-dir test_refs
```

### Test apply-signature:
```bash
# Create test FASTA
echo -e ">chr1\nACGTACGTACGTACGTACGTACGTACGT" > test.fa

# Apply signature
basebuddy apply-signature \
  test.fa \
  --bundled-type sbs \
  --signature-name SBS1 \
  --num-mutations 5 \
  --output test_mutated.fa
```

---

## Implementation Notes

### Code Organization
1. Both commands follow existing CLI patterns
2. Use existing runner functions (no duplication)
3. Proper error handling with BaseBuddy exceptions
4. Consistent parameter naming

### Future Enhancements
1. Add progress bar for downloads
2. Support for compressed signature files
3. Batch signature application
4. Signature visualization tools
5. Integration with cloud storage (S3, GCS)

---

## Impact

These additions make BaseBuddy more complete by exposing functionality that was previously hidden:

1. **download-reference** - Simplifies reference genome acquisition
2. **apply-signature** - Enables signature-based mutation studies
3. Together they enable new workflows for cancer genomics research

The commands integrate seamlessly with existing BaseBuddy features, allowing users to create complete simulation pipelines from reference download through read generation and quality control.
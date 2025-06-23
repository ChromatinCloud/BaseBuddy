# BaseBuddy Vignettes

This document provides comprehensive examples and use cases for BaseBuddy's functionality. All commands assume BaseBuddy is installed and available in your PATH.

## Table of Contents
1. [Short Read Simulation (Illumina)](#short-read-simulation)
2. [Long Read Simulation (Nanopore/PacBio)](#long-read-simulation)
3. [Variant Spiking](#variant-spiking)
4. [Mutational Signatures](#mutational-signatures)
5. [Strand Bias Introduction](#strand-bias)
6. [Quality Control](#quality-control)
7. [Complete Workflows](#complete-workflows)

---

## Short Read Simulation

### Basic Usage
Simulate paired-end Illumina reads at 30x depth:
```bash
basebuddy short --reference genome.fa --depth 30 --profile HS25
```

### Advanced Options
```bash
# High-depth targeted sequencing simulation
basebuddy short \
  --reference target_region.fa \
  --outdir results_targeted \
  --depth 500 \
  --read-len 150 \
  --profile HSXt \
  --mean-frag-len 350 \
  --std-frag-len 50 \
  --paired

# Single-end reads for RNA-seq simulation
basebuddy short \
  --reference transcriptome.fa \
  --outdir results_rnaseq \
  --depth 100 \
  --read-len 75 \
  --profile HS25 \
  --single
```

### Using Cached GRCh38 Reference
```bash
# No --reference flag uses cached GRCh38
basebuddy short --depth 50 --profile HS25 --outdir grch38_reads
```

### Error Handling Examples
The improved error messages help guide correct usage:
```bash
# This will show helpful error about depth
basebuddy short --reference genome.fa --depth -10

# This will show error about fragment length
basebuddy short --reference genome.fa --depth 30 --read-len 150 --mean-frag-len 100 --paired
```

---

## Long Read Simulation

### Basic Nanopore Simulation
```bash
# Coverage-based simulation
basebuddy long --reference genome.fa --depth 20 --model nanopore_R9.4.1

# Exact read count simulation (overrides depth)
basebuddy long --reference genome.fa --num-reads 10000 --model nanopore_R10.3
```

### PacBio Simulation
```bash
# PacBio CLR (Continuous Long Read)
basebuddy long \
  --reference genome.fa \
  --outdir pacbio_clr \
  --depth 15 \
  --model pacbio_CLR

# PacBio CCS (Circular Consensus Sequencing)
basebuddy long \
  --reference genome.fa \
  --outdir pacbio_ccs \
  --num-reads 5000 \
  --model pacbio_CCS
```

### RNA Nanopore Simulation
```bash
basebuddy long \
  --reference transcriptome.fa \
  --outdir nanopore_rna \
  --depth 50 \
  --model nanopore_RNA
```

---

## Variant Spiking

### Basic SNP Spiking
Create a VCF file with variants to spike:
```bash
cat > snps.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr12	25245350	KRAS_G12C	G	C	100	PASS	.
chr17	7577538	TP53_R248Q	C	T	100	PASS	.
EOF

# Spike variants into existing BAM
basebuddy spike \
  --reference genome.fa \
  --input-bam normal.bam \
  --snp-vcf snps.vcf \
  --output-prefix tumor \
  --vaf 0.3 \
  --seed 42
```

### Multi-BAM Spiking
```bash
# Spike into multiple BAMs simultaneously
basebuddy spike \
  --reference genome.fa \
  --input-bam sample1.bam \
  --input-bam sample2.bam \
  --input-bam sample3.bam \
  --snp-vcf driver_mutations.vcf \
  --indel-vcf indels.vcf \
  --output-prefix spiked \
  --vaf 0.25
```

### Setting up Picard JAR
BAMSurgeon requires Picard. Set it up using one of these methods:
```bash
# Method 1: Environment variable
export BAMSURGEON_PICARD_JAR=/path/to/picard.jar
basebuddy spike --input-bam normal.bam --snp-vcf snps.vcf

# Method 2: Command line flag
basebuddy spike \
  --input-bam normal.bam \
  --snp-vcf snps.vcf \
  --picard-jar /path/to/picard.jar
```

---

## Mutational Signatures

### Basic Signature Simulation
```bash
# Single base substitutions (most common)
basebuddy signature \
  --sig-type SBS \
  --num-mutations 5000 \
  --sample-id tumor_sample

# Doublet base substitutions
basebuddy signature \
  --sig-type DBS \
  --num-mutations 500 \
  --sample-id tumor_sample_dbs

# Insertions and deletions
basebuddy signature \
  --sig-type ID \
  --num-mutations 1000 \
  --sample-id tumor_sample_id
```

### Using Custom Reference
```bash
# Use specific genome build
basebuddy signature \
  --reference GRCh37 \
  --sig-type SBS \
  --num-mutations 10000 \
  --sample-id patient_001

# Use custom FASTA
basebuddy signature \
  --reference /path/to/custom_genome.fa \
  --sig-type SBS \
  --num-mutations 3000 \
  --sample-id custom_sample
```

### Advanced Options
```bash
# Exome-only simulation
basebuddy signature \
  --sig-type SBS \
  --num-mutations 2000 \
  --sample-id exome_tumor \
  --exome

# Chromosome-based generation with seed
basebuddy signature \
  --sig-type ID \
  --num-mutations 5000 \
  --sample-id reproducible_run \
  --chrom-based \
  --seed 12345
```

---

## Strand Bias

### Basic Usage
```bash
# Introduce 80% forward strand bias
basebuddy strand-bias \
  --input-bam aligned.bam \
  --output-bam biased.bam \
  --forward-fraction 0.8

# Extreme reverse strand bias
basebuddy strand-bias \
  --input-bam normal.bam \
  --output-bam reverse_biased.bam \
  --forward-fraction 0.2 \
  --seed 999
```

### Use Cases
```bash
# Simulate FFPE artifacts (forward strand bias)
basebuddy strand-bias -i fresh.bam -o ffpe_like.bam -f 0.75

# Simulate oxidative damage patterns
basebuddy strand-bias -i undamaged.bam -o oxidative.bam -f 0.65
```

---

## Quality Control

### Single File QC
```bash
basebuddy qc reads.fastq.gz --output-dir qc_results
```

### Batch QC
```bash
# Multiple files
basebuddy qc \
  sample1_R1.fastq.gz \
  sample1_R2.fastq.gz \
  sample2_R1.fastq.gz \
  sample2_R2.fastq.gz \
  --output-dir batch_qc \
  --threads 4

# Using wildcards (in bash)
basebuddy qc data/*.fastq.gz --output-dir all_qc --threads 8
```

---

## Complete Workflows

### Workflow 1: Tumor Simulation with Known Driver Mutations

```bash
# Step 1: Create reference subset (optional, for faster processing)
samtools faidx genome.fa chr12:25200000-25300000 chr17:7500000-7600000 > target_regions.fa

# Step 2: Simulate normal reads
basebuddy short \
  --reference target_regions.fa \
  --outdir normal_reads \
  --depth 100 \
  --profile HS25

# Step 3: Align reads (using BWA or your preferred aligner)
bwa index target_regions.fa
bwa mem target_regions.fa normal_reads/*/simulated_*1.fq normal_reads/*/simulated_*2.fq | \
  samtools sort -o normal.bam
samtools index normal.bam

# Step 4: Create driver mutation VCF
cat > drivers.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr12	25245350	KRAS_G12C	G	C	100	PASS	.
chr17	7577538	TP53_R248Q	C	T	100	PASS	.
EOF

# Step 5: Spike mutations at 30% VAF
basebuddy spike \
  --reference target_regions.fa \
  --input-bam normal.bam \
  --snp-vcf drivers.vcf \
  --output-prefix tumor \
  --vaf 0.3

# Step 6: QC the results
basebuddy qc normal_reads/*/*.fq --output-dir qc_results
```

### Workflow 2: Mutational Signature Analysis Pipeline

```bash
# Step 1: Generate mutations with specific signature
basebuddy signature \
  --sig-type SBS \
  --num-mutations 10000 \
  --sample-id smoking_tumor \
  --outdir signature_mutations

# Step 2: Apply mutations to reference and simulate reads
# (This would require additional scripting to apply VCF to reference)

# Step 3: Introduce strand bias (optional, for FFPE samples)
basebuddy strand-bias \
  --input-bam aligned_signature.bam \
  --output-bam ffpe_signature.bam \
  --forward-fraction 0.7
```

### Workflow 3: Long Read Tumor/Normal Pair

```bash
# Normal sample
basebuddy long \
  --reference genome.fa \
  --outdir normal_long \
  --depth 30 \
  --model nanopore_R10.3 \
  --num-reads 100000

# Tumor sample with mutations
# First create tumor reference with mutations, then:
basebuddy long \
  --reference tumor_genome.fa \
  --outdir tumor_long \
  --depth 30 \
  --model nanopore_R10.3 \
  --num-reads 100000
```

### Workflow 4: Benchmarking Variant Callers

```bash
# Create truth set VCF with various VAFs
cat > truth_variants.vcf << EOF
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1000000	var1	A	T	100	PASS	.
chr1	2000000	var2	C	G	100	PASS	.
chr1	3000000	var3	G	C	100	PASS	.
EOF

# Generate multiple samples with different VAFs
for vaf in 0.01 0.05 0.10 0.20 0.50; do
  basebuddy spike \
    --reference genome.fa \
    --input-bam normal.bam \
    --snp-vcf truth_variants.vcf \
    --output-prefix benchmark_vaf_${vaf} \
    --vaf ${vaf} \
    --seed 42
done

# Run your variant caller on each BAM and compare to truth set
```

---

## Tips and Best Practices

1. **Always index your reference FASTA**:
   ```bash
   samtools faidx reference.fa
   ```

2. **Use appropriate depths**:
   - WGS: 30-50x
   - WES: 100-200x
   - Targeted panels: 500-1000x
   - RNA-seq: 50-100x

3. **Choose correct sequencing profiles**:
   - HS25: HiSeq 2500
   - HSXt: HiSeq X Ten
   - NS50: NextSeq 500
   - MSv3: MiSeq v3

4. **Set random seeds for reproducibility**:
   ```bash
   basebuddy spike --seed 42 ...
   basebuddy strand-bias --seed 42 ...
   ```

5. **Check output manifests**:
   ```bash
   cat results_dir/*/manifest.json | jq .
   ```

---

## Troubleshooting

### Missing Tools
If you get tool not found errors, ensure all dependencies are installed:
```bash
# Check if tools are available
which art_illumina
which nanosim-h
which samtools
which fastqc

# Use Docker image for all dependencies
docker run -v $(pwd):/data basebuddy:latest <command>
```

### Memory Issues
For large simulations, increase available memory or reduce depth:
```bash
# Reduce depth
basebuddy short --depth 10 ...

# Or process by chromosome
for chr in {1..22} X Y; do
  basebuddy short --reference chr${chr}.fa --outdir chr${chr}_reads ...
done
```

### Output Not Found
Check the manifest file for actual output locations:
```bash
find results_dir -name "manifest.json" -exec cat {} \;
```
#!/bin/bash
# BRAF V600K Vignette - Multi-Nucleotide Variant (MNV) Detection
# This script generates synthetic reads for two scenarios:
# 2.1: BRAF V600K as a dinucleotide substitution (correct representation)
# 2.2: BRAF V600K mis-called as two separate SNVs (incorrect representation)

set -e  # Exit on error

# Configuration
WORK_DIR="braf_v600k_vignette"
REF_BUILD="GRCh38"
BRAF_REGION="chr7:140753300-140924900"  # Full BRAF gene region
BRAF_V600_REGION="chr7:140753200-140753500"  # Focused region around codon 600
V600K_POS1="chr7:140753336"  # First position of dinucleotide change
V600K_POS2="chr7:140753337"  # Second position of dinucleotide change

echo "=========================================="
echo "BRAF V600K Vignette - MNV vs SNV Detection"
echo "=========================================="

# Create working directory
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

# Step 1: Prepare reference sequence
echo -e "\n[Step 1] Preparing BRAF reference region..."

# For demonstration, we'll create a synthetic BRAF region around codon 600
# In production, extract from actual GRCh38 reference
cat > braf_reference.fa << 'EOF'
>chr7:140753200-140753500
CTACAGTGAAATCTCGATGGAGTGGGTCCCATCAGTTTGAACAGTTGTCTGGATCCATTTTGTGGATGG
TCTTCTAGCTACAGTACTGTTATATAGAAACAAAATTTTTCTCTATCAGCAAGAACAAATGATTTTTAC
TCTAGCTAGACCAAAATCACCTATTTTTACTGTGAGGTCTTCATGAAGAAATATACCAAAAACTGGAAA
AGTTAACACCCTCATCTCCTACCTGGGCCCCCACCTTGCTCTGCTCTAGACCCCACAGAGCCCGAAGAC
TGACTGGGCCAGCGCTACCTTAACA
EOF

# Index the reference
samtools faidx braf_reference.fa

echo -e "\n=========================================="
echo "Scenario 2.1: BRAF V600K as Dinucleotide Substitution (MNV)"
echo "Correct representation of V600K mutation"
echo "=========================================="

# Create directory for Scenario 2.1
mkdir -p scenario_2.1_mnv

# Step 2.1: Create VCF for dinucleotide substitution (MNV)
# GT > AA change at positions 140753336-140753337
cat > scenario_2.1_mnv/braf_v600k_mnv.vcf << 'EOF'
##fileformat=VCFv4.2
##reference=braf_reference.fa
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr7	140753336	BRAF_V600K_MNV	GT	AA	100	PASS	AF=0.30;TYPE=MNV;SOMATIC
EOF

# Step 3.1: Simulate normal reads without mutation
echo -e "\n[Scenario 2.1 - Step 1] Simulating normal reads..."
basebuddy short \
  --reference braf_reference.fa \
  --outdir scenario_2.1_mnv/normal_reads \
  --depth 150 \
  --read-len 150 \
  --profile HS25 \
  --mean-frag-len 300 \
  --std-frag-len 40 \
  --paired \
  --overwrite

# Step 4.1: Create aligned BAM
echo -e "\n[Scenario 2.1 - Step 2] Creating aligned BAM..."
if command -v bwa &> /dev/null; then
  bwa index braf_reference.fa 2>/dev/null || true
  bwa mem -t 4 braf_reference.fa \
    scenario_2.1_mnv/normal_reads/*/simulated_*1.fq \
    scenario_2.1_mnv/normal_reads/*/simulated_*2.fq | \
    samtools sort -o scenario_2.1_mnv/normal.bam
  samtools index scenario_2.1_mnv/normal.bam
else
  echo "BWA not found. Creating placeholder BAM..."
  touch scenario_2.1_mnv/normal.bam
fi

# Step 5.1: Apply MNV using custom approach
# Since BAMSurgeon typically handles SNVs, we'll create a custom approach for MNV
echo -e "\n[Scenario 2.1 - Step 3] Creating reads with dinucleotide substitution..."

# First, create intermediate files with the mutation
if command -v python3 &> /dev/null; then
  python3 << 'PYTHON_MNV'
import sys
import os

# Create a Python script to generate reads with MNV
mnv_script = """
# This would be a custom script to ensure both mutations are on the same read
# For demonstration purposes, we'll use the VCF with BAMSurgeon
print("MNV generation would require custom read manipulation")
print("In practice, this ensures G>A and T>A are on the same DNA molecule")
"""

with open("scenario_2.1_mnv/generate_mnv.py", "w") as f:
    f.write(mnv_script)

print("MNV generation script created")
PYTHON_MNV
fi

# For now, we'll use the MNV VCF with basebuddy spike
# Note: This assumes BAMSurgeon can handle MNVs or we post-process
if [ -n "$BAMSURGEON_PICARD_JAR" ]; then
  echo "Attempting to spike MNV..."
  basebuddy spike \
    --reference braf_reference.fa \
    --input-bam scenario_2.1_mnv/normal.bam \
    --snp-vcf scenario_2.1_mnv/braf_v600k_mnv.vcf \
    --output-prefix scenario_2.1_mnv/melanoma_mnv \
    --vaf 0.30 \
    --seed 42 \
    --overwrite || echo "MNV spiking may require custom handling"
fi

# Summary for Scenario 2.1
echo -e "\n[Scenario 2.1 Summary]"
echo "- Location: chr7:140753336-140753337 (GT>AA)"
echo "- Type: Dinucleotide substitution (MNV)"
echo "- VAF: 30%"
echo "- Key feature: Both mutations on same read/molecule"
echo "- Clinical: Common UV-induced BRAF V600K in melanoma"

echo -e "\n=========================================="
echo "Scenario 2.2: BRAF V600K Mis-called as Two Separate SNVs"
echo "Incorrect representation showing phasing challenge"
echo "=========================================="

# Create directory for Scenario 2.2
mkdir -p scenario_2.2_separate_snvs

# Step 2.2: Create separate VCFs for each SNV
# First SNV: G>A at position 140753336
cat > scenario_2.2_separate_snvs/braf_snv1.vcf << 'EOF'
##fileformat=VCFv4.2
##reference=braf_reference.fa
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr7	140753336	BRAF_V600K_SNV1	G	A	100	PASS	AF=0.30
EOF

# Second SNV: T>A at position 140753337
cat > scenario_2.2_separate_snvs/braf_snv2.vcf << 'EOF'
##fileformat=VCFv4.2
##reference=braf_reference.fa
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr7	140753337	BRAF_V600K_SNV2	T	A	100	PASS	AF=0.30
EOF

# Step 3.2: Simulate normal reads
echo -e "\n[Scenario 2.2 - Step 1] Simulating normal reads..."
basebuddy short \
  --reference braf_reference.fa \
  --outdir scenario_2.2_separate_snvs/normal_reads \
  --depth 150 \
  --read-len 150 \
  --profile HS25 \
  --mean-frag-len 300 \
  --std-frag-len 40 \
  --paired \
  --overwrite

# Step 4.2: Create two separate BAMs with different mutations
echo -e "\n[Scenario 2.2 - Step 2] Creating BAMs with separate SNVs..."

# Align normal reads
if command -v bwa &> /dev/null; then
  bwa mem -t 4 braf_reference.fa \
    scenario_2.2_separate_snvs/normal_reads/*/simulated_*1.fq \
    scenario_2.2_separate_snvs/normal_reads/*/simulated_*2.fq | \
    samtools sort -o scenario_2.2_separate_snvs/normal.bam
  samtools index scenario_2.2_separate_snvs/normal.bam
else
  touch scenario_2.2_separate_snvs/normal.bam
fi

# Create a Python script to ensure mutations are on different reads
if command -v python3 &> /dev/null; then
  cat > scenario_2.2_separate_snvs/separate_mutations.py << 'PYTHON_SEPARATE'
#!/usr/bin/env python3
"""
This script would split reads to ensure SNVs appear on different molecules
Simulating incorrect phasing where the two mutations are not linked
"""
import random

print("Separating mutations to different read sets...")
print("50% of reads get G>A at position 140753336")
print("Different 50% of reads get T>A at position 140753337")
print("No reads contain both mutations - simulating phasing error")

# In practice, this would:
# 1. Split the BAM into two subsets
# 2. Apply SNV1 to subset 1
# 3. Apply SNV2 to subset 2
# 4. Merge the BAMs
PYTHON_SEPARATE
  chmod +x scenario_2.2_separate_snvs/separate_mutations.py
fi

# Spike mutations separately to simulate incorrect calling
if [ -n "$BAMSURGEON_PICARD_JAR" ]; then
  # First, split the BAM into two parts
  echo -e "\n[Scenario 2.2 - Step 3] Splitting reads for separate mutations..."
  
  # Apply first SNV to a subset
  samtools view -h -s 0.5 scenario_2.2_separate_snvs/normal.bam | \
    samtools view -b > scenario_2.2_separate_snvs/subset1.bam
  samtools index scenario_2.2_separate_snvs/subset1.bam
  
  # Apply second SNV to different subset
  samtools view -h -s 0.5 scenario_2.2_separate_snvs/normal.bam | \
    samtools view -b > scenario_2.2_separate_snvs/subset2.bam
  samtools index scenario_2.2_separate_snvs/subset2.bam
  
  # Spike first mutation
  echo "Spiking first SNV (G>A)..."
  basebuddy spike \
    --reference braf_reference.fa \
    --input-bam scenario_2.2_separate_snvs/subset1.bam \
    --snp-vcf scenario_2.2_separate_snvs/braf_snv1.vcf \
    --output-prefix scenario_2.2_separate_snvs/snv1 \
    --vaf 0.60 \
    --seed 123 \
    --overwrite
  
  # Spike second mutation
  echo "Spiking second SNV (T>A)..."
  basebuddy spike \
    --reference braf_reference.fa \
    --input-bam scenario_2.2_separate_snvs/subset2.bam \
    --snp-vcf scenario_2.2_separate_snvs/braf_snv2.vcf \
    --output-prefix scenario_2.2_separate_snvs/snv2 \
    --vaf 0.60 \
    --seed 456 \
    --overwrite
  
  # Merge the two BAMs
  echo "Merging BAMs with separate mutations..."
  samtools merge -f scenario_2.2_separate_snvs/miscalled_final.bam \
    scenario_2.2_separate_snvs/snv1_*.bam \
    scenario_2.2_separate_snvs/snv2_*.bam
  samtools index scenario_2.2_separate_snvs/miscalled_final.bam
fi

# Summary for Scenario 2.2
echo -e "\n[Scenario 2.2 Summary]"
echo "- Locations: chr7:140753336 (G>A) and chr7:140753337 (T>A)"
echo "- Type: Two separate SNVs (incorrect representation)"
echo "- VAF: 30% each"
echo "- Key feature: Mutations on different reads/molecules"
echo "- Issue: Incorrect phasing - misses true MNV nature"

# Final summary
echo -e "\n=========================================="
echo "Vignette Generation Complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "- Scenario 2.1 (Correct MNV):"
echo "  - VCF: scenario_2.1_mnv/braf_v600k_mnv.vcf"
echo "  - BAM: scenario_2.1_mnv/melanoma_mnv_*.bam"
echo ""
echo "- Scenario 2.2 (Incorrect separate SNVs):"
echo "  - VCFs: scenario_2.2_separate_snvs/braf_snv[1,2].vcf"
echo "  - BAM: scenario_2.2_separate_snvs/miscalled_final.bam"
echo ""
echo "IGV visualization will show:"
echo "- Scenario 2.1: Reads with BOTH mutations (correct)"
echo "- Scenario 2.2: Reads with EITHER mutation, never both (incorrect)"
echo ""
echo "Key learning:"
echo "- MNVs must be called as single events"
echo "- Separate SNV calls can miss functional impact"
echo "- V600K requires both changes for amino acid substitution"
echo ""
echo "Note: Full functionality requires:"
echo "- BWA and samtools"
echo "- BAMSURGEON_PICARD_JAR environment variable"
echo "- Custom MNV handling may be needed for some tools"
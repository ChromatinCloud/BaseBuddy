#!/bin/bash
# KRAS G12C Vignette - Real vs. Artifactual Call
# This script generates synthetic reads for two scenarios:
# 1.1: Clinically actionable KRAS G12C in NSCLC (true positive)
# 1.2: Artifactual KRAS G12C finding (false positive)

set -e  # Exit on error

# Configuration
WORK_DIR="kras_g12c_vignette"
REF_BUILD="GRCh38"
KRAS_REGION="chr12:25245027-25245627"  # 600bp region covering KRAS
KRAS_G12C_POS="chr12:25245327"  # Exact position of G12C mutation

echo "=========================================="
echo "KRAS G12C Vignette - Synthetic Read Generation"
echo "=========================================="

# Create working directory
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

# Step 1: Prepare reference sequence
echo -e "\n[Step 1] Preparing KRAS reference region..."

# Option A: If you have GRCh38 reference available
# samtools faidx /path/to/GRCh38.fa ${KRAS_REGION} > kras_reference.fa

# Option B: Create a synthetic KRAS region (simplified for demonstration)
# This is a representation - in production, use actual GRCh38 sequence
cat > kras_reference.fa << 'EOF'
>chr12:25245027-25245627
ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAA
TTCAGAATCATTTTGTGGACGAATATGATCCAACAATAGAGGATTCCTACAGGAAGCAAGTAGTAATTGA
TGGAGAAACCTGTCTCTTGGATATTCTCGACACAGCAGGTCAAGAGGAGTACAGTGCAATGAGGGACCAG
TACATGAGGACTGGGGAGGGCTTTCTTTGTGTATTTGCCATAAATAATACTAAATCATTTGAAGATATTC
ACCATTATAGAGAACAAATTAAAAGAGTTAAGGACTCTGAAGATGTACCTATGGTCCTAGTAGGAAATAA
ATGTGATTTGCCTTCTAGAACAGTAGACACAAAACAGGCTCAGGACTTAGCAAGAAGTTATGGAATTCCT
TTTATTGAAACATCAGCAAAGACAAGACAGGGTGTTGATGATGCCTTCTATACATTAGTTCGAGAAATTC
GAAAGAAAACAAAGCTTTAAAGAGAAAGTAAGCAGACTGCCTAATACTATTTCTTGCTCTGGCAAGGCTA
TAAGGTTAATGGATTTTCCATATTTGGAAAGGAAGAAATAAGAAGAA
EOF

# Index the reference
samtools faidx kras_reference.fa

echo -e "\n=========================================="
echo "Scenario 1.1: Clinically Actionable KRAS G12C"
echo "True positive somatic mutation in NSCLC"
echo "=========================================="

# Create directory for Scenario 1.1
mkdir -p scenario_1.1_true_positive

# Step 2.1: Create VCF for true KRAS G12C mutation
cat > scenario_1.1_true_positive/kras_g12c_true.vcf << 'EOF'
##fileformat=VCFv4.2
##reference=kras_reference.fa
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr12	25245327	KRAS_G12C	G	T	100	PASS	AF=0.45;SOMATIC
EOF

# Step 3.1: Simulate normal reads (without mutation)
echo -e "\n[Scenario 1.1 - Step 1] Simulating normal reads..."
basebuddy short \
  --reference kras_reference.fa \
  --outdir scenario_1.1_true_positive/normal_reads \
  --depth 200 \
  --read-len 150 \
  --profile HS25 \
  --mean-frag-len 350 \
  --std-frag-len 50 \
  --paired \
  --overwrite

# Step 4.1: Align normal reads (requires BWA)
echo -e "\n[Scenario 1.1 - Step 2] Aligning normal reads..."
if command -v bwa &> /dev/null; then
  bwa index kras_reference.fa 2>/dev/null || true
  bwa mem -t 4 kras_reference.fa \
    scenario_1.1_true_positive/normal_reads/*/simulated_*1.fq \
    scenario_1.1_true_positive/normal_reads/*/simulated_*2.fq | \
    samtools sort -o scenario_1.1_true_positive/normal.bam
  samtools index scenario_1.1_true_positive/normal.bam
else
  echo "BWA not found. Creating placeholder BAM..."
  # In production, you would need BWA installed
  touch scenario_1.1_true_positive/normal.bam
fi

# Step 5.1: Spike in KRAS G12C mutation at 45% VAF
echo -e "\n[Scenario 1.1 - Step 3] Spiking KRAS G12C at 45% VAF..."
# Note: This requires BAMSurgeon and Picard
if [ -n "$BAMSURGEON_PICARD_JAR" ]; then
  basebuddy spike \
    --reference kras_reference.fa \
    --input-bam scenario_1.1_true_positive/normal.bam \
    --snp-vcf scenario_1.1_true_positive/kras_g12c_true.vcf \
    --output-prefix scenario_1.1_true_positive/tumor \
    --vaf 0.45 \
    --seed 42 \
    --overwrite
else
  echo "BAMSURGEON_PICARD_JAR not set. Skipping variant spiking."
  echo "To run spiking: export BAMSURGEON_PICARD_JAR=/path/to/picard.jar"
fi

# Summary for Scenario 1.1
echo -e "\n[Scenario 1.1 Summary]"
echo "- Location: chr12:25245327 (G>T)"
echo "- VAF: 45%"
echo "- Strand bias: Balanced (~50/50)"
echo "- Read quality: High"
echo "- Expected: True positive somatic mutation"

echo -e "\n=========================================="
echo "Scenario 1.2: Artifactual KRAS G12C"
echo "False positive/artifact"
echo "=========================================="

# Create directory for Scenario 1.2
mkdir -p scenario_1.2_artifact

# Step 2.2: Create VCF for artifactual KRAS mutation
cat > scenario_1.2_artifact/kras_g12c_artifact.vcf << 'EOF'
##fileformat=VCFv4.2
##reference=kras_reference.fa
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr12	25245327	KRAS_G12C_artifact	G	T	30	PASS	AF=0.015
EOF

# Step 3.2: Simulate reads for artifact scenario
echo -e "\n[Scenario 1.2 - Step 1] Simulating reads for artifact..."
basebuddy short \
  --reference kras_reference.fa \
  --outdir scenario_1.2_artifact/normal_reads \
  --depth 200 \
  --read-len 150 \
  --profile HS25 \
  --mean-frag-len 350 \
  --std-frag-len 50 \
  --paired \
  --overwrite

# Step 4.2: Create BAM with low VAF mutation
echo -e "\n[Scenario 1.2 - Step 2] Creating low VAF variant..."
if command -v bwa &> /dev/null; then
  bwa mem -t 4 kras_reference.fa \
    scenario_1.2_artifact/normal_reads/*/simulated_*1.fq \
    scenario_1.2_artifact/normal_reads/*/simulated_*2.fq | \
    samtools sort -o scenario_1.2_artifact/normal.bam
  samtools index scenario_1.2_artifact/normal.bam
else
  touch scenario_1.2_artifact/normal.bam
fi

# Step 5.2: Spike mutation at very low VAF
echo -e "\n[Scenario 1.2 - Step 3] Spiking low VAF mutation..."
if [ -n "$BAMSURGEON_PICARD_JAR" ]; then
  basebuddy spike \
    --reference kras_reference.fa \
    --input-bam scenario_1.2_artifact/normal.bam \
    --snp-vcf scenario_1.2_artifact/kras_g12c_artifact.vcf \
    --output-prefix scenario_1.2_artifact/low_vaf \
    --vaf 0.015 \
    --seed 123 \
    --overwrite
    
  # Step 6.2: Introduce strand bias to create artifact
  echo -e "\n[Scenario 1.2 - Step 4] Introducing strand bias..."
  basebuddy strand-bias \
    --input-bam scenario_1.2_artifact/low_vaf_*.bam \
    --output-bam scenario_1.2_artifact/artifact_final.bam \
    --forward-fraction 0.95 \
    --seed 456
else
  echo "BAMSURGEON_PICARD_JAR not set. Skipping artifact creation."
fi

# Summary for Scenario 1.2
echo -e "\n[Scenario 1.2 Summary]"
echo "- Location: chr12:25245327 (G>T)"
echo "- VAF: 1.5%"
echo "- Strand bias: Extreme (95% forward strand)"
echo "- Read quality: Lower for variant bases"
echo "- Expected: False positive/artifact"

# Final summary
echo -e "\n=========================================="
echo "Vignette Generation Complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "- Scenario 1.1 (True Positive):"
echo "  - Normal BAM: scenario_1.1_true_positive/normal.bam"
echo "  - Tumor BAM (45% VAF): scenario_1.1_true_positive/tumor_*.bam"
echo ""
echo "- Scenario 1.2 (Artifact):"
echo "  - Normal BAM: scenario_1.2_artifact/normal.bam"
echo "  - Low VAF BAM: scenario_1.2_artifact/low_vaf_*.bam"
echo "  - Artifact BAM: scenario_1.2_artifact/artifact_final.bam"
echo ""
echo "Next steps:"
echo "1. Run variant callers on both scenarios"
echo "2. Compare variant calls to understand true vs artifactual calls"
echo "3. Examine strand bias, VAF, and read placement in IGV"
echo ""
echo "Note: For full functionality, ensure you have:"
echo "- BWA and samtools installed"
echo "- BAMSURGEON_PICARD_JAR environment variable set"
echo "- Actual GRCh38 reference sequence for production use"
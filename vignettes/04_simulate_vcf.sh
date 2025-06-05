#!/usr/bin/env bash
set -euo pipefail

echo "### Vignette 04: Simulate Reads Directly from VCF ###"

mkdir -p vignette_04_output

# Prepare a reference FASTA
cat > vignette_04_output/ref_simvcf.fa <<EOF
>chrY_vcf
GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA
GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA
EOF
docker run --rm -v "${PWD}/vignette_04_output:/data" -u "$(id -u):$(id -g)" staphb/samtools:latest samtools faidx /data/ref_simvcf.fa

# Prepare a VCF with germline and somatic variants
cat > vignette_04_output/variants_for_sim.vcf <<EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chrY_vcf,length=160>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	GERMLINE_SAMPLE	TUMOR_SAMPLE
chrY_vcf	15	germ_snp_het	A	T	100	PASS	.	GT	0/1	0/0
chrY_vcf	30	germ_snp_hom	C	G	100	PASS	.	GT	1/1	0/0
chrY_vcf	45	somatic_snp	G	C	100	PASS	.	GT:AF	0/0:0.0	0/1:0.3
chrY_vcf	60	somatic_indel	TACA	T	100	PASS	.	GT:AF	0/0:0.0	0/1:0.15
EOF

# Run basebuddy simvcf (hypothetical subcommand)
# This command would simulate reads based on the reference and incorporate variants from the VCF.
# For TUMOR_SAMPLE, it would use AF; for GERMLINE_SAMPLE, GT.
docker run --rm \
  -v "${PWD}/vignette_04_output:/data" \
  basebuddy simvcf /data/ref_simvcf.fa /data/variants_for_sim.vcf \
    --sample TUMOR_SAMPLE \
    --depth 30 \
    --readlen 100 \
    --outdir /data/sim_reads_from_vcf_tumor \
    --seed 456

docker run --rm \
  -v "${PWD}/vignette_04_output:/data" \
  basebuddy simvcf /data/ref_simvcf.fa /data/variants_for_sim.vcf \
    --sample GERMLINE_SAMPLE \
    --depth 20 \
    --readlen 100 \
    --outdir /data/sim_reads_from_vcf_germline \
    --seed 789


# Verify the output FASTQs
echo "Verifying output FASTQ files..."
if [[ -s vignette_04_output/sim_reads_from_vcf_tumor/reads_R1.fq.gz && \
      -s vignette_04_output/sim_reads_from_vcf_tumor/reads_R2.fq.gz && \
      -s vignette_04_output/sim_reads_from_vcf_germline/reads_R1.fq.gz && \
      -s vignette_04_output/sim_reads_from_vcf_germline/reads_R2.fq.gz ]]; then
  echo "Vignette 04: Simulate from VCF example SUCCEEDED."
  echo "Output FASTQs are in vignette_04_output/sim_reads_from_vcf_tumor/ and sim_reads_from_vcf_germline/"
else
  echo "ERROR (Vignette 04): Output FASTQ files from VCF simulation missing or empty." >&2
  exit 1
fi

echo "Cleaning up vignette_04_output..."
# rm -rf vignette_04_output # Uncomment to clean up
echo "### Vignette 04 Complete ###"
echo ""


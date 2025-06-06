#!/usr/bin/env bash
set -euo pipefail

# Create a 100-bp FASTA and index it
cat > toy.fa <<EOF
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
EOF
samtools faidx toy.fa

# Align an empty read to produce a BAM (so BAMSurgeon can modify it):
echo -e "@HD\tVN:1.0\tSO:unsorted" > toy.sam
echo -e "r1\t0\tchr1\t1\t60\t100M\t*\t0\t0\tACGT...ACGT\tFFFFFFFFF" >> toy.sam
samtools view -bS toy.sam > toy.bam
samtools sort toy.bam -o toy.sorted.bam
mv toy.sorted.bam toy.bam
samtools index toy.bam

# Make a simple VCF with one SNV at chr1:10
cat > variants.vcf <<EOF
##fileformat=VCFv4.2
#CHROM POS ID REF ALT QUAL FILTER INFO
chr1 10 . A G . . .
EOF

# Run spike
docker run --rm \
  -v "${PWD}/toy.fa:/data/toy.fa" \
  -v "${PWD}/toy.fa.fai:/data/toy.fa.fai" \
  -v "${PWD}/toy.bam:/data/toy.bam" \
  -v "${PWD}/toy.bam.bai:/data/toy.bam.bai" \
  -v "${PWD}/variants.vcf:/data/variants.vcf" \
  -v "${PWD}:/data" \
  basebuddy spike /data/toy.fa /data/toy.bam /data/variants.vcf \
    --vaf 0.1 --out-bam /data/spiked.bam --seed 42

# Verify the spiked.bam exists and can be indexed
if samtools quickcheck spiked.bam; then
  echo "Spike-variants example succeeded."
else
  echo "ERROR: spiked.bam failed quickcheck." >&2
  exit 1
fi

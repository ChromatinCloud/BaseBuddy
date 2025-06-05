#!/usr/bin/env bash
set -euo pipefail

echo "### Vignette 05: Simulate Reads with Structural Variants ###"

mkdir -p vignette_05_output

# Prepare a larger reference FASTA for SVs
cat > vignette_05_output/ref_sv.fa <<EOF
>chrS
GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA
GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA
GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA
GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA
GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA
EOF
docker run --rm -v "${PWD}/vignette_05_output:/data" -u "$(id -u):$(id -g)" staphb/samtools:latest samtools faidx /data/ref_sv.fa

# Prepare a BEDPE file describing a deletion and a duplication
# chr1 pos1 pos2 chr1 pos3 pos4 name score strand1 strand2 type
cat > vignette_05_output/structural_variants.bedpe <<EOF
chrS	50	51	chrS	100	101	DEL001	.	+	-	DEL
chrS	150	151	chrS	200	201	DUP001	.	+	+	DUP:TANDEM
EOF

# Run basebuddy simsv (hypothetical subcommand for SV simulation)
# This could also be an option within 'short' or 'long' e.g., --sv-bedpe
docker run --rm \
  -v "${PWD}/vignette_05_output:/data" \
  basebuddy simsv /data/ref_sv.fa /data/structural_variants.bedpe \
    --depth 40 \
    --readlen 150 \
    --insert-mean 500 \
    --insert-sd 75 \
    --outdir /data/sim_reads_with_sv \
    --seed 101

# Verify the output FASTQs (existence is a basic check; actual SV validation is complex)
echo "Verifying output FASTQ files..."
if [[ -s vignette_05_output/sim_reads_with_sv/reads_R1.fq.gz && \
      -s vignette_05_output/sim_reads_with_sv/reads_R2.fq.gz ]]; then
  echo "Vignette 05: Simulate Structural Variants example SUCCEEDED."
  echo "Output FASTQs are in vignette_05_output/sim_reads_with_sv/"
  echo "Note: Verifying SV presence in reads requires alignment and SV calling."
else
  echo "ERROR (Vignette 05): Output FASTQ files from SV simulation missing or empty." >&2
  exit 1
fi

echo "Cleaning up vignette_05_output..."
# rm -rf vignette_05_output # Uncomment to clean up
echo "### Vignette 05 Complete ###"
echo ""


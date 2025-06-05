#!/usr/bin/env bash
set -euo pipefail

echo "### Vignette 02: Advanced Long-Read Simulation ###"

# Prepare a tiny FASTA (can reuse if previous vignette created it, or make specific)
mkdir -p vignette_02_output
cat > vignette_02_output/long_ref.fa <<EOF
>chrM_long
AAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGG
AAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGG
AAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGG
AAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGGAAACCCTTTGGG
EOF

# Run advanced long-read simulation
# Assuming new options: --length-mean, --length-sd, --error-model for long reads
docker run --rm \
  -v "${PWD}/vignette_02_output:/data" \
  basebuddy long /data/long_ref.fa \
    --depth 15 \
    --model PacBio_HiFi \
    --length-mean 100 \
    --length-sd 20 \
    --outdir /data/sim_reads_advanced_long

# Verify the output FASTQ
echo "Verifying output FASTQ file..."
if [[ -s vignette_02_output/sim_reads_advanced_long/reads.fq.gz ]]; then
  echo "Vignette 02: Advanced Long-read simulation example SUCCEEDED."
  echo "Output FASTQ is in vignette_02_output/sim_reads_advanced_long/"
else
  echo "ERROR (Vignette 02): Output FASTQ file missing or empty." >&2
  exit 1
fi

echo "Cleaning up vignette_02_output..."
# rm -rf vignette_02_output # Uncomment to clean up
echo "### Vignette 02 Complete ###"
echo ""


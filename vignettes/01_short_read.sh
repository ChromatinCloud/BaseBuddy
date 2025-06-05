#!/usr/bin/env bash
set -euo pipefail

# 1a. Prepare a tiny FASTA
cat > mini.fa <<EOF
>chr1
ACGTACGTACGTACGT
EOF

# 1b. Run short-read simulation
docker run --rm \
  -v "${PWD}/mini.fa:/data/mini.fa" \
  -v "${PWD}:/data" \
  basebuddy short /data/mini.fa --depth 3 --readlen 50 --outdir /data/out_short

# 1c. Verify the output FASTQs
if [[ -s out_short/reads1.fq && -s out_short/reads2.fq ]]; then
  echo "Short-read example succeeded."
else
  echo "ERROR: fastq files missing or empty" >&2
  exit 1
fi

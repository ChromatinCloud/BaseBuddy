#!/usr/bin/env bash
set -euo pipefail

echo "1️⃣ Rebuilding Docker image..."
DOCKER_BUILDKIT=1 docker build --progress=plain -t basebuddy:latest .

echo "2️⃣ Checking version..."
docker run --rm basebuddy version

echo "3️⃣ CLI help..."
docker run --rm basebuddy --help

echo "4️⃣ Subcommand help..."
for subcmd in short long spike signature strand-bias; do
  echo "----- basebuddy $subcmd --help -----"
  docker run --rm basebuddy $subcmd --help || true
done

# Optional: real smoke test for each command if files exist
ROOT="$(pwd)"

echo "5️⃣ Testing short-read sim (ART)..."
if [ -f "$ROOT/mini.fa" ]; then
  docker run --rm \
    -v "$ROOT/mini.fa:/data/mini.fa" \
    basebuddy short /data/mini.fa --depth 3 --outdir /data/results_short || echo "short-read sim failed"
else
  echo "SKIP: mini.fa not found"
fi

echo "6️⃣ Testing long-read sim (NanoSim-h)..."
if [ -f "$ROOT/mini.fa" ]; then
  docker run --rm \
    -v "$ROOT/mini.fa:/data/mini.fa" \
    basebuddy long /data/mini.fa --depth 1 --model nanopore_R9.4.1 --outdir /data/results_long || echo "long-read sim failed"
else
  echo "SKIP: mini.fa not found"
fi

echo "7️⃣ Testing spike-in (BAMSurgeon)..."
if [ -f "$ROOT/mini.fa" ] && [ -f "$ROOT/spike_in.bam" ] && [ -f "$ROOT/spike_in.bam.bai" ] && [ -f "$ROOT/variants.vcf" ]; then
  docker run --rm \
    -v "$ROOT/mini.fa:/data/mini.fa" \
    -v "$ROOT/spike_in.bam:/data/spike_in.bam" \
    -v "$ROOT/spike_in.bam.bai:/data/spike_in.bam.bai" \
    -v "$ROOT/variants.vcf:/data/variants.vcf" \
    basebuddy spike /data/mini.fa /data/spike_in.bam /data/variants.vcf \
      --vaf 0.05 \
      --out-bam /data/results_spiked.bam \
      --seed 123 || echo "spike-in failed"
else
  echo "SKIP: spike-in BAM or VCF missing"
fi

echo "8️⃣ Testing mutational signature (SigProfilerSimulator)..."
if [ -f "$ROOT/mini.fa" ]; then
  docker run --rm \
    -v "$ROOT/mini.fa:/data/mini.fa" \
    -v "$ROOT:/data" \
    basebuddy signature /data/mini.fa \
      --outdir /data/results_sig \
      --sig-type SBS \
      --num-mutations 20 \
      --sample-id testSig || echo "signature sim failed"
else
  echo "SKIP: mini.fa not found"
fi

echo "Smoke test complete."


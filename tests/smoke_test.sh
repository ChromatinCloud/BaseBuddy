#!/usr/bin/env bash
set -euo pipefail

# --- Configuration ---
# Set rebuild_image to 1 to force a rebuild of the Docker image.
# Set rebuild_image to 0 to skip the rebuild step and use the existing image.
rebuild_image=0

# --- Conditional Docker Image Rebuild ---
if [ "$rebuild_image" -eq 1 ]; then
    echo "1️⃣ Rebuilding Docker image..."
    DOCKER_BUILDKIT=1 docker build --progress=plain -t basebuddy:latest .
else
    echo "1️⃣ Skipping Docker image rebuild (rebuild_image=0)."
fi

echo "2️⃣ Checking version..."
docker run --rm basebuddy version

echo "3️⃣ CLI help..."
docker run --rm basebuddy --help

echo "4️⃣ Subcommand help..."
for subcmd in short long spike signature strand-bias; do # Assuming 'apply-signature' might be added here later or tested separately
  echo "----- basebuddy $subcmd --help -----"
  docker run --rm basebuddy $subcmd --help || true
done

echo "4.5️⃣ Testing local signature loading (apply-signature command)..."
# This test assumes 'apply-signature' CLI command will be added and will use
# the apply_signature_to_fasta runner.
# smoke_test.fa is expected at /opt/app/src/basebuddy/data/smoke_test_data/smoke_test.fa in the image
docker run --rm basebuddy:latest basebuddy apply-signature \
  /opt/app/src/basebuddy/data/smoke_test_data/smoke_test.fa \
  mutated_smoke.fa \
  SBS1 \
  10 \
  --output-root /tmp/apply_sig_smoke_test_output \
  --overwrite || echo "Smoke test for apply-signature with SBS1 FAILED (Note: CLI command might not be implemented yet)"
echo "Local signature loading test complete."

# Optional: real smoke test for each command if files exist
ROOT="$(pwd)" # This will be the repository root when running 'bash tests/smoke_test.sh'

echo "5️⃣ Testing short-read sim (ART)..."
if [ -f "$ROOT/mini.fa" ]; then
  docker run --rm \
    -v "$ROOT/mini.fa:/data/mini.fa" \
    -v "$ROOT/results_short_smoke:/data/results_short" \
    basebuddy short /data/mini.fa --depth 3 --outdir /data/results_short --overwrite || echo "short-read sim failed"
else
  echo "SKIP: mini.fa not found"
fi

echo "6️⃣ Testing long-read sim (NanoSim-h)..."
if [ -f "$ROOT/mini.fa" ]; then
  docker run --rm \
    -v "$ROOT/mini.fa:/data/mini.fa" \
    -v "$ROOT/results_long_smoke:/data/results_long" \
    basebuddy long /data/mini.fa --depth 1 --model nanopore_R9.4.1 --outdir /data/results_long --overwrite || echo "long-read sim failed"
else
  echo "SKIP: mini.fa not found"
fi

echo "7️⃣ Testing spike-in (addsnv.py based - was BAMSurgeon)..." # Assuming spike_variants runner uses addsnv.py
# This test needs review based on current spike_variants runner capabilities (VCF vs variant list)
# For now, keeping it similar to original, but it might fail if VCF input is removed from CLI for spike.
if [ -f "$ROOT/mini.fa" ] && [ -f "$ROOT/spike_in.bam" ] && [ -f "$ROOT/spike_in.bam.bai" ] && [ -f "$ROOT/variants.vcf" ]; then
  docker run --rm \
    -v "$ROOT/mini.fa:/data/mini.fa" \
    -v "$ROOT/spike_in.bam:/data/spike_in.bam" \
    -v "$ROOT/spike_in.bam.bai:/data/spike_in.bam.bai" \
    -v "$ROOT/variants.vcf:/data/variants.vcf" \
    -v "$ROOT/results_spike_smoke:/data/results_spike" \
    basebuddy spike --reference-fasta /data/mini.fa \
      --input-bams /data/spike_in.bam \
      --variants-vcf /data/variants.vcf \
      --output-prefix-for-bam results_spiked \
      --output-root-dir /data/results_spike \
      --overwrite || echo "spike-in failed (Note: CLI might have changed for variant input)"
else
  echo "SKIP: spike_in.bam, variants.vcf or mini.fa not found for spike-in test"
fi

echo "8️⃣ Testing mutational signature simulation (SigProfilerSimulator)..."
if [ -f "$ROOT/mini.fa" ]; then
  docker run --rm \
    -v "$ROOT/mini.fa:/data/mini.fa" \
    -v "$ROOT/results_sig_smoke:/data/results_sig" \
    basebuddy signature /data/mini.fa \
      --outdir /data/results_sig \
      --sig-type SBS \
      --num-mutations 20 \
      --sample-id testSigSmoke \
      --overwrite || echo "signature sim failed"
else
  echo "SKIP: mini.fa not found for signature sim test"
fi

echo "Smoke test complete."
# BaseBuddy

**Version:** 0.1.0

BaseBuddy is a lightweight toolkit to simulate and manipulate sequencing data:
- **Short‐read simulation** (Illumina/ART)
- **Long‐read simulation** (Nanopore/NanoSim‐h)
- **Spike‐in variants** (BAMSurgeon)
- **Mutational signature simulation** (SigProfilerSimulator)
- **Introduce strand bias** into existing BAM files

---

## Table of Contents

1. [Installation](#installation)  
   1.1 [Python (Conda or venv)](#python‐conda‐or‐venv)  
   1.2 [Docker](#docker)  
2. [Getting a Reference FASTA](#getting‐a‐reference‐fasta)  
3. [Usage Examples](#usage)  
   3.1 [Short‐Read Simulation](#short‐read‐simulation)  
   3.2 [Long‐Read Simulation](#long‐read‐simulation)  
   3.3 [Spike‐In Variants](#spike‐in‐variants)  
   3.4 [Mutational Signature Simulation](#mutational‐signature‐simulation)  
   3.5 [Strand‐Bias Introduction](#strand‐bias‐introduction)  
4. [Edge Cases & Troubleshooting](#edge‐cases‐troubleshooting)  
   4.1 [Missing Dependencies](#missing‐dependencies)  
   4.2 [Reference FASTA Issues](#reference‐fasta‐issues)  
   4.3 [Memory & Disk Considerations](#memory‐disk)  
   4.4 [Docker File‐Sharing Errors (macOS)](#docker‐filesharing‐errors)  
5. [FAQ](#faq)  
6. [License](#license)

---

## Installation

### Python (Conda or `venv`)

1. **Clone the repository** (or download release tarball):
   ```bash
   git clone https://github.com/yourusername/BaseBuddy.git
   cd BaseBuddy

    2.  Create a Conda environment (recommended) or a venv:
Conda

# If you have mamba
mamba env create -f environment.yml
mamba activate basebuddy
pip install -e .

venv (Python 3.10+)

python3 -m venv .venv
source .venv/bin/activate
# Install minimal dependencies (typer). You will need ART, NanoSim‐h, SAMtools, BAMSurgeon, SigProfilerSimulator on your path.
pip install typer
pip install -e .


    3.  Verify installation:

basebuddy version
# Expect: "BaseBuddy 0.1.0"



Docker
    1.  Build the Docker image (from repo root, which contains Dockerfile):

DOCKER_BUILDKIT=1 docker build --progress=plain -t basebuddy:latest .


    2.  Verify:

docker run --rm basebuddy version
# Expect: "BaseBuddy 0.1.0"


    3.  (Optional) You may push it to a registry:

docker tag basebuddy:latest yourrepo/basebuddy:latest
docker push yourrepo/basebuddy:latest



⸻

Getting a Reference FASTA

Many tools (ART, NanoSim‐h, BAMSurgeon, SigProfilerSimulator) require a reference FASTA. You can:
    1.  Download GRCh38 from UCSC or Ensembl and index it locally:

# Example (NCBI FTP)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz
mv GRCh38_latest_genomic.fna grch38.fa

# Create index for SAMtools, BAMSurgeon, etc.
samtools faidx grch38.fa


    2.  Store it in a known directory (e.g. ~/refs/GRCh38/).
By default, BaseBuddy expects you to supply the path to your FASTA.
If no FASTA is found and you have a web connection, you can do:

# (Illustrative only; BaseBuddy does not auto‐download by default)
mkdir -p ~/refs/GRCh38
cd ~/refs/GRCh38
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
ln -s Homo_sapiens.GRCh38.dna.primary_assembly.fa grch38.fa
samtools faidx grch38.fa


    3.  Custom FASTA: If you wish to simulate reads only around a small locus (e.g., FGFR2 region), create a subset FASTA using samtools:

# Extract ±50 kb around FGFR2 (chr10:123,000,000–123,200,000)
samtools faidx grch38.fa chr10:122950000-123250000 > fgfr2_locus.fa
samtools faidx fgfr2_locus.fa

Then use fgfr2_locus.fa in place of reference.

⸻

Usage

Common Options
    •   --outdir (or --out-bam) can be absolute or relative.
    •   All commands check that inputs exist before proceeding.
    •   Make sure tools like art_illumina, nanosim-h, samtools, addsnv.py, and SigProfilerSimulator are on your PATH.

⸻

3.1. Short‐Read Simulation

# Example: simulate 50× coverage, 100 bp reads, HIS.eq profile, output to ./out_short
basebuddy short /path/to/grch38.fa --depth 50 --readlen 100 --profile HS25 --outdir ./out_short

    •   Edge Cases:
    •   If reference.fa is missing or not indexed, you get FileNotFoundError: FASTA reference not found: /path/to/ref.fa.
    •   If art_illumina is not on PATH, you get RuntimeError: art_illumina not found in PATH.
    •   Both --depth and --readlen must be positive integers; otherwise, you see an error and exit early.
    •   Output:
    •   Two FASTQ files: reads1.fq and reads2.fq, within ./out_short/.
    •   The default fragment size is 200±10 bp, with paired‐end flag (-p).

⸻

3.2. Long‐Read Simulation

# Example: simulate 20× coverage Nanopore reads (R9.4.1 model), output to ./out_long
basebuddy long /path/to/grch38.fa --depth 20 --model nanopore_R9.4.1 --outdir ./out_long

    •   Edge Cases:
    •   If nanosim-h is missing, you’ll see RuntimeError: nanosim-h not found in PATH.
    •   If --depth ≤ 0, you see an error and no simulation runs.
    •   If --model is invalid or NanoSim‐h cannot find its model files, it prints NanoSim’s error.
    •   Output:
    •   By default, a folder like ./out_long/longreads/ containing a FASTQ (or multiple FASTQs) depending on NanoSim version.

⸻

3.3. Spike‐In Variants

# Assumptions: 
#  - You have a sorted & indexed BAM: spike_in.bam + spike_in.bai
#  - You have variants.vcf
#  - Reference FASTA (with .fai): grch38.fa + grch38.fa.fai

basebuddy spike /path/to/grch38.fa spike_in.bam variants.vcf \
  --vaf 0.1 \
  --out-bam spiked_10pct.bam \
  --seed 42

    •   Edge Cases:
    •   If spike_in.bam or variants.vcf is missing → FileNotFoundError.
    •   If variants.vcf has no header or malformed entries, BAMSurgeon will error out.
    •   --vaf must be >0 and <1; invalid values cause an early exit.
    •   If the input BAM isn’t indexed, BaseBuddy automatically runs samtools index spike_in.bam.
    •   Output:
    •   spiked_10pct.bam plus its index spiked_10pct.bam.bai.

⸻

3.4. Mutational Signature Simulation

# Standard (built‐in) genome:
basebuddy signature GRCh38 --sig-type SBS --num-mutations 500 --sample-id tumorA --outdir ./out_sig

# Custom FASTA:
basebuddy signature /path/to/fgfr2_locus.fa \
  --sig-type SBS \
  --num-mutations 200 \
  --sample-id FGFR2_test \
  --outdir ./out_sig_fgfr2

    •   Edge Cases:
    •   If reference does not exist → FileNotFoundError.
    •   If SigProfilerSimulator is not installed → RuntimeError: SigProfilerSimulator not installed…
    •   If --sig-type is not in {SBS, DBS, ID}, you see an error about valid types.
    •   If --num-mutations ≤ 0, you see an error.
    •   Output:
    •   A directory structure like ./out_sig/SBS/tumorA/ containing:
    •   tumorA_signatures.vcf (the simulated VCF)
    •   Mutation counts and summary files.

⸻

3.5. Strand‐Bias Introduction

# Assume you just produced `spiked_10pct.bam` and its index
basebuddy strand-bias spiked_10pct.bam \
  --out-bam biased_10pct.bam \
  --forward-fraction 0.8 \
  --seed 123

    •   Edge Cases:
    •   If spiked_10pct.bam or its .bai is missing → FileNotFoundError or auto‐index.
    •   If --forward-fraction < 0 or > 1 → ValueError.
    •   Insufficient disk space to create intermediate BAMs (plus/minus) may cause CalledProcessError.
    •   Output:
    •   biased_10pct.bam + biased_10pct.bai.
    •   Temporary folder strandtemp/ created under the output directory, then removed.

⸻

4. Edge Cases & Troubleshooting

4.1. Missing Dependencies
    •   art_illumina:
    •   Error: RuntimeError: art_illumina not found in PATH.
    •   Solution: Install via Bioconda (conda install -c bioconda art) or from ART GitHub.
    •   nanosim-h:
    •   Error: RuntimeError: nanosim-h not found in PATH.
    •   Solution: Install via Bioconda (conda install -c bioconda nanosim-h) or see NanoSim GitHub.
    •   samtools:
    •   Required by multiple steps. Install via Bioconda (conda install -c bioconda samtools).
    •   addsnv.py (BAMSurgeon):
    •   If missing, install via Bioconda (conda install -c bioconda bamsurgeon) or clone BAMSurgeon repo and add to PATH.
    •   SigProfilerSimulator:
    •   Install via pip install SigProfilerSimulator.
    •   If you see ImportError: cannot import name 'SigProfilerSimulatorFunc', ensure your version matches the code’s expected interface. We imported SigProfilerSimulator directly and provided both input_type and input_data.

4.2. Reference FASTA Issues
    •   No .fai index:
    •   Error: FileNotFoundError: FASTA reference not found: /path/to/ref.fa or “no SQ lines” when running ART.
    •   Solution: Run samtools faidx /path/to/ref.fa before using BaseBuddy.
    •   Malformed FASTA:
    •   If your FASTA has nonstandard headers or missing bases, ART/NanoSim‐h may fail. Validate with:

samtools faidx /path/to/ref.fa



4.3. Memory & Disk Considerations
    •   simulate_long (NanoSim‐h) can produce multi‐gigabyte FASTQ files if depth is high. Monitor disk space:

df -h /your/output/directory


    •   spike_variants: Merging and indexing can require extra disk. Temporary index files (.bai) appear in the same folder as the BAM.

4.4. Docker File‐Sharing Errors (macOS)

If your mounted host directory is not in Docker’s File Sharing list, you see:

docker: Error response from daemon: Mounts denied: 
The path /Users/… is not shared from the host. You can configure shared paths from Docker -> Preferences → Resources → File Sharing.

    •   Solution:
    1.  Open Docker Desktop → Preferences → Resources → File Sharing.
    2.  Add /Users/lauferva/Desktop/Professional/Projects/2025/GOAL/BaseBuddy (or parent) to the list.
    3.  Save & restart Docker.






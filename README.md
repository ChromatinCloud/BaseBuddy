## BaseBuddy

**Version:** 0.1.0

BaseBuddy is a lightweight toolkit to simulate and manipulate sequencing data:
- **Short‐read simulation** (Illumina/ART)
- **Long‐read simulation** (Nanopore/NanoSim‐h)
- **Spike‐in variants** (SNPs and Indels via BAMSurgeon)
- **Mutational signature simulation** (SigProfilerSimulator)
- **Introduce strand bias** into existing BAM files
- **Read Quality Control** (via FastQC integration)

---

## Table of Contents

1. [Installation](#installation)  
   1.1 [Python (Conda or venv)](#python‐conda‐or‐venv)  
   1.2 [Docker](#docker)  
2. [Directory Structure](#directory-structure)
3. [Getting a Reference FASTA](#getting‐a‐reference‐fasta)  
4. [Usage Examples](#usage)  
   3.1 [Short‐Read Simulation](#short‐read‐simulation)  
   3.2 [Long‐Read Simulation](#long‐read‐simulation)  
   3.3 [Spike‐In Variants](#spike‐in‐variants)  
   3.4 [Mutational Signature Simulation](#mutational‐signature‐simulation)  
   3.5 [Strand‐Bias Introduction](#strand‐bias‐introduction)  
   3.6 [Read Quality Control](#read‐quality‐control)
4. [Edge Cases & Troubleshooting](#edge‐cases‐troubleshooting)  
   4.1 [Missing Dependencies](#missing‐dependencies)  
   4.2 [Reference FASTA Issues](#reference‐fasta‐issues)  
   4.3 [Memory & Disk Considerations](#memory‐disk)  
   4.4 [Docker File‐Sharing Errors (macOS)](#docker‐filesharing‐errors)  
5. [FAQ](#faq)  
6. [License](#license)

---

## Directory Structure

```
BaseBuddy/
├── src/              # Source code
├── tests/            # Test suite
├── docs/             # Documentation
│   ├── installation/ # Installation guides
│   ├── parameters/   # Parameter guides (BRAF, KRAS, etc.)
│   └── usage/        # Usage tutorials
├── scripts/          # Utility scripts
│   ├── setup/        # Installation and setup scripts
│   └── fixes/        # Common issue fixes
├── examples/         # Example workflows and demos
├── data/             # Data files
├── logs/             # Log files
└── vignettes/        # Detailed use cases
```

### Quick Start
- **Run GUI**: `python run_gui.py`
- **Run CLI**: `basebuddy --help`
- **View examples**: See `examples/` directory
- **Read docs**: See `docs/` directory

---

## Installation

## Python (Conda or `venv`)

1. **Clone the repository** (or download release tarball):
   ```bash
   git clone https://github.com/yourusername/BaseBuddy.git
   cd BaseBuddy

    2.  Create a Conda environment (recommended) or a venv:
Conda

## If you have mamba
# Note: If environment.yml is missing or fails, please refer to
# docs/dependencies.md for manual dependency installation instructions.
# The full list includes Python packages and external CLI tools.
mamba env create -f environment.yml
mamba activate basebuddy
pip install -e .

venv (Python 3.10+)

python3 -m venv .venv
source .venv/bin/activate
## Install minimal dependencies (typer). You will need ART, NanoSim‐h, SAMtools, BAMSurgeon, SigProfilerSimulator on your path.
pip install typer
pip install -e .


    3.  Verify installation:

basebuddy version
## Expect: "BaseBuddy 0.1.0"



Docker
    1.  Build the Docker image (from repo root, which contains Dockerfile):

DOCKER_BUILDKIT=1 docker build --progress=plain -t basebuddy:latest .


    2.  Verify:

docker run --rm basebuddy version
## Expect: "BaseBuddy 0.1.0"


    3.  (Optional) You may push it to a registry:

docker tag basebuddy:latest yourrepo/basebuddy:latest
docker push yourrepo/basebuddy:latest



⸻

Getting a Reference FASTA

Many tools (ART, NanoSim‐h, BAMSurgeon, SigProfilerSimulator) require a reference FASTA. You can:
    1.  Download GRCh38 from UCSC or Ensembl and index it locally:

## Example (NCBI FTP)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz
mv GRCh38_latest_genomic.fna grch38.fa

## Create index for SAMtools, BAMSurgeon, etc.
samtools faidx grch38.fa


    2.  Store it in a known directory (e.g. ~/refs/GRCh38/).
By default, BaseBuddy expects you to supply the path to your FASTA.
If no FASTA is found and you have a web connection, you can do:

## (Illustrative only; BaseBuddy does not auto‐download by default)
mkdir -p ~/refs/GRCh38
cd ~/refs/GRCh38
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
ln -s Homo_sapiens.GRCh38.dna.primary_assembly.fa grch38.fa
samtools faidx grch38.fa


    3.  Custom FASTA: If you wish to simulate reads only around a small locus (e.g., FGFR2 region), create a subset FASTA using samtools:

## Extract ±50 kb around FGFR2 (chr10:123,000,000–123,200,000)
samtools faidx grch38.fa chr10:122950000-123250000 > fgfr2_locus.fa
samtools faidx fgfr2_locus.fa

Then use fgfr2_locus.fa in place of reference.

---

## Graphical User Interface (GUI)

BaseBuddy also offers a graphical user interface (GUI) for easier access to its features, including:
*   Short and Long Read Simulation
*   Variant Spiking (SNPs and Indels)
*   Germline Variant Simulation (FASTA modification + Read Simulation)
*   Read Quality Control (FastQC)
*   Applying Signatures to FASTA

To run the GUI, ensure all dependencies (including `customtkinter`) are installed, then execute:
```bash
python -m src.basebuddy.gui.main_app
```

---

⸻

Usage

Common Options
    •   --outdir (or --out-bam) can be absolute or relative.
    •   All commands check that inputs exist before proceeding.
    •   Make sure tools like art_illumina, nanosim-h, samtools, addsnv.py, and SigProfilerSimulator are on your PATH.

⸻

3.1. Short‐Read Simulation

## Example: simulate 50× coverage, 100 bp reads, HIS.eq profile, output to ./out_short
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

## Example: simulate 20× coverage Nanopore reads (R9.4.1 model), output to ./out_long
basebuddy long /path/to/grch38.fa --depth 20 --model nanopore_R9.4.1 --outdir ./out_long

    •   Edge Cases:
    •   If nanosim-h is missing, you’ll see RuntimeError: nanosim-h not found in PATH.
    •   If --depth ≤ 0, you see an error and no simulation runs.
    •   If --model is invalid or NanoSim‐h cannot find its model files, it prints NanoSim’s error.
    •   Output:
    •   By default, a folder like ./out_long/longreads/ containing a FASTQ (or multiple FASTQs) depending on NanoSim version.

⸻

3.3. Spike‐In Variants (SNPs and Indels)

The `spike` command allows you to introduce SNVs (Single Nucleotide Variants) and/or Indels (Insertions/Deletions) into one or more existing BAM files using BAMSurgeon.

**Assumptions:**
*   You have one or more sorted and indexed input BAM files.
*   You have a reference FASTA file (e.g., `grch38.fa`) indexed with `samtools faidx`.
*   You have VCF files specifying the SNPs and/or Indels to spike.
*   BAMSurgeon (including `addsnv.py` and `addindel.py`) and its dependencies (like `samtools` and a Picard JAR) are installed and accessible.

**Example:**

To spike SNPs from `snps.vcf` and Indels from `indels.vcf` into `input1.bam` and `input2.bam`:

```bash
basebuddy spike \
  --input-bam input1.bam \
  --input-bam input2.bam \
  --snp-vcf snps.vcf \
  --indel-vcf indels.vcf \
  --reference /path/to/grch38.fa \
  --output-prefix spiked_output \
  --vaf 0.25 \
  --seed 123 \
  --picard-jar /path/to/picard.jar
```

**Key Options:**
*   `--input-bam` / `-i`: Path to an input BAM file. This option can be used multiple times for multiple input BAMs. (Required)
*   `--snp-vcf`: VCF file containing SNPs to spike. (At least one of `--snp-vcf` or `--indel-vcf` is required)
*   `--indel-vcf`: VCF file containing Indels to spike.
*   `--reference` / `-r`: Path to the reference FASTA file.
*   `--output-prefix` / `-p`: Prefix for output files. Final BAMs will be named like `{output_prefix}_{input_bam_stem}_final_sorted.bam` and placed in a run-specific subdirectory within the directory of the prefix (e.g., if prefix is `results/spiked_run`, output is in `results/spike_variants_run_timestamp/...`).
*   `--vaf`: Target Variant Allele Frequency for the spiked variants (default: 0.05).
*   `--seed`: Random seed for reproducibility.
*   `--picard-jar`: Path to `picard.jar`. BAMSurgeon requires this. Alternatively, set the `BAMSURGEON_PICARD_JAR` environment variable.
*   `--overwrite`: Overwrite the output directory if it already exists.

**Edge Cases & Output:**
*   If input BAMs or VCF files are missing, a `FileNotFoundError` or `FileError` will occur.
*   If `addsnv.py`, `addindel.py`, or `samtools` are not found in `$PATH`, a `ConfigurationError` is raised.
*   The command creates a run-specific output directory (e.g., `spike_variants_YYYYMMDD_HHMMSS/`) where all outputs are stored, including:
    *   Final BAM files (one for each input BAM, e.g., `spiked_output_input1_final_sorted.bam`).
    *   Index files (`.bai`) for each output BAM.
    *   Log VCFs from `addsnv.py` and `addindel.py` (if variants were processed).
    *   An IGV session file for each processed BAM.
    *   A main `manifest.json` file summarizing the run.
*   Input BAMs are automatically indexed if `.bai` files are missing (controlled by `--auto-index-bam`).

⸻

3.4. Mutational Signature Simulation

## Standard (built‐in) genome:
basebuddy signature GRCh38 --sig-type SBS --num-mutations 500 --sample-id tumorA --outdir ./out_sig

## Custom FASTA:
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

---

#### Applying Signatures with Bundled Data (API/GUI Functionality)

BaseBuddy also includes an internal capability to apply mutational signatures directly to FASTA files using bundled canonical signature definitions. This is primarily accessed via its Python API (e.g., for GUI operations) and allows for more direct control over the mutation process.

*   **Bundled Signatures**: BaseBuddy aims to pre-package a set of canonical COSMIC v3.3 signatures (GRCh37 context) for Single Base Substitutions (SBS), Doublet Base Substitutions (DBS), and Insertions/Deletions (ID). This means you do not need to download these specific signature sets from external sources like Synapse or COSMIC for this functionality when using these specific internal tools. For example, "SBS1", "SBS5" are available. *(The actual set of fully bundled master signature files like `sbs_grch37_cosmic_v3.3.tsv` is progressively being added to `src/basebuddy/data/signatures/`.)*
*   **Functionality**: This feature (e.g., via the `apply_signature_to_fasta` API function) allows you to take an input FASTA file, select one of these bundled signatures (or provide a path to a custom signature matrix in a compatible TSV format), specify the number of mutations, and generate a new FASTA file with the applied mutations.
*   *(Note: CLI commands for this specific internal signature application pathway will be detailed once available.)*

---

⸻

3.5. Strand‐Bias Introduction

## Assume you just produced `spiked_10pct.bam` and its index
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

3.6 Read Quality Control

BaseBuddy can run FastQC on your FASTQ files to generate standard quality control reports.

**Example:**

To run FastQC on two FASTQ files and save reports to `./qc_results/my_qc_run/`:

```bash
basebuddy qc reads_1.fastq.gz reads_2.fastq.gz --output-dir ./qc_results --run-name my_qc_run
```

**Key Options:**
*   `fastq_files...`: One or more input FASTQ files.
*   `--output-dir` / `-o`: Specify the main directory where QC results will be stored. A run-specific subdirectory will be created inside this.
*   `--run-name`: Optional name for the QC run, used for the subdirectory.
*   `--threads` / `-t`: Number of threads FastQC should use.
*   `--overwrite`: Overwrite the output subdirectory if it exists.

**Output:**
*   A run-specific output directory (e.g., `qc_results/my_qc_run_YYYYMMDD_HHMMSS/`).
*   Inside this, FastQC will create its standard output structure for each input file (e.g., `reads_1_fastqc/fastqc_report.html`).
*   A `manifest_fastqc_run.json` file summarizing the inputs and report locations.

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

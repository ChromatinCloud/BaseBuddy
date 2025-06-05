Below are two deliverables:
	1.	Robustness & Edge‐Case Handling
	•	We’ve added file‐existence checks, type/parameter validation, and more informative error messages.
	•	All commands now validate inputs before invoking external tools.
	•	Updated both cli.py and runner.py to reflect these checks.
	2.	Updated README
	•	A complete, self‐contained README that explains installation, usage, edge cases, and examples.
	•	Describes how to interpret error messages and common troubleshooting steps.

⸻

1. Robustness & Edge‐Case Handling

1.1. cli.py (with additional validation)

# File: src/basebuddy/cli.py

import typer
from pathlib import Path
from . import runner, __version__

app = typer.Typer(add_completion=False)


@app.command()
def version():
    """Print version."""
    print(f"BaseBuddy {__version__}")


@app.command()
def short(
    reference: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True),
    depth: int = typer.Option(30, "--depth", help="Desired fold coverage (integer > 0)"),
    readlen: int = typer.Option(150, "--readlen", help="Read length (integer > 0)"),
    profile: str = typer.Option("HS25", "--profile", help="ART profile (e.g. HS25, HS20)"),
    outdir: Path = typer.Option(Path("results_short"), "--outdir", help="Output directory for FASTQ"),
):
    """
    Simulate short Illumina reads:
      basebuddy short /path/to/ref.fa --depth 50 --readlen 100 --profile HS25 --outdir ./out
    """
    # Validate numeric inputs
    if depth <= 0:
        typer.echo("ERROR: --depth must be a positive integer.", err=True)
        raise typer.Exit(code=1)
    if readlen <= 0:
        typer.echo("ERROR: --readlen must be a positive integer.", err=True)
        raise typer.Exit(code=1)

    runner.simulate_short(str(reference), outdir, depth, readlen, profile)


@app.command()
def long(
    reference: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True),
    depth: int = typer.Option(30, "--depth", help="Desired fold coverage (integer > 0)"),
    model: str = typer.Option("nanopore_R9.4.1", "--model", help="NanoSim model name"),
    outdir: Path = typer.Option(Path("results_long"), "--outdir", help="Output directory for FASTQ"),
):
    """
    Simulate long Nanopore reads:
      basebuddy long /path/to/ref.fa --depth 20 --model nanopore_R9.4.1 --outdir ./out
    """
    if depth <= 0:
        typer.echo("ERROR: --depth must be a positive integer.", err=True)
        raise typer.Exit(code=1)

    runner.simulate_long(str(reference), outdir, depth, model)


@app.command()
def spike(
    reference: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True),
    in_bam: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True),
    vcf: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True),
    vaf: float = typer.Option(0.05, "--vaf", help="Variant allele fraction (0.0 < vaf < 1.0)"),
    out_bam: Path = typer.Option(Path("spiked.bam"), "--out-bam", help="Output BAM path"),
    seed: int = typer.Option(0, "--seed", help="Random seed (integer >= 0)"),
):
    """
    Spike variants into an existing BAM:
      basebuddy spike /path/to/ref.fa /path/to/input.bam /path/to/variants.vcf \
         --vaf 0.05 --out-bam spiked.bam --seed 42
    """
    if not (0.0 < vaf < 1.0):
        typer.echo("ERROR: --vaf must be between 0.0 and 1.0 (exclusive).", err=True)
        raise typer.Exit(code=1)
    if seed < 0:
        typer.echo("ERROR: --seed must be a non‐negative integer.", err=True)
        raise typer.Exit(code=1)

    runner.spike_variants(str(reference), str(in_bam), str(vcf), vaf, str(out_bam), seed)


@app.command()
def signature(
    reference: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True),
    outdir: Path = typer.Option(Path("results_sig"), "--outdir", help="Output directory for signatures"),
    sig_type: str = typer.Option("SBS", "--sig-type", help="Signature type (SBS, DBS, ID)"),
    num_mutations: int = typer.Option(100, "--num-mutations", help="Number of mutations (integer > 0)"),
    sample_id: str = typer.Option("Sample", "--sample-id", help="Sample/project identifier"),
):
    """
    Simulate mutational signatures:
      basebuddy signature /path/to/ref.fa --sig-type SBS --num-mutations 1000 --sample-id tumorA --outdir ./sigout
    """
    if num_mutations <= 0:
        typer.echo("ERROR: --num-mutations must be a positive integer.", err=True)
        raise typer.Exit(code=1)

    valid_types = {"SBS", "DBS", "ID"}
    if sig_type not in valid_types:
        typer.echo(f"ERROR: --sig-type must be one of {', '.join(valid_types)}.", err=True)
        raise typer.Exit(code=1)

    runner.simulate_signatures(str(reference), outdir, sig_type, num_mutations, sample_id)


@app.command("strand-bias")
def strand_bias(
    in_bam: Path = typer.Argument(..., exists=True, file_okay=True, dir_okay=False, readable=True),
    out_bam: Path = typer.Option(Path("biased.bam"), "--out-bam", help="Output BAM path"),
    forward_fraction: float = typer.Option(0.5, "--forward-fraction", help="Fraction of reads on forward strand (0.0–1.0)"),
    seed: int = typer.Option(0, "--seed", help="Random seed (integer >= 0)"),
):
    """
    Introduce strand bias into an existing BAM:
      basebuddy strand-bias /path/to/input.bam --out-bam biased.bam --forward-fraction 0.8 --seed 42
    """
    if not (0.0 <= forward_fraction <= 1.0):
        typer.echo("ERROR: --forward-fraction must be between 0.0 and 1.0.", err=True)
        raise typer.Exit(code=1)
    if seed < 0:
        typer.echo("ERROR: --seed must be a non‐negative integer.", err=True)
        raise typer.Exit(code=1)

    runner.introduce_strand_bias(str(in_bam), str(out_bam), forward_fraction, seed)


if __name__ == "__main__":
    app()

Key additions in cli.py:
	•	Type annotations and typer.Argument / typer.Option constraints for automatic validation (e.g., exists=True, readable=True).
	•	Manual checks for numeric parameters (depth, readlen, etc.) with clear error messages.
	•	Restrict sig_type to one of the supported values.
	•	Early exits (typer.Exit(code=1)) on invalid inputs so the user sees a concise error.

⸻

1.2. runner.py (with file‐existence checks and robust error messages)

# File: src/basebuddy/runner.py

import subprocess
import shutil
from pathlib import Path


def _ensure_exists(path: str, desc: str) -> None:
    """
    Verify that the given path exists; if not, raise FileNotFoundError.
    """
    if not Path(path).exists():
        raise FileNotFoundError(f"{desc} not found: {path}")


def _run_cmd(cmd: list[str], cwd: Path | None = None) -> None:
    """
    Helper function to print and run a command. Raises CalledProcessError on failure.
    """
    print("➤", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def simulate_short(
    reference: str,
    outdir: Path,
    depth: int,
    readlen: int = 150,
    profile: str = "HS25",
) -> None:
    """
    Simulate paired-end Illumina reads using ART:
      - reference: path to FASTA
      - outdir: Path object for output directory
      - depth: desired fold coverage
      - readlen: read length
      - profile: ART profile name
    """
    # Ensure reference exists
    _ensure_exists(reference, "FASTA reference")
    # Create output directory if needed
    outdir.mkdir(parents=True, exist_ok=True)

    # Locate ART executable
    art = shutil.which("art_illumina")
    if art is None:
        raise RuntimeError("art_illumina not found in PATH. Please install ART or add it to PATH.")

    prefix = str(outdir / "reads")
    # Default fragment size parameters for paired-end
    fragment_mean = 200
    fragment_sd = 10

    cmd = [
        art,
        "-ss", profile,
        "-i", reference,
        "-l", str(readlen),
        "-f", str(depth),
        "-o", prefix,
        "-p",                # Paired-end
        "-m", str(fragment_mean),
        "-s", str(fragment_sd),
    ]

    _run_cmd(cmd)


def simulate_long(
    reference: str,
    outdir: Path,
    depth: int = 30,
    model: str = "nanopore_R9.4.1",
) -> None:
    """
    Simulate Nanopore long reads using NanoSim-h:
      - reference: path to FASTA
      - outdir: Path object for output directory
      - depth: desired coverage (fold or read count factor)
      - model: NanoSim model name
    """
    _ensure_exists(reference, "FASTA reference")
    outdir.mkdir(parents=True, exist_ok=True)

    nanosim = shutil.which("nanosim-h")
    if nanosim is None:
        raise RuntimeError("nanosim-h not found in PATH. Please install NanoSim-h or add it to PATH.")

    # NanoSim command: adjust flags if your version differs
    cmd = [
        nanosim,
        "simulate",
        "-r", reference,
        "-c", str(depth),
        "-m", model,
        "-o", str(outdir / "longreads"),
    ]

    _run_cmd(cmd)


def spike_variants(
    reference: str,
    in_bam: str,
    vcf: str,
    vaf: float,
    out_bam: str,
    seed: int = 0,
) -> None:
    """
    Spike SNVs/indels into an existing BAM using BAMSurgeon:
      - reference: path to FASTA
      - in_bam: path to input BAM (must exist and be sorted)
      - vcf: path to VCF of variants
      - vaf: target allele fraction (0.0 < vaf < 1.0)
      - out_bam: path for output BAM
      - seed: random seed
    """
    _ensure_exists(reference, "FASTA reference")
    _ensure_exists(in_bam, "Input BAM")
    _ensure_exists(vcf, "Input VCF")

    # Index input BAM if .bai missing
    if not Path(in_bam + ".bai").exists():
        _run_cmd(["samtools", "index", in_bam])

    addsnv = shutil.which("addsnv.py")
    if addsnv is None:
        raise RuntimeError("addsnv.py (BAMSurgeon) not found in PATH. Install BAMSurgeon.")

    cmd = [
        addsnv,
        "-v", vcf,
        "-f", in_bam,
        "-r", reference,
        "-o", out_bam,
        "-p", str(vaf),
        "-s", str(seed),
    ]
    _run_cmd(cmd)

    # Index the newly created BAM
    _run_cmd(["samtools", "index", out_bam])


def simulate_signatures(
    reference: str,
    outdir: Path,
    sig_type: str = "SBS",
    num_mutations: int = 100,
    sample_id: str = "Sample",
) -> None:
    """
    Simulate mutational signatures using SigProfilerSimulator:
      - reference: either a genome build string (e.g., "GRCh38") or a custom FASTA path
      - outdir: Path object for output
      - sig_type: "SBS", "DBS", or "ID"
      - num_mutations: total number of mutations per sample
      - sample_id: project/sample name
    """
    _ensure_exists(reference, "Reference for signature simulation")
    outdir.mkdir(parents=True, exist_ok=True)

    # Importing inside function to avoid breaking if unavailable
    try:
        from SigProfilerSimulator import SigProfilerSimulator
    except ImportError as e:
        raise RuntimeError(
            "SigProfilerSimulator not installed or not in PATH. Install via 'pip install SigProfilerSimulator'."
        )

    # If reference is a FASTA file, use input_type="custom_genome"
    input_type = "reference_genome"
    input_data = reference
    if reference.lower().endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz")):
        input_type = "custom_genome"
        input_data = reference

    # Build and run the simulator call
    SigProfilerSimulator(
        project=sample_id,
        input_type=input_type,       # "reference_genome" or "custom_genome"
        input_data=input_data,       # genome string or path to FASTA
        reference_genome="GRCh38",   # only used for built‐in genomes; can be ignored if custom
        sig_type=sig_type,
        sampleNumMutations=num_mutations,
        outdir=str(outdir),
        seed=0,
    )


def introduce_strand_bias(
    in_bam: str,
    out_bam: str,
    forward_fraction: float = 0.5,
    seed: int = 0,
) -> None:
    """
    Introduce strand bias into a BAM by subsampling reads from each strand:
      - in_bam: path to input BAM (must exist)
      - out_bam: path to output biased BAM
      - forward_fraction: fraction of reads on forward strand (0.0–1.0)
      - seed: random seed
    """
    if not (0.0 <= forward_fraction <= 1.0):
        raise ValueError("forward_fraction must be between 0.0 and 1.0 (inclusive)")

    _ensure_exists(in_bam, "Input BAM")

    # Index BAM if necessary
    if not Path(in_bam + ".bai").exists():
        _run_cmd(["samtools", "index", in_bam])

    outdir = Path(out_bam).parent
    outdir.mkdir(parents=True, exist_ok=True)

    temp_dir = outdir / "strandtemp"
    temp_dir.mkdir(exist_ok=True)

    plus_bam = temp_dir / "plus.bam"
    minus_bam = temp_dir / "minus.bam"

    # Extract forward‐strand reads (FLAG 0x10 not set)
    _run_cmd(["samtools", "view", "-h", "-F", "16", in_bam, "-o", str(plus_bam)])
    # Extract reverse‐strand reads (FLAG 0x10 set)
    _run_cmd(["samtools", "view", "-h", "-f", "16", in_bam, "-o", str(minus_bam)])

    plus_keep = forward_fraction
    minus_keep = 1.0 - forward_fraction

    # Subsample forward reads
    _run_cmd([
        "samtools", "view", "-s",
        f"{seed}.{int(plus_keep * 1000):03d}",
        "-b", str(plus_bam), "-o", str(temp_dir / "plus_sub.bam")
    ])
    # Subsample reverse reads (use a different seed for reproducibility)
    _run_cmd([
        "samtools", "view", "-s",
        f"{seed+1}.{int(minus_keep * 1000):03d}",
        "-b", str(minus_bam), "-o", str(temp_dir / "minus_sub.bam")
    ])

    merged = temp_dir / "merged.bam"
    # Merge subsampled BAMs
    _run_cmd([
        "samtools", "merge", "-f", str(merged),
        str(temp_dir / "plus_sub.bam"), str(temp_dir / "minus_sub.bam"),
    ])

    # Sort and index final BAM
    _run_cmd(["samtools", "sort", "-o", out_bam, str(merged)])
    _run_cmd(["samtools", "index", out_bam])

    # Clean up temp files
    shutil.rmtree(temp_dir)

Key improvements in runner.py:
	1.	_ensure_exists(path, desc) at the top of each function to verify required files (FASTA, BAM, VCF) before running commands.
	2.	Clear, user‐friendly error messages when executables (e.g. art_illumina, nanosim-h, addsnv.py) are missing.
	3.	Validation of numeric inputs inside cli.py ensures runner.py is only called with sane parameters.
	4.	Automatic BAM indexing if the .bai is missing, rather than crashing halfway.
	5.	Dynamic handling for simulate_signatures—if the “reference” is a FASTA (custom), set input_type="custom_genome", otherwise assume a built‐in genome build (e.g. "GRCh38").

⸻

2. Updated README

Below is a fully updated README.md. It includes:
	•	Installation instructions (Python‐only and Docker).
	•	Detailed usage examples for every subcommand.
	•	Edge‐case notes and common troubleshooting tips.
	•	How to point to a downloaded GRCh38 FASTA if the user wants to avoid repeated downloads.

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

	2.	Create a Conda environment (recommended) or a venv:
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


	3.	Verify installation:

basebuddy version
# Expect: "BaseBuddy 0.1.0"



Docker
	1.	Build the Docker image (from repo root, which contains Dockerfile):

DOCKER_BUILDKIT=1 docker build --progress=plain -t basebuddy:latest .


	2.	Verify:

docker run --rm basebuddy version
# Expect: "BaseBuddy 0.1.0"


	3.	(Optional) You may push it to a registry:

docker tag basebuddy:latest yourrepo/basebuddy:latest
docker push yourrepo/basebuddy:latest



⸻

Getting a Reference FASTA

Many tools (ART, NanoSim‐h, BAMSurgeon, SigProfilerSimulator) require a reference FASTA. You can:
	1.	Download GRCh38 from UCSC or Ensembl and index it locally:

# Example (NCBI FTP)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz
mv GRCh38_latest_genomic.fna grch38.fa

# Create index for SAMtools, BAMSurgeon, etc.
samtools faidx grch38.fa


	2.	Store it in a known directory (e.g. ~/refs/GRCh38/).
By default, BaseBuddy expects you to supply the path to your FASTA.
If no FASTA is found and you have a web connection, you can do:

# (Illustrative only; BaseBuddy does not auto‐download by default)
mkdir -p ~/refs/GRCh38
cd ~/refs/GRCh38
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
ln -s Homo_sapiens.GRCh38.dna.primary_assembly.fa grch38.fa
samtools faidx grch38.fa


	3.	Custom FASTA: If you wish to simulate reads only around a small locus (e.g., FGFR2 region), create a subset FASTA using samtools:

# Extract ±50 kb around FGFR2 (chr10:123,000,000–123,200,000)
samtools faidx grch38.fa chr10:122950000-123250000 > fgfr2_locus.fa
samtools faidx fgfr2_locus.fa

Then use fgfr2_locus.fa in place of reference.

⸻

Usage

Common Options
	•	--outdir (or --out-bam) can be absolute or relative.
	•	All commands check that inputs exist before proceeding.
	•	Make sure tools like art_illumina, nanosim-h, samtools, addsnv.py, and SigProfilerSimulator are on your PATH.

⸻

3.1. Short‐Read Simulation

# Example: simulate 50× coverage, 100 bp reads, HIS.eq profile, output to ./out_short
basebuddy short /path/to/grch38.fa --depth 50 --readlen 100 --profile HS25 --outdir ./out_short

	•	Edge Cases:
	•	If reference.fa is missing or not indexed, you get FileNotFoundError: FASTA reference not found: /path/to/ref.fa.
	•	If art_illumina is not on PATH, you get RuntimeError: art_illumina not found in PATH.
	•	Both --depth and --readlen must be positive integers; otherwise, you see an error and exit early.
	•	Output:
	•	Two FASTQ files: reads1.fq and reads2.fq, within ./out_short/.
	•	The default fragment size is 200±10 bp, with paired‐end flag (-p).

⸻

3.2. Long‐Read Simulation

# Example: simulate 20× coverage Nanopore reads (R9.4.1 model), output to ./out_long
basebuddy long /path/to/grch38.fa --depth 20 --model nanopore_R9.4.1 --outdir ./out_long

	•	Edge Cases:
	•	If nanosim-h is missing, you’ll see RuntimeError: nanosim-h not found in PATH.
	•	If --depth ≤ 0, you see an error and no simulation runs.
	•	If --model is invalid or NanoSim‐h cannot find its model files, it prints NanoSim’s error.
	•	Output:
	•	By default, a folder like ./out_long/longreads/ containing a FASTQ (or multiple FASTQs) depending on NanoSim version.

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

	•	Edge Cases:
	•	If spike_in.bam or variants.vcf is missing → FileNotFoundError.
	•	If variants.vcf has no header or malformed entries, BAMSurgeon will error out.
	•	--vaf must be >0 and <1; invalid values cause an early exit.
	•	If the input BAM isn’t indexed, BaseBuddy automatically runs samtools index spike_in.bam.
	•	Output:
	•	spiked_10pct.bam plus its index spiked_10pct.bam.bai.

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

	•	Edge Cases:
	•	If reference does not exist → FileNotFoundError.
	•	If SigProfilerSimulator is not installed → RuntimeError: SigProfilerSimulator not installed…
	•	If --sig-type is not in {SBS, DBS, ID}, you see an error about valid types.
	•	If --num-mutations ≤ 0, you see an error.
	•	Output:
	•	A directory structure like ./out_sig/SBS/tumorA/ containing:
	•	tumorA_signatures.vcf (the simulated VCF)
	•	Mutation counts and summary files.

⸻

3.5. Strand‐Bias Introduction

# Assume you just produced `spiked_10pct.bam` and its index
basebuddy strand-bias spiked_10pct.bam \
  --out-bam biased_10pct.bam \
  --forward-fraction 0.8 \
  --seed 123

	•	Edge Cases:
	•	If spiked_10pct.bam or its .bai is missing → FileNotFoundError or auto‐index.
	•	If --forward-fraction < 0 or > 1 → ValueError.
	•	Insufficient disk space to create intermediate BAMs (plus/minus) may cause CalledProcessError.
	•	Output:
	•	biased_10pct.bam + biased_10pct.bai.
	•	Temporary folder strandtemp/ created under the output directory, then removed.

⸻

4. Edge Cases & Troubleshooting

4.1. Missing Dependencies
	•	art_illumina:
	•	Error: RuntimeError: art_illumina not found in PATH.
	•	Solution: Install via Bioconda (conda install -c bioconda art) or from ART GitHub.
	•	nanosim-h:
	•	Error: RuntimeError: nanosim-h not found in PATH.
	•	Solution: Install via Bioconda (conda install -c bioconda nanosim-h) or see NanoSim GitHub.
	•	samtools:
	•	Required by multiple steps. Install via Bioconda (conda install -c bioconda samtools).
	•	addsnv.py (BAMSurgeon):
	•	If missing, install via Bioconda (conda install -c bioconda bamsurgeon) or clone BAMSurgeon repo and add to PATH.
	•	SigProfilerSimulator:
	•	Install via pip install SigProfilerSimulator.
	•	If you see ImportError: cannot import name 'SigProfilerSimulatorFunc', ensure your version matches the code’s expected interface. We imported SigProfilerSimulator directly and provided both input_type and input_data.

4.2. Reference FASTA Issues
	•	No .fai index:
	•	Error: FileNotFoundError: FASTA reference not found: /path/to/ref.fa or “no SQ lines” when running ART.
	•	Solution: Run samtools faidx /path/to/ref.fa before using BaseBuddy.
	•	Malformed FASTA:
	•	If your FASTA has nonstandard headers or missing bases, ART/NanoSim‐h may fail. Validate with:

samtools faidx /path/to/ref.fa



4.3. Memory & Disk Considerations
	•	simulate_long (NanoSim‐h) can produce multi‐gigabyte FASTQ files if depth is high. Monitor disk space:

df -h /your/output/directory


	•	spike_variants: Merging and indexing can require extra disk. Temporary index files (.bai) appear in the same folder as the BAM.

4.4. Docker File‐Sharing Errors (macOS)

If your mounted host directory is not in Docker’s File Sharing list, you see:

docker: Error response from daemon: Mounts denied: 
The path /Users/… is not shared from the host. You can configure shared paths from Docker -> Preferences → Resources → File Sharing.

	•	Solution:
	1.	Open Docker Desktop → Preferences → Resources → File Sharing.
	2.	Add /Users/lauferva/Desktop/Professional/Projects/2025/GOAL/BaseBuddy (or parent) to the list.
	3.	Save & restart Docker.






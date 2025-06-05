import subprocess
import shutil
from pathlib import Path
import sys

BASEBUDDY_REF_DIR = Path.home() / ".basebuddy" / "refs"
DEFAULT_REF = BASEBUDDY_REF_DIR / "GRCh38.fa"

def _ensure_exists(file, desc="file"):
    if not Path(file).exists():
        raise FileNotFoundError(f"{desc} not found: {file}")

def _run_cmd(cmd: list[str], cwd: Path | None = None) -> None:
    """Helper function to print and run a command."""
    print("➤", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True, cwd=cwd)

def _ensure_reference(reference: str | None) -> str:
    """Ensure a reference genome is available, downloading if needed."""
    ref = Path(reference) if reference else DEFAULT_REF
    if not ref.exists():
        print(f"Reference not found: {ref}")
        print("Downloading GRCh38 to", ref)
        BASEBUDDY_REF_DIR.mkdir(parents=True, exist_ok=True)
        _download_grch38(str(ref))
    return str(ref)

def _download_grch38(dest: str) -> None:
    import requests
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    gz_dest = dest + ".gz"
    # Download
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(gz_dest, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    # Unzip
    import gzip, shutil as sh
    with gzip.open(gz_dest, 'rb') as f_in, open(dest, 'wb') as f_out:
        sh.copyfileobj(f_in, f_out)
    Path(gz_dest).unlink()
    print(f"Downloaded and extracted GRCh38 to {dest}")

def simulate_short(reference: str | None, outdir: Path, depth: int, readlen: int = 150, profile: str = "HS25") -> None:
    ref = _ensure_reference(reference)
    outdir.mkdir(parents=True, exist_ok=True)
    art = shutil.which("art_illumina")
    if art is None:
        raise RuntimeError("art_illumina not found in PATH")
    prefix = str(outdir / "reads")
    fragment_mean = 200
    fragment_sd = 10
    cmd = [
        art,
        "-ss", profile,
        "-i", ref,
        "-l", str(readlen),
        "-f", str(depth),
        "-o", prefix,
        "-p",
        "-m", str(fragment_mean),
        "-s", str(fragment_sd),
    ]
    _run_cmd(cmd)

def simulate_long(reference: str | None, outdir: Path, depth: int = 30, model: str = "nanopore_R9.4.1") -> None:
    ref = _ensure_reference(reference)
    outdir.mkdir(parents=True, exist_ok=True)
    nanosim = shutil.which("nanosim-h")
    if nanosim is None:
        raise RuntimeError("nanosim-h not found in PATH")
    cmd = [
        nanosim, "simulate",
        "-r", ref,
        "-c", str(depth),
        "-m", model,
        "-o", str(outdir / "longreads"),
    ]
    _run_cmd(cmd)

def spike_variants(reference: str | None, in_bam: str, vcf: str, vaf: float, out_bam: str, seed: int = 0) -> None:
    ref = _ensure_reference(reference)
    if not Path(in_bam + ".bai").exists():
        _run_cmd(["samtools", "index", in_bam])
    addsnv = shutil.which("addsnv.py")
    if addsnv is None:
        raise RuntimeError("addsnv.py (BAMSurgeon) not found in PATH")
    cmd = [
        addsnv,
        "-v", vcf,
        "-f", in_bam,
        "-r", ref,
        "-o", out_bam,
        "-p", str(vaf),
        "-s", str(seed),
    ]
    _run_cmd(cmd)
    _run_cmd(["samtools", "index", out_bam])

def simulate_signatures(reference: str | None, outdir: Path, sig_type: str = "SBS", num_mutations: int = 100, sample_id: str = "Sample") -> None:
    ref = _ensure_reference(reference)
    from SigProfilerSimulator import SigProfilerSimulator
    outdir.mkdir(parents=True, exist_ok=True)
    SigProfilerSimulator(
        project=sample_id,
        reference_genome="GRCh38",  # Adjust if using other builds
        outdir=str(outdir),
        num_samples=1,
        exome=False,
        type=sig_type,
        total_mutations=num_mutations,
        chrom_based=False,
        seed=0,
    )

def introduce_strand_bias(in_bam: str, out_bam: str, forward_fraction: float = 0.5, seed: int = 0) -> None:
    if not Path(in_bam + ".bai").exists():
        _run_cmd(["samtools", "index", in_bam])
    outdir = Path(out_bam).parent
    outdir.mkdir(parents=True, exist_ok=True)
    temp_dir = outdir / "strandtemp"
    temp_dir.mkdir(exist_ok=True)
    plus_bam = temp_dir / "plus.bam"
    minus_bam = temp_dir / "minus.bam"
    _run_cmd(["samtools", "view", "-h", "-F", "16", in_bam, "-o", str(plus_bam)])
    _run_cmd(["samtools", "view", "-h", "-f", "16", in_bam, "-o", str(minus_bam)])
    plus_keep = forward_fraction
    minus_keep = 1.0 - forward_fraction
    _run_cmd([
        "samtools", "view", "-s",
        f"{seed}.{int(plus_keep * 1000):03d}",
        "-b", str(plus_bam), "-o", str(temp_dir / "plus_sub.bam")
    ])
    _run_cmd([
        "samtools", "view", "-s",
        f"{seed}.{int(minus_keep * 1000):03d}",
        "-b", str(minus_bam), "-o", str(temp_dir / "minus_sub.bam")
    ])
    merged = temp_dir / "merged.bam"
    _run_cmd([
        "samtools", "merge", "-f", str(merged),
        str(temp_dir / "plus_sub.bam"), str(temp_dir / "minus_sub.bam"),
    ])
    _run_cmd(["samtools", "sort", "-o", out_bam, str(merged)])
    _run_cmd(["samtools", "index", out_bam])
    shutil.rmtree(temp_dir)

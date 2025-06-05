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
    reference: Path = None,
    depth: int = 30,
    readlen: int = 150,
    profile: str = "HS25",
    outdir: Path = Path("results_short"),
):
    """
    Simulate short Illumina reads:
      basebuddy short --depth 50 --readlen 100 --profile HS25 --outdir ./out
    (Reference genome optional; will use cached GRCh38 if not specified)
    """
    runner.simulate_short(str(reference) if reference else None, outdir, depth, readlen, profile)

@app.command()
def long(
    reference: Path = None,
    depth: int = 30,
    model: str = "nanopore_R9.4.1",
    outdir: Path = Path("results_long"),
):
    """
    Simulate long Nanopore reads:
      basebuddy long --depth 20 --model nanopore_R9.4.1 --outdir ./out
    (Reference genome optional; will use cached GRCh38 if not specified)
    """
    runner.simulate_long(str(reference) if reference else None, outdir, depth, model)

@app.command()
def spike(
    reference: Path = None,
    in_bam: Path = None,
    vcf: Path = None,
    vaf: float = 0.05,
    out_bam: Path = Path("spiked.bam"),
    seed: int = 0,
):
    """
    Spike variants into an existing BAM:
      basebuddy spike --in-bam /path/to/input.bam --vcf /path/to/variants.vcf --vaf 0.05 --out-bam spiked.bam --seed 42
    (Reference genome optional; will use cached GRCh38 if not specified)
    """
    runner.spike_variants(
        str(reference) if reference else None,
        str(in_bam),
        str(vcf),
        vaf,
        str(out_bam),
        seed,
    )

@app.command()
def signature(
    reference: Path = None,
    outdir: Path = Path("results_sig"),
    sig_type: str = "SBS",
    num_mutations: int = 100,
    sample_id: str = "Sample",
):
    """
    Simulate mutational signatures:
      basebuddy signature --sig-type SBS --num-mutations 1000 --sample-id tumorA --outdir ./sigout
    (Reference genome optional; will use cached GRCh38 if not specified)
    """
    runner.simulate_signatures(
        str(reference) if reference else None,
        outdir,
        sig_type,
        num_mutations,
        sample_id,
    )

@app.command("strand-bias")
def strand_bias(
    in_bam: Path,
    out_bam: Path = Path("biased.bam"),
    forward_fraction: float = 0.5,
    seed: int = 0,
):
    """
    Introduce strand bias into an existing BAM:
      basebuddy strand-bias /path/to/input.bam --out-bam biased.bam --forward-fraction 0.8 --seed 42
    """
    runner.introduce_strand_bias(str(in_bam), str(out_bam), forward_fraction, seed)

if __name__ == "__main__":
    app()

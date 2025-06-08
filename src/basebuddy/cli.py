import typer
from pathlib import Path
from . import runner, __version__
from basebuddy import utils as bb_utils # Add if not present, or adjust alias

app = typer.Typer(add_completion=False)

@app.command()
def version():
    """Print version."""
    print(f"BaseBuddy {__version__}")

@app.command()
def short(
    reference: Path = typer.Option(None, help="Path to the reference FASTA file. Uses cached GRCh38 if not specified."),
    outdir: Path = typer.Option(Path("results_short"), "--outdir", "-o", help="Output directory for simulation results."),
    depth: int = typer.Option(30, "-d", "--depth", help="Sequencing depth."),
    readlen: int = typer.Option(150, "-l", "--read-len", help="Length of the reads."),
    profile: str = typer.Option("HS25", "-p", "--profile", help="ART sequencing profile (e.g., HS25, MSv3)."),
    mean_fragment_length: int = typer.Option(200, "--mean-frag-len", help="Mean fragment length for paired-end reads."),
    std_dev_fragment_length: int = typer.Option(10, "--std-frag-len", help="Standard deviation of fragment length for paired-end reads."),
    is_paired_end: bool = typer.Option(True, "--paired/--single", help="Simulate paired-end or single-end reads."),
    art_platform: str = typer.Option("illumina", "--platform", help="ART sequencing platform (e.g., illumina, 454)."),
    overwrite_output: bool = typer.Option(False, "--overwrite", help="Overwrite output directory if it exists."),
    auto_index_fasta: bool = typer.Option(True, "--auto-index/--no-auto-index", help="Automatically index FASTA if .fai is missing.")
):
    """
    Simulate short reads using ART (e.g., Illumina).
    Example: basebuddy short --reference ref.fa --outdir ./out_short --depth 50 --readlen 100 --profile HS25
    """
    cli_params = {
        "reference_fasta": str(reference) if reference else None,
        "output_root_dir": str(outdir),
        "depth": depth,
        "read_length": readlen,
        "art_profile": profile,
        "mean_fragment_length": mean_fragment_length,
        "std_dev_fragment_length": std_dev_fragment_length,
        "is_paired_end": is_paired_end,
        "art_platform": art_platform,
        "overwrite_output": overwrite_output,
        "auto_index_fasta": auto_index_fasta,
    }

    try:
        results = runner.simulate_short(
            output_root_dir=outdir,
            reference_fasta=str(reference) if reference else None,
            depth=depth,
            read_length=readlen, # Runner uses read_length
            art_profile=profile,
            command_params=cli_params,
            mean_fragment_length=mean_fragment_length,
            std_dev_fragment_length=std_dev_fragment_length,
            is_paired_end=is_paired_end,
            overwrite_output=overwrite_output,
            art_platform=art_platform,
            auto_index_fasta=auto_index_fasta
        )
        typer.secho(f"Short read simulation completed successfully for run: {results.get('run_name')}", fg=typer.colors.GREEN)
        typer.echo(f"Output directory: {results.get('output_directory')}")
        if results.get('manifest_path'):
            typer.echo(f"Manifest file: {results.get('manifest_path')}")

    except bb_utils.BaseBuddyConfigError as e:
        typer.secho(f"Configuration Error: {e}", fg=typer.colors.RED, err=True)
        if e.details:
            typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except bb_utils.BaseBuddyToolError as e:
        typer.secho(f"Tool Execution Error: {e}", fg=typer.colors.RED, err=True)
        if e.details:
            typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except bb_utils.BaseBuddyFileError as e:
        typer.secho(f"File Error: {e}", fg=typer.colors.RED, err=True)
        if e.details:
            typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.secho(f"An unexpected error occurred: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

@app.command()
def long(
    reference: Path = typer.Option(None, help="Path to the reference FASTA file. Uses cached GRCh38 if not specified."),
    outdir: Path = typer.Option(Path("results_long"), "--outdir", "-o", help="Output directory for simulation results."),
    depth: int = typer.Option(30, "-d", "--depth", help="Sequencing depth (approximate coverage)."),
    model: str = typer.Option("nanopore_R9.4.1", "-m", "--model", help="NanoSim-h model profile (e.g., nanopore_R9.4.1, pacbio_sequel)."),
    num_reads: Optional[int] = typer.Option(None, "-N", "--num-reads", help="Exact number of reads to simulate (overrides depth if specified)."),
    overwrite_output: bool = typer.Option(False, "--overwrite", help="Overwrite output directory if it exists."),
    auto_index_fasta: bool = typer.Option(True, "--auto-index/--no-auto-index", help="Automatically index FASTA if .fai is missing."),
):
    """
    Simulate long reads using NanoSim-h (e.g., Nanopore, PacBio).
    Example: basebuddy long --reference ref.fa --outdir ./out_long --depth 20 --model nanopore_R9.4.1
    """
    cli_params = {
        "reference_fasta": str(reference) if reference else None,
        "output_root_dir": str(outdir),
        "depth": depth,
        "model": model,
        "num_reads": num_reads, # This will be passed to runner via command_params
        "overwrite_output": overwrite_output,
        "auto_index_fasta": auto_index_fasta,
    }

    try:
        results = runner.simulate_long(
            output_root_dir=outdir,
            reference_fasta=str(reference) if reference else None,
            depth=depth,
            model=model,
            command_params=cli_params,
            overwrite_output=overwrite_output,
            auto_index_fasta=auto_index_fasta
        )
        typer.secho(f"Long read simulation completed successfully for run: {results.get('run_name')}", fg=typer.colors.GREEN)
        typer.echo(f"Output directory: {results.get('output_directory')}")
        if results.get('manifest_path'):
            typer.echo(f"Manifest file: {results.get('manifest_path')}")

    except bb_utils.BaseBuddyConfigError as e:
        typer.secho(f"Configuration Error: {e}", fg=typer.colors.RED, err=True)
        if e.details:
            typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except bb_utils.BaseBuddyToolError as e:
        typer.secho(f"Tool Execution Error: {e}", fg=typer.colors.RED, err=True)
        if e.details:
            typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except bb_utils.BaseBuddyFileError as e:
        typer.secho(f"File Error: {e}", fg=typer.colors.RED, err=True)
        if e.details:
            typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.secho(f"An unexpected error occurred: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

@app.command()
def spike(
    ctx: typer.Context, # Added context
    reference: Optional[Path] = typer.Option(None, "--reference", "-r", help="Path to the reference FASTA. Uses cached GRCh38 if not specified."),
    in_bam: List[Path] = typer.Option(..., "--input-bam", "--in-bam", "-i", help="Input BAM file(s) to spike variants into. Required."),
    snp_vcf_file: Optional[Path] = typer.Option(None, "--snp-vcf", help="VCF file with SNVs to spike."),
    indel_vcf_file: Optional[Path] = typer.Option(None, "--indel-vcf", help="VCF file with Indels to spike."),
    output_prefix: str = typer.Option("spiked_vars", "--output-prefix", "-p", help="Prefix for output BAM files and run directory."),
    vaf: float = typer.Option(0.05, "--vaf", help="Variant Allele Fraction for spiked variants."),
    seed: int = typer.Option(0, "--seed", help="Random seed for spiking."),
    overwrite_output: bool = typer.Option(False, "--overwrite", help="Overwrite output directory if it exists."),
    auto_index_input_bam: bool = typer.Option(True, "--auto-index-bam/--no-auto-index-bam", help="Automatically index input BAMs if needed."),
    auto_index_fasta: bool = typer.Option(True, "--auto-index-fasta/--no-auto-index-fasta", help="Automatically index FASTA if needed."),
    picard_jar_path: Optional[Path] = typer.Option(None, "--picard-jar", help="Path to picard.jar, required by BAMSurgeon. Can also be set via BAMSURGEON_PICARD_JAR env var.") # New option
):
    """
    Spike SNVs and/or Indels into existing BAM file(s) using VCF files.
    Example: basebuddy spike -i original.bam --snp-vcf snps.vcf --indel-vcf indels.vcf -p spiked_output
    """
    if not snp_vcf_file and not indel_vcf_file:
        typer.secho("Error: At least one VCF file (--snp-vcf or --indel-vcf) must be provided.", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    if not in_bam:
        typer.secho("Error: At least one input BAM file must be provided with --input-bam.", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    variants_list_for_runner: List[Dict[str, Any]] = []
    if snp_vcf_file:
        typer.secho(f"Note: --snp-vcf ('{snp_vcf_file}') provided. VCF processing logic in runner is pending.", fg=typer.colors.YELLOW, err=True)
    if indel_vcf_file:
        typer.secho(f"Note: --indel-vcf ('{indel_vcf_file}') provided. VCF processing logic in runner is pending.", fg=typer.colors.YELLOW, err=True)

    input_bams_str = [str(p) for p in in_bam]

    output_prefix_path = Path(output_prefix)
    output_root_for_runner = output_prefix_path.parent
    output_name_for_runner = output_prefix_path.name

    command_params_for_runner = {
        "reference_fasta": str(reference) if reference else None,
        "input_bams": input_bams_str,
        "num_input_bams": len(input_bams_str),
        "snp_vcf_file": str(snp_vcf_file) if snp_vcf_file else None,
        "indel_vcf_file": str(indel_vcf_file) if indel_vcf_file else None,
        "num_variants_to_spike_requested": len(variants_list_for_runner),
        "output_prefix_for_bam": output_name_for_runner,
        "vaf": vaf,
        "seed": seed,
        "overwrite_output": overwrite_output,
        "auto_index_input_bam": auto_index_input_bam,
        "auto_index_fasta": auto_index_fasta,
        "picard_jar": str(picard_jar_path) if picard_jar_path else None, # Add this
        "output_root_dir": str(output_root_for_runner)
    }

    try:
        # This call is to the *current* runner.spike_variants which expects `variants_list`.
        # This will be changed when `runner.spike_variants` is refactored in the next step.
        results = runner.spike_variants(
            output_root_dir=output_root_for_runner,
            reference_fasta=str(reference) if reference else None,
            input_bams=input_bams_str,
            variants_list=variants_list_for_runner,
            output_prefix_for_bam=output_name_for_runner,
            run_name=None,
            command_params=command_params_for_runner,
            overwrite_output=overwrite_output,
            auto_index_input_bam=auto_index_input_bam,
            auto_index_fasta=auto_index_fasta
        )
        typer.secho(f"Variant spiking call (structure pending full VCF support in runner) completed for run: {results.get('run_name')}", fg=typer.colors.GREEN)
        typer.echo(f"Output directory: {results.get('output_directory')}")
        if results.get('manifest_path'):
            typer.echo(f"Manifest: {results.get('manifest_path')}")
        if results.get("errors_per_bam") and results["errors_per_bam"]:
            typer.secho("Some BAMs encountered errors during processing:", fg=typer.colors.YELLOW)
            for error_info in results["errors_per_bam"]:
                typer.secho(f"  BAM: {error_info['source_bam']}, Error: {error_info['error']}", fg=typer.colors.YELLOW)

    except bb_utils.BaseBuddyConfigError as e:
        typer.secho(f"Configuration Error: {e}", fg=typer.colors.RED, err=True)
        if e.details: typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except bb_utils.BaseBuddyToolError as e:
        typer.secho(f"Tool Execution Error: {e}", fg=typer.colors.RED, err=True)
        if e.details: typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except bb_utils.BaseBuddyFileError as e:
        typer.secho(f"File Error: {e}", fg=typer.colors.RED, err=True)
        if e.details: typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.secho(f"An unexpected error occurred: {type(e).__name__} - {e}", fg=typer.colors.RED, err=True)
        # import traceback
        # typer.secho(traceback.format_exc(), fg=typer.colors.YELLOW, err=True) # Uncomment for debug
        raise typer.Exit(code=1)

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
    in_bam: Path = typer.Option(..., "--input-bam", "-i", help="Input BAM file to process.", exists=True, file_okay=True, dir_okay=False, readable=True),
    out_bam: Path = typer.Option(Path("biased.bam"), "--output-bam", "-o", help="Output BAM file with introduced bias."),
    forward_fraction: float = typer.Option(0.5, "--forward-fraction", "-f", min=0.0, max=1.0, help="Fraction of reads to keep from the forward strand."),
    seed: int = typer.Option(0, "--seed", "-s", help="Random seed for subsampling."),
):
    """
    Introduce strand bias into an existing BAM file by subsampling reads based on strand.
    Example: basebuddy strand-bias -i original.bam -o biased.bam -f 0.8
    """
    try:
        # Ensure output directory exists for out_bam
        out_bam.parent.mkdir(parents=True, exist_ok=True)

        runner.introduce_strand_bias(str(in_bam), str(out_bam), forward_fraction, seed)
        typer.secho(f"Strand bias introduction complete. Output: {out_bam}", fg=typer.colors.GREEN)
        # Check for index (BAI) file. Common locations: file.bam.bai or file.bai
        bai_path_1 = out_bam.with_suffix(out_bam.suffix + ".bai")
        bai_path_2 = out_bam.with_suffix(".bai") # Handles cases like out.sorted.bam -> out.sorted.bai
        if bai_path_1.exists() or bai_path_2.exists():
            typer.echo("Output BAM has been indexed.")
        else:
            typer.secho("Warning: Output BAM index not found. Indexing may be required for some tools.", fg=typer.colors.YELLOW)

    except bb_utils.BaseBuddyConfigError as e:
        typer.secho(f"Configuration Error: {e}", fg=typer.colors.RED, err=True)
        if e.details: typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except bb_utils.BaseBuddyToolError as e:
        typer.secho(f"Tool Execution Error: {e}", fg=typer.colors.RED, err=True)
        if e.details: typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except bb_utils.BaseBuddyFileError as e:
        typer.secho(f"File Error: {e}", fg=typer.colors.RED, err=True)
        if e.details: typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.secho(f"An unexpected error occurred: {type(e).__name__} - {e}", fg=typer.colors.RED, err=True)
        # import traceback
        # typer.secho(traceback.format_exc(), fg=typer.colors.YELLOW, err=True) # For debugging
        raise typer.Exit(code=1)

@app.command("qc")
def qc(
    fastq_files: List[Path] = typer.Argument(..., help="One or more FASTQ files to process."),
    output_dir: Path = typer.Option(Path("results_qc"), "--output-dir", "-o", help="Output directory for FastQC results.", writable=True, file_okay=False, dir_okay=True),
    run_name: Optional[str] = typer.Option(None, "--run-name", help="Optional name for this specific QC run."),
    threads: int = typer.Option(1, "--threads", "-t", help="Number of threads for FastQC to use per file processing."),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite output directory if it exists.")
):
    """
    Run FastQC on one or more FASTQ files to generate quality control reports.
    Example: basebuddy qc reads1.fastq.gz reads2.fastq.gz -o ./qc_output
    """
    if not fastq_files:
        typer.secho("Error: At least one FASTQ file must be provided.", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        typer.secho(f"Error creating output directory {output_dir}: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    fastq_file_paths_str = [str(f.resolve()) for f in fastq_files]

    command_params_for_runner = {
        "threads": threads
    }

    try:
        results = runner.run_fastq_qc(
            fastq_files=fastq_file_paths_str,
            output_dir_str=str(output_dir.resolve()),
            run_name=run_name,
            command_params=command_params_for_runner,
            overwrite_output=overwrite
        )

        typer.secho(f"FastQC run '{results.get('run_name')}' completed successfully.", fg=typer.colors.GREEN)
        typer.echo(f"Output directory: {results.get('output_directory')}")
        if results.get('manifest_path'):
            typer.echo(f"Manifest file: {results.get('manifest_path')}")

        reports = results.get("qc_reports", [])
        if reports:
            typer.echo("\nGenerated Reports:")
            for report_info in reports:
                input_fq = report_info.get('input_fastq')
                html_path = report_info.get('full_html_path')
                status = "OK" if html_path and Path(html_path).exists() else f"ERROR ({report_info.get('error', 'unknown')})"
                typer.echo(f"  - {input_fq}: {html_path if html_path else 'N/A'} (Status: {status})")
        else:
            typer.echo("No QC reports were generated (or found). Check logs for details.")


    except bb_utils.BaseBuddyConfigError as e:
        typer.secho(f"Configuration Error: {e}", fg=typer.colors.RED, err=True)
        if hasattr(e, 'details') and e.details: typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except bb_utils.BaseBuddyToolError as e:
        typer.secho(f"Tool Execution Error: {e}", fg=typer.colors.RED, err=True)
        if hasattr(e, 'details') and e.details: typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except bb_utils.BaseBuddyFileError as e:
        typer.secho(f"File Error: {e}", fg=typer.colors.RED, err=True)
        if hasattr(e, 'details') and e.details: typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.secho(f"An unexpected error occurred: {type(e).__name__} - {e}", fg=typer.colors.RED, err=True)
        # import traceback
        # typer.secho(traceback.format_exc(), fg=typer.colors.YELLOW, err=True) # For debugging
        raise typer.Exit(code=1)

if __name__ == "__main__":
    app()
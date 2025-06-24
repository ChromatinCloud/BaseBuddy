import typer
from pathlib import Path
from typing import Optional, List, Dict, Any
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
    auto_index_fasta: bool = typer.Option(True, "--auto-index/--no-auto-index", help="Automatically index FASTA if .fai is missing."),
    auto_align: bool = typer.Option(False, "--auto-align", help="Automatically align reads to reference and output BAM file."),
    aligner: str = typer.Option("bwa", "--aligner", help="Aligner to use for auto-align (currently only bwa is supported)."),
    genome_build: str = typer.Option(None, "--build", "-b", help="Genome build (e.g., hg19, hg38, GRCh37, GRCh38). For documentation only.")
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
        "auto_align": auto_align,
        "aligner": aligner,
        "genome_build": genome_build,
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
            auto_index_fasta=auto_index_fasta,
            auto_align=auto_align,
            aligner=aligner,
            genome_build=genome_build
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

    # Validate VCF files exist if provided
    if snp_vcf_file and not snp_vcf_file.exists():
        typer.secho(f"Error: SNP VCF file not found: {snp_vcf_file}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    if indel_vcf_file and not indel_vcf_file.exists():
        typer.secho(f"Error: Indel VCF file not found: {indel_vcf_file}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

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
        # Call the runner with VCF paths
        results = runner.spike_variants(
            output_root_dir=output_root_for_runner,
            reference_fasta=str(reference) if reference else None,
            input_bams=input_bams_str,
            output_prefix_for_bam=output_name_for_runner,
            run_name=None,
            command_params=command_params_for_runner,
            snp_vcf_path=str(snp_vcf_file) if snp_vcf_file else None,
            indel_vcf_path=str(indel_vcf_file) if indel_vcf_file else None,
            overwrite_output=overwrite_output,
            auto_index_input_bam=auto_index_input_bam,
            auto_index_fasta=auto_index_fasta
        )
        typer.secho(f"Variant spiking completed successfully for run: {results.get('run_name')}", fg=typer.colors.GREEN)
        typer.echo(f"Output directory: {results.get('output_directory')}")
        if results.get('manifest_path'):
            typer.echo(f"Manifest: {results.get('manifest_path')}")
        
        # Report on spiked variants
        if results.get('igv_session_xml_path'):
            typer.echo(f"IGV session file: {results.get('igv_session_xml_path')}")
        if results.get('output_bams'):
            typer.echo(f"Generated {len(results.get('output_bams', []))} spiked BAM file(s)")
        
        if results.get("errors_per_bam") and results["errors_per_bam"]:
            typer.secho("Some BAMs encountered errors during processing:", fg=typer.colors.YELLOW)
            for error_info in results["errors_per_bam"]:
                typer.secho(f"  BAM: {error_info['source_bam']}, Error: {error_info['error']}", fg=typer.colors.YELLOW)

    except bb_utils.BaseBuddyConfigError as e:
        typer.secho(f"Configuration Error: {e}", fg=typer.colors.RED, err=True)
        if e.details: 
            typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        if "picard" in str(e).lower():
            typer.secho("\nHint: Set BAMSURGEON_PICARD_JAR environment variable or use --picard-jar flag", fg=typer.colors.YELLOW)
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
    reference: Optional[Path] = typer.Option(None, "--reference", "-r", help="Path to reference FASTA. Uses cached GRCh38 if not specified."),
    outdir: Path = typer.Option(Path("results_sig"), "--outdir", "-o", help="Output directory for signature simulation results."),
    sig_type: str = typer.Option("SBS", "--sig-type", "-t", help="Signature type: SBS (single base substitution), DBS (doublet), ID (indel)."),
    num_mutations: int = typer.Option(1000, "--num-mutations", "-n", help="Number of mutations to simulate."),
    sample_id: str = typer.Option("Sample", "--sample-id", "-s", help="Sample identifier for output files."),
    exome: bool = typer.Option(False, "--exome", help="Simulate on exome regions only (requires genome build)."),
    chrom_based: bool = typer.Option(False, "--chrom-based", help="Generate mutations chromosome by chromosome."),
    seed: int = typer.Option(0, "--seed", help="Random seed for reproducibility."),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite output directory if it exists.")
):
    """
    Simulate mutational signatures using SigProfilerSimulator.
    
    Examples:
      basebuddy signature --sig-type SBS --num-mutations 5000 --sample-id tumor_sample
      basebuddy signature --reference hg38.fa --sig-type ID --num-mutations 100 --exome
    """
    # Validate signature type
    valid_sig_types = ["SBS", "DBS", "ID"]
    if sig_type not in valid_sig_types:
        typer.secho(f"Error: Invalid signature type '{sig_type}'. Must be one of: {', '.join(valid_sig_types)}", 
                   fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    
    # Validate num_mutations
    if num_mutations <= 0:
        typer.secho(f"Error: Number of mutations must be positive, got {num_mutations}", 
                   fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    
    command_params = {
        "sig_type": sig_type,
        "num_mutations": num_mutations,
        "sample_id": sample_id,
        "exome": exome,
        "chrom_based": chrom_based,
        "seed": seed,
        "overwrite_output": overwrite
    }
    
    try:
        results = runner.simulate_signatures(
            reference_fasta=str(reference) if reference else None,
            output_root_dir=outdir,
            sig_type=sig_type,
            num_mutations=num_mutations,
            sample_id=sample_id,
            run_name=None,
            command_params=command_params,
            overwrite_output=overwrite
        )
        
        typer.secho(f"Mutational signature simulation completed for run: {results.get('run_name')}", 
                   fg=typer.colors.GREEN)
        typer.echo(f"Output directory: {results.get('output_directory')}")
        if results.get('manifest_path'):
            typer.echo(f"Manifest file: {results.get('manifest_path')}")
        
        # Report generated files
        output_files = results.get('output_files', [])
        if output_files:
            typer.echo("\nGenerated files:")
            for file_info in output_files:
                typer.echo(f"  - {file_info.get('name')}: {file_info.get('path')}")
        else:
            typer.secho("Warning: No output files detected. Check logs for details.", fg=typer.colors.YELLOW)
            
    except bb_utils.BaseBuddyConfigError as e:
        typer.secho(f"Configuration Error: {e}", fg=typer.colors.RED, err=True)
        if e.details:
            typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        if "sigprofilersimulator" in str(e).lower():
            typer.secho("\nHint: Install SigProfilerSimulator with: pip install SigProfilerSimulator", 
                       fg=typer.colors.YELLOW)
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
        typer.secho(f"An unexpected error occurred: {type(e).__name__} - {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

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

@app.command("download-reference")
def download_reference(
    url: str = typer.Argument(..., help="URL of the reference file to download"),
    output_dir: Path = typer.Option(Path("references"), "--output-dir", "-o", help="Output directory for downloaded reference"),
    filename: Optional[str] = typer.Option(None, "--filename", "-f", help="Filename to save as (defaults to URL filename)"),
    checksum: Optional[str] = typer.Option(None, "--checksum", help="Expected checksum (SHA256) for verification"),
    checksum_type: str = typer.Option("sha256", "--checksum-type", help="Checksum algorithm (sha256, md5)"),
    run_name: Optional[str] = typer.Option(None, "--run-name", help="Custom run name"),
    timeout: int = typer.Option(10800, "--timeout", help="Download timeout in seconds (default: 3 hours)"),
    overwrite: bool = typer.Option(False, "--overwrite", help="Overwrite existing files")
):
    """
    Download reference genomes or other reference files.
    
    Examples:
      basebuddy download-reference https://example.com/GRCh38.fa.gz --checksum abc123...
      basebuddy download-reference ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    """
    # Import here to avoid circular imports
    from basebuddy.bb_runners import download_reference_runner
    
    # Determine filename
    if filename is None:
        filename = url.split('/')[-1]
        if not filename:
            typer.secho("Error: Could not determine filename from URL. Please specify --filename", 
                       fg=typer.colors.RED, err=True)
            raise typer.Exit(code=1)
    
    # Validate checksum if provided
    if checksum and checksum_type not in ["sha256", "md5", "sha1"]:
        typer.secho(f"Error: Unsupported checksum type '{checksum_type}'", 
                   fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    
    command_params = {
        "url": url,
        "filename": filename,
        "checksum": checksum,
        "checksum_type": checksum_type,
        "timeout": timeout,
        "overwrite": overwrite
    }
    
    try:
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run download
        download_reference_runner(
            output_root_dir=output_dir,
            run_name=run_name or bb_utils.generate_unique_run_name("download_ref"),
            command_params=command_params,
            download_url=url,
            destination_filename=filename,
            expected_checksum=checksum or "skip",  # bb_runners expects a value
            checksum_algorithm=checksum_type,
            timeout_download=float(timeout),
            overwrite_output=overwrite
        )
        
        typer.secho(f"\nSuccessfully downloaded: {filename}", fg=typer.colors.GREEN)
        typer.echo(f"Location: {output_dir / filename}")
        
    except bb_utils.BaseBuddyFileError as e:
        typer.secho(f"File Error: {e}", fg=typer.colors.RED, err=True)
        if e.details:
            typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except bb_utils.BaseBuddyToolError as e:
        typer.secho(f"Download Error: {e}", fg=typer.colors.RED, err=True)
        if e.details:
            typer.secho(f"Details: {e.details}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.secho(f"Unexpected error: {type(e).__name__} - {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)


@app.command("apply-signature")
def apply_signature(
    reference: Path = typer.Argument(..., help="Input reference FASTA file"),
    signature_file: Optional[Path] = typer.Option(None, "--signature-file", "-s", help="Path to signature TSV file"),
    signature_name: Optional[str] = typer.Option(None, "--signature-name", "-n", help="Name of signature to apply (e.g., SBS1)"),
    output: Path = typer.Option(Path("mutated.fa"), "--output", "-o", help="Output FASTA file"),
    num_mutations: int = typer.Option(1000, "--num-mutations", "-m", help="Number of mutations to apply"),
    seed: int = typer.Option(42, "--seed", help="Random seed for reproducibility"),
    bundled_type: Optional[str] = typer.Option(None, "--bundled-type", help="Use bundled signatures (sbs, dbs, id)")
):
    """
    Apply mutational signatures to a reference FASTA file.
    
    Examples:
      # Use bundled signatures
      basebuddy apply-signature genome.fa --bundled-type sbs --signature-name SBS1 --num-mutations 5000
      
      # Use custom signature file
      basebuddy apply-signature genome.fa --signature-file custom_sigs.tsv --signature-name MySignature
    """
    # Validate inputs
    if not reference.exists():
        typer.secho(f"Error: Reference file not found: {reference}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    
    if not signature_file and not bundled_type:
        typer.secho("Error: Must specify either --signature-file or --bundled-type", 
                   fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    
    if not signature_name:
        typer.secho("Error: Must specify --signature-name", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    
    try:
        # Determine signature file path
        if bundled_type:
            # Use bundled signature files
            from pathlib import Path as PathlibPath
            data_dir = PathlibPath(__file__).parent / "data"
            
            bundled_files = {
                "sbs": "sbs_grch37_cosmic_v3.3.tsv",
                "dbs": "dbs_grch37_cosmic_v3.3.tsv",
                "id": "id_grch37_cosmic_v3.3.tsv"
            }
            
            if bundled_type.lower() not in bundled_files:
                typer.secho(f"Error: Invalid bundled type '{bundled_type}'. Choose from: sbs, dbs, id",
                           fg=typer.colors.RED, err=True)
                raise typer.Exit(code=1)
            
            signature_file = data_dir / bundled_files[bundled_type.lower()]
            if not signature_file.exists():
                typer.secho(f"Error: Bundled signature file not found: {signature_file}",
                           fg=typer.colors.RED, err=True)
                raise typer.Exit(code=1)
        
        # Import signature utilities
        from basebuddy.signature_utils import parse_signature_matrix_tsv, SignatureFormatError
        
        # Parse signature file
        typer.echo(f"Loading signatures from: {signature_file}")
        signatures = parse_signature_matrix_tsv(signature_file)
        
        if signature_name not in signatures:
            available = list(signatures.keys())
            typer.secho(f"Error: Signature '{signature_name}' not found in file", 
                       fg=typer.colors.RED, err=True)
            typer.echo(f"Available signatures: {', '.join(available[:10])}{'...' if len(available) > 10 else ''}")
            raise typer.Exit(code=1)
        
        # Apply signature to FASTA
        typer.echo(f"Applying {num_mutations} mutations from signature '{signature_name}'...")
        
        # Create a temporary run directory
        output_dir = output.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        command_params = {
            "signature_source": "bundled" if bundled_type else "custom",
            "signature_file": str(signature_file),
            "signature_name": signature_name,
            "num_mutations": num_mutations,
            "seed": seed
        }
        
        # The runner expects signature_id_or_path to be either a signature ID or file path
        # For bundled signatures, just pass the signature name
        # For custom files, we need to pass the file path and let the runner parse it
        if bundled_type:
            signature_id_or_path = signature_name
        else:
            signature_id_or_path = str(signature_file)
        
        result = runner.apply_signature_to_fasta(
            output_root_dir=output_dir,
            run_name=None,  # Will be auto-generated
            command_params=command_params,
            input_fasta_path_str=str(reference),
            output_fasta_name=output.name,
            signature_id_or_path=signature_id_or_path,
            num_mutations=num_mutations,
            target_regions=None,
            overwrite_output=True,  # We already checked if file exists
            auto_index_input_fasta=True
        )
        
        # Move the output file to the desired location if needed
        generated_file = Path(result["output_directory"]) / output.name
        if generated_file != output and generated_file.exists():
            import shutil
            shutil.move(str(generated_file), str(output))
        
        typer.secho(f"\nSuccessfully created mutated FASTA: {output}", fg=typer.colors.GREEN)
        typer.echo(f"Applied {num_mutations} mutations using signature '{signature_name}'")
        
        # Index the output if it's a FASTA
        if output.suffix.lower() in [".fa", ".fasta", ".fna"]:
            typer.echo("Indexing output FASTA...")
            import subprocess
            try:
                subprocess.run(["samtools", "faidx", str(output)], check=True, capture_output=True)
                typer.echo("Indexing complete.")
            except (subprocess.CalledProcessError, FileNotFoundError):
                typer.secho("Warning: Could not index output FASTA (samtools not found)", 
                           fg=typer.colors.YELLOW)
        
    except SignatureFormatError as e:
        typer.secho(f"Signature Format Error: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        typer.secho(f"File Error: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        typer.secho(f"Unexpected error: {type(e).__name__} - {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
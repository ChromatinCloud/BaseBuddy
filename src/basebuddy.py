# File 3: basebuddy.py (Main CLI)

import argparse
import logging
import sys
from pathlib import Path
from bb_utils import BaseBuddyError # Import base custom error
from bb_runners import (
    simulate_short_reads_runner,
    spike_variants_runner,
    download_reference_runner
)

__version__ = "0.2.0" # Example version

# --- Global Logger Setup ---
def setup_logging(log_level_str: str = "INFO", log_file: Optional[str] = None):
    """Configures global logging."""
    numeric_level = getattr(logging, log_level_str.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level_str}")

    # Use a more informative format
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    date_format = "%Y-%m-%d %H:%M:%S"
    
    handlers: List[logging.Handler] = [logging.StreamHandler(sys.stderr)] # Log to stderr by default
    if log_file:
        try:
            file_handler = logging.FileHandler(log_file, mode='a') # Append mode
            handlers.append(file_handler)
        except IOError as e:
            print(f"Warning: Could not open log file {log_file}. Logging to stderr only. Error: {e}", file=sys.stderr)

    logging.basicConfig(level=numeric_level, format=log_format, datefmt=date_format, handlers=handlers)
    # Suppress overly verbose logs from common libraries if needed
    # logging.getLogger("shutil").setLevel(logging.WARNING)


def create_main_parser() -> argparse.ArgumentParser:
    """Creates the main ArgumentParser for BaseBuddy CLI."""
    parser = argparse.ArgumentParser(
        description=f"BaseBuddy CLI - Sequencing Data Simulation Toolkit (v{__version__}).",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Example: basebuddy short --reference ref.fa --depth 30 --readlen 100 --profile HS25 --outdir ./sim_out"
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging verbosity (default: INFO)."
    )
    parser.add_argument(
        "--log-file", type=str, default=None,
        help="Path to a file for logging output (appends if file exists)."
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Allow overwriting of existing output files/directories where applicable."
    )


    subparsers = parser.add_subparsers(dest="command", title="Available commands",
                                       help="Sub-command help", required=True)
    
    # --- Subparser for "short" read simulation (ART) ---
    parser_short = subparsers.add_parser(
        "short",
        help="Simulate short reads using ART.",
        description="Simulates short sequencing reads from a reference genome using the ART toolkit."
    )
    parser_short.add_argument(
        "--reference", type=str, required=True,
        help="Path to the reference FASTA file (must be indexed, or use samtools faidx)."
    )
    parser_short.add_argument(
        "--depth", type=int, required=True, help="Desired mean sequencing depth."
    )
    parser_short.add_argument(
        "--readlen", type=int, required=True, help="Length of simulated reads (e.g., 100, 150)."
    )
    parser_short.add_argument(
        "--profile", type=str, required=True,
        help="ART sequencing system profile (e.g., HS25 for Illumina HiSeq 2500)."
    )
    parser_short.add_argument(
        "--outdir", type=str, required=True,
        help="Output directory for FASTQ files and any other ART outputs."
    )
    parser_short.add_argument(
        "--art-platform", type=str, default="illumina", choices=["illumina"], # Extend if ART supports more you wrap
        help="ART simulation platform (default: illumina)."
    )
    parser_short.add_argument(
        "--fragmean", type=int, default=200,
        help="Mean fragment length for paired-end reads (default: 200bp)."
    )
    parser_short.add_argument(
        "--fragstd", type=int, default=10,
        help="Standard deviation of fragment length for paired-end reads (default: 10bp)."
    )
    parser_short.add_argument(
        "--single-end", action="store_true",
        help="Generate single-end reads (default is paired-end)."
    )
    parser_short.add_argument(
        "--timeout", type=float, default=3600.0,
        help="Timeout in seconds for the ART simulation process (default: 3600s / 1hr)."
    )
    parser_short.set_defaults(func=handle_short_reads_command)

    # --- Subparser for "spike" variants ---
    parser_spike = subparsers.add_parser(
        "spike",
        help="Spike variants into an existing BAM file.",
        description="Introduces variants from a VCF file into reads in a BAM file using 'addsnv.py' (conceptual)."
    )
    parser_spike.add_argument(
        "--reference", type=str, required=True, help="Path to the reference FASTA file."
    )
    parser_spike.add_argument(
        "--in_bam", type=str, required=True, help="Path to the input BAM file (must be indexed)."
    )
    parser_spike.add_argument(
        "--vcf", type=str, required=True, help="Path to VCF file containing variants to spike."
    )
    parser_spike.add_argument(
        "--out_prefix", type=str, required=True,
        help="Prefix for the output BAM file (e.g., 'spiked_reads/sample1'). '.bam' will be appended."
    )
    parser_spike.add_argument(
        "--no-auto-index-input", action="store_false", dest="auto_index_input_bam",
        help="Disable automatic indexing of input BAM if .bai is missing (will raise error instead)."
    )
    parser_spike.add_argument(
        "--timeout", type=float, default=7200.0,
        help="Timeout in seconds for the variant spiking process (default: 7200s / 2hr)."
    )
    # Add other necessary parameters for addsnv.py here
    parser_spike.set_defaults(func=handle_spike_variants_command)

    # --- Subparser for "download-ref" ---
    parser_download = subparsers.add_parser(
        "download-ref",
        help="Download a reference genome and verify its checksum.",
        description="Downloads a file (e.g., reference genome), verifies its checksum, and indexes if FASTA."
    )
    parser_download.add_argument(
        "--url", type=str, required=True, help="URL of the file to download."
    )
    parser_download.add_argument(
        "--destination", type=str, required=True, help="Local path to save the downloaded file."
    )
    parser_download.add_argument(
        "--checksum", type=str, required=True, help="Expected checksum (e.g., SHA256) of the file."
    )
    parser_download.add_argument(
        "--algo", type=str, default="sha256", choices=["sha256", "md5", "sha1"], # Add more as needed
        help="Checksum algorithm (default: sha256)."
    )
    parser_download.add_argument(
        "--timeout", type=float, default=10800.0,
        help="Timeout for download process (default: 10800s / 3hr)."
    )
    parser_download.set_defaults(func=handle_download_ref_command)

    return parser

# --- Command Handler Functions ---
def handle_short_reads_command(args: argparse.Namespace, global_args: argparse.Namespace):
    simulate_short_reads_runner(
        reference_fasta=args.reference,
        depth=args.depth,
        read_length=args.readlen,
        art_profile=args.profile,
        output_directory=args.outdir,
        mean_fragment_length=args.fragmean,
        std_dev_fragment_length=args.fragstd,
        is_paired_end=not args.single_end,
        overwrite_output=global_args.overwrite,
        art_platform=args.art_platform,
        timeout=args.timeout
    )

def handle_spike_variants_command(args: argparse.Namespace, global_args: argparse.Namespace):
    spike_variants_runner(
        reference_fasta=args.reference,
        input_bam=args.in_bam,
        vcf_file=args.vcf,
        output_bam_prefix=args.out_prefix,
        overwrite_output=global_args.overwrite,
        auto_index_input_bam=args.auto_index_input_bam,
        timeout=args.timeout
        # Pass other specific args for addsnv.py from 'args' if defined
    )

def handle_download_ref_command(args: argparse.Namespace, global_args: argparse.Namespace):
    # global_args.overwrite might not be directly applicable here unless we check if destination exists
    if Path(args.destination).exists() and not global_args.overwrite:
        raise BaseBuddyFileError(
            f"Destination file '{args.destination}' already exists. Use --overwrite to replace it."
        )
    download_reference_runner(
        download_url=args.url,
        destination_path_str=args.destination,
        expected_checksum=args.checksum,
        checksum_algorithm=args.algo,
        timeout_download=args.timeout
    )

def cli_entry_point():
    """Main entry point for the BaseBuddy CLI."""
    parser = create_main_parser()
    # Separate global args from command-specific args if needed for clarity or pre-processing
    # For now, parse all args together.
    # To parse global args first: parser.parse_known_args()
    all_args = parser.parse_args()

    # Setup logging based on global args
    setup_logging(log_level_str=all_args.log_level, log_file=all_args.log_file)
    logger = logging.getLogger(__name__) # Get logger for this main script
    logger.info(f"BaseBuddy CLI started. Command: {all_args.command}. Version: {__version__}")
    logger.debug(f"Full arguments: {all_args}")


    try:
        # Pass all_args to handlers; they can pick what they need (or pass specific sub-parser args and global_args)
        # For better separation, one might pass all_args.subcommand_args and all_args (for global)
        all_args.func(all_args, all_args) # Passing all_args twice to simulate (command_args, global_args)
                                          # A cleaner way is to use parse_known_args if global args need separate handling
    except BaseBuddyError as e:
        logger.error(f"A BaseBuddy operation failed: {e}", exc_info=False) # exc_info=False because BaseBuddyError formats its details
        print(f"ERROR: {e}", file=sys.stderr) # User-facing error
        sys.exit(1)
    except KeyboardInterrupt:
        logger.warning("Operation cancelled by user (KeyboardInterrupt).")
        print("\nOperation cancelled by user.", file=sys.stderr)
        sys.exit(130) # Standard exit code for Ctrl+C
    except Exception as e:
        logger.critical("An unexpected critical error occurred.", exc_info=True) # exc_info=True for full traceback for unexpected errors
        print(f"CRITICAL UNEXPECTED ERROR: {e}. Check logs for more details.", file=sys.stderr)
        sys.exit(2)
    
    logger.info("BaseBuddy CLI command finished successfully.")

if __name__ == "__main__":
    cli_entry_point()
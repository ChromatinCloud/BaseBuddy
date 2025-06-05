# File 3: basebuddy.py (Main CLI - Updated for FASTA auto-indexing option)

import argparse
import logging
import sys
import json # For list-outputs
from pathlib import Path
from typing import Dict, Any, Optional, List # For list-outputs
import copy # For creating manifest parameters

from bb_utils import BaseBuddyError, BaseBuddyFileError, generate_unique_run_name, read_run_manifest # New imports
from bb_runners import (
    simulate_short_reads_runner,
    spike_variants_runner,
    download_reference_runner
)

__version__ = "0.3.0" # Example version update (or increment if you prefer, e.g., 0.3.1)

# --- Global Logger Setup ( 그대로 유지 ) ---
def setup_logging(log_level_str: str = "INFO", log_file: Optional[str] = None):
    numeric_level = getattr(logging, log_level_str.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level_str}")
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    date_format = "%Y-%m-%d %H:%M:%S"
    handlers: List[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if log_file:
        try:
            file_handler = logging.FileHandler(log_file, mode='a')
            handlers.append(file_handler)
        except IOError as e:
            print(f"Warning: Could not open log file {log_file}. Logging to stderr only. Error: {e}", file=sys.stderr)
    logging.basicConfig(level=numeric_level, format=log_format, datefmt=date_format, handlers=handlers)

def _get_common_output_args(parser: argparse.ArgumentParser):
    """Adds common output arguments to a subparser."""
    group = parser.add_argument_group("Output Options")
    group.add_argument(
        "--run-name", type=str, default=None,
        help="Specify a name for this run. Outputs will be in 'output_root/run_name'. "
             "If not provided, a unique name will be generated (e.g., command_YYYYMMDD_HHMMSS)."
    )
    return group

# --- New Helper for Common Reference Arguments ---
def _get_common_reference_args(parser: argparse.ArgumentParser):
    """Adds common reference and FASTA indexing arguments to a subparser."""
    group = parser.add_argument_group("Reference Options")
    group.add_argument(
        "--reference", type=str, required=True,
        help="Path to the reference FASTA file."
    )
    group.add_argument(
        "--no-auto-index-fasta", action="store_false", dest="auto_index_fasta", default=True,
        help="Disable automatic indexing of the user-provided FASTA reference if its .fai index is missing. "
             "If disabled and index is missing, an error will be raised."
    )
    return group

def _extract_manifest_params(args: argparse.Namespace) -> Dict[str, Any]:
    """Extracts and cleans parameters for manifest from argparse namespace."""
    params_for_manifest = vars(copy.deepcopy(args))
    # Remove sensitive or irrelevant info for manifest if needed
    params_for_manifest.pop('func', None)
    params_for_manifest.pop('log_level', None)
    params_for_manifest.pop('log_file', None)
    params_for_manifest.pop('output_root', None)
    # auto_index_fasta and auto_index_input_bam are operational flags, record them
    # params_for_manifest.pop('auto_index_fasta', None) # Keep it to record the setting
    # params_for_manifest.pop('auto_index_input_bam', None) # Keep it
    return params_for_manifest


def create_main_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=f"BaseBuddy CLI - Sequencing Data Simulation Toolkit (v{__version__}).",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Example: basebuddy --output-root ./my_sims short --reference ref.fa --depth 30 ..."
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    # Global args
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging verbosity (default: INFO)."
    )
    parser.add_argument(
        "--log-file", type=str, default=None,
        help="Path to a file for logging output (appends if file exists)."
    )
    parser.add_argument(
        "--output-root", type=str, default="./basebuddy_outputs",
        help="Root directory for all generated outputs (default: ./basebuddy_outputs)."
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Allow overwriting of existing output files/subdirectories within a run directory if it exists. "
             "Also allows overwriting of the run directory itself if it's not empty (use with caution)."
    )

    subparsers = parser.add_subparsers(dest="command", title="Available commands",
                                       help="Sub-command help", required=True)
    
    # --- Subparser for "short" read simulation ---
    parser_short = subparsers.add_parser("short", help="Simulate short reads using ART.")
    _get_common_reference_args(parser_short) # Use helper for --reference and --no-auto-index-fasta
    # Specific args for 'short'
    parser_short.add_argument("--depth", type=int, required=True, help="Desired mean sequencing depth.")
    parser_short.add_argument("--readlen", type=int, required=True, help="Length of simulated reads.")
    parser_short.add_argument("--profile", type=str, required=True, help="ART sequencing system profile.")
    parser_short.add_argument("--art-platform", type=str, default="illumina", choices=["illumina"], help="ART simulation platform.")
    parser_short.add_argument("--fragmean", type=int, default=200, help="Mean fragment length.")
    parser_short.add_argument("--fragstd", type=int, default=10, help="Std dev of fragment length.")
    parser_short.add_argument("--single-end", action="store_true", help="Generate single-end reads.")
    parser_short.add_argument("--timeout", type=float, default=3600.0, help="Timeout for ART simulation (seconds).")
    _get_common_output_args(parser_short) # Use helper for --run-name
    parser_short.set_defaults(func=handle_short_reads_command)

    # --- Subparser for "spike" variants ---
    parser_spike = subparsers.add_parser("spike", help="Spike variants into a BAM file.")
    _get_common_reference_args(parser_spike) # Use helper for --reference and --no-auto-index-fasta
    # Specific args for 'spike'
    parser_spike.add_argument("--in_bam", type=str, required=True, help="Path to the input BAM file.")
    parser_spike.add_argument("--vcf", type=str, required=True, help="Path to VCF file with variants.")
    parser_spike.add_argument(
        "--out_prefix_name", type=str, required=True,
        help="Name (prefix) for the output BAM file within the run directory (e.g., 'spiked_output'). '.bam' will be appended."
    )
    # Note: The original had "--no-auto-index-input", dest="auto_index_input_bam".
    # Assuming this should be default=True behavior (auto_index_input_bam=True unless flag is given)
    parser_spike.add_argument(
        "--no-auto-index-input-bam", action="store_false", dest="auto_index_input_bam", default=True,
        help="Disable automatic indexing of input BAM if .bai is missing (will raise error instead)."
    )
    parser_spike.add_argument("--timeout", type=float, default=7200.0, help="Timeout for variant spiking (seconds).")
    _get_common_output_args(parser_spike) # Use helper for --run-name
    parser_spike.set_defaults(func=handle_spike_variants_command)

    # --- Subparser for "download-ref" ---
    parser_download = subparsers.add_parser("download-ref", help="Download a reference genome.")
    parser_download.add_argument("--url", type=str, required=True, help="URL of the file to download.")
    parser_download.add_argument("--filename", type=str, required=True, help="Filename to save the downloaded file as (within the run directory).")
    parser_download.add_argument("--checksum", type=str, required=True, help="Expected checksum (e.g., SHA256).")
    parser_download.add_argument("--algo", type=str, default="sha256", choices=["sha256", "md5", "sha1"], help="Checksum algorithm.")
    parser_download.add_argument("--timeout", type=float, default=10800.0, help="Timeout for download (seconds).")
    _get_common_output_args(parser_download) # Use helper for --run-name
    parser_download.set_defaults(func=handle_download_ref_command)

    # --- Subparser for "list-outputs" (or "ls") ---
    parser_ls = subparsers.add_parser(
        "list-outputs", aliases=["ls"],
        help="List available output runs and their contents.",
        description="Scans the output root directory for completed runs and lists their details from manifest files."
    )
    parser_ls.add_argument(
        "run_name_pattern", nargs="?", type=str, default=None, 
        help="Optional: Name or pattern of specific run(s) to list. If not provided, lists all runs. "
             "Simple glob patterns like '*' or 'prefix_*' are supported."
    )
    parser_ls.add_argument(
        "--show-all-files", action="store_true",
        help="If listing a specific run, show all files in its directory, not just those in manifest."
    )
    parser_ls.set_defaults(func=handle_list_outputs_command)

    return parser

# --- Command Handler Functions ---
def handle_short_reads_command(args: argparse.Namespace, global_args: argparse.Namespace):
    output_root = Path(global_args.output_root)
    run_name = args.run_name if args.run_name else generate_unique_run_name("short")
    manifest_params = _extract_manifest_params(args)
    
    simulate_short_reads_runner(
        output_root_dir=output_root,
        run_name=run_name,
        command_params=manifest_params,
        reference_fasta=args.reference,
        depth=args.depth,
        read_length=args.readlen,
        art_profile=args.profile,
        mean_fragment_length=args.fragmean,
        std_dev_fragment_length=args.fragstd,
        is_paired_end=not args.single_end,
        # overwrite_output from global_args is now included in manifest_params if needed by runner's _prepare_run_output_directory
        # or directly through global_args.overwrite in the runner. For consistency, let's ensure it's accessible.
        # command_params should contain 'overwrite': global_args.overwrite
        art_platform=args.art_platform,
        timeout=args.timeout,
        auto_index_fasta=args.auto_index_fasta # Pass the new flag
    )

def handle_spike_variants_command(args: argparse.Namespace, global_args: argparse.Namespace):
    output_root = Path(global_args.output_root)
    run_name = args.run_name if args.run_name else generate_unique_run_name("spike")
    manifest_params = _extract_manifest_params(args)

    spike_variants_runner(
        output_root_dir=output_root,
        run_name=run_name,
        command_params=manifest_params, # This will contain in_bam, vcf, etc.
        reference_fasta=args.reference, # Explicitly passed
        # input_bam=args.in_bam, # Now expected to be in command_params if runner needs it
        # vcf_file=args.vcf,     # Now expected to be in command_params if runner needs it
        output_bam_prefix_rel=args.out_prefix_name,
        # overwrite_output=global_args.overwrite,
        auto_index_input_bam=args.auto_index_input_bam,
        auto_index_fasta=args.auto_index_fasta, # Pass the new flag
        timeout=args.timeout
    )

def handle_download_ref_command(args: argparse.Namespace, global_args: argparse.Namespace):
    output_root = Path(global_args.output_root)
    run_name = args.run_name if args.run_name else generate_unique_run_name(f"download_{Path(args.filename).stem}")
    manifest_params = _extract_manifest_params(args)
    
    # Pass global_args.overwrite to the runner for it to handle specifically for downloaded file
    download_reference_runner(
        output_root_dir=output_root,
        run_name=run_name,
        command_params=manifest_params,
        download_url=args.url,
        destination_filename=args.filename,
        expected_checksum=args.checksum,
        checksum_algorithm=args.algo,
        timeout_download=args.timeout,
        overwrite_output=global_args.overwrite # Pass global overwrite
    )

def handle_list_outputs_command(args: argparse.Namespace, global_args: argparse.Namespace):
    logger = logging.getLogger(__name__) 
    output_root = Path(global_args.output_root).resolve()
    if not output_root.is_dir():
        print(f"Output root directory '{output_root}' does not exist or is not a directory.", file=sys.stderr)
        logger.error(f"Output root directory '{output_root}' not found for list-outputs.")
        return

    print(f"--- Listing Outputs in: {output_root} ---")
    
    found_runs = False
    run_dirs_to_scan: List[Path] = []

    if args.run_name_pattern:
        for item in output_root.glob(args.run_name_pattern):
            if item.is_dir():
                run_dirs_to_scan.append(item)
        if not run_dirs_to_scan:
            print(f"No run directories found matching pattern '{args.run_name_pattern}' in {output_root}.")
            return
    else:
        for item in output_root.iterdir():
            if item.is_dir() and (item / "manifest.json").is_file():
                run_dirs_to_scan.append(item)
        if not run_dirs_to_scan:
            print(f"No completed runs found in {output_root} (looking for 'manifest.json' in subdirectories).")
            return
            
    run_dirs_to_scan.sort()

    for run_dir in run_dirs_to_scan:
        found_runs = True
        print(f"\nRun: {run_dir.name}")
        print(f"  Full Path: {run_dir}")
        manifest_path = run_dir / "manifest.json"
        manifest_data = read_run_manifest(manifest_path)

        if manifest_data:
            print(f"  Command: {manifest_data.get('command', 'N/A')}")
            print(f"  Timestamp (UTC): {manifest_data.get('timestamp_utc', 'N/A')}")
            print(f"  Status: {manifest_data.get('status', 'N/A')}")
            if manifest_data.get('reference_genome_path'):
                 print(f"  Reference Genome: {manifest_data.get('reference_genome_path')}")
            print("  Manifest Outputs:")
            for out_file in manifest_data.get("outputs", []):
                rel_path = out_file.get('path', 'N/A')
                abs_path = run_dir / rel_path
                exists_str = "(exists)" if abs_path.exists() else "(MISSING!)"
                print(f"    - {out_file.get('name', 'Unnamed')} ({out_file.get('type', 'Unknown')}): {rel_path} {exists_str}")
        else:
            print(f"  Manifest file ('manifest.json') missing or unreadable in {run_dir.name}.")
            if args.show_all_files: 
                 print("  Directory Contents (raw listing):")
                 try:
                     for item in run_dir.iterdir():
                         print(f"    - {item.name} {'(dir)' if item.is_dir() else '(file)'}")
                 except OSError as e:
                     print(f"    Could not list directory contents: {e}")
        
        if args.show_all_files and manifest_data: 
            print("  Other Files in Directory (not in manifest):")
            manifested_files = {item['path'] for item in manifest_data.get("outputs", [])}
            other_files_found = False
            try:
                for item in run_dir.iterdir():
                    if item.name not in manifested_files and item.name != "manifest.json":
                        other_files_found = True
                        print(f"    - {item.name} {'(dir)' if item.is_dir() else '(file)'}")
                if not other_files_found:
                    print("    (No other files)")
            except OSError as e:
                print(f"    Could not list directory contents: {e}")
                
    if not found_runs and not args.run_name_pattern: 
        print(f"No run directories with 'manifest.json' found in '{output_root}'.")

    print("\n--- End of Listing ---")


def cli_entry_point():
    parser = create_main_parser()
    all_args = parser.parse_args()

    setup_logging(log_level_str=all_args.log_level, log_file=all_args.log_file)
    logger_main = logging.getLogger(__name__) 
    logger_main.info(f"BaseBuddy CLI started. Command: {all_args.command}. Version: {__version__}")
    logger_main.debug(f"Full arguments: {all_args}")

    # Make global_args explicitly available if needed, though all_args contains everything
    # For the func(args, global_args) pattern, we can pass all_args for both if handlers expect it.
    # A more robust way could be to separate global_args from command_specific_args if they truly diverge.
    # current_command_args = all_args 
    # global_command_args = all_args # simplified for now
    
    try:
        # The design where command_params (from _extract_manifest_params(all_args)) and specific args (like args.reference)
        # are both passed to runners can be streamlined.
        # For now, handlers pass necessary args explicitly.
        all_args.func(all_args, all_args) # Passing all_args for both 'command-specific' and 'global' for simplicity
    except BaseBuddyError as e:
        logger_main.error(f"A BaseBuddy operation failed: {e}", exc_info=False)
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        logger_main.warning("Operation cancelled by user (KeyboardInterrupt).")
        print("\nOperation cancelled by user.", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        logger_main.critical("An unexpected critical error occurred.", exc_info=True)
        print(f"CRITICAL UNEXPECTED ERROR: {e}. Check logs for more details.", file=sys.stderr)
        sys.exit(2)
    
    logger_main.info("BaseBuddy CLI command finished successfully.")

if __name__ == "__main__":
    cli_entry_point()
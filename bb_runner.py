# File 2: bb_runners.py

from pathlib import Path
import logging
from typing import List, Dict, Any
from bb_utils import (
    run_external_cmd,
    BaseBuddyInputError,
    BaseBuddyFileError,
    BaseBuddyConfigError,
    BaseBuddyToolError,
    ensure_file_exists,
    ensure_directory_exists,
    check_fasta_indexed,
    check_bam_indexed,
    find_tool_path,
    verify_file_checksum
)

logger = logging.getLogger(__name__)

# Define known ART profiles if applicable, for validation
KNOWN_ART_PROFILES = {
    "illumina": ["HS25", "HSXn", "HSXt", "MSv1", "MSv3", "MiS", "NS50"],
    # Add other platforms if ART supports them and you want to list profiles
}


def simulate_short_reads_runner(
    reference_fasta: str,
    depth: int,
    read_length: int,
    art_profile: str,
    output_directory: str,
    mean_fragment_length: int = 200,
    std_dev_fragment_length: int = 10,
    is_paired_end: bool = True,
    overwrite_output: bool = False,
    art_platform: str = "illumina", # e.g., illumina, 454
    timeout: float = 3600.0 # 1 hour timeout for ART
) -> None:
    """
    Runner function for simulating short reads using ART.
    Validates inputs, prepares environment, runs ART, and checks output.
    """
    logger.info("Initiating short read simulation process...")

    # --- Input Parameter Validation ---
    if depth <= 0:
        raise BaseBuddyInputError(f"Sequencing depth must be a positive integer, got {depth}.")
    if read_length <= 0:
        raise BaseBuddyInputError(f"Read length must be a positive integer, got {read_length}.")
    if mean_fragment_length <= 0:
        raise BaseBuddyInputError(f"Mean fragment length must be positive, got {mean_fragment_length}.")
    if std_dev_fragment_length < 0: # Can be 0 for fixed size
        raise BaseBuddyInputError(f"Fragment length standard deviation cannot be negative, got {std_dev_fragment_length}.")
    
    if art_platform.lower() in KNOWN_ART_PROFILES:
        if art_profile not in KNOWN_ART_PROFILES[art_platform.lower()]:
            logger.warning(f"ART profile '{art_profile}' not in known list for platform '{art_platform}'. "
                           f"Known: {KNOWN_ART_PROFILES[art_platform.lower()]}")
            # Depending on strictness, this could be an error or just a warning.
    elif not art_profile:
         raise BaseBuddyInputError("An ART sequencing profile (e.g., HS25) must be specified.")


    # --- Pre-flight Checks: Tools and Files ---
    art_exe_name = f"art_{art_platform}" # e.g., art_illumina, art_454
    art_exe_path = find_tool_path(art_exe_name)
    samtools_exe_path = find_tool_path("samtools")

    ref_path = ensure_file_exists(reference_fasta, "Reference FASTA")
    check_fasta_indexed(ref_path, samtools_exe_path)

    out_dir_path = Path(output_directory)
    if out_dir_path.exists() and not overwrite_output and any(out_dir_path.iterdir()):
        raise BaseBuddyFileError(
            f"Output directory '{out_dir_path}' is not empty. Use --overwrite or specify a different directory."
        )
    ensure_directory_exists(out_dir_path, "Output directory", create=True)
    
    output_prefix = out_dir_path / f"simulated_{art_platform}_reads"

    # --- Construct and Run ART Command ---
    # Consult ART documentation for precise command structure for the chosen platform.
    # This is a generic example, primarily for art_illumina.
    cmd = [
        art_exe_path,
        "-ss", art_profile,      # Sequencing system profile
        "-i", str(ref_path),     # Input reference FASTA
        "-l", str(read_length),  # Read length
        "-f", str(depth),        # Fold coverage
        "-o", str(output_prefix),# Output file prefix
        "-m", str(mean_fragment_length),
        "-s", str(std_dev_fragment_length),
        # Potentially other flags like --noALN if alignment not needed
    ]
    if is_paired_end:
        cmd.append("-p")  # Paired-end flag

    logger.info(f"Preparing to run {art_exe_name}...")
    run_external_cmd(cmd, timeout_seconds=timeout, stream_output=True) # Stream output for long processes

    # --- Post-flight: Verify ART Output ---
    logger.info(f"{art_exe_name} simulation completed. Verifying output files...")
    expected_files = []
    if is_paired_end:
        expected_files.append(output_prefix.with_name(output_prefix.name + "1.fq"))
        expected_files.append(output_prefix.with_name(output_prefix.name + "2.fq"))
    else:
        expected_files.append(output_prefix.with_suffix(".fq"))
    # ART might also produce .aln files
    # expected_files.append(output_prefix.with_suffix(".aln"))


    for f_path in expected_files:
        ensure_file_exists(f_path, f"Expected ART output file")
    
    logger.info(f"Short read simulation successful. Output is in '{out_dir_path}'.")


def spike_variants_runner(
    reference_fasta: str,
    input_bam: str,
    vcf_file: str,
    output_bam_prefix: str,
    overwrite_output: bool = False,
    auto_index_input_bam: bool = True,
    # Specific parameters for `addsnv.py` or similar tool would go here:
    # e.g., target_vaf: Optional[float] = None,
    timeout: float = 7200.0 # 2 hours timeout for potentially complex BAM operations
) -> None:
    """
    Runner function for spiking variants into a BAM file using a conceptual 'addsnv.py'.
    """
    logger.info("Initiating variant spiking process...")

    # --- Input Parameter Validation & Pre-flight Checks ---
    # Example: if target_vaf is not None and not (0.0 < target_vaf <= 1.0):
    #     raise BaseBuddyInputError(f"Target VAF must be between 0.0 and 1.0, got {target_vaf}")

    addsnv_exe_path = find_tool_path("addsnv.py") # Assuming addsnv.py is in PATH and executable
    samtools_exe_path = find_tool_path("samtools")

    ref_path = ensure_file_exists(reference_fasta, "Reference FASTA")
    check_fasta_indexed(ref_path, samtools_exe_path)

    in_bam_path = ensure_file_exists(input_bam, "Input BAM")
    check_bam_indexed(in_bam_path, samtools_exe_path, auto_index_if_missing=auto_index_input_bam)

    vcf_path = ensure_file_exists(vcf_file, "Input VCF")

    output_final_bam = Path(output_bam_prefix + ".bam")
    ensure_directory_exists(output_final_bam.parent, "Output directory", create=True)

    if output_final_bam.exists() and not overwrite_output:
        raise BaseBuddyFileError(
            f"Output BAM file '{output_final_bam}' already exists. Use --overwrite or specify a different prefix."
        )

    # --- Construct and Run Variant Spiking Command ---
    # This is highly dependent on `addsnv.py`'s actual interface.
    # This example assumes `addsnv.py` takes input files and an output file path.
    # It might output to stdout, requiring piping to `samtools sort`.
    
    # Let's use a temporary name for the direct output of addsnv.py, in case sorting is needed
    temp_addsnv_out_bam = output_final_bam.with_name(output_final_bam.stem + "_temp_addsnv" + ".bam")

    cmd_addsnv = [
        "python", addsnv_exe_path, # If addsnv.py is a Python script not directly executable
        # Or just [addsnv_exe_path] if it has a shebang and execute permissions
        "--reference", str(ref_path),
        "--in_bam", str(in_bam_path),
        "--vcf", str(vcf_path),
        "--out_bam", str(temp_addsnv_out_bam) # Assuming it writes directly
        # Add other addsnv.py specific args here, e.g. --vaf, --threads
    ]

    logger.info("Preparing to run variant spiker (addsnv.py)...")
    run_external_cmd(cmd_addsnv, timeout_seconds=timeout, stream_output=True)
    ensure_file_exists(temp_addsnv_out_bam, "Temporary output BAM from addsnv.py")

    # --- Post-processing: Sort and Index BAM (if necessary) ---
    # Many variant manipulators output unsorted BAMs or BAMs that need re-indexing.
    # Assuming addsnv.py might not produce a sorted/indexed BAM.
    logger.info(f"Sorting BAM: {temp_addsnv_out_bam} -> {output_final_bam}")
    cmd_sort = [
        samtools_exe_path, "sort",
        str(temp_addsnv_out_bam),
        "-o", str(output_final_bam),
        # Add threads, memory options for samtools sort if needed
        # "-@","4", "-m", "4G"
    ]
    run_external_cmd(cmd_sort, timeout_seconds=timeout/2) # Give reasonable time for sort
    
    # Clean up temporary unsorted BAM
    try:
        temp_addsnv_out_bam.unlink()
        logger.debug(f"Removed temporary BAM: {temp_addsnv_out_bam}")
    except OSError as e:
        logger.warning(f"Could not remove temporary BAM {temp_addsnv_out_bam}: {e}")

    logger.info(f"Indexing final BAM: {output_final_bam}")
    check_bam_indexed(output_final_bam, samtools_exe_path, auto_index_if_missing=True) # Index the final sorted BAM

    ensure_file_exists(output_final_bam, "Final output BAM from variant spiking")
    logger.info(f"Variant spiking successful. Final output BAM: '{output_final_bam}'.")


def download_reference_runner(
    download_url: str,
    destination_path_str: str,
    expected_checksum: str,
    checksum_algorithm: str = "sha256",
    timeout_download: float = 10800.0 # 3 hours for large downloads
) -> None:
    """
    Downloads a reference file, verifies its checksum, and indexes it if FASTA.
    NOTE: Actual robust download logic (with retries, progress) is complex.
          This uses a conceptual `wget` or `curl` call.
    """
    logger.info(f"Preparing to download reference from {download_url} to {destination_path_str}.")
    destination_path = Path(destination_path_str)
    ensure_directory_exists(destination_path.parent, "Download destination directory", create=True)

    # --- Download ---
    # Using curl as an example for robust download. wget is another option.
    # For true production, consider using a Python library like 'requests' with streaming and retries.
    curl_exe_path = find_tool_path("curl")
    cmd_download = [
        curl_exe_path,
        "-L", # Follow redirects
        "-o", str(destination_path), # Output file
        download_url
    ]
    logger.info(f"Starting download using curl...")
    try:
        run_external_cmd(cmd_download, timeout_seconds=timeout_download, stream_output=True)
    except BaseBuddyToolError as e:
        # If download failed, try to clean up partial file
        if destination_path.exists():
            try: destination_path.unlink()
            except OSError: logger.warning(f"Could not remove partial download: {destination_path}")
        raise # Re-raise the original error

    ensure_file_exists(destination_path, "Downloaded reference file")
    logger.info("Download completed.")

    # --- Verify Checksum ---
    verify_file_checksum(destination_path, expected_checksum, checksum_algorithm)

    # --- Index if FASTA ---
    if destination_path.suffix.lower() in [".fa", ".fasta", ".fna"]:
        logger.info("Attempting to index downloaded FASTA file...")
        samtools_exe_path = find_tool_path("samtools")
        try:
            run_external_cmd([samtools_exe_path, "faidx", str(destination_path)])
            check_fasta_indexed(destination_path, samtools_exe_path) # Verify .fai was created
            logger.info(f"FASTA file '{destination_path.name}' indexed successfully.")
        except (BaseBuddyToolError, BaseBuddyFileError) as e:
            logger.error(f"Failed to index downloaded FASTA: {e}")
            # Decide if this is a critical failure or a warning
            raise BaseBuddyError(f"Downloaded reference {destination_path.name}, but failed to index it.", details=str(e))
    
    logger.info(f"Reference '{destination_path.name}' processed successfully.")
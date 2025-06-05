# File 2: bb_runners.py (Updated)

from pathlib import Path
import logging
from typing import List, Dict, Any, Optional, cast # Added cast
import copy # For deepcopying parameters for manifest

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
    verify_file_checksum,
    generate_igv_session_xml, # New import
    write_run_manifest        # New import
)

logger = logging.getLogger(__name__)

KNOWN_ART_PROFILES = {
    "illumina": ["HS25", "HSXn", "HSXt", "MSv1", "MSv3", "MiS", "NS50"],
}

def _prepare_run_output_directory(
    output_root_dir: Path,
    run_name: str,
    overwrite_output: bool
) -> Path:
    """Creates and returns the specific output directory for the current run."""
    run_output_dir = output_root_dir / run_name
    if run_output_dir.exists() and not overwrite_output and any(run_output_dir.iterdir()):
        raise BaseBuddyFileError(
            f"Run output directory '{run_output_dir}' is not empty and overwrite is not permitted. "
            "Use --overwrite, --run-name with a new name, or clear the directory."
        )
    ensure_directory_exists(run_output_dir, "Run output directory", create=True)
    logger.info(f"Using run output directory: {run_output_dir.resolve()}")
    return run_output_dir


def simulate_short_reads_runner(
    # New args for output management
    output_root_dir: Path,
    run_name: str,
    command_params: Dict[str, Any], # For manifest
    # Original args
    reference_fasta: str,
    depth: int,
    read_length: int,
    art_profile: str,
    # output_directory: str, # This is now derived from output_root_dir and run_name
    mean_fragment_length: int = 200,
    std_dev_fragment_length: int = 10,
    is_paired_end: bool = True,
    overwrite_output: bool = False, # Passed via command_params or directly
    art_platform: str = "illumina",
    timeout: float = 3600.0
) -> None:
    logger.info(f"Initiating short read simulation for run: {run_name}")
    
    run_output_dir = _prepare_run_output_directory(output_root_dir, run_name, overwrite_output)
    manifest_params = copy.deepcopy(command_params) # Store params before modification

    # --- Input Parameter Validation (remains mostly the same) ---
    if depth <= 0:
        raise BaseBuddyInputError(f"Sequencing depth must be a positive integer, got {depth}.")
    # ... (other validations as before) ...

    # --- Pre-flight Checks: Tools and Files ---
    art_exe_name = f"art_{art_platform}"
    art_exe_path = find_tool_path(art_exe_name)
    samtools_exe_path = find_tool_path("samtools")

    ref_path = ensure_file_exists(reference_fasta, "Reference FASTA").resolve() # Resolve for absolute path
    check_fasta_indexed(ref_path, samtools_exe_path)
    
    # Define output prefix relative to run_output_dir
    art_output_prefix = run_output_dir / f"simulated_{art_platform}_reads"

    # --- Construct and Run ART Command ---
    cmd = [
        art_exe_path, "-ss", art_profile, "-i", str(ref_path),
        "-l", str(read_length), "-f", str(depth), "-o", str(art_output_prefix),
        "-m", str(mean_fragment_length), "-s", str(std_dev_fragment_length)
    ]
    if is_paired_end:
        cmd.append("-p")

    logger.info(f"Preparing to run {art_exe_name}...")
    run_external_cmd(cmd, timeout_seconds=timeout, stream_output=True, cwd=run_output_dir) # cwd might be useful for ART

    # --- Post-flight: Verify ART Output & Manifest ---
    logger.info(f"{art_exe_name} simulation completed. Verifying output files...")
    output_files_manifest: List[Dict[str,str]] = []
    
    # Paths for manifest should be relative to run_output_dir
    if is_paired_end:
        r1_path = art_output_prefix.with_name(art_output_prefix.name + "1.fq")
        r2_path = art_output_prefix.with_name(art_output_prefix.name + "2.fq")
        ensure_file_exists(r1_path, "Expected ART output R1 FASTQ")
        ensure_file_exists(r2_path, "Expected ART output R2 FASTQ")
        output_files_manifest.append({"name": "Simulated Reads (R1)", "path": r1_path.name, "type": "FASTQ"})
        output_files_manifest.append({"name": "Simulated Reads (R2)", "path": r2_path.name, "type": "FASTQ"})
    else:
        r_path = art_output_prefix.with_suffix(".fq")
        ensure_file_exists(r_path, "Expected ART output FASTQ")
        output_files_manifest.append({"name": "Simulated Reads", "path": r_path.name, "type": "FASTQ"})

    # ART may produce other files like .aln; add them if needed
    aln_path = art_output_prefix.with_suffix(".aln")
    if aln_path.exists():
         output_files_manifest.append({"name": "ART Alignment ALN", "path": aln_path.name, "type": "ALN"})
    
    manifest_path = run_output_dir / "manifest.json"
    write_run_manifest(
        manifest_path, run_name, "short", manifest_params,
        output_files_manifest, reference_genome_path=str(ref_path)
    )
    output_files_manifest.append({"name": "Run Manifest", "path": manifest_path.name, "type": "MANIFEST"})

    logger.info(f"Short read simulation successful for run '{run_name}'. Output in '{run_output_dir}'.")
    print(f"\n--- Outputs for run: {run_name} ---")
    print(f"Run directory: {run_output_dir.resolve()}")
    for item in output_files_manifest:
        print(f"  {item['name']} ({item['type']}): {item['path']}")
    print("--- End of Summary ---")


def spike_variants_runner(
    # New args
    output_root_dir: Path,
    run_name: str,
    command_params: Dict[str, Any],
    # Original args
    reference_fasta: str,
    input_bam: str,
    vcf_file: str,
    output_bam_prefix_rel: str, # Relative prefix within the run_output_dir
    overwrite_output: bool = False,
    auto_index_input_bam: bool = True,
    timeout: float = 7200.0
) -> None:
    logger.info(f"Initiating variant spiking for run: {run_name}")
    
    run_output_dir = _prepare_run_output_directory(output_root_dir, run_name, overwrite_output)
    manifest_params = copy.deepcopy(command_params)

    # --- Pre-flight Checks ---
    addsnv_exe_path = find_tool_path("addsnv.py")
    samtools_exe_path = find_tool_path("samtools")

    ref_path = ensure_file_exists(reference_fasta, "Reference FASTA").resolve()
    check_fasta_indexed(ref_path, samtools_exe_path)

    in_bam_path = ensure_file_exists(input_bam, "Input BAM").resolve()
    check_bam_indexed(in_bam_path, samtools_exe_path, auto_index_if_missing=auto_index_input_bam)

    vcf_path = ensure_file_exists(vcf_file, "Input VCF").resolve()

    # Output paths are relative to run_output_dir
    output_final_bam_name = output_bam_prefix_rel + ".bam"
    output_final_bam_path = run_output_dir / output_final_bam_name
    
    # If output_final_bam_path already exists (e.g. due to a previous overwrite_output=True run that failed mid-way for sort)
    # we might need to handle it depending on exact overwrite logic.
    # For now, the _prepare_run_output_directory handles the top-level run_output_dir.
    # If output_final_bam_path specifically exists and overwrite is true, we allow it.
    if output_final_bam_path.exists() and not overwrite_output:
         raise BaseBuddyFileError(
            f"Target output BAM file '{output_final_bam_path}' already exists and overwrite is not permitted."
        )

    temp_addsnv_out_bam_name = output_bam_prefix_rel + "_temp_addsnv.bam"
    temp_addsnv_out_bam_path = run_output_dir / temp_addsnv_out_bam_name

    # --- Construct and Run Variant Spiking Command ---
    cmd_addsnv = [
        "python", addsnv_exe_path,
        "--reference", str(ref_path), "--in_bam", str(in_bam_path),
        "--vcf", str(vcf_path), "--out_bam", str(temp_addsnv_out_bam_path)
    ]
    logger.info("Preparing to run variant spiker (addsnv.py)...")
    run_external_cmd(cmd_addsnv, timeout_seconds=timeout, stream_output=True, cwd=run_output_dir)
    ensure_file_exists(temp_addsnv_out_bam_path, "Temporary output BAM from addsnv.py")

    # --- Post-processing: Sort and Index BAM ---
    logger.info(f"Sorting BAM: {temp_addsnv_out_bam_path.name} -> {output_final_bam_path.name}")
    cmd_sort = [
        samtools_exe_path, "sort", str(temp_addsnv_out_bam_path),
        "-o", str(output_final_bam_path)
    ]
    run_external_cmd(cmd_sort, timeout_seconds=timeout/2, cwd=run_output_dir)
    
    try:
        temp_addsnv_out_bam_path.unlink()
        logger.debug(f"Removed temporary BAM: {temp_addsnv_out_bam_path.name}")
    except OSError as e:
        logger.warning(f"Could not remove temporary BAM {temp_addsnv_out_bam_path.name}: {e}")

    check_bam_indexed(output_final_bam_path, samtools_exe_path, auto_index_if_missing=True)
    ensure_file_exists(output_final_bam_path, "Final output BAM from variant spiking")

    # --- IGV Session and Manifest ---
    output_files_manifest: List[Dict[str,str]] = []
    output_files_manifest.append({"name": "Spiked Variants BAM", "path": output_final_bam_path.name, "type": "BAM"})
    
    bai_path_1 = output_final_bam_path.with_suffix(output_final_bam_path.suffix + ".bai").name
    bai_path_2 = output_final_bam_path.with_suffix(".bai").name
    if (run_output_dir / bai_path_1).exists():
        output_files_manifest.append({"name": "Spiked BAM Index", "path": bai_path_1, "type": "BAI"})
    elif (run_output_dir / bai_path_2).exists():
        output_files_manifest.append({"name": "Spiked BAM Index", "path": bai_path_2, "type": "BAI"})

    # Generate IGV session
    igv_session_file_name = "igv_session_spike.xml"
    igv_session_file_path = run_output_dir / igv_session_file_name
    igv_tracks = [
        {"name": "Spiked BAM", "path": output_final_bam_path.name, "format": "bam", "type": "alignment"},
        # The original VCF could also be added, ensure its path is relative or accessible
        {"name": "Input VCF", "path": str(vcf_path.resolve()), "format": "vcf", "type": "variant"} # Making VCF path absolute for robustness
    ]
    generate_igv_session_xml(igv_session_file_path, str(ref_path), igv_tracks)
    output_files_manifest.append({"name": "IGV Session (Spike)", "path": igv_session_file_name, "type": "IGV_SESSION"})

    manifest_path = run_output_dir / "manifest.json"
    # Ensure parameters for manifest are serializable (e.g. Path objects to str)
    # command_params already handled this in the CLI part by creating manifest_params from args
    write_run_manifest(
        manifest_path, run_name, "spike", manifest_params,
        output_files_manifest, reference_genome_path=str(ref_path)
    )
    output_files_manifest.append({"name": "Run Manifest", "path": manifest_path.name, "type": "MANIFEST"})

    logger.info(f"Variant spiking successful for run '{run_name}'. Final BAM: '{output_final_bam_path.name}'.")
    print(f"\n--- Outputs for run: {run_name} ---")
    print(f"Run directory: {run_output_dir.resolve()}")
    for item in output_files_manifest:
        print(f"  {item['name']} ({item['type']}): {item['path']}")
    print("--- End of Summary ---")


def download_reference_runner(
    # New args
    output_root_dir: Path,
    run_name: str,
    command_params: Dict[str, Any],
    # Original args
    download_url: str,
    destination_filename: str, # Filename, not full path
    expected_checksum: str,
    checksum_algorithm: str = "sha256",
    timeout_download: float = 10800.0,
    overwrite_output: bool = False
) -> None:
    logger.info(f"Preparing to download reference for run '{run_name}' from {download_url} to {destination_filename}.")
    
    run_output_dir = _prepare_run_output_directory(output_root_dir, run_name, overwrite_output)
    manifest_params = copy.deepcopy(command_params)

    destination_path = run_output_dir / destination_filename # Full path to save the file

    if destination_path.exists() and not overwrite_output:
        # If file exists and overwrite is false, we can choose to use existing or error.
        # For a download, it's safer to error or require specific re-verification.
        # Or, if checksum matches, use existing. For simplicity, error out here.
        raise BaseBuddyFileError(
             f"Destination file '{destination_path}' already exists in run directory and overwrite is false."
        )

    # --- Download ---
    curl_exe_path = find_tool_path("curl")
    cmd_download = [curl_exe_path, "-L", "-o", str(destination_path), download_url]
    logger.info(f"Starting download using curl...")
    try:
        run_external_cmd(cmd_download, timeout_seconds=timeout_download, stream_output=True, cwd=run_output_dir)
    except BaseBuddyToolError as e:
        if destination_path.exists():
            try: destination_path.unlink()
            except OSError: logger.warning(f"Could not remove partial download: {destination_path}")
        raise

    ensure_file_exists(destination_path, "Downloaded reference file")
    logger.info("Download completed.")

    # --- Verify Checksum & Manifest ---
    verify_file_checksum(destination_path, expected_checksum, checksum_algorithm)
    output_files_manifest: List[Dict[str,str]] = []
    output_files_manifest.append({"name": "Downloaded File", "path": destination_path.name, "type": "REFERENCE_FILE"})

    # --- Index if FASTA ---
    ref_genome_path_for_manifest = str(destination_path.resolve()) # Store absolute path for potential later use
    if destination_path.suffix.lower() in [".fa", ".fasta", ".fna"]:
        logger.info("Attempting to index downloaded FASTA file...")
        samtools_exe_path = find_tool_path("samtools")
        try:
            run_external_cmd([samtools_exe_path, "faidx", str(destination_path)], cwd=run_output_dir)
            fai_path = destination_path.with_suffix(destination_path.suffix + ".fai")
            check_fasta_indexed(destination_path, samtools_exe_path) # Verify .fai was created
            output_files_manifest.append({"name": "FASTA Index", "path": fai_path.name, "type": "FAI"})
            logger.info(f"FASTA file '{destination_path.name}' indexed successfully.")
        except (BaseBuddyToolError, BaseBuddyFileError) as e:
            logger.error(f"Failed to index downloaded FASTA: {e}")
            write_run_manifest( # Write manifest even on partial failure
                run_output_dir / "manifest.json", run_name, "download-ref", manifest_params,
                output_files_manifest, reference_genome_path=ref_genome_path_for_manifest, status="failed_indexing"
            )
            raise BaseBuddyError(f"Downloaded reference {destination_path.name}, but failed to index it.", details=str(e))
    
    manifest_path = run_output_dir / "manifest.json"
    write_run_manifest(
        manifest_path, run_name, "download-ref", manifest_params,
        output_files_manifest, reference_genome_path=ref_genome_path_for_manifest
    )
    output_files_manifest.append({"name": "Run Manifest", "path": manifest_path.name, "type": "MANIFEST"})

    logger.info(f"Reference '{destination_path.name}' for run '{run_name}' processed successfully.")
    print(f"\n--- Outputs for run: {run_name} ---")
    print(f"Run directory: {run_output_dir.resolve()}")
    for item in output_files_manifest:
        print(f"  {item['name']} ({item['type']}): {item['path']}")
    print("--- End of Summary ---")



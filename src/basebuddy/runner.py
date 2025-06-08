import subprocess
import shutil
from pathlib import Path
import sys
import logging
from typing import Dict, Any, Optional, List
import copy
import tempfile # For temporary VCF
import os # For os.fsencode if needed with tempfile paths, though str(Path) is usually fine
import traceback # For detailed error logging
import pysam # For FastaFile in apply_signature_to_fasta

from basebuddy import utils as bb_utils
# Assuming signature_utils is in the same directory or PYTHONPATH is set up
# For relative import if runner.py and signature_utils.py are in src/basebuddy/
from . import signature_utils as sig_utils


logger = logging.getLogger(__name__)

BASEBUDDY_REF_DIR = Path.home() / ".basebuddy" / "refs"
DEFAULT_REF = BASEBUDDY_REF_DIR / "GRCh38.fa"

# These _ensure_exists, _ensure_reference, _download_grch38 functions
# are kept for other runners until they are refactored.
# _run_cmd was replaced by subprocess.run in non-refactored runners, will be bb_utils.run_external_cmd
def _ensure_exists(file, desc="file"):
    if not Path(file).exists():
        raise FileNotFoundError(f"{desc} not found: {file}")

def _ensure_reference(reference: str | None) -> str:
    ref = Path(reference) if reference else DEFAULT_REF
    if not ref.exists():
        print(f"Reference not found: {ref}, placeholder for download logic.")
        raise FileNotFoundError(f"Reference {ref} not found and auto-download not implemented here.")
    return str(ref)

def _download_grch38(dest: str) -> None:
    logger.warning("Legacy _download_grch38 called. This should be handled by download_reference_runner.")
    (Path(dest).parent).mkdir(parents=True, exist_ok=True)
    Path(dest).touch()
    logger.info(f"Placeholder: Created dummy GRCh38 at {dest}")


def simulate_short(
    output_root_dir: Path,
    reference_fasta: str,
    depth: int,
    read_length: int,
    art_profile: str,
    run_name: Optional[str] = None,
    command_params: Optional[Dict[str, Any]] = None,
    mean_fragment_length: int = 200,
    std_dev_fragment_length: int = 10,
    is_paired_end: bool = True,
    overwrite_output: bool = False,
    art_platform: str = "illumina",
    timeout: float = 3600.0,
    auto_index_fasta: bool = True,
    variants_list: Optional[List[Dict[str, Any]]] = None,
    genomic_ranges: Optional[List[str]] = None
) -> Dict[str, Any]:

    effective_run_name = run_name if run_name else bb_utils.generate_unique_run_name(f"simulate_short_{art_platform}")
    logger.info(f"Initiating short read simulation for run: {effective_run_name}")

    temp_mutated_fasta_path: Optional[Path] = None
    temp_subset_fasta_path: Optional[Path] = None

    try:
        # --- Start of the code block that needed to be indented ---
        manifest_params = copy.deepcopy(command_params) if command_params is not None else {}
        manifest_params.update({
            "run_name": effective_run_name, "reference_fasta": reference_fasta, "depth": depth,
            "read_length": read_length, "art_profile": art_profile,
            "mean_fragment_length": mean_fragment_length, "std_dev_fragment_length": std_dev_fragment_length,
            "is_paired_end": is_paired_end, "art_platform": art_platform,
            "auto_index_fasta": auto_index_fasta, "timeout": timeout,
            "output_root_dir": str(output_root_dir), "overwrite_output": overwrite_output
        })
        id_prefix = manifest_params.get("id_prefix", f"simulated_{art_platform}_reads")
        manifest_params["id_prefix"] = id_prefix

        run_output_dir = bb_utils.prepare_run_output_dir(Path(output_root_dir), effective_run_name, overwrite_output)

        # Parameter validation
        if depth <= 0: raise bb_utils.BaseBuddyInputError(f"Sequencing depth must be a positive integer, got {depth}.")
        if read_length <= 0: raise bb_utils.BaseBuddyInputError(f"Read length must be a positive integer, got {read_length}.")
        if is_paired_end:
            if mean_fragment_length <= 0: raise bb_utils.BaseBuddyInputError(f"Mean fragment length must be positive for paired-end reads, got {mean_fragment_length}.")
            if std_dev_fragment_length < 0: raise bb_utils.BaseBuddyInputError(f"Standard deviation of fragment length cannot be negative, got {std_dev_fragment_length}.")

        # Tool paths
        art_exe_name = f"art_{art_platform}"
        # Command object for ART, this will also check for art_exe_name's existence via find_tool_path
        # BaseBuddyConfigError raised by Command constructor if art_exe_name not found will propagate up.
        art_cmd_builder = bb_utils.Command(art_exe_name)

        # samtools_exe_path might still be needed if samtools is used directly later for other things
        # For now, assume it's still needed for check_fasta_indexed etc.
        samtools_exe_path = bb_utils.find_tool_path("samtools")

        # Reference FASTA handling
        original_reference_path_obj = bb_utils.ensure_file_exists(reference_fasta, "Reference FASTA").resolve()
        bb_utils.check_fasta_indexed(original_reference_path_obj, samtools_exe_path, auto_index_if_missing=auto_index_fasta)
        manifest_params["original_reference_fasta"] = str(original_reference_path_obj)

        reference_after_variants = original_reference_path_obj # Start with original

        # 1. Apply variants if provided
        if variants_list and len(variants_list) > 0:
            logger.info(f"Applying {len(variants_list)} variants to reference FASTA {original_reference_path_obj}.")
            temp_mutated_fasta_path = run_output_dir / f"temp_mutated_ref_{effective_run_name}.fa"
            bb_utils.apply_variants_to_fasta(original_reference_path_obj, variants_list, temp_mutated_fasta_path)
            bb_utils.check_fasta_indexed(temp_mutated_fasta_path, samtools_exe_path, auto_index_if_missing=True)
            reference_after_variants = temp_mutated_fasta_path # This becomes the source for potential subsetting
            manifest_params["variants_applied_to_reference"] = True
            manifest_params["num_variants_applied"] = len(variants_list)
            manifest_params["mutated_reference_source"] = str(temp_mutated_fasta_path.name)
        else:
            manifest_params["variants_applied_to_reference"] = False

        current_reference_for_art = reference_after_variants # This will be used by ART, possibly after subsetting

        # 2. Extract genomic ranges if provided
        if genomic_ranges and len(genomic_ranges) > 0:
            logger.info(f"Extracting {len(genomic_ranges)} genomic ranges from {reference_after_variants}.")
            temp_subset_fasta_path = run_output_dir / f"temp_subset_ref_{effective_run_name}.fa"
            bb_utils.extract_ranges_to_fasta(reference_after_variants, genomic_ranges, temp_subset_fasta_path, samtools_path=samtools_exe_path)
            bb_utils.check_fasta_indexed(temp_subset_fasta_path, samtools_exe_path, auto_index_if_missing=True)
            current_reference_for_art = temp_subset_fasta_path # ART will use the subset
            manifest_params["genomic_ranges_applied"] = True
            manifest_params["num_genomic_ranges"] = len(genomic_ranges)
            manifest_params["genomic_ranges_requested"] = genomic_ranges
            manifest_params["subset_reference_fasta_path_for_art"] = str(temp_subset_fasta_path.name)
        else:
            manifest_params["genomic_ranges_applied"] = False

        art_output_prefix_path = run_output_dir / id_prefix

        # ART command construction using Command class
        art_cmd_builder.add_option("-ss", art_profile)
        art_cmd_builder.add_option("-i", str(current_reference_for_art))
        art_cmd_builder.add_option("-l", str(read_length))
        art_cmd_builder.add_option("-f", str(depth))
        art_cmd_builder.add_option("-o", str(art_output_prefix_path))

        if is_paired_end:
            art_cmd_builder.add_flag("-p")
            art_cmd_builder.add_option("-m", str(mean_fragment_length))
            art_cmd_builder.add_option("-s", str(std_dev_fragment_length))

        art_cmd_builder.add_option("-qs", manifest_params.get("quality_shift"))
        if is_paired_end: # qs2 only if paired-end
            art_cmd_builder.add_option("-qs2", manifest_params.get("quality_shift2"))
        art_cmd_builder.add_option("-rs", manifest_params.get("random_seed"))
        art_cmd_builder.add_flag("-na", condition=manifest_params.get("no_aln_output", False))

        art_cmd_parts = art_cmd_builder.get_command_parts()
        logger.info(f"Preparing to run {art_exe_name} with command: {' '.join(art_cmd_parts)}")
        bb_utils.run_external_cmd(art_cmd_parts, timeout_seconds=timeout, stream_output=True, cwd=run_output_dir)

        logger.info(f"{art_exe_name} simulation completed. Verifying output files...")
        output_files_manifest: List[Dict[str,str]] = []

        if is_paired_end:
            r1_path = art_output_prefix_path.with_name(art_output_prefix_path.name + "1.fq")
            r2_path = art_output_prefix_path.with_name(art_output_prefix_path.name + "2.fq")
            bb_utils.ensure_file_exists(r1_path, "Expected ART output R1 FASTQ")
            bb_utils.ensure_file_exists(r2_path, "Expected ART output R2 FASTQ")
            output_files_manifest.extend([
                {"name": "Simulated Reads (R1)", "path": r1_path.name, "type": "FASTQ"},
                {"name": "Simulated Reads (R2)", "path": r2_path.name, "type": "FASTQ"}
            ])
        else:
            r_path = art_output_prefix_path.with_suffix(".fq")
            bb_utils.ensure_file_exists(r_path, "Expected ART output FASTQ")
            output_files_manifest.append({"name": "Simulated Reads", "path": r_path.name, "type": "FASTQ"})

        if not manifest_params.get("no_aln_output", False):
            aln_path = art_output_prefix_path.with_suffix(".aln")
            if aln_path.exists():
                bb_utils.ensure_file_exists(aln_path, "ART Alignment ALN")
                output_files_manifest.append({"name": "ART Alignment ALN", "path": aln_path.name, "type": "ALN_ART"})

        manifest_path = run_output_dir / "manifest.json"
        bb_utils.write_run_manifest(
            manifest_path=manifest_path, run_name=effective_run_name,
            command_name=f"simulate_short_reads_{art_platform}", parameters=manifest_params,
            output_files=output_files_manifest, reference_genome_path=str(original_reference_path_obj) # Always report original ref
        )
        logger.info(f"Short read simulation successful for run '{effective_run_name}'. Output in '{run_output_dir}'.")

        return_dict = {
            "run_name": effective_run_name, "output_directory": str(run_output_dir.resolve()),
            "output_files": output_files_manifest, "manifest_path": str(manifest_path.resolve()),
            "reference_fasta_used": str(original_reference_path_obj) # Report original ref here too
        }
        if temp_mutated_fasta_path:
            return_dict["mutated_reference_info"] = {
                "mutated_fasta_path": str(temp_mutated_fasta_path.resolve()), # Path to the intermediate mutated FASTA
                "num_variants_applied": len(variants_list) if variants_list else 0
            }
        if temp_subset_fasta_path:
            return_dict["subset_reference_info"] = {
                "subset_fasta_path_used_for_art": str(temp_subset_fasta_path.resolve()), # Path to the subset FASTA actually used
                "genomic_ranges_requested": genomic_ranges
            }
        return return_dict
    
    finally:
        # Cleanup temp_mutated_fasta_path (variants applied)
        if temp_mutated_fasta_path and temp_mutated_fasta_path.exists():
            logger.debug(f"Cleaning up temporary mutated FASTA: {temp_mutated_fasta_path}")
            try:
                fai_path_mut = temp_mutated_fasta_path.with_suffix(temp_mutated_fasta_path.suffix + ".fai")
                if fai_path_mut.exists():
                    fai_path_mut.unlink(missing_ok=True)
                temp_mutated_fasta_path.unlink(missing_ok=True)
            except OSError as e:
                logger.warning(f"Could not delete temporary mutated FASTA file {temp_mutated_fasta_path} or its index: {e}")

        # Cleanup temp_subset_fasta_path (ranges extracted)
        if temp_subset_fasta_path and temp_subset_fasta_path.exists():
            logger.debug(f"Cleaning up temporary subset FASTA: {temp_subset_fasta_path}")
            try:
                fai_path_sub = temp_subset_fasta_path.with_suffix(temp_subset_fasta_path.suffix + ".fai")
                if fai_path_sub.exists():
                    fai_path_sub.unlink(missing_ok=True)
                temp_subset_fasta_path.unlink(missing_ok=True)
            except OSError as e:
                logger.warning(f"Could not delete temporary subset FASTA file {temp_subset_fasta_path} or its index: {e}")


def spike_variants(
    output_root_dir: Path,
    reference_fasta: str,
    input_bams: List[str], # List of path strings
    output_prefix_for_bam: str,
    run_name: Optional[str] = None,
    command_params: Optional[Dict[str, Any]] = None, # Expected to contain vaf, seed, picard_jar path
    snp_vcf_path: Optional[str] = None,
    indel_vcf_path: Optional[str] = None,
    overwrite_output: bool = False,
    auto_index_input_bam: bool = True, # For input_bams
    auto_index_fasta: bool = True # For reference_fasta
) -> Dict[str, Any]:

    effective_run_name = run_name if run_name else bb_utils.generate_unique_run_name("spike_variants_run")
    logger.info(f"Initiating variant spiking for run: {effective_run_name} on {len(input_bams)} BAMs.")

    if command_params is None:
        command_params = {}

    run_output_dir = bb_utils.prepare_run_output_dir(Path(output_root_dir), effective_run_name, overwrite_output)

    # Prepare manifest parameters (deep copy to avoid modifying caller's dict)
    manifest_params = copy.deepcopy(command_params)
    manifest_params.update({
        "run_name": effective_run_name,
        "reference_fasta": reference_fasta,
        "input_bams_original_paths": [str(Path(b).name) for b in input_bams], # Store basenames for manifest
        "num_input_bams": len(input_bams),
        "output_prefix_for_bam": output_prefix_for_bam,
        "snp_vcf_path": snp_vcf_path,
        "indel_vcf_path": indel_vcf_path,
        "overwrite_output": overwrite_output,
        "auto_index_input_bam": auto_index_input_bam,
        "auto_index_fasta": auto_index_fasta,
        "output_root_dir": str(output_root_dir)
        # vaf and seed are expected to be in command_params from the CLI
    })

    final_picard_jar_path_for_tools: Optional[str] = None
    try:
        addsnv_exe_path = bb_utils.find_tool_path("addsnv.py")
        addindel_exe_path = bb_utils.find_tool_path("addindel.py")
        samtools_exe_path = bb_utils.find_tool_path("samtools")

        picard_jar_cli_arg = command_params.get("picard_jar")
        picard_jar_env_var = os.environ.get('BAMSURGEON_PICARD_JAR')

        if picard_jar_cli_arg:
            p_path = Path(picard_jar_cli_arg)
            if not p_path.is_file(): raise bb_utils.BaseBuddyFileError(f"Specified Picard JAR path does not exist or is not a file: {picard_jar_cli_arg}")
            final_picard_jar_path_for_tools = str(p_path.resolve())
        elif picard_jar_env_var:
            p_path = Path(picard_jar_env_var)
            if not p_path.is_file(): raise bb_utils.BaseBuddyFileError(f"BAMSURGEON_PICARD_JAR environment variable points to non-existent file: {picard_jar_env_var}")
            final_picard_jar_path_for_tools = str(p_path.resolve())

        if not final_picard_jar_path_for_tools:
            # BaseBuddy should ensure this is provided to avoid bamsurgeon failing cryptically.
            raise bb_utils.BaseBuddyConfigError(
                "Picard JAR path is required by BAMSurgeon tools (addsnv.py, addindel.py). "
                "Please provide it via the '--picard-jar' option in 'command_params' (e.g. from CLI a --picard-jar option to spike) "
                "or set the BAMSURGEON_PICARD_JAR environment variable."
            )

    except bb_utils.BaseBuddyConfigError as e:
        logger.error(f"Tool configuration error for spiking: {e.details if hasattr(e, 'details') and e.details else str(e)}")
        raise

    ref_path_obj = bb_utils.ensure_file_exists(reference_fasta, "Reference FASTA").resolve()
    bb_utils.check_fasta_indexed(ref_path_obj, samtools_exe_path, auto_index_if_missing=auto_index_fasta)

    temp_snp_bed_file_path: Optional[Path] = None
    temp_indel_bed_file_path: Optional[Path] = None
    # VAF from CLI (via command_params), applied to all variants from VCFs for now
    cli_vaf = manifest_params.get("vaf", 0.5)

    try: # Wrap VCF parsing and temp file generation
        if snp_vcf_path:
            logger.info(f"Processing SNP VCF: {snp_vcf_path}")
            actual_snp_vcf_path = bb_utils.ensure_file_exists(snp_vcf_path, "SNP VCF").resolve()
            temp_snp_bed_file_path = run_output_dir / f"temp_snps_for_addsnv_{effective_run_name}.bed"
            snp_bed_lines = []
            with pysam.VariantFile(str(actual_snp_vcf_path)) as vcf:
                for rec in vcf.fetch():
                    if rec.alts and len(rec.ref) == 1 and len(rec.alts[0]) == 1: # Basic SNV
                        # addsnv.py BED format: chrom, 1-based start, 1-based end, vaf, altbase
                        snp_bed_lines.append(f"{rec.chrom}\t{rec.pos + 1}\t{rec.pos + 1}\t{cli_vaf}\t{rec.alts[0]}\n")
            if snp_bed_lines:
                with open(temp_snp_bed_file_path, "w") as f: f.writelines(snp_bed_lines)
                logger.info(f"Created temporary SNP BED: {temp_snp_bed_file_path} ({len(snp_bed_lines)} SNVs).")
                manifest_params["temp_snp_bed_file"] = str(temp_snp_bed_file_path.name)
            else:
                logger.info(f"No valid SNVs found in {snp_vcf_path} to process for addsnv.py."); temp_snp_bed_file_path = None

        if indel_vcf_path:
            logger.info(f"Processing Indel VCF: {indel_vcf_path}")
            actual_indel_vcf_path = bb_utils.ensure_file_exists(indel_vcf_path, "Indel VCF").resolve()
            temp_indel_bed_file_path = run_output_dir / f"temp_indels_for_addindel_{effective_run_name}.bed"
            indel_bed_lines = []
            with pysam.VariantFile(str(actual_indel_vcf_path)) as vcf:
                for rec in vcf.fetch():
                    if not rec.alts: continue
                    vcf_pos_1based, ref_len, alt_len = rec.pos + 1, len(rec.ref), len(rec.alts[0])
                    # addindel.py BED format: chrom, start, end, vaf, type ('INS' or 'DEL'), ins_sequence (for INS)
                    if ref_len == 1 and alt_len > 1: # Insertion: e.g. REF=A, ALT=AGGG -> INS=GGG
                        ins_seq = rec.alts[0][1:]
                        indel_bed_lines.append(f"{rec.chrom}\t{vcf_pos_1based}\t{vcf_pos_1based}\t{cli_vaf}\tINS\t{ins_seq}\n")
                    elif alt_len == 1 and ref_len > 1: # Deletion: e.g. REF=AGGG, ALT=A -> DEL=GGG
                        del_start_1based = vcf_pos_1based + 1
                        del_end_1based = vcf_pos_1based + (ref_len - alt_len)
                        indel_bed_lines.append(f"{rec.chrom}\t{del_start_1based}\t{del_end_1based}\t{cli_vaf}\tDEL\n")
            if indel_bed_lines:
                with open(temp_indel_bed_file_path, "w") as f: f.writelines(indel_bed_lines)
                logger.info(f"Created temporary Indel BED: {temp_indel_bed_file_path} ({len(indel_bed_lines)} Indels).")
                manifest_params["temp_indel_bed_file"] = str(temp_indel_bed_file_path.name)
            else:
                logger.info(f"No valid Indels found in {indel_vcf_path} to process for addindel.py."); temp_indel_bed_file_path = None
    except Exception as e:
        if temp_snp_bed_file_path and temp_snp_bed_file_path.exists(): temp_snp_bed_file_path.unlink(missing_ok=True)
        if temp_indel_bed_file_path and temp_indel_bed_file_path.exists(): temp_indel_bed_file_path.unlink(missing_ok=True)
        logger.error(f"Error during VCF parsing or temp BED file creation: {e}")
        raise bb_utils.BaseBuddyFileError(f"Failed during VCF parsing for spike command.", details=str(e))

    all_results_per_bam: List[Dict[str, Any]] = []
    all_errors_per_bam: List[Dict[str, str]] = []
    overall_output_files_manifest: List[Dict[str,str]] = []
    processed_bam_count = 0
    seed_for_tool = manifest_params.get("seed", 0)

    for input_bam_idx, input_bam_str_path in enumerate(input_bams):
        input_bam_path_obj = bb_utils.ensure_file_exists(input_bam_str_path, f"Input BAM {input_bam_str_path}").resolve()
        bb_utils.check_bam_indexed(input_bam_path_obj, samtools_exe_path, auto_index_input_bam)
        input_bam_stem = input_bam_path_obj.stem
        logger.info(f"Processing input BAM ({input_bam_idx+1}/{len(input_bams)}): {input_bam_path_obj.name}")

        current_bam_to_process: Path = input_bam_path_obj
        made_changes_to_this_bam = False
        intermediate_files_for_this_input_bam: List[Path] = []

        snp_log_vcf_path = run_output_dir / f"{output_prefix_for_bam}_{input_bam_stem}_addsnv_log.vcf"
        indel_log_vcf_path = run_output_dir / f"{output_prefix_for_bam}_{input_bam_stem}_addindel_log.vcf"

        try:
            if temp_snp_bed_file_path and temp_snp_bed_file_path.exists():
                bam_after_snps_path = run_output_dir / f"{output_prefix_for_bam}_{input_bam_stem}_TEMP_with_snps.bam"
                intermediate_files_for_this_input_bam.append(bam_after_snps_path)
                logger.info(f"Adding SNPs to {current_bam_to_process.name} -> {bam_after_snps_path.name}")

                cmd_builder = bb_utils.Command("python").add_option(None, addsnv_exe_path)
                cmd_builder.add_option("-v", str(temp_snp_bed_file_path))
                cmd_builder.add_option("-f", str(current_bam_to_process))
                cmd_builder.add_option("-r", str(ref_path_obj))
                cmd_builder.add_option("-o", str(bam_after_snps_path))
                cmd_builder.add_option("-m", str(cli_vaf))
                cmd_builder.add_option("--seed", str(seed_for_tool))
                cmd_builder.add_option("--tmpdir", str(run_output_dir / f"addsnv_tmp_{input_bam_stem}"))
                cmd_builder.add_option("--picardjar", final_picard_jar_path_for_tools)
                cmd_builder.add_option("--vcf", str(snp_log_vcf_path))

                bb_utils.run_external_cmd(cmd_builder.get_command_parts(), stream_output=True, cwd=run_output_dir)
                current_bam_to_process = bam_after_snps_path
                made_changes_to_this_bam = True

            if temp_indel_bed_file_path and temp_indel_bed_file_path.exists():
                name_suffix = "_TEMP_with_indels.bam" if not made_changes_to_this_bam else "_TEMP_with_snps_indels.bam"
                bam_after_indels_path = run_output_dir / f"{output_prefix_for_bam}_{input_bam_stem}{name_suffix}"
                intermediate_files_for_this_input_bam.append(bam_after_indels_path)
                logger.info(f"Adding Indels to {current_bam_to_process.name} -> {bam_after_indels_path.name}")

                cmd_builder = bb_utils.Command("python").add_option(None, addindel_exe_path)
                cmd_builder.add_option("-v", str(temp_indel_bed_file_path))
                cmd_builder.add_option("-f", str(current_bam_to_process))
                cmd_builder.add_option("-r", str(ref_path_obj))
                cmd_builder.add_option("-o", str(bam_after_indels_path))
                cmd_builder.add_option("-m", str(cli_vaf)) # VAF for indels
                cmd_builder.add_option("--seed", str(seed_for_tool))
                cmd_builder.add_option("--tmpdir", str(run_output_dir / f"addindel_tmp_{input_bam_stem}"))
                cmd_builder.add_option("--picardjar", final_picard_jar_path_for_tools)
                cmd_builder.add_option("--vcf", str(indel_log_vcf_path))

                bb_utils.run_external_cmd(cmd_builder.get_command_parts(), stream_output=True, cwd=run_output_dir)
                current_bam_to_process = bam_after_indels_path
                made_changes_to_this_bam = True

            if not made_changes_to_this_bam:
                logger.info(f"No variants specified or found valid to apply to {input_bam_path_obj.name}.")
                all_results_per_bam.append({
                    "source_bam": str(input_bam_path_obj.resolve()), "status": "skipped_no_variants_to_apply",
                    "output_bam": str(input_bam_path_obj.resolve())
                })
                overall_output_files_manifest.append({"name": f"Unchanged BAM ({input_bam_stem})", "path": input_bam_path_obj.name, "type": "BAM"})
                # Copy original BAM and index to output if no changes, or let user handle? For now, just note it.
                continue

            final_sorted_bam_path = run_output_dir / f"{output_prefix_for_bam}_{input_bam_stem}_final_sorted.bam"
            logger.info(f"Sorting BAM: {current_bam_to_process.name} -> {final_sorted_bam_path.name}")
            cmd_sort = bb_utils.Command("samtools").add_option(None, "sort").add_option("-o", str(final_sorted_bam_path)).add_option(None, str(current_bam_to_process))
            bb_utils.run_external_cmd(cmd_sort.get_command_parts(), cwd=run_output_dir)

            bb_utils.check_bam_indexed(final_sorted_bam_path, samtools_exe_path, auto_index_if_missing=True)
            final_bam_obj = bb_utils.ensure_file_exists(final_sorted_bam_path, "Final sorted BAM")
            final_bai_obj = Path(str(final_bam_obj) + ".bai")
            if not final_bai_obj.exists(): final_bai_obj = final_bam_obj.with_suffix(".bai")

            igv_session_path = run_output_dir / f"igv_session_spike_{input_bam_stem}.xml"
            igv_tracks = [{"name": f"Spiked BAM ({input_bam_stem})", "path": final_bam_obj.name, "type": "alignment"}]
            if snp_log_vcf_path.exists(): igv_tracks.append({"name": f"SNV Log ({input_bam_stem})", "path": snp_log_vcf_path.name, "type": "variant"})
            if indel_log_vcf_path.exists(): igv_tracks.append({"name": f"Indel Log ({input_bam_stem})", "path": indel_log_vcf_path.name, "type": "variant"})
            bb_utils.generate_igv_session_xml(igv_session_path, str(ref_path_obj), igv_tracks)

            current_bam_manifest_entries = [
                {"name": f"Final Spiked BAM ({input_bam_stem})", "path": final_bam_obj.name, "type": "BAM"},
                {"name": f"BAM Index ({input_bam_stem})", "path": final_bai_obj.name if final_bai_obj.exists() else "index_not_found.bai", "type": "BAI"},
                {"name": f"IGV Session ({input_bam_stem})", "path": igv_session_path.name, "type": "IGV_SESSION"}
            ]
            if snp_log_vcf_path.exists(): current_bam_manifest_entries.append({"name": f"AddSNV Log VCF ({input_bam_stem})", "path": snp_log_vcf_path.name, "type": "VCF_LOG"})
            if indel_log_vcf_path.exists(): current_bam_manifest_entries.append({"name": f"AddIndel Log VCF ({input_bam_stem})", "path": indel_log_vcf_path.name, "type": "VCF_LOG"})
            overall_output_files_manifest.extend(current_bam_manifest_entries)

            all_results_per_bam.append({
                "source_bam": str(input_bam_path_obj.resolve()), "status": "success",
                "output_bam": str(final_bam_obj.resolve()),
                "output_bam_index": str(final_bai_obj.resolve()) if final_bai_obj.exists() else None,
                "igv_session_file": str(igv_session_path.resolve())
            })
            processed_bam_count += 1
            logger.info(f"Successfully processed {input_bam_path_obj.name} -> {final_bam_obj.name}")

        except Exception as e_bam:
            logger.error(f"Error processing BAM {input_bam_path_obj.name}: {type(e_bam).__name__} - {e_bam}")
            details = e_bam.details if hasattr(e_bam, 'details') and e_bam.details else traceback.format_exc()
            all_errors_per_bam.append({"source_bam": str(input_bam_path_obj.resolve()), "error": str(e_bam), "details": details})
        finally:
            logger.debug(f"Cleaning up intermediate files for {input_bam_path_obj.name}: {intermediate_files_for_this_input_bam}")
            for f_path in intermediate_files_for_this_input_bam:
                if f_path.exists(): f_path.unlink(missing_ok=True)
                # Attempt to remove .bai for temp bams too
                bai_path_temp = f_path.with_suffix(f_path.suffix + ".bai")
                if bai_path_temp.exists(): bai_path_temp.unlink(missing_ok=True)
                bai_path_temp_alt = f_path.with_suffix(".bai")
                if bai_path_temp_alt.exists(): bai_path_temp_alt.unlink(missing_ok=True)

    try:
        if temp_snp_bed_file_path and temp_snp_bed_file_path.exists():
            logger.debug(f"Cleaning up temporary SNP BED file: {temp_snp_bed_file_path}")
            temp_snp_bed_file_path.unlink(missing_ok=True)
        if temp_indel_bed_file_path and temp_indel_bed_file_path.exists():
            logger.debug(f"Cleaning up temporary Indel BED file: {temp_indel_bed_file_path}")
            temp_indel_bed_file_path.unlink(missing_ok=True)
    except Exception as e_cleanup:
        logger.warning(f"Error during final cleanup of temporary BED files: {e_cleanup}")

    manifest_params["num_bams_processed_successfully"] = processed_bam_count
    manifest_params["num_bams_failed"] = len(all_errors_per_bam)
    status_detail = f"Completed spiking. Processed BAMs: {processed_bam_count}, Failed: {len(all_errors_per_bam)}."

    no_variants_applied_at_all = True
    if temp_snp_bed_file_path and temp_snp_bed_file_path.exists() and processed_bam_count > 0 : no_variants_applied_at_all = False
    if temp_indel_bed_file_path and temp_indel_bed_file_path.exists() and processed_bam_count > 0 : no_variants_applied_at_all = False

    if (snp_vcf_path or indel_vcf_path) and no_variants_applied_at_all and processed_bam_count == 0 and len(input_bams) > 0 :
         status_detail += " No valid SNPs or Indels found in provided VCFs to apply, or all BAM processing failed early."
    elif not snp_vcf_path and not indel_vcf_path:
        status_detail += " No SNP or Indel VCFs provided."

    manifest_params["status_detail"] = status_detail

    main_manifest_path = run_output_dir / "manifest_spike_run.json"
    bb_utils.write_run_manifest(
        manifest_path=main_manifest_path, run_name=effective_run_name,
        command_name="spike_variants", parameters=manifest_params,
        output_files=overall_output_files_manifest, reference_genome_path=str(ref_path_obj)
    )

    return {
        "run_name": effective_run_name, "output_directory": str(run_output_dir.resolve()),
        "results_per_bam": all_results_per_bam, "errors_per_bam": all_errors_per_bam,
        "manifest_path": str(main_manifest_path.resolve()), "reference_fasta_used": str(ref_path_obj)
    }

def simulate_long(
    output_root_dir: Path,
    reference_fasta: str,
    depth: int,
    model: str,
    run_name: Optional[str] = None,
    command_params: Optional[Dict[str, Any]] = None,
    overwrite_output: bool = False,
    auto_index_fasta: bool = True,
    timeout: float = 3600.0,
    variants_list: Optional[List[Dict[str, Any]]] = None,
    genomic_ranges: Optional[List[str]] = None
) -> Dict[str, Any]:
    effective_run_name = run_name if run_name else bb_utils.generate_unique_run_name(f"simulate_long_nanosim")
    logger.info(f"Initiating long read simulation for run: {effective_run_name}")

    temp_mutated_fasta_path: Optional[Path] = None
    temp_subset_fasta_path: Optional[Path] = None

    try:
        manifest_params = copy.deepcopy(command_params) if command_params is not None else {}
        manifest_params.update({
            "run_name": effective_run_name, "reference_fasta": reference_fasta, "depth": depth,
            "model": model, "overwrite_output": overwrite_output,
            "auto_index_fasta": auto_index_fasta, "timeout": timeout,
            "output_root_dir": str(output_root_dir)
        })
        nanosim_internal_out_prefix = "nanosim_reads"

        run_output_dir = bb_utils.prepare_run_output_dir(Path(output_root_dir), effective_run_name, overwrite_output)

        # Parameter validation
        if depth <= 0 and not (manifest_params.get("num_reads") and int(manifest_params["num_reads"]) > 0) : # Check depth if num_reads is not primary
             raise bb_utils.BaseBuddyInputError(f"Sequencing depth must be a positive integer if num_reads is not specified, got {depth}.")
        if not model: raise bb_utils.BaseBuddyInputError("NanoSim model must be specified.")

        # Tool Paths
        # Command object for NanoSim-h, also checks for existence.
        # BaseBuddyConfigError raised by Command constructor if nanosim-h not found will propagate.
        nanosim_cmd_builder = bb_utils.Command("nanosim-h")

        samtools_exe_path = bb_utils.find_tool_path("samtools") # Assuming still needed

        # Reference FASTA handling
        original_reference_path_obj = bb_utils.ensure_file_exists(reference_fasta, "Reference FASTA").resolve()
        bb_utils.check_fasta_indexed(original_reference_path_obj, samtools_exe_path, auto_index_if_missing=auto_index_fasta)
        manifest_params["original_reference_fasta"] = str(original_reference_path_obj)

        reference_after_variants = original_reference_path_obj # Start with original

        # 1. Apply variants if provided
        if variants_list and len(variants_list) > 0:
            logger.info(f"Applying {len(variants_list)} variants to reference FASTA {original_reference_path_obj} for NanoSim-h.")
            temp_mutated_fasta_path = run_output_dir / f"temp_mutated_ref_nanosim_{effective_run_name}.fa"
            bb_utils.apply_variants_to_fasta(original_reference_path_obj, variants_list, temp_mutated_fasta_path)
            bb_utils.check_fasta_indexed(temp_mutated_fasta_path, samtools_exe_path, auto_index_if_missing=True)
            reference_after_variants = temp_mutated_fasta_path
            manifest_params["variants_applied_to_reference"] = True
            manifest_params["num_variants_applied"] = len(variants_list)
            manifest_params["mutated_reference_source"] = str(temp_mutated_fasta_path.name)
        else:
            manifest_params["variants_applied_to_reference"] = False

        current_reference_for_nanosim = reference_after_variants # Will be used by NanoSim, possibly after subsetting

        # 2. Extract genomic ranges if provided
        if genomic_ranges and len(genomic_ranges) > 0:
            logger.info(f"Extracting {len(genomic_ranges)} genomic ranges from {reference_after_variants} for NanoSim-h.")
            temp_subset_fasta_path = run_output_dir / f"temp_subset_ref_nanosim_{effective_run_name}.fa"
            bb_utils.extract_ranges_to_fasta(reference_after_variants, genomic_ranges, temp_subset_fasta_path, samtools_path=samtools_exe_path)
            bb_utils.check_fasta_indexed(temp_subset_fasta_path, samtools_exe_path, auto_index_if_missing=True)
            current_reference_for_nanosim = temp_subset_fasta_path # NanoSim will use the subset
            manifest_params["genomic_ranges_applied"] = True
            manifest_params["num_genomic_ranges"] = len(genomic_ranges)
            manifest_params["genomic_ranges_requested"] = genomic_ranges
            manifest_params["subset_reference_fasta_path_for_nanosim"] = str(temp_subset_fasta_path.name)
        else:
            manifest_params["genomic_ranges_applied"] = False

        nanosim_output_prefix_path = run_output_dir / nanosim_internal_out_prefix

        # NanoSim-h command construction using Command class
        nanosim_cmd_builder.add_option(None, "simulate") # "simulate" is a sub-command for nanosim-h
        nanosim_cmd_builder.add_option("-r", str(current_reference_for_nanosim))
        nanosim_cmd_builder.add_option("-m", model)
        nanosim_cmd_builder.add_option("-o", str(nanosim_output_prefix_path))

        num_reads_param = manifest_params.get("num_reads")
        # final_cmd_parts_for_nanosim: Optional[List[str]] = None # Not needed, builder holds parts

        if num_reads_param is not None:
            try:
                num_reads_val = int(num_reads_param)
                if num_reads_val > 0:
                    nanosim_cmd_builder.add_option("-N", str(num_reads_val))
                    if "depth" in manifest_params:
                         logger.info("Using num_reads for NanoSim; 'depth' parameter will be ignored by NanoSim if also set.")
                elif depth > 0:
                    nanosim_cmd_builder.add_option("-c", str(depth))
                    logger.warning(f"num_reads parameter '{num_reads_param}' is not a positive integer. Falling back to depth '{depth}'.")
                else:
                    raise bb_utils.BaseBuddyInputError("Either 'depth' or 'num_reads' must be a positive integer for NanoSim when num_reads is specified but invalid.")
            except ValueError:
                if depth > 0:
                    nanosim_cmd_builder.add_option("-c", str(depth))
                    logger.warning(f"num_reads parameter '{num_reads_param}' is not a valid integer. Falling back to depth '{depth}'.")
                else:
                    raise bb_utils.BaseBuddyInputError("num_reads is not a valid integer and depth is not positive for NanoSim.")
        elif depth > 0:
            nanosim_cmd_builder.add_option("-c", str(depth))
        else:
             raise bb_utils.BaseBuddyInputError("Either 'depth' or 'num_reads' must be a positive integer for NanoSim.")

        final_cmd_parts_for_nanosim = nanosim_cmd_builder.get_command_parts()
        logger.info(f"Running NanoSim-h with command: {' '.join(final_cmd_parts_for_nanosim)}")
        bb_utils.run_external_cmd(final_cmd_parts_for_nanosim, timeout_seconds=timeout, stream_output=True, cwd=run_output_dir)

        output_files_manifest: List[Dict[str,str]] = []
        simulated_fastq = nanosim_output_prefix_path.with_name(nanosim_output_prefix_path.name + "_aligned_reads.fastq")
        if simulated_fastq.exists():
            bb_utils.ensure_file_exists(simulated_fastq, "Simulated long reads FASTQ")
            output_files_manifest.append({"name": "Simulated Long Reads (FASTQ)", "path": simulated_fastq.name, "type": "FASTQ"})
        else:
            simulated_fasta = nanosim_output_prefix_path.with_name(nanosim_output_prefix_path.name + "_aligned_reads.fasta")
            if simulated_fasta.exists():
                bb_utils.ensure_file_exists(simulated_fasta, "Simulated long reads FASTA")
                output_files_manifest.append({"name": "Simulated Long Reads (FASTA)", "path": simulated_fasta.name, "type": "FASTA"})
            else:
                logger.warning(f"Could not find primary NanoSim output FASTQ/FASTA: {simulated_fastq.name} or {simulated_fasta.name}")

        log_file = nanosim_output_prefix_path.with_name(nanosim_output_prefix_path.name + "_log.txt")
        if log_file.exists():
            output_files_manifest.append({"name": "NanoSim Log", "path": log_file.name, "type": "LOG"})

        manifest_path = run_output_dir / "manifest.json"
        bb_utils.write_run_manifest(
            manifest_path=manifest_path, run_name=effective_run_name, command_name="simulate_long_nanosim",
            parameters=manifest_params, output_files=output_files_manifest,
            reference_genome_path=str(original_reference_path_obj) # Always report original
        )
        logger.info(f"Long read simulation successful for run '{effective_run_name}'. Output in '{run_output_dir}'.")

        return_dict = {
            "run_name": effective_run_name,
            "output_directory": str(run_output_dir.resolve()),
            "output_files": output_files_manifest,
            "manifest_path": str(manifest_path.resolve()),
            "reference_fasta_used": str(original_reference_path_obj) # Always report original
        }
        if temp_mutated_fasta_path:
            return_dict["mutated_reference_info"] = {
                "mutated_fasta_path": str(temp_mutated_fasta_path.resolve()),
                "num_variants_applied": len(variants_list) if variants_list else 0
            }
        if temp_subset_fasta_path:
            return_dict["subset_reference_info"] = {
                "subset_fasta_path_used_for_nanosim": str(temp_subset_fasta_path.resolve()),
                "genomic_ranges_requested": genomic_ranges
            }
        return return_dict

    finally:
        # Cleanup temp_mutated_fasta_path (variants applied)
        if temp_mutated_fasta_path and temp_mutated_fasta_path.exists():
            logger.debug(f"Cleaning up temporary mutated FASTA for NanoSim: {temp_mutated_fasta_path}")
            try:
                fai_path_mut = temp_mutated_fasta_path.with_suffix(temp_mutated_fasta_path.suffix + ".fai")
                if fai_path_mut.exists():
                    fai_path_mut.unlink(missing_ok=True)
                temp_mutated_fasta_path.unlink(missing_ok=True)
            except OSError as e:
                logger.warning(f"Could not delete temporary mutated FASTA file {temp_mutated_fasta_path} or its index: {e}")

        # Cleanup temp_subset_fasta_path (ranges extracted)
        if temp_subset_fasta_path and temp_subset_fasta_path.exists():
            logger.debug(f"Cleaning up temporary subset FASTA for NanoSim: {temp_subset_fasta_path}")
            try:
                fai_path_sub = temp_subset_fasta_path.with_suffix(temp_subset_fasta_path.suffix + ".fai")
                if fai_path_sub.exists():
                    fai_path_sub.unlink(missing_ok=True)
                temp_subset_fasta_path.unlink(missing_ok=True)
            except OSError as e:
                logger.warning(f"Could not delete temporary subset FASTA file {temp_subset_fasta_path} or its index: {e}")


def apply_signature_to_fasta(
    output_root_dir: Path,
    run_name: Optional[str], # If None, will be auto-generated
    command_params: Dict[str, Any], # For manifest, should contain user-facing params
    input_fasta_path_str: str,
    output_fasta_name: str, # e.g., "mutated_by_sbs1.fa"
    signature_id_or_path: str, # e.g., "SBS1" or "path/to/custom_sig.tsv"
    num_mutations: int,
    target_regions: Optional[List[str]] = None, # Placeholder, NOT IMPLEMENTED in this version
    overwrite_output: bool = False,
    auto_index_input_fasta: bool = True # For the input FASTA
) -> Dict[str, Any]:
    """
    Applies a specified SBS signature to an input FASTA file.
    The specified number of mutations will be distributed across the contigs.
    """
    _logger = logging.getLogger(__name__) # Local logger for this function

    effective_run_name = run_name if run_name else bb_utils.generate_unique_run_name(
        f"applysig_{Path(output_fasta_name).stem}"
    )
    _logger.info(f"Initiating 'apply_signature_to_fasta' for run: {effective_run_name}")

    run_output_dir = bb_utils.prepare_run_output_dir(
        Path(output_root_dir), effective_run_name, overwrite_output
    )

    actual_command_params = copy.deepcopy(command_params) # Base for manifest
    actual_command_params.update({ # Add/override with specific operational params
        "run_name": effective_run_name,
        "input_fasta_path": input_fasta_path_str,
        "output_fasta_name": output_fasta_name,
        "signature_id_or_path": signature_id_or_path,
        "num_mutations_requested": num_mutations,
        "target_regions": target_regions if target_regions else "whole_genome", # For manifest
        "overwrite_output": overwrite_output,
        "auto_index_input_fasta": auto_index_input_fasta
    })

    if num_mutations <= 0:
        raise bb_utils.BaseBuddyInputError("Number of mutations must be a positive integer.")
    if target_regions:
        _logger.warning("target_regions parameter is provided but not yet implemented in this version. Mutations will be applied genome-wide.")

    input_fasta_path = bb_utils.ensure_file_exists(input_fasta_path_str, "Input FASTA for signature application")
    samtools_exe_path = bb_utils.find_tool_path("samtools")
    bb_utils.check_fasta_indexed(input_fasta_path, samtools_exe_path, auto_index_if_missing=auto_index_input_fasta)

    output_modified_fasta_path = run_output_dir / output_fasta_name
    if output_modified_fasta_path.exists() and not overwrite_output:
        raise bb_utils.BaseBuddyFileError(
            f"Output FASTA file '{output_modified_fasta_path}' already exists and overwrite is False."
        )

    # Load Signature Profile
    bundled_signatures_dir = Path(sig_utils.__file__).resolve().parent / "data" / "signatures"
    _logger.debug(f"Bundled signatures directory: {bundled_signatures_dir}")

    sig_file_to_parse: Optional[Path] = None
    active_signature_name: str = "" # Must be set
    signature_profile_to_apply: Dict[str, float] = {}

    # Determine if signature_id_or_path is a bundled ID or a direct path
    # And load the appropriate signature profile.
    # This logic assumes SBS for now, as apply_sbs_mutations_to_sequence is SBS-specific.

    potential_sig_id = signature_id_or_path.upper() # Normalize for prefix checking
    sig_type_prefix = None
    if potential_sig_id.startswith("SBS"):
        sig_type_prefix = "sbs"
    elif potential_sig_id.startswith("DBS"):
        sig_type_prefix = "dbs"
    elif potential_sig_id.startswith("ID"):
        sig_type_prefix = "id"

    if sig_type_prefix and not Path(signature_id_or_path).is_file(): # Check if it's NOT a file path AND has a known prefix
        _logger.info(f"Treating '{signature_id_or_path}' as a bundled signature ID of type '{sig_type_prefix}'.")
        master_file_path = sig_utils.get_bundled_signature_master_file_path(sig_type_prefix, bundled_signatures_dir)
        if not master_file_path or not master_file_path.exists():
            # Fallback to old method for sbs1.tsv / sbs5.tsv if master files are not yet populated
            if sig_type_prefix == "sbs" and (signature_id_or_path.lower() == "sbs1" or signature_id_or_path.lower() == "sbs5"):
                 _logger.warning(f"Master file for '{sig_type_prefix}' not found at '{master_file_path}'. Attempting to load '{signature_id_or_path}.tsv' directly.")
                 sig_file_to_parse = bundled_signatures_dir / f"{signature_id_or_path.lower()}.tsv"
                 if not sig_file_to_parse.exists():
                     raise bb_utils.BaseBuddyInputError(f"Fallback bundled signature file '{sig_file_to_parse}' also not found.")
                 all_signatures_in_file = sig_utils.parse_signature_matrix_tsv(sig_file_to_parse)
                 active_signature_name = signature_id_or_path
            else:
                raise bb_utils.BaseBuddyInputError(f"Bundled master signature file for type '{sig_type_prefix}' not found at '{master_file_path}'.")
        else:
            all_signatures_in_file = sig_utils.parse_signature_matrix_tsv(master_file_path)
            active_signature_name = signature_id_or_path # The specific ID like "SBS1"

        if active_signature_name not in all_signatures_in_file:
            raise bb_utils.BaseBuddyInputError(f"Signature ID '{active_signature_name}' not found in master file '{sig_file_to_parse or master_file_path}'. Available: {list(all_signatures_in_file.keys())}")
        signature_profile_to_apply = all_signatures_in_file[active_signature_name]

    else: # Treat as a direct file path to a custom signature matrix
        _logger.info(f"Treating '{signature_id_or_path}' as a direct path to a custom signature matrix file.")
        custom_file_path = bb_utils.ensure_file_exists(signature_id_or_path, "Custom signature matrix file")
        all_signatures_in_file = sig_utils.parse_signature_matrix_tsv(custom_file_path)

        # Determine which signature to use from the custom file
        file_stem_name = Path(signature_id_or_path).stem
        if file_stem_name in all_signatures_in_file:
            active_signature_name = file_stem_name
        elif len(all_signatures_in_file) == 1:
            active_signature_name = list(all_signatures_in_file.keys())[0]
            _logger.warning(f"Using the only signature '{active_signature_name}' found in custom file '{custom_file_path}'.")
        else:
            raise bb_utils.BaseBuddyInputError(
                f"Could not determine which signature to use from custom file '{custom_file_path}'. "
                f"File stem '{file_stem_name}' not found and multiple signatures exist: {list(all_signatures_in_file.keys())}. "
                "Hint: Name your custom file like 'MySignature.tsv' if 'MySignature' is the column header, or ensure only one signature is in the TSV."
            )
        signature_profile_to_apply = all_signatures_in_file[active_signature_name]

    actual_command_params["active_signature_name_used"] = active_signature_name
    _logger.info(f"Successfully loaded signature profile for '{active_signature_name}'.")

    # Process FASTA and Apply Mutations
    _logger.info(f"Applying up to {num_mutations} mutations using signature '{active_signature_name}' to '{input_fasta_path.name}' -> '{output_modified_fasta_path.name}'.")

    total_applied_mutations = 0
    try:
        with pysam.FastaFile(str(input_fasta_path)) as fasta_in, \
             open(output_modified_fasta_path, 'w') as fasta_out:

            if not fasta_in.references:
                raise bb_utils.BaseBuddyFileError(f"Input FASTA file '{input_fasta_path}' contains no sequences/references.")

            for contig_name in fasta_in.references:
                if total_applied_mutations >= num_mutations:
                    original_seq = fasta_in.fetch(contig_name)
                    fasta_out.write(f">{contig_name}\n{original_seq}\n")
                    continue

                original_seq = fasta_in.fetch(contig_name)
                _logger.debug(f"Processing contig '{contig_name}' (length: {len(original_seq)}).")

                # Currently, application logic is SBS-specific.
                # Future: if sig_type_prefix == "dbs", call dbs_context_and_apply, etc.
                if not active_signature_name.lower().startswith("sbs"): # Basic check
                    _logger.warning(f"Signature '{active_signature_name}' does not seem to be an SBS signature. SBS-specific mutation application will be attempted.")

                available_contexts = sig_utils.get_sbs_mutation_contexts(original_seq)
                if not available_contexts:
                    _logger.warning(f"No mutable SBS contexts found in contig '{contig_name}'. Writing original sequence for this contig.")
                    fasta_out.write(f">{contig_name}\n{original_seq}\n")
                    continue

                mutations_to_attempt_on_contig = num_mutations - total_applied_mutations

                _logger.debug(f"Attempting to apply {mutations_to_attempt_on_contig} mutations to '{contig_name}'. "
                             f"Available contexts: {len(available_contexts)}")

                mutated_seq, applied_on_contig = sig_utils.apply_sbs_mutations_to_sequence(
                    original_seq,
                    mutations_to_attempt_on_contig,
                    signature_profile_to_apply,
                    available_contexts
                )
                fasta_out.write(f">{contig_name}\n{mutated_seq}\n")
                total_applied_mutations += applied_on_contig
                _logger.info(f"Applied {applied_on_contig} mutations to contig '{contig_name}'. Total applied so far: {total_applied_mutations}.")

    except Exception as e:
        _logger.error(f"Error during FASTA processing or mutation application: {e}")
        if output_modified_fasta_path.exists():
            try: output_modified_fasta_path.unlink()
            except OSError: _logger.warning(f"Could not remove incomplete output FASTA: {output_modified_fasta_path}")
        raise bb_utils.BaseBuddyToolError(f"Failed to apply signature to FASTA: {e}", command=["apply_signature_to_fastaInternal"], stderr=str(e))


    _logger.info(f"Finished applying mutations. Total applied: {total_applied_mutations} / Requested: {num_mutations}.")

    bb_utils.check_fasta_indexed(output_modified_fasta_path, samtools_exe_path, auto_index_if_missing=True)
    output_index_path = output_modified_fasta_path.with_suffix(output_modified_fasta_path.suffix + ".fai")

    actual_command_params["num_mutations_applied"] = total_applied_mutations
    output_files_manifest = [
        {"name": "Mutated FASTA by signature", "path": output_modified_fasta_path.name, "type": "FASTA"},
    ]
    if output_index_path.exists():
        output_files_manifest.append({"name": "Mutated FASTA Index", "path": output_index_path.name, "type": "FASTA_INDEX"})

    manifest_path = run_output_dir / "manifest.json"
    bb_utils.write_run_manifest(
        manifest_path=manifest_path,
        run_name=effective_run_name,
        command_name="apply_signature_to_fasta",
        parameters=actual_command_params,
        output_files=output_files_manifest,
        reference_genome_path=str(input_fasta_path.resolve())
    )

    return {
        "run_name": effective_run_name,
        "output_directory": str(run_output_dir.resolve()),
        "output_modified_fasta_path": str(output_modified_fasta_path.resolve()),
        "output_modified_fasta_index_path": str(output_index_path.resolve()) if output_index_path.exists() else None,
        "manifest_path": str(manifest_path.resolve()),
        "num_mutations_applied": total_applied_mutations,
        "num_mutations_requested": num_mutations,
        "signature_used": active_signature_name,
        "original_input_fasta": str(input_fasta_path.resolve())
    }


def simulate_signatures(
    output_root_dir: Path,
    reference_input: str,
    sig_type: str,
    num_mutations: int,
    sample_id: str = "Sample1",
    run_name: Optional[str] = None,
    command_params: Optional[Dict[str, Any]] = None,
    exome: bool = False,
    overwrite_output: bool = False,
    auto_index_fasta: bool = True
) -> Dict[str, Any]:

    effective_run_name = run_name if run_name else bb_utils.generate_unique_run_name(f"simulate_signatures_{sig_type}")
    logger.info(f"Initiating mutational signature simulation for run: {effective_run_name}")

    manifest_params = copy.deepcopy(command_params) if command_params is not None else {}
    manifest_params.update({
        "run_name": effective_run_name, "reference_input": reference_input, "sig_type": sig_type,
        "num_mutations": num_mutations, "sample_id": sample_id, "exome": exome,
        "overwrite_output": overwrite_output, "auto_index_fasta": auto_index_fasta,
        "output_root_dir": str(output_root_dir)
    })

    run_output_dir = bb_utils.prepare_run_output_dir(Path(output_root_dir), effective_run_name, overwrite_output)

    genome_build_arg: Optional[str] = None
    custom_genome_path_arg: Optional[str] = None
    ref_display_for_manifest = reference_input
    ref_path_for_igv: Optional[str] = None # Not used by this runner for IGV session

    if Path(reference_input).is_file() or \
       any(reference_input.endswith(ext) for ext in [".fa", ".fasta", ".fna"]) or \
       "/" in reference_input or "\\" in reference_input:

        logger.info(f"Treating reference_input '{reference_input}' as a custom FASTA path.")
        ref_path_obj = bb_utils.ensure_file_exists(reference_input, "Reference FASTA").resolve()
        if auto_index_fasta:
            samtools_exe_path = bb_utils.find_tool_path("samtools")
            bb_utils.check_fasta_indexed(ref_path_obj, samtools_exe_path, auto_index_if_missing=True)
        custom_genome_path_arg = str(ref_path_obj)
        ref_display_for_manifest = str(ref_path_obj)
        genome_build_name = Path(reference_input).stem.split('.')[0]
        if genome_build_name in ["GRCh37", "GRCh38", "hg19", "hg38", "mm10"]:
             genome_build_arg = genome_build_name
        else:
            genome_build_arg = manifest_params.get("genome_build_str_for_custom", "GRCh38")
        logger.info(f"Using genome_build='{genome_build_arg}' and custom_genome='{custom_genome_path_arg}' for SigProfilerSimulator.")
    else:
        logger.info(f"Treating reference_input '{reference_input}' as a genome build string.")
        genome_build_arg = reference_input
        ref_display_for_manifest = genome_build_arg # Store the string (e.g. "GRCh38")

    try:
        from SigProfilerSimulator import SigProfilerSimulator
        logger.info(f"Running SigProfilerSimulator for sample '{sample_id}' in '{run_output_dir}'...")
        sps = SigProfilerSimulator(
            project=sample_id,
            genome_build=genome_build_arg,
            custom_genome=custom_genome_path_arg,
            outdir=str(run_output_dir),
            num_samples=1,
            exome=exome,
            type=[sig_type],
            total_mutations=[num_mutations],
            chrom_based=manifest_params.get("chrom_based", False),
            seed=manifest_params.get("seed", 0)
        )
        # Assuming SigProfilerSimulator logs its own progress adequately
        logger.info("SigProfilerSimulator finished successfully.")

    except ImportError:
        logger.error("SigProfilerSimulator library not found. Please install it.")
        raise bb_utils.BaseBuddyConfigError(
            "SigProfilerSimulator library not found. Please install it using: pip install SigProfilerSimulator"
        )
    except Exception as e:
        logger.error(f"SigProfilerSimulator failed during execution: {e}")
        raise bb_utils.BaseBuddyToolError(
            message=f"SigProfilerSimulator tool encountered an error during execution: {str(e)}",
            command=["SigProfilerSimulator (library call)"],
            stderr=str(e)
        )

    output_files_manifest: List[Dict[str,str]] = []
    actual_sps_output_dir = run_output_dir / sample_id / sig_type
    vcf_filename = f"{sample_id}.{sig_type}.all" # This is the VCF with all mutations
    output_vcf_path = actual_sps_output_dir / vcf_filename

    if output_vcf_path.exists():
        bb_utils.ensure_file_exists(output_vcf_path, "Simulated VCF from SigProfilerSimulator")
        relative_vcf_path = output_vcf_path.relative_to(run_output_dir)
        output_files_manifest.append({"name": f"Simulated VCF ({sig_type})", "path": str(relative_vcf_path), "type": "VCF"})
    else:
        logger.warning(f"Primary output VCF from SigProfilerSimulator not found at expected path: {output_vcf_path}")
        if actual_sps_output_dir.exists():
            logger.warning(f"Contents of {actual_sps_output_dir}: {list(actual_sps_output_dir.iterdir())}")
        else:
            logger.warning(f"SigProfilerSimulator output directory {actual_sps_output_dir} does not exist.")

    manifest_path = run_output_dir / "manifest.json"
    bb_utils.write_run_manifest(
        manifest_path=manifest_path, run_name=effective_run_name, command_name="simulate_signatures",
        parameters=manifest_params, output_files=output_files_manifest,
        reference_genome_path=ref_display_for_manifest
    )
    logger.info(f"Mutational signature simulation successful for run '{effective_run_name}'. Output in '{run_output_dir}'.")

    return {
        "run_name": effective_run_name,
        "output_directory": str(run_output_dir.resolve()),
        "output_files": output_files_manifest,
        "manifest_path": str(manifest_path.resolve()),
        "reference_used": ref_display_for_manifest,
        "actual_sps_output_dir": str(actual_sps_output_dir.resolve())
    }


def introduce_strand_bias(in_bam: str, out_bam: str, forward_fraction: float = 0.5, seed: int = 0) -> None:
    logger.info(f"Introducing strand bias to '{in_bam}', outputting to '{out_bam}'.")
    in_bam_path = bb_utils.ensure_file_exists(in_bam, "Input BAM for strand bias")

    # Ensure output directory for out_bam exists
    out_bam_path_obj = Path(out_bam)
    bb_utils.ensure_directory_exists(out_bam_path_obj.parent, "Output directory for biased BAM", create=True)

    # Indexing input BAM if necessary
    # Check common BAI locations: next to BAM or with .bai appended to full name
    bai_path_1 = in_bam_path.with_suffix(in_bam_path.suffix + ".bai")
    bai_path_2 = Path(str(in_bam_path) + ".bai") # Handles cases like file.sorted.bam -> file.sorted.bam.bai

    if not bai_path_1.exists() and not bai_path_2.exists():
        logger.info(f"Indexing input BAM for strand bias: {in_bam_path}")
        index_cmd_builder = bb_utils.Command("samtools")
        index_cmd_builder.add_option(None, "index")
        index_cmd_builder.add_option(None, str(in_bam_path))
        bb_utils.run_external_cmd(index_cmd_builder.get_command_parts(), cwd=in_bam_path.parent, stream_output=True)

    outdir = out_bam_path_obj.parent # This is already ensured by ensure_directory_exists
    temp_dir = outdir / f"strandtemp_{Path(in_bam).stem}_{seed}" # More unique temp dir

    try:
        bb_utils.ensure_directory_exists(temp_dir, "Temporary directory for strand bias", create=True)

        plus_bam = temp_dir / "plus.bam"
        minus_bam = temp_dir / "minus.bam"

        logger.debug(f"Separating forward strand reads to: {plus_bam}")
        view_plus_cmd = bb_utils.Command("samtools")
        view_plus_cmd.add_option(None, "view")
        view_plus_cmd.add_flag("-h")
        view_plus_cmd.add_option("-F", "16") # Select reads NOT reverse-complemented (i.e., forward strand)
        view_plus_cmd.add_option(None, str(in_bam_path))
        view_plus_cmd.add_option("-o", str(plus_bam))
        bb_utils.run_external_cmd(view_plus_cmd.get_command_parts(), cwd=temp_dir, stream_output=True)

        logger.debug(f"Separating reverse strand reads to: {minus_bam}")
        view_minus_cmd = bb_utils.Command("samtools")
        view_minus_cmd.add_option(None, "view")
        view_minus_cmd.add_flag("-h")
        view_minus_cmd.add_option("-f", "16") # Select reads that ARE reverse-complemented
        view_minus_cmd.add_option(None, str(in_bam_path))
        view_minus_cmd.add_option("-o", str(minus_bam))
        bb_utils.run_external_cmd(view_minus_cmd.get_command_parts(), cwd=temp_dir, stream_output=True)

        plus_keep_fraction_str = f"{seed}.{int(forward_fraction * 1000):03d}"
        minus_keep_fraction_str = f"{seed}.{int((1.0 - forward_fraction) * 1000):03d}"

        plus_sub_bam = temp_dir / "plus_sub.bam"
        minus_sub_bam = temp_dir / "minus_sub.bam"

        logger.debug(f"Subsampling forward strand reads (fraction: {forward_fraction}, seed string: {plus_keep_fraction_str})")
        subsample_plus_cmd = bb_utils.Command("samtools")
        subsample_plus_cmd.add_option(None, "view")
        subsample_plus_cmd.add_option("-s", plus_keep_fraction_str)
        subsample_plus_cmd.add_flag("-b") # Output BAM
        subsample_plus_cmd.add_option(None, str(plus_bam))
        subsample_plus_cmd.add_option("-o", str(plus_sub_bam))
        bb_utils.run_external_cmd(subsample_plus_cmd.get_command_parts(), cwd=temp_dir, stream_output=True)

        logger.debug(f"Subsampling reverse strand reads (fraction: {1.0-forward_fraction}, seed string: {minus_keep_fraction_str})")
        subsample_minus_cmd = bb_utils.Command("samtools")
        subsample_minus_cmd.add_option(None, "view")
        subsample_minus_cmd.add_option("-s", minus_keep_fraction_str)
        subsample_minus_cmd.add_flag("-b") # Output BAM
        subsample_minus_cmd.add_option(None, str(minus_bam))
        subsample_minus_cmd.add_option("-o", str(minus_sub_bam))
        bb_utils.run_external_cmd(subsample_minus_cmd.get_command_parts(), cwd=temp_dir, stream_output=True)

        merged_bam = temp_dir / "merged.bam"
        logger.debug(f"Merging subsampled strand BAMs to: {merged_bam}")
        merge_cmd = bb_utils.Command("samtools")
        merge_cmd.add_option(None, "merge")
        merge_cmd.add_flag("-f") # Overwrite if exists (should not with unique temp_dir)
        merge_cmd.add_option(None, str(merged_bam))
        merge_cmd.add_option(None, str(plus_sub_bam))
        merge_cmd.add_option(None, str(minus_sub_bam))
        bb_utils.run_external_cmd(merge_cmd.get_command_parts(), cwd=temp_dir, stream_output=True)

        logger.info(f"Sorting merged BAM to final output: {out_bam_path_obj}")
        sort_cmd = bb_utils.Command("samtools")
        sort_cmd.add_option(None, "sort")
        sort_cmd.add_option(None, str(merged_bam))
        sort_cmd.add_option("-o", str(out_bam_path_obj))
        bb_utils.run_external_cmd(sort_cmd.get_command_parts(), cwd=temp_dir, stream_output=True) # cwd can be outdir too

        logger.info(f"Indexing final biased BAM: {out_bam_path_obj}")
        final_index_cmd = bb_utils.Command("samtools")
        final_index_cmd.add_option(None, "index")
        final_index_cmd.add_option(None, str(out_bam_path_obj))
        bb_utils.run_external_cmd(final_index_cmd.get_command_parts(), cwd=outdir, stream_output=True)

        logger.info("Strand bias introduction process completed successfully.")

    finally:
        if temp_dir.exists():
            logger.debug(f"Cleaning up temporary directory: {temp_dir}")
            shutil.rmtree(temp_dir)


def run_germline_simulation_workflow(
    output_root_dir: Path,
    reference_fasta_path: str,
    germline_vcf_path: str,
    read_sim_type: str, # "short" or "long"
    read_sim_params: Dict[str, Any], # Params for simulate_short/simulate_long
    run_name: Optional[str] = None,
    command_params: Optional[Dict[str, Any]] = None, # For manifest
    overwrite_output: bool = False,
    auto_index_fasta: bool = True, # For intermediate and final FASTAs
    keep_intermediate_files: bool = False # New param to control cleanup
) -> Dict[str, Any]:

    effective_run_name = run_name if run_name else bb_utils.generate_unique_run_name(f"germline_sim_{read_sim_type}")
    logger.info(f"Initiating germline simulation workflow for run: {effective_run_name}")

    if command_params is None: command_params = {}
    run_output_dir = bb_utils.prepare_run_output_dir(Path(output_root_dir), effective_run_name, overwrite_output)

    manifest_params = copy.deepcopy(command_params)
    manifest_params.update({
        "run_name": effective_run_name, "reference_fasta_original": reference_fasta_path,
        "germline_vcf_path": germline_vcf_path, "read_sim_type": read_sim_type,
        "read_sim_params_input": copy.deepcopy(read_sim_params), # Store initial sim params
        "overwrite_output": overwrite_output, "auto_index_fasta": auto_index_fasta,
        "keep_intermediate_files": keep_intermediate_files,
        "output_root_dir": str(output_root_dir)
    })

    try:
        perl_exe_path = bb_utils.find_tool_path("perl")
        simug_script_path = bb_utils.find_tool_path("simuG.pl")
        samtools_exe_path = bb_utils.find_tool_path("samtools")
    except bb_utils.BaseBuddyConfigError as e:
        logger.error(f"A required tool for germline simulation was not found: {e.details if hasattr(e, 'details') else str(e)}")
        raise

    original_ref_path_obj = bb_utils.ensure_file_exists(reference_fasta_path, "Original reference FASTA").resolve()
    bb_utils.check_fasta_indexed(original_ref_path_obj, samtools_exe_path, auto_index_if_missing=auto_index_fasta)

    temp_snp_vcf_for_simug: Optional[Path] = None
    temp_indel_vcf_for_simug: Optional[Path] = None
    intermediate_simug_outputs_to_clean: List[Path] = []

    final_modified_fasta_path: Optional[Path] = None
    simulated_reads_results: Dict[str, Any] = {"status": "pending_read_simulation"}


    try:
        # VCF Splitting
        actual_germline_vcf_path = bb_utils.ensure_file_exists(germline_vcf_path, "Germline VCF").resolve()
        temp_snp_vcf_for_simug = run_output_dir / f"temp_snps_for_simug_{effective_run_name}.vcf"
        temp_indel_vcf_for_simug = run_output_dir / f"temp_indels_for_simug_{effective_run_name}.vcf"
        intermediate_simug_outputs_to_clean.extend([temp_snp_vcf_for_simug, temp_indel_vcf_for_simug])

        snp_vcf_writer, indel_vcf_writer = None, None
        snp_count, indel_count = 0, 0

        with pysam.VariantFile(str(actual_germline_vcf_path)) as vcf_in:
            # Create writers with header from input VCF
            snp_vcf_writer = pysam.VariantFile(str(temp_snp_vcf_for_simug), "w", header=vcf_in.header)
            indel_vcf_writer = pysam.VariantFile(str(temp_indel_vcf_for_simug), "w", header=vcf_in.header)

            for rec in vcf_in.fetch():
                if not rec.alts or len(rec.alts) > 1: # Skip no-alt or multi-allelic for simplicity with SimuG
                    logger.warning(f"Skipping multi-allelic or no-alt record in germline VCF: {rec.chrom}:{rec.pos+1}")
                    continue

                is_snp = len(rec.ref) == 1 and len(rec.alts[0]) == 1
                is_indel = (len(rec.ref) == 1 and len(rec.alts[0]) > 1) or \
                           (len(rec.ref) > 1 and len(rec.alts[0]) == 1)

                if is_snp:
                    snp_vcf_writer.write(rec); snp_count += 1
                elif is_indel:
                    indel_vcf_writer.write(rec); indel_count += 1
                # else: logger.debug(f"Skipping complex/other variant type: {rec.chrom}:{rec.pos+1}")

        if snp_vcf_writer: snp_vcf_writer.close()
        if indel_vcf_writer: indel_vcf_writer.close()

        logger.info(f"Split germline VCF: {snp_count} SNPs written to {temp_snp_vcf_for_simug.name}, {indel_count} Indels to {temp_indel_vcf_for_simug.name}")
        manifest_params["temp_snp_vcf_for_simug"] = temp_snp_vcf_for_simug.name if snp_count > 0 else None
        manifest_params["temp_indel_vcf_for_simug"] = temp_indel_vcf_for_simug.name if indel_count > 0 else None

        # Sequential SimuG Execution
        current_fasta_for_simug = original_ref_path_obj

        if snp_count > 0:
            simug_snp_prefix_str = str(run_output_dir / f"{effective_run_name}_snp_modified")
            cmd_builder_snp = bb_utils.Command(perl_exe_path).add_option(None, simug_script_path)
            cmd_builder_snp.add_option("-refseq", str(current_fasta_for_simug))
            cmd_builder_snp.add_option("-snp_vcf", str(temp_snp_vcf_for_simug))
            cmd_builder_snp.add_option("-prefix", simug_snp_prefix_str)

            logger.info(f"Running SimuG for SNPs: {' '.join(cmd_builder_snp.get_command_parts())}")
            bb_utils.run_external_cmd(cmd_builder_snp.get_command_parts(), cwd=run_output_dir)

            current_fasta_for_simug = Path(simug_snp_prefix_str + ".simulated.genome.fa")
            bb_utils.ensure_file_exists(current_fasta_for_simug, "SimuG SNP-modified FASTA")
            intermediate_simug_outputs_to_clean.append(current_fasta_for_simug)
            # SimuG creates other files based on prefix; add them to cleanup
            intermediate_simug_outputs_to_clean.append(Path(simug_snp_prefix_str + ".variant.locations.txt"))
            intermediate_simug_outputs_to_clean.append(Path(simug_snp_prefix_str + ".output.vcf"))
            logger.info(f"FASTA after SNP modification: {current_fasta_for_simug.name}")
            manifest_params["simug_snp_modified_fasta"] = current_fasta_for_simug.name

        if indel_count > 0:
            # Use a different prefix for the next SimuG run to avoid overwriting its own logs if same prefix used
            current_prefix_name_part = "_snp_indel_modified" if snp_count > 0 else "_indel_modified"
            simug_indel_prefix_str = str(run_output_dir / f"{effective_run_name}{current_prefix_name_part}")

            cmd_builder_indel = bb_utils.Command(perl_exe_path).add_option(None, simug_script_path)
            cmd_builder_indel.add_option("-refseq", str(current_fasta_for_simug))
            cmd_builder_indel.add_option("-indel_vcf", str(temp_indel_vcf_for_simug))
            cmd_builder_indel.add_option("-prefix", simug_indel_prefix_str)

            logger.info(f"Running SimuG for Indels: {' '.join(cmd_builder_indel.get_command_parts())}")
            bb_utils.run_external_cmd(cmd_builder_indel.get_command_parts(), cwd=run_output_dir)

            current_fasta_for_simug = Path(simug_indel_prefix_str + ".simulated.genome.fa")
            bb_utils.ensure_file_exists(current_fasta_for_simug, "SimuG Indel-modified FASTA")
            # If SNPs were applied, the FASTA from that step is now intermediate
            if snp_count > 0 and current_fasta_for_simug != intermediate_simug_outputs_to_clean[-1]:
                 # This condition might be true if simug_indel_prefix_str is different from simug_snp_prefix_str
                 pass # The previous current_fasta_for_simug is already in intermediate_simug_outputs_to_clean
            intermediate_simug_outputs_to_clean.append(current_fasta_for_simug) # Add the latest one
            intermediate_simug_outputs_to_clean.append(Path(simug_indel_prefix_str + ".variant.locations.txt"))
            intermediate_simug_outputs_to_clean.append(Path(simug_indel_prefix_str + ".output.vcf"))
            logger.info(f"FASTA after Indel modification: {current_fasta_for_simug.name}")
            manifest_params["simug_indel_modified_fasta"] = current_fasta_for_simug.name

        final_modified_fasta_path = run_output_dir / f"{effective_run_name}_final_germline_genome.fa"
        # current_fasta_for_simug should point to the latest version (snp+indel, or just snp, or just indel, or original if no variants)
        if current_fasta_for_simug != original_ref_path_obj and current_fasta_for_simug.exists():
            logger.info(f"Renaming SimuG output {current_fasta_for_simug.name} to {final_modified_fasta_path.name}")
            current_fasta_for_simug.rename(final_modified_fasta_path)
            # Remove the last SimuG output from intermediate_simug_outputs_to_clean as it's now the final one
            if final_modified_fasta_path in intermediate_simug_outputs_to_clean:
                intermediate_simug_outputs_to_clean.remove(final_modified_fasta_path)
        elif current_fasta_for_simug == original_ref_path_obj : # No variants applied by SimuG
            logger.info("No SNPs or Indels applied by SimuG. Using original reference for read simulation.")
            shutil.copy(original_ref_path_obj, final_modified_fasta_path)
            logger.info(f"Copied original reference to {final_modified_fasta_path.name} as no germline variants were applied/found.")
        else: # Should not happen if ensure_file_exists worked
            raise bb_utils.BaseBuddyFileError("Modified FASTA from SimuG not found unexpectedly.")

        bb_utils.check_fasta_indexed(final_modified_fasta_path, samtools_exe_path, auto_index_if_missing=auto_index_fasta)
        manifest_params["final_modified_germline_fasta"] = final_modified_fasta_path.name

        # -- START MODIFICATION AREA --
        logger.info(f"Starting read simulation ({read_sim_type}) using modified FASTA: {final_modified_fasta_path}")
        manifest_params["read_simulation_status"] = "in_progress"

        sub_run_output_root = run_output_dir
        sub_run_name = f"{effective_run_name}_{read_sim_type}_reads"

        sub_command_params_for_manifest = copy.deepcopy(read_sim_params)
        sub_command_params_for_manifest.update({
            "run_name": sub_run_name,
            "reference_fasta_used_for_simulation": str(final_modified_fasta_path.name),
            "original_germline_vcf": germline_vcf_path,
            "parent_workflow_run_name": effective_run_name,
        })

        try:
            if read_sim_type.lower() == "short":
                sim_results = simulate_short(
                    output_root_dir=sub_run_output_root,
                    run_name=sub_run_name,
                    command_params=sub_command_params_for_manifest,
                    reference_fasta=str(final_modified_fasta_path),
                    depth=read_sim_params.get("depth", 50),
                    read_length=read_sim_params.get("read_length", 150),
                    art_profile=read_sim_params.get("art_profile", "HS25"),
                    mean_fragment_length=read_sim_params.get("mean_fragment_length", 400),
                    std_dev_fragment_length=read_sim_params.get("std_dev_fragment_length", 50),
                    is_paired_end=read_sim_params.get("is_paired_end", True),
                    art_platform=read_sim_params.get("art_platform", "illumina"),
                    overwrite_output=overwrite_output,
                    auto_index_fasta=auto_index_fasta
                )
                simulated_reads_results = sim_results
            elif read_sim_type.lower() == "long":
                sim_results = simulate_long(
                    output_root_dir=sub_run_output_root,
                    run_name=sub_run_name,
                    command_params=sub_command_params_for_manifest,
                    reference_fasta=str(final_modified_fasta_path),
                    depth=read_sim_params.get("depth", 30),
                    model=read_sim_params.get("model", "nanopore_R9.4.1"),
                    overwrite_output=overwrite_output,
                    auto_index_fasta=auto_index_fasta
                )
                simulated_reads_results = sim_results
            else:
                err_msg = f"Unsupported read_sim_type: {read_sim_type}. Must be 'short' or 'long'."
                logger.error(err_msg)
                raise bb_utils.BaseBuddyInputError(err_msg)

            manifest_params["read_simulation_status"] = "completed"
            manifest_params["read_simulation_outputs"] = simulated_reads_results.get("output_files", [])
            manifest_params["read_simulation_run_name"] = sub_run_name
            logger.info(f"Read simulation completed. Results: {simulated_reads_results.get('output_directory')}")

        except Exception as e_sim:
            logger.error(f"Read simulation failed: {type(e_sim).__name__} - {e_sim}")
            logger.error(traceback.format_exc())
            manifest_params["read_simulation_status"] = "failed"
            manifest_params["read_simulation_error"] = str(e_sim)
            raise bb_utils.BaseBuddyToolError(f"Read simulation sub-step failed: {e_sim}", command=[read_sim_type, "simulation"])
        # -- END MODIFICATION AREA --

    except Exception as e:
        logger.error(f"Error in germline simulation workflow: {type(e).__name__} - {e}")
        logger.error(traceback.format_exc())
        # Write manifest even on failure, with error status
        manifest_params["workflow_status"] = "failed"
        manifest_params["workflow_error"] = str(e)
        bb_utils.write_run_manifest(
            manifest_path = run_output_dir / "manifest_germline_run.json",
            run_name = effective_run_name, command_name = "run_germline_simulation_workflow",
            parameters = manifest_params, output_files = [], # No guaranteed output files on error
            reference_genome_path = str(original_ref_path_obj)
        )
        raise # Re-raise the exception to be caught by CLI handler
    finally:
        if not keep_intermediate_files:
            logger.info("Cleaning up intermediate files for germline workflow...")
            for f_path in intermediate_simug_outputs_to_clean:
                if f_path.exists() and f_path != final_modified_fasta_path : # Don't delete final renamed FASTA
                    logger.debug(f"Deleting intermediate: {f_path.name}")
                    f_path.unlink(missing_ok=True)
                    # Attempt to delete associated SimuG log/VCF files if prefixes match intermediate FASTA names
                    for simug_ext in [".variant.locations.txt", ".output.vcf", ".simulated.genome.fa.fai"]:
                        simug_extra_file = f_path.with_name(f_path.name.replace(".simulated.genome.fa", "") + simug_ext)
                        if simug_extra_file.exists():
                            logger.debug(f"Deleting SimuG intermediate: {simug_extra_file.name}")
                            simug_extra_file.unlink(missing_ok=True)
        else:
            logger.info("Skipping cleanup of intermediate files as per 'keep_intermediate_files=True'.")


    manifest_params["workflow_status"] = "completed"
    # Collect all key output files for the main manifest
    main_output_files = [
        {"name": "Final Modified Germline FASTA", "path": final_modified_fasta_path.name, "type": "FASTA"},
        {"name": "FASTA Index for Modified Genome", "path": final_modified_fasta_path.name + ".fai", "type": "FAI"}
    ]
    if "output_files" in simulated_reads_results:
        sub_run_name = manifest_params.get("read_simulation_run_name", f"{effective_run_name}_{read_sim_type}_reads")
        for f_info in simulated_reads_results["output_files"]:
            main_output_files.append({
                "name": f"{read_sim_type.capitalize()} Reads: {f_info['name']}",
                "path": str(Path(sub_run_name) / f_info['path']),
                "type": f_info['type']
            })
    if "manifest_path" in simulated_reads_results:
         sub_run_name = manifest_params.get("read_simulation_run_name", f"{effective_run_name}_{read_sim_type}_reads")
         main_output_files.append({
                "name": f"{read_sim_type.capitalize()} Reads Simulation Manifest",
                "path": str(Path(sub_run_name) / Path(simulated_reads_results["manifest_path"]).name),
                "type": "MANIFEST_SUB"
            })

    bb_utils.write_run_manifest(
        manifest_path = run_output_dir / "manifest_germline_run.json",
        run_name = effective_run_name, command_name = "run_germline_simulation_workflow",
        parameters = manifest_params,
        output_files = main_output_files,
        reference_genome_path = str(original_ref_path_obj)
    )

    return {
        "run_name": effective_run_name,
        "output_directory": str(run_output_dir.resolve()),
        "final_modified_fasta_path": str(final_modified_fasta_path.resolve()) if final_modified_fasta_path else None,
        "simulated_reads_results": simulated_reads_results,
        "manifest_path": str((run_output_dir / "manifest_germline_run.json").resolve())
    }


def run_fastq_qc(
    fastq_files: List[str],
    output_dir_str: str, # This is the root for this QC run's outputs
    run_name: Optional[str] = None,
    command_params: Optional[Dict[str, Any]] = None, # For manifest, threads for FastQC etc.
    overwrite_output: bool = False
) -> Dict[str, Any]:

    effective_run_name = run_name if run_name else bb_utils.generate_unique_run_name("fastq_qc_run")
    logger.info(f"Initiating FastQC run: {effective_run_name} for {len(fastq_files)} file(s).")

    if command_params is None: command_params = {}

    qc_run_main_output_dir = bb_utils.prepare_run_output_dir(Path(output_dir_str), effective_run_name, overwrite_output)

    manifest_params = copy.deepcopy(command_params)
    manifest_params.update({
        "run_name": effective_run_name,
        "input_fastq_files": [str(Path(f).name) for f in fastq_files],
        "num_fastq_files": len(fastq_files),
        "output_qc_run_dir": str(qc_run_main_output_dir.resolve()),
        "overwrite_output": overwrite_output
    })

    try:
        # fastqc_exe_path is found here and used by Command. If just "fastqc" is passed to Command, it finds it.
        # Using bb_utils.Command("fastqc") directly is cleaner.
        pass # find_tool_path for "fastqc" will be implicitly called by Command("fastqc")
    except bb_utils.BaseBuddyConfigError as e: # This try-except is likely not needed here if Command handles it
        logger.error(f"A required tool for FastQC was not found: {e.details if hasattr(e, 'details') else str(e)}")
        raise

    qc_reports_generated = []
    all_output_files_for_manifest: List[Dict[str,str]] = []

    num_threads_for_fastqc = command_params.get("threads", 1)

    for fastq_file_path_str in fastq_files:
        fastq_path_obj = bb_utils.ensure_file_exists(fastq_file_path_str, f"Input FASTQ {fastq_file_path_str}").resolve()
        fastq_basename = fastq_path_obj.name

        logger.info(f"Running FastQC on: {fastq_basename}")

        # fastqc_cmd_builder will find 'fastqc' using find_tool_path internally
        fastqc_cmd_builder = bb_utils.Command("fastqc")
        fastqc_cmd_builder.add_option("--outdir", str(qc_run_main_output_dir))
        fastqc_cmd_builder.add_option("--threads", str(num_threads_for_fastqc))
        fastqc_cmd_builder.add_option(None, str(fastq_path_obj))

        try:
            bb_utils.run_external_cmd(fastqc_cmd_builder.get_command_parts(), stream_output=True, cwd=qc_run_main_output_dir)

            stripped_name = fastq_basename
            for ext in ['.gz', '.bz2', '.txt', '.sam', '.bam', '.fastq', '.fq']:
                if stripped_name.lower().endswith(ext):
                    stripped_name = stripped_name[:-len(ext)]
            for ext in ['.fastq', '.fq']:
                 if stripped_name.lower().endswith(ext):
                    stripped_name = stripped_name[:-len(ext)]

            fastqc_output_subdir_name = stripped_name + "_fastqc"
            report_html_path = qc_run_main_output_dir / fastqc_output_subdir_name / "fastqc_report.html"
            report_zip_path = qc_run_main_output_dir / (stripped_name + "_fastqc.zip")

            if report_html_path.exists():
                logger.info(f"FastQC report generated: {report_html_path}")
                relative_html_path = report_html_path.relative_to(qc_run_main_output_dir)
                qc_reports_generated.append({
                    "input_fastq": fastq_basename,
                    "html_report_path": str(relative_html_path),
                    "full_html_path": str(report_html_path.resolve())
                })
                all_output_files_for_manifest.append({
                    "name": f"FastQC HTML Report for {fastq_basename}",
                    "path": str(relative_html_path),
                    "type": "FASTQC_HTML_REPORT"
                })
                if report_zip_path.exists():
                    relative_zip_path = report_zip_path.relative_to(qc_run_main_output_dir)
                    all_output_files_for_manifest.append({
                        "name": f"FastQC ZIP Data for {fastq_basename}",
                        "path": str(relative_zip_path),
                        "type": "FASTQC_ZIP_DATA"
                    })
            else:
                logger.warning(f"FastQC HTML report not found at expected path: {report_html_path}")
                qc_reports_generated.append({
                    "input_fastq": fastq_basename,
                    "html_report_path": None,
                    "error": "Report HTML not found post-execution."
                })

        except bb_utils.BaseBuddyToolError as e_tool:
            logger.error(f"FastQC failed for {fastq_basename}: {e_tool}")
            qc_reports_generated.append({ "input_fastq": fastq_basename, "html_report_path": None, "error": str(e_tool.details if hasattr(e_tool, 'details') else e_tool)})
        except Exception as e_unknown:
            logger.error(f"An unexpected error occurred processing {fastq_basename} with FastQC: {e_unknown}")
            qc_reports_generated.append({ "input_fastq": fastq_basename, "html_report_path": None, "error": f"Unexpected error: {str(e_unknown)}"})

    manifest_params["qc_reports_generated_details"] = qc_reports_generated
    manifest_params["status_detail"] = f"FastQC processing completed for {len(fastq_files)} file(s)."

    main_manifest_path = qc_run_main_output_dir / "manifest_fastqc_run.json"
    bb_utils.write_run_manifest(
        manifest_path=main_manifest_path,
        run_name=effective_run_name,
        command_name="run_fastq_qc",
        parameters=manifest_params,
        output_files=all_output_files_for_manifest
    )

    return {
        "run_name": effective_run_name,
        "output_directory": str(qc_run_main_output_dir.resolve()),
        "qc_reports": qc_reports_generated,
        "manifest_path": str(main_manifest_path.resolve())
    }

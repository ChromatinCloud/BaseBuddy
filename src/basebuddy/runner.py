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
        art_exe_path = bb_utils.find_tool_path(art_exe_name)
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

        # ART command construction
        cmd = [art_exe_path, "-ss", art_profile, "-i", str(current_reference_for_art),
               "-l", str(read_length), "-f", str(depth), "-o", str(art_output_prefix_path)]
        if is_paired_end:
            cmd.extend(["-p", "-m", str(mean_fragment_length), "-s", str(std_dev_fragment_length)])
        if manifest_params.get("quality_shift") is not None: cmd.extend(["-qs", str(manifest_params["quality_shift"])])
        if manifest_params.get("quality_shift2") is not None and is_paired_end: cmd.extend(["-qs2", str(manifest_params["quality_shift2"])])
        if manifest_params.get("random_seed") is not None: cmd.extend(["-rs", str(manifest_params["random_seed"])])
        if manifest_params.get("no_aln_output", False): cmd.append("-na")

        logger.info(f"Preparing to run {art_exe_name}...")
        bb_utils.run_external_cmd(cmd, timeout_seconds=timeout, stream_output=True, cwd=run_output_dir)

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
    input_bams: List[str],
    variants_list: List[Dict[str, Any]],
    output_prefix_for_bam: str,
    run_name: Optional[str] = None,
    command_params: Optional[Dict[str, Any]] = None,
    overwrite_output: bool = False,
    auto_index_input_bam: bool = True,
    auto_index_fasta: bool = True,
    timeout: float = 7200.0  # Timeout per BAM, effectively
) -> Dict[str, Any]:

    effective_run_name = run_name if run_name else bb_utils.generate_unique_run_name("spike_variants_run")
    logger.info(f"Initiating variant spiking for run: {effective_run_name} on {len(input_bams)} BAMs.")

    run_output_dir = bb_utils.prepare_run_output_dir(Path(output_root_dir), effective_run_name, overwrite_output)

    base_manifest_params = copy.deepcopy(command_params) if command_params is not None else {}
    base_manifest_params.pop("vcf", None)
    base_manifest_params.pop("vcf_file", None)
    base_manifest_params.pop("input_bam", None) # Remove single bam key if present
    base_manifest_params.pop("input_bams", None) # Remove if somehow already there

    base_manifest_params.update({
        "run_name": effective_run_name, "reference_fasta": reference_fasta,
        "input_bams": [str(Path(b).name) for b in input_bams], # Store list of basenames
        "num_input_bams": len(input_bams),
        "num_variants_to_spike": len(variants_list),
        "output_prefix_for_bam": output_prefix_for_bam, "overwrite_output": overwrite_output,
        "auto_index_input_bam": auto_index_input_bam, "auto_index_fasta": auto_index_fasta,
        "timeout_per_bam": timeout, "output_root_dir": str(output_root_dir)
    })
    vaf = base_manifest_params.get("vaf", 0.5)
    seed = base_manifest_params.get("seed", 0)

    temp_vcf_file_path_str = ""
    all_results_per_bam: List[Dict[str, Any]] = []
    all_errors_per_bam: List[Dict[str, str]] = []
    overall_output_files_manifest: List[Dict[str,str]] = []

    try:
        # --- Temporary VCF Generation (once for the run) ---
        temp_vcf_path = run_output_dir / "temp_variants_for_run.vcf"
        temp_vcf_file_path_str = str(temp_vcf_path)

        with open(temp_vcf_path, "w") as f_vcf:
            f_vcf.write("##fileformat=VCFv4.2\n")
            try:
                fai_file_for_contigs = Path(reference_fasta + ".fai")
                if fai_file_for_contigs.exists():
                    with open(fai_file_for_contigs, "r") as fai:
                        for line in fai:
                            parts = line.split("\t")
                            f_vcf.write(f"##contig=<ID={parts[0]},length={parts[1]}>\n")
            except Exception as e:
                logger.warning(f"Could not read .fai file to generate contig headers for temp VCF: {e}")

            f_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for var in variants_list:
                chrom = var.get("chromosome")
                pos = var.get("position")
                ref = var.get("ref_allele")
                alt = var.get("alt_allele")
                if not all([chrom, pos, ref, alt]):
                    logger.warning(f"Skipping invalid variant entry for temp VCF: {var}")
                    continue
                f_vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
        logger.info(f"Common temporary VCF for spiking created at {temp_vcf_path}")

        addsnv_exe_path = bb_utils.find_tool_path("addsnv.py")
        samtools_exe_path = bb_utils.find_tool_path("samtools")

        ref_path_obj = bb_utils.ensure_file_exists(reference_fasta, "Reference FASTA").resolve()
        bb_utils.check_fasta_indexed(ref_path_obj, samtools_exe_path, auto_index_if_missing=auto_index_fasta)

        for input_bam_str in input_bams:
            input_bam_stem = Path(input_bam_str).stem
            logger.info(f"Processing BAM: {input_bam_str} (stem: {input_bam_stem})")

            output_files_for_this_bam: List[Dict[str,str]] = []

            try:
                in_bam_path_obj = bb_utils.ensure_file_exists(input_bam_str, f"Input BAM {input_bam_str}").resolve()
                bb_utils.check_bam_indexed(in_bam_path_obj, samtools_exe_path, auto_index_if_missing=auto_index_input_bam)

                final_bam_name = f"{output_prefix_for_bam}_{input_bam_stem}.bam"
                final_bam_path = run_output_dir / final_bam_name

                if final_bam_path.exists() and not overwrite_output:
                    raise bb_utils.BaseBuddyFileError(
                        f"Target output BAM file '{final_bam_path}' for input '{input_bam_str}' already exists and overwrite is not permitted."
                    )

                temp_addsnv_out_bam_name = f"{output_prefix_for_bam}_{input_bam_stem}_temp_addsnv.bam"
                temp_addsnv_out_bam_path = run_output_dir / temp_addsnv_out_bam_name

                cmd_addsnv = [
                    "python", addsnv_exe_path,
                    "--reference", str(ref_path_obj), "--in_bam", str(in_bam_path_obj),
                    "--vcf", str(temp_vcf_path), "--out_bam", str(temp_addsnv_out_bam_path),
                    "-p", str(vaf), "-s", str(seed)
                ]
                logger.info(f"Running addsnv.py for {input_bam_str}...")
                bb_utils.run_external_cmd(cmd_addsnv, timeout_seconds=timeout*0.8, stream_output=True, cwd=run_output_dir)
                bb_utils.ensure_file_exists(temp_addsnv_out_bam_path, f"Temporary output BAM from addsnv.py for {input_bam_str}")

                logger.info(f"Sorting BAM: {temp_addsnv_out_bam_path.name} -> {final_bam_path.name}")
                cmd_sort = [samtools_exe_path, "sort", str(temp_addsnv_out_bam_path), "-o", str(final_bam_path)]
                bb_utils.run_external_cmd(cmd_sort, timeout_seconds=timeout*0.1, cwd=run_output_dir)

                try:
                    temp_addsnv_out_bam_path.unlink()
                    logger.debug(f"Removed temporary addsnv BAM: {temp_addsnv_out_bam_path.name}")
                except OSError as e:
                    logger.warning(f"Could not remove temporary addsnv BAM {temp_addsnv_out_bam_path.name}: {e}")

                bb_utils.check_bam_indexed(final_bam_path, samtools_exe_path, auto_index_if_missing=True)
                final_bam_path_obj = bb_utils.ensure_file_exists(final_bam_path, f"Final output BAM for {input_bam_str}")

                final_bai_path_str_loop = ""
                bai_path_1_loop = final_bam_path.with_suffix(final_bam_path.suffix + ".bai")
                bai_path_2_loop = final_bam_path.with_suffix(".bai")
                if bai_path_1_loop.exists(): final_bai_path_str_loop = bai_path_1_loop.name
                elif bai_path_2_loop.exists(): final_bai_path_str_loop = bai_path_2_loop.name

                output_files_for_this_bam.append(
                    {"name": f"Spiked Variants BAM ({input_bam_stem})", "path": final_bam_path_obj.name, "type": "BAM"}
                )
                if final_bai_path_str_loop:
                    output_files_for_this_bam.append(
                        {"name": f"Spiked BAM Index ({input_bam_stem})", "path": final_bai_path_str_loop, "type": "BAI"}
                    )

                igv_session_file_name = f"igv_session_spike_{input_bam_stem}.xml"
                igv_session_file_path = run_output_dir / igv_session_file_name
                igv_tracks = [
                    {"name": f"Spiked BAM ({input_bam_stem})", "path": final_bam_path_obj.name, "type": "alignment"},
                    {"name": "Spiked Variants (common VCF)", "path": str(temp_vcf_path.name), "type": "variant"} # Common VCF
                ]
                bb_utils.generate_igv_session_xml(igv_session_file_path, str(ref_path_obj), igv_tracks)
                output_files_for_this_bam.append(
                    {"name": f"IGV Session ({input_bam_stem})", "path": igv_session_file_name, "type": "IGV_SESSION"}
                )

                overall_output_files_manifest.extend(output_files_for_this_bam)

                all_results_per_bam.append({
                    "source_bam": input_bam_str,
                    "output_bam": str(final_bam_path_obj.resolve()),
                    "output_bam_index": str((run_output_dir / final_bai_path_str_loop).resolve()) if final_bai_path_str_loop else None,
                    "igv_session_file": str(igv_session_file_path.resolve()),
                    "status": "success"
                })
                logger.info(f"Variant spiking successful for BAM '{input_bam_str}'. Final BAM: '{final_bam_path_obj.name}'.")

            except Exception as e_bam:
                logger.error(f"Error processing BAM {input_bam_str}: {e_bam}")
                all_errors_per_bam.append({"source_bam": input_bam_str, "error": str(e_bam), "details": traceback.format_exc() if isinstance(e_bam, BaseException) else "" })
                # Continue to the next BAM

        # Single manifest for the entire run
        manifest_path = run_output_dir / "manifest_spike_run.json"
        bb_utils.write_run_manifest(
            manifest_path=manifest_path, run_name=effective_run_name,
            command_name="spike_variants_batch",
            parameters=base_manifest_params, output_files=overall_output_files_manifest,
            reference_genome_path=str(ref_path_obj)
        )

        return {
            "run_name": effective_run_name,
            "output_directory": str(run_output_dir.resolve()),
            "results_per_bam": all_results_per_bam,
            "errors_per_bam": all_errors_per_bam,
            "manifest_path": str(manifest_path.resolve()),
            "reference_fasta_used": str(ref_path_obj)
        }
    finally:
        if temp_vcf_file_path_str and Path(temp_vcf_file_path_str).exists():
            try:
                Path(temp_vcf_file_path_str).unlink()
                logger.debug(f"Removed common temporary VCF: {temp_vcf_file_path_str}")
            except OSError as e:
                logger.warning(f"Could not remove common temporary VCF {temp_vcf_file_path_str}: {e}")

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
        nanosim_exe_path = bb_utils.find_tool_path("nanosim-h")
        samtools_exe_path = bb_utils.find_tool_path("samtools")

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

        # NanoSim-h command construction
        cmd_base = [nanosim_exe_path, "simulate", "-r", str(current_reference_for_nanosim), "-m", model, "-o", str(nanosim_output_prefix_path)]

        # Prefer num_reads if provided and valid, otherwise use depth.
        num_reads_param = manifest_params.get("num_reads")
        if num_reads_param is not None:
            try:
                num_reads_val = int(num_reads_param)
                if num_reads_val > 0:
                    cmd = cmd_base + ["-N", str(num_reads_val)]
                    if "depth" in manifest_params: # Ensure depth is not also defining coverage
                         logger.info("Using num_reads for NanoSim; 'depth' parameter will be ignored by NanoSim if also set.")
                elif depth > 0: # num_reads invalid, fallback to depth
                    cmd = cmd_base + ["-c", str(depth)]
                    logger.warning(f"num_reads parameter '{num_reads_param}' is not a positive integer. Falling back to depth '{depth}'.")
                else: # Both num_reads and depth are invalid
                    raise bb_utils.BaseBuddyInputError("Either 'depth' or 'num_reads' must be a positive integer for NanoSim.")
            except ValueError: # num_reads not an int
                if depth > 0:
                    cmd = cmd_base + ["-c", str(depth)]
                    logger.warning(f"num_reads parameter '{num_reads_param}' is not a valid integer. Falling back to depth '{depth}'.")
                else:
                    raise bb_utils.BaseBuddyInputError("num_reads is not a valid integer and depth is not positive for NanoSim.")
        elif depth > 0: # num_reads not provided, use depth
            cmd = cmd_base + ["-c", str(depth)]
        else: # Neither num_reads nor depth is valid
             raise bb_utils.BaseBuddyInputError("Either 'depth' or 'num_reads' must be a positive integer for NanoSim.")


        logger.info(f"Running NanoSim-h with command: {' '.join(cmd)}")
        bb_utils.run_external_cmd(cmd, timeout_seconds=timeout, stream_output=True, cwd=run_output_dir)

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

    from SigProfilerSimulator import SigProfilerSimulator

    try:
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
    except Exception as e:
        logger.error(f"SigProfilerSimulator failed: {e}")
        raise bb_utils.BaseBuddyToolError(message=f"SigProfilerSimulator tool failed: {str(e)}", command=["SigProfilerSimulator", "...details..."])

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
    _ensure_exists(in_bam)
    if not Path(in_bam + ".bai").exists():
        subprocess.run(["samtools", "index", in_bam], check=True)
    outdir = Path(out_bam).parent
    outdir.mkdir(parents=True, exist_ok=True)
    temp_dir = outdir / "strandtemp"
    temp_dir.mkdir(exist_ok=True)
    plus_bam = temp_dir / "plus.bam"
    minus_bam = temp_dir / "minus.bam"
    subprocess.run(["samtools", "view", "-h", "-F", "16", in_bam, "-o", str(plus_bam)], check=True)
    subprocess.run(["samtools", "view", "-h", "-f", "16", in_bam, "-o", str(minus_bam)], check=True)
    plus_keep = forward_fraction
    minus_keep = 1.0 - forward_fraction
    subprocess.run([
        "samtools", "view", "-s",
        f"{seed}.{int(plus_keep * 1000):03d}",
        "-b", str(plus_bam), "-o", str(temp_dir / "plus_sub.bam")
    ], check=True)
    subprocess.run([
        "samtools", "view", "-s",
        f"{seed}.{int(minus_keep * 1000):03d}",
        "-b", str(minus_bam), "-o", str(temp_dir / "minus_sub.bam")
    ], check=True)
    merged = temp_dir / "merged.bam"
    subprocess.run([
        "samtools", "merge", "-f", str(merged),
        str(temp_dir / "plus_sub.bam"), str(temp_dir / "minus_sub.bam"),
    ], check=True)
    subprocess.run(["samtools", "sort", "-o", out_bam, str(merged)], check=True)
    subprocess.run(["samtools", "index", out_bam], check=True)
    shutil.rmtree(temp_dir)

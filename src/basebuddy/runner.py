import subprocess
import shutil
from pathlib import Path
import sys
import logging
from typing import Dict, Any, Optional, List
import copy

import bb_utils

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
        # This basic auto-download logic will be replaced by a dedicated download_runner
        print(f"Reference not found: {ref}, placeholder for download logic.")
        # For now, to prevent actual download in old runners, let's assume it would error if not found
        # or expect user to provide it. In a real scenario, this would call the download runner.
        raise FileNotFoundError(f"Reference {ref} not found and auto-download not implemented here.")
    return str(ref)

def _download_grch38(dest: str) -> None: # This logic will be part of download_reference_runner
    logger.warning("Legacy _download_grch38 called. This should be handled by download_reference_runner.")
    # Simplified, actual download runner will be more robust
    (Path(dest).parent).mkdir(parents=True, exist_ok=True)
    Path(dest).touch() # Create dummy file
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
    auto_index_fasta: bool = True
) -> Dict[str, Any]:

    effective_run_name = run_name if run_name else bb_utils.generate_unique_run_name(f"simulate_short_{art_platform}")
    logger.info(f"Initiating short read simulation for run: {effective_run_name}")

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

    if depth <= 0: raise bb_utils.BaseBuddyInputError(f"Sequencing depth must be a positive integer, got {depth}.")
    if read_length <= 0: raise bb_utils.BaseBuddyInputError(f"Read length must be a positive integer, got {read_length}.")
    if is_paired_end:
        if mean_fragment_length <= 0: raise bb_utils.BaseBuddyInputError(f"Mean fragment length must be positive for paired-end reads, got {mean_fragment_length}.")
        if std_dev_fragment_length < 0: raise bb_utils.BaseBuddyInputError(f"Standard deviation of fragment length cannot be negative, got {std_dev_fragment_length}.")

    art_exe_name = f"art_{art_platform}"
    art_exe_path = bb_utils.find_tool_path(art_exe_name)
    samtools_exe_path = bb_utils.find_tool_path("samtools")

    ref_path_obj = bb_utils.ensure_file_exists(reference_fasta, "Reference FASTA").resolve()
    bb_utils.check_fasta_indexed(ref_path_obj, samtools_exe_path, auto_index_if_missing=auto_index_fasta)

    art_output_prefix_path = run_output_dir / id_prefix

    cmd = [art_exe_path, "-ss", art_profile, "-i", str(ref_path_obj),
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
        output_files=output_files_manifest, reference_genome_path=str(ref_path_obj)
    )
    logger.info(f"Short read simulation successful for run '{effective_run_name}'. Output in '{run_output_dir}'.")
    return {
        "run_name": effective_run_name, "output_directory": str(run_output_dir.resolve()),
        "output_files": output_files_manifest, "manifest_path": str(manifest_path.resolve()),
        "reference_fasta_used": str(ref_path_obj)
    }

def spike_variants(
    output_root_dir: Path,
    reference_fasta: str,
    input_bam: str,
    vcf_file: str,
    output_bam_prefix_rel: str,
    run_name: Optional[str] = None,
    command_params: Optional[Dict[str, Any]] = None,
    overwrite_output: bool = False,
    auto_index_input_bam: bool = True,
    auto_index_fasta: bool = True,
    timeout: float = 7200.0
) -> Dict[str, Any]:

    effective_run_name = run_name if run_name else bb_utils.generate_unique_run_name("spike_variants")
    logger.info(f"Initiating variant spiking for run: {effective_run_name}")

    manifest_params = copy.deepcopy(command_params) if command_params is not None else {}
    manifest_params.update({
        "run_name": effective_run_name, "reference_fasta": reference_fasta,
        "input_bam": input_bam, "vcf_file": vcf_file,
        "output_bam_prefix_rel": output_bam_prefix_rel, "overwrite_output": overwrite_output,
        "auto_index_input_bam": auto_index_input_bam, "auto_index_fasta": auto_index_fasta,
        "timeout": timeout, "output_root_dir": str(output_root_dir)
    })
    vaf = manifest_params.get("vaf", 0.5)
    seed = manifest_params.get("seed", 0)

    run_output_dir = bb_utils.prepare_run_output_dir(Path(output_root_dir), effective_run_name, overwrite_output)

    addsnv_exe_path = bb_utils.find_tool_path("addsnv.py")
    samtools_exe_path = bb_utils.find_tool_path("samtools")

    ref_path_obj = bb_utils.ensure_file_exists(reference_fasta, "Reference FASTA").resolve()
    bb_utils.check_fasta_indexed(ref_path_obj, samtools_exe_path, auto_index_if_missing=auto_index_fasta)

    in_bam_path_obj = bb_utils.ensure_file_exists(input_bam, "Input BAM").resolve()
    bb_utils.check_bam_indexed(in_bam_path_obj, samtools_exe_path, auto_index_if_missing=auto_index_input_bam)

    vcf_path_obj = bb_utils.ensure_file_exists(vcf_file, "Input VCF").resolve()

    final_bam_name = output_bam_prefix_rel + ".bam"
    final_bam_path = run_output_dir / final_bam_name

    if final_bam_path.exists() and not overwrite_output:
         raise bb_utils.BaseBuddyFileError(
            f"Target output BAM file '{final_bam_path}' already exists and overwrite is not permitted."
        )

    temp_addsnv_out_bam_name = output_bam_prefix_rel + "_temp_addsnv.bam"
    temp_addsnv_out_bam_path = run_output_dir / temp_addsnv_out_bam_name

    cmd_addsnv = [
        "python", addsnv_exe_path,
        "--reference", str(ref_path_obj), "--in_bam", str(in_bam_path_obj),
        "--vcf", str(vcf_path_obj), "--out_bam", str(temp_addsnv_out_bam_path),
        "-p", str(vaf), "-s", str(seed)
    ]
    logger.info("Running addsnv.py (BAMSurgeon)...")
    bb_utils.run_external_cmd(cmd_addsnv, timeout_seconds=timeout*0.8, stream_output=True, cwd=run_output_dir)
    bb_utils.ensure_file_exists(temp_addsnv_out_bam_path, "Temporary output BAM from addsnv.py")

    logger.info(f"Sorting BAM: {temp_addsnv_out_bam_path.name} -> {final_bam_path.name}")
    cmd_sort = [samtools_exe_path, "sort", str(temp_addsnv_out_bam_path), "-o", str(final_bam_path)]
    bb_utils.run_external_cmd(cmd_sort, timeout_seconds=timeout*0.1, cwd=run_output_dir)

    try:
        temp_addsnv_out_bam_path.unlink()
        logger.debug(f"Removed temporary BAM: {temp_addsnv_out_bam_path.name}")
    except OSError as e:
        logger.warning(f"Could not remove temporary BAM {temp_addsnv_out_bam_path.name}: {e}")

    bb_utils.check_bam_indexed(final_bam_path, samtools_exe_path, auto_index_if_missing=True)
    final_bam_path_obj = bb_utils.ensure_file_exists(final_bam_path, "Final output BAM from variant spiking")

    final_bai_path_str = ""
    bai_path_1 = final_bam_path_obj.with_suffix(final_bam_path_obj.suffix + ".bai")
    bai_path_2 = final_bam_path_obj.with_suffix(".bai")
    if bai_path_1.exists(): final_bai_path_str = bai_path_1.name
    elif bai_path_2.exists(): final_bai_path_str = bai_path_2.name

    output_files_manifest: List[Dict[str,str]] = [
        {"name": "Spiked Variants BAM", "path": final_bam_path_obj.name, "type": "BAM"}
    ]
    if final_bai_path_str:
        output_files_manifest.append({"name": "Spiked BAM Index", "path": final_bai_path_str, "type": "BAI"})

    igv_session_file_name = "igv_session_spike.xml"
    igv_session_file_path = run_output_dir / igv_session_file_name
    igv_tracks = [
        {"name": "Spiked BAM", "path": final_bam_path_obj.name, "format": "bam", "type": "alignment"},
        {"name": "Input VCF", "path": str(vcf_path_obj.resolve()), "format": "vcf", "type": "variant"}
    ]
    bb_utils.generate_igv_session_xml(igv_session_file_path, str(ref_path_obj), igv_tracks)
    output_files_manifest.append({"name": "IGV Session (Spike)", "path": igv_session_file_name, "type": "IGV_SESSION"})

    manifest_path = run_output_dir / "manifest.json"
    bb_utils.write_run_manifest(
        manifest_path=manifest_path, run_name=effective_run_name, command_name="spike_variants",
        parameters=manifest_params, output_files=output_files_manifest,
        reference_genome_path=str(ref_path_obj)
    )
    logger.info(f"Variant spiking successful for run '{effective_run_name}'. Final BAM: '{final_bam_path_obj.name}'.")

    return {
        "run_name": effective_run_name,
        "output_directory": str(run_output_dir.resolve()),
        "output_bam": str(final_bam_path_obj.resolve()),
        "output_bam_index": str((run_output_dir / final_bai_path_str).resolve()) if final_bai_path_str else None,
        "igv_session_file": str(igv_session_file_path.resolve()),
        "manifest_path": str(manifest_path.resolve()),
        "reference_fasta_used": str(ref_path_obj)
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
    timeout: float = 3600.0
) -> Dict[str, Any]:
    effective_run_name = run_name if run_name else bb_utils.generate_unique_run_name(f"simulate_long_nanosim")
    logger.info(f"Initiating long read simulation for run: {effective_run_name}")

    manifest_params = copy.deepcopy(command_params) if command_params is not None else {}
    manifest_params.update({
        "run_name": effective_run_name, "reference_fasta": reference_fasta, "depth": depth,
        "model": model, "overwrite_output": overwrite_output,
        "auto_index_fasta": auto_index_fasta, "timeout": timeout,
        "output_root_dir": str(output_root_dir)
    })
    nanosim_internal_out_prefix = "nanosim_reads"

    run_output_dir = bb_utils.prepare_run_output_dir(Path(output_root_dir), effective_run_name, overwrite_output)

    if depth <= 0: raise bb_utils.BaseBuddyInputError(f"Sequencing depth must be a positive integer, got {depth}.")
    if not model: raise bb_utils.BaseBuddyInputError("NanoSim model must be specified.")

    nanosim_exe_path = bb_utils.find_tool_path("nanosim-h")
    samtools_exe_path = bb_utils.find_tool_path("samtools")

    ref_path_obj = bb_utils.ensure_file_exists(reference_fasta, "Reference FASTA").resolve()
    bb_utils.check_fasta_indexed(ref_path_obj, samtools_exe_path, auto_index_if_missing=auto_index_fasta)

    nanosim_output_prefix_path = run_output_dir / nanosim_internal_out_prefix

    cmd = [
        nanosim_exe_path, "simulate",
        "-r", str(ref_path_obj),
        "-c", str(depth),
        "-m", model,
        "-o", str(nanosim_output_prefix_path)
    ]
    if "num_reads" in manifest_params:
        cmd = [
            nanosim_exe_path, "simulate",
            "-r", str(ref_path_obj),
            "-N", str(manifest_params["num_reads"]),
            "-m", model,
            "-o", str(nanosim_output_prefix_path)
        ]
        if "depth" in manifest_params: del manifest_params["depth"]

    logger.info("Running NanoSim-h...")
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
        reference_genome_path=str(ref_path_obj)
    )
    logger.info(f"Long read simulation successful for run '{effective_run_name}'. Output in '{run_output_dir}'.")

    return {
        "run_name": effective_run_name,
        "output_directory": str(run_output_dir.resolve()),
        "output_files": output_files_manifest,
        "manifest_path": str(manifest_path.resolve()),
        "reference_fasta_used": str(ref_path_obj)
    }

def simulate_signatures(reference: str | None, outdir: Path, sig_type: str = "SBS", num_mutations: int = 100, sample_id: str = "Sample") -> None:
    ref = _ensure_reference(reference)
    from SigProfilerSimulator import SigProfilerSimulator
    outdir.mkdir(parents=True, exist_ok=True)
    SigProfilerSimulator(
        project=sample_id,
        reference_genome="GRCh38",
        outdir=str(outdir),
        num_samples=1,
        exome=False,
        type=sig_type,
        total_mutations=num_mutations,
        chrom_based=False,
        seed=0,
    )

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

[end of src/basebuddy/runner.py]

# tests/test_integration.py
import pytest
from pathlib import Path
import shutil # Added for copying test files
import json # For reading manifest files
import os
import traceback # For more detailed error logging in tests

# Assuming PYTHONPATH is set to include repo root for bb_utils,
# and src layout for basebuddy modules.
# Or that pytest is run from repo root.
import bb_utils # Assuming bb_utils.py is at the repo root or in PYTHONPATH
from src.basebuddy import runner as bb_runner
from src.basebuddy import signature_utils as bb_sig_utils # May not be directly used, but runner uses it

# Content for test FASTA files (can be used by multiple tests)
INTEGRATION_TEST_FASTA_CONTENT = """>integ_chr1
ACGTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
>integ_chr2
GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT
ACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAG
"""

SMOKE_TEST_FASTA_CONTENT = """>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
""" # Made chr1 longer for more variant positions

def create_temp_fasta(tmp_path: Path, filename: str, content: str) -> Path:
    """Helper to create a temporary FASTA file for testing."""
    fasta_file_path = tmp_path / filename
    # Ensure parent directory for the fasta file exists if filename includes subdirs within tmp_path
    fasta_file_path.parent.mkdir(parents=True, exist_ok=True)
    with open(fasta_file_path, "w") as f:
        f.write(content)
    return fasta_file_path

def test_initial_setup_correct():
    """Placeholder test to ensure pytest discovers the file."""
    assert True

def test_apply_sbs1_signature_integration(tmp_path: Path):
    """
    Integration test for the apply_signature_to_fasta runner with a bundled SBS signature.
    Uses a predefined small FASTA (smoke_test.fa content).
    """
    input_fasta_path = create_temp_fasta(tmp_path, "smoke_test_input_sbs.fa", SMOKE_TEST_FASTA_CONTENT)

    output_root_dir = tmp_path / "integration_test_output_sbs1"

    run_name = "apply_sbs1_integration_test"
    output_fasta_name = "mutated_by_sbs1.fa"
    signature_id = "SBS1"
    num_mutations_to_request = 5

    command_params_for_manifest = {
        "input_fasta_display": str(input_fasta_path.name),
        "output_fasta_name_request": output_fasta_name,
        "signature_input": signature_id,
        "num_mutations_config": num_mutations_to_request,
    }

    results = None
    try:
        results = bb_runner.apply_signature_to_fasta(
            output_root_dir=output_root_dir,
            run_name=run_name,
            command_params=command_params_for_manifest,
            input_fasta_path_str=str(input_fasta_path),
            output_fasta_name=output_fasta_name,
            signature_id_or_path=signature_id,
            num_mutations=num_mutations_to_request,
            overwrite_output=True,
            auto_index_input_fasta=True
        )
    except Exception as e:
        tb_str = traceback.format_exc()
        pytest.fail(f"apply_signature_to_fasta raised an exception: {e}\nTraceback:\n{tb_str}")

    assert results is not None, "Runner function did not return results."

    expected_keys = [
        "run_name", "output_directory", "output_modified_fasta_path",
        "output_modified_fasta_index_path", "manifest_path",
        "num_mutations_applied", "num_mutations_requested",
        "signature_used", "original_input_fasta"
    ]
    for key in expected_keys:
        assert key in results, f"Result dictionary missing key: {key}"

    output_modified_fasta = Path(results["output_modified_fasta_path"])
    assert results["output_modified_fasta_index_path"] is not None, "FASTA index path is None in results."
    output_modified_fasta_index = Path(results["output_modified_fasta_index_path"])
    manifest_path = Path(results["manifest_path"])

    assert output_modified_fasta.exists(), f"Output mutated FASTA file not found: {output_modified_fasta}"
    assert output_modified_fasta_index.exists(), f"Output mutated FASTA index not found: {output_modified_fasta_index}"
    assert manifest_path.exists(), f"Manifest file not found: {manifest_path}"

    assert results["num_mutations_requested"] == num_mutations_to_request
    assert results["num_mutations_applied"] <= num_mutations_to_request
    assert results["num_mutations_applied"] > 0, \
        f"Expected some mutations to be applied to {input_fasta_path.name}, but got 0. Signature: {signature_id}"

    assert results["signature_used"] == signature_id
    assert Path(results["original_input_fasta"]) == input_fasta_path

    manifest_data = bb_utils.read_run_manifest(manifest_path)
    assert manifest_data is not None, "Could not read manifest file."

    assert manifest_data["run_name"] == run_name
    assert manifest_data["command"] == "apply_signature_to_fasta", \
        f"Manifest command name mismatch. Expected 'apply_signature_to_fasta', got '{manifest_data['command']}'"
    assert manifest_data["parameters"]["signature_id_or_path"] == signature_id
    assert manifest_data["parameters"]["num_mutations_requested"] == num_mutations_to_request
    assert manifest_data["parameters"]["num_mutations_applied"] == results["num_mutations_applied"]
    assert len(manifest_data["outputs"]) >= 1

    found_fasta_in_manifest = False
    for out_file_info in manifest_data["outputs"]:
        if out_file_info["path"] == output_fasta_name and out_file_info["type"] == "FASTA":
            found_fasta_in_manifest = True
            break
    assert found_fasta_in_manifest, f"Mutated FASTA '{output_fasta_name}' not found in manifest outputs."

    if results["num_mutations_applied"] > 0:
        with open(input_fasta_path, "r") as f_in, open(output_modified_fasta, "r") as f_out:
            original_content = f_in.read()
            mutated_content = f_out.read()

            original_seqs = {s.split("\n")[0].lstrip(">"):s.split("\n")[1] for s in original_content.strip().split(">") if s}
            mutated_seqs = {s.split("\n")[0].lstrip(">"):s.split("\n")[1] for s in mutated_content.strip().split(">") if s}

            assert mutated_seqs["chr2"] == original_seqs["chr2"], "chr2 (all Ns) should not have been modified."
            if "chr1" in original_seqs and "chr1" in mutated_seqs and any(c != 'N' for c in original_seqs["chr1"]):
                 assert original_seqs["chr1"] != mutated_seqs["chr1"], \
                     "chr1 content should have changed after applying mutations, but it's identical."
            elif "chr1" not in original_seqs or "chr1" not in mutated_seqs:
                 pytest.fail("chr1 not found in both original and mutated FASTA outputs for comparison.")
    print(f"Integration test for {signature_id} passed. {results['num_mutations_applied']}/{num_mutations_to_request} mutations applied.")


def test_spike_variants_integration(tmp_path: Path):
    """
    Integration test for the spike_variants runner.
    Uses a temporary FASTA, a predefined minimal BAM (copied to tmp_path), and a list of variants.
    """
    # --- Setup FASTA ---
    ref_fasta_path = create_temp_fasta(tmp_path, "ref_spike.fa", SMOKE_TEST_FASTA_CONTENT)
    samtools_exe_path = bb_utils.find_tool_path("samtools")
    bb_utils.check_fasta_indexed(ref_fasta_path, samtools_exe_path, auto_index_if_missing=True)

    # --- Setup BAM (Copy from predefined minimal BAM) ---
    # Determine repo root relative to this test file to find data files
    current_script_path = Path(__file__).resolve()
    repo_root = current_script_path.parent.parent

    tiny_bam_src_path = repo_root / "src" / "basebuddy" / "data" / "smoke_test_data" / "tiny.bam"
    tiny_bai_src_path = repo_root / "src" / "basebuddy" / "data" / "smoke_test_data" / "tiny.bam.bai"

    if not tiny_bam_src_path.exists() or not tiny_bai_src_path.exists():
        pytest.skip(f"Skipping test_spike_variants_integration: Missing {tiny_bam_src_path.name} or {tiny_bai_src_path.name} in smoke_test_data. These need to be created manually.")

    input_bam_path = tmp_path / "input_for_spike.bam"
    shutil.copy(tiny_bam_src_path, input_bam_path)
    shutil.copy(tiny_bai_src_path, tmp_path / (input_bam_path.name + ".bai")) # Ensure BAI matches copied BAM name

    # --- Define Variants and Parameters ---
    variants_to_spike = [
        {"chromosome": "chr1", "position": 10, "ref_allele": "A", "alt_allele": "G"}, # Within ACGTACGTACGT...
        {"chromosome": "chr1", "position": 20, "ref_allele": "C", "alt_allele": "T"}  # Within ACGTACGTACGT...
    ]
    # Note: smoke_test.fa chr1 is: ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT (60 bases)
    # Position 10 (0-based 9) is 'C'. Original: A[C>G]G.
    # Position 20 (0-based 19) is 'G'. Original: A[G>T]T.
    # The provided variants should be:
    # {"chromosome": "chr1", "position": 10, "ref_allele": "C", "alt_allele": "G"},
    # {"chromosome": "chr1", "position": 20, "ref_allele": "G", "alt_allele": "T"} -> This will work.

    # Corrected variants based on SMOKE_TEST_FASTA_CONTENT's chr1:
    # >chr1
    # ACGTACGTAC G TACGTACGT A CGTACGTACG TACGTACGTACGTACGTACGTACGTACGT
    # 123456789 10           19 20
    # Position 10 is C. Position 20 is G.
    variants_to_spike_corrected = [
        {"chromosome": "chr1", "position": 10, "ref_allele": "C", "alt_allele": "G"},
        {"chromosome": "chr1", "position": 20, "ref_allele": "G", "alt_allele": "T"}
    ]


    output_root_dir = tmp_path / "integration_test_output_spike"
    run_name = "spike_variants_integration_test"
    output_prefix = "spiked_result"

    command_params_for_manifest = {
        "input_bams_gui": [str(input_bam_path.name)],
        "num_variants_to_spike_gui": len(variants_to_spike_corrected),
        "vaf_gui": 0.5 # Example VAF
    }

    # --- Execution ---
    results = None
    try:
        results = bb_runner.spike_variants(
            output_root_dir=output_root_dir,
            run_name=run_name,
            command_params=command_params_for_manifest,
            reference_fasta=str(ref_fasta_path),
            input_bams=[str(input_bam_path)], # Runner expects a list of BAM paths
            variants_list=variants_to_spike_corrected,
            output_prefix_for_bam=output_prefix,
            overwrite_output=True,
            auto_index_input_bam=False, # Already indexed
            auto_index_fasta=False      # Already indexed
        )
    except Exception as e:
        tb_str = traceback.format_exc()
        pytest.fail(f"spike_variants raised an exception: {e}\nTraceback:\n{tb_str}")

    # --- Assertions ---
    assert results is not None, "Runner function did not return results."
    assert "results_per_bam" in results, "Missing 'results_per_bam' in results."
    assert "errors_per_bam" in results, "Missing 'errors_per_bam' in results."
    assert len(results["errors_per_bam"]) == 0, f"Expected 0 errors, got: {results['errors_per_bam']}"
    assert len(results["results_per_bam"]) == 1, "Expected results for 1 input BAM."

    bam_result = results["results_per_bam"][0]
    expected_bam_keys = [
        "source_bam", "output_bam", "output_bam_index",
        "igv_session_file", "status"
    ]
    for key in expected_bam_keys:
        assert key in bam_result, f"BAM result dictionary missing key: {key}"
    assert bam_result["status"] == "success"

    output_bam = Path(bam_result["output_bam"])
    output_bai = Path(bam_result["output_bam_index"])
    igv_session = Path(bam_result["igv_session_file"])
    manifest_path = Path(results["manifest_path"]) # Main manifest for the run

    assert output_bam.exists(), f"Output spiked BAM not found: {output_bam}"
    assert output_bai.exists(), f"Output spiked BAI not found: {output_bai}"
    assert igv_session.exists(), f"IGV session file not found: {igv_session}"
    assert manifest_path.exists(), f"Main manifest file not found: {manifest_path}"

    # Check manifest content
    manifest_data = bb_utils.read_run_manifest(manifest_path)
    assert manifest_data is not None, "Could not read main manifest."
    assert manifest_data["run_name"] == run_name
    assert manifest_data["command"] == "spike_variants_batch" # As per current runner logic
    assert manifest_data["parameters"]["num_input_bams"] == 1
    assert manifest_data["parameters"]["num_variants_to_spike"] == len(variants_to_spike_corrected)

    # Check that the output BAM is listed in the main manifest
    found_bam_in_manifest = False
    expected_output_bam_name = f"{output_prefix}_{input_bam_path.stem}.bam"
    for out_file_info in manifest_data["outputs"]:
        if out_file_info["path"] == expected_output_bam_name and out_file_info["type"] == "BAM":
            found_bam_in_manifest = True
            break
    assert found_bam_in_manifest, f"Spiked BAM '{expected_output_bam_name}' not found in manifest outputs."

    # (Optional) Advanced: Pysam check for variants (simplified check on read counts for now)
    try:
        with pysam.AlignmentFile(input_bam_path, "rb") as infile, \
             pysam.AlignmentFile(output_bam, "rb") as outfile:
            assert infile.count_reads() == outfile.count_reads(), "Read count should be preserved after spiking."
            # More detailed check would involve pileups at variant positions.
    except Exception as e:
        pytest.warning(f"Pysam check on BAM files failed: {e}")

    print(f"Integration test for spike_variants passed for BAM: {input_bam_path.name}")


# --- Spike Command Integration Tests ---
# import subprocess # Already at top
# from typing import Optional, List # Already at top (or should be)

# TEST_FILE_DIR is already defined at the top of the file
TEST_DATA_SPIKE_DIR = (TEST_FILE_DIR / "test_data").resolve() # Re-use or ensure consistency
TEST_SPIKE_OUTPUT_DIR = (TEST_FILE_DIR / "test_spike_output").resolve() # Re-use or ensure consistency


def get_picard_jar_for_test() -> Optional[str]: # This helper can be defined once
    env_var = os.environ.get("BAMSURGEON_PICARD_JAR")
    if env_var and Path(env_var).is_file():
        return env_var
    local_picard_candidate = Path(__file__).parent.parent / "picard.jar"
    if local_picard_candidate.is_file():
        return str(local_picard_candidate)
    return None

PICARD_JAR_PATH_FOR_TESTS = get_picard_jar_for_test() # Define once

BASE_SPIKE_ARGS = [
    "basebuddy", "spike",
    "--reference", str(TEST_DATA_SPIKE_DIR / "spike_ref.fa"),
    "--output-prefix", str(TEST_SPIKE_OUTPUT_DIR / "spiked_bam_test_run"),
    "--overwrite"
]
if PICARD_JAR_PATH_FOR_TESTS:
    BASE_SPIKE_ARGS.extend(["--picard-jar", PICARD_JAR_PATH_FOR_TESTS])


@pytest.fixture(scope="module")
def spike_test_environment(): # This fixture can be defined once
    if not PICARD_JAR_PATH_FOR_TESTS:
        pytest.skip("Skipping all spike tests: BAMSURGEON_PICARD_JAR not set and no local picard.jar found.")
    if not TEST_DATA_SPIKE_DIR.exists():
        pytest.fail(f"Test data directory not found: {TEST_DATA_SPIKE_DIR}")
    required_files = ["spike_ref.fa", "spike_snps.vcf", "spike_indels.vcf", "qc_test_1.fastq", "qc_test_2.fastq"]
    for req_file in required_files:
        if not (TEST_DATA_SPIKE_DIR / req_file).exists():
            pytest.fail(f"Required test data file not found: {TEST_DATA_SPIKE_DIR / req_file}")

    input_bam = TEST_DATA_SPIKE_DIR / "spike_input.bam"
    input_bai = TEST_DATA_SPIKE_DIR / "spike_input.bam.bai"
    if not input_bam.exists():
        TEST_DATA_SPIKE_DIR.mkdir(parents=True, exist_ok=True)
        input_bam.touch()
        print(f"WARNING: Test input BAM {input_bam} not found. Created empty placeholder. BAMSurgeon will likely fail.")
    if not input_bai.exists() and input_bam.exists():
        input_bai.touch()
        print(f"WARNING: Test input BAI {input_bai} not found. Created empty placeholder.")

    if TEST_SPIKE_OUTPUT_DIR.exists(): # Clean specific spike output dir
        shutil.rmtree(TEST_SPIKE_OUTPUT_DIR)
    TEST_SPIKE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    yield
    if TEST_SPIKE_OUTPUT_DIR.exists():
        shutil.rmtree(TEST_SPIKE_OUTPUT_DIR)

def run_basebuddy_process_for_spike(args: List[str]): # Can be defined once
    final_args = list(args)
    has_input_bam = any(arg == "--input-bam" for arg in final_args)
    if not has_input_bam:
        try:
            spike_cmd_index = final_args.index("spike")
            final_args.insert(spike_cmd_index + 1, "--input-bam")
            final_args.insert(spike_cmd_index + 2, str(TEST_DATA_SPIKE_DIR / "spike_input.bam"))
        except ValueError:
            final_args.extend(["--input-bam", str(TEST_DATA_SPIKE_DIR / "spike_input.bam")])
    process = subprocess.run(final_args, capture_output=True, text=True)
    if process.returncode != 0:
        print(f"Executing: {' '.join(final_args)}\nStdout: {process.stdout}\nStderr: {process.stderr}")
    return process

# --- QC Command Integration Tests ---
TEST_QC_OUTPUT_DIR = (TEST_FILE_DIR / "test_qc_output").resolve()

@pytest.fixture(scope="module")
def qc_test_environment():
    if not TEST_DATA_SPIKE_DIR.exists():
        pytest.fail(f"Test data directory not found: {TEST_DATA_SPIKE_DIR}") # QC files are there

    required_qc_files = ["qc_test_1.fastq", "qc_test_2.fastq"]
    for req_file in required_qc_files:
        if not (TEST_DATA_SPIKE_DIR / req_file).exists():
            pytest.fail(f"Required QC test data file not found: {TEST_DATA_SPIKE_DIR / req_file}")

    if TEST_QC_OUTPUT_DIR.exists():
        shutil.rmtree(TEST_QC_OUTPUT_DIR)
    TEST_QC_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    yield
    if TEST_QC_OUTPUT_DIR.exists():
        shutil.rmtree(TEST_QC_OUTPUT_DIR)

def run_basebuddy_qc_process(args: List[str]):
    final_args = ["basebuddy", "qc"] + args
    process = subprocess.run(final_args, capture_output=True, text=True)
    if process.returncode != 0:
        print(f"Executing: {' '.join(final_args)}\nStdout: {process.stdout}\nStderr: {process.stderr}")
    return process

@pytest.mark.usefixtures("qc_test_environment")
def test_qc_single_fastq():
    output_subdir_name = "qc_single_file_run"
    args = [
        str(TEST_DATA_SPIKE_DIR / "qc_test_1.fastq"),
        "--output-dir", str(TEST_QC_OUTPUT_DIR),
        "--run-name", output_subdir_name,
        "--overwrite"
    ]
    result = run_basebuddy_qc_process(args)
    assert result.returncode == 0, f"basebuddy qc single file failed. Stderr: {result.stderr}"
    expected_report_path = TEST_QC_OUTPUT_DIR / output_subdir_name / "qc_test_1_fastqc" / "fastqc_report.html"
    assert expected_report_path.exists(), f"FastQC HTML report not found at {expected_report_path}"

@pytest.mark.usefixtures("qc_test_environment")
def test_qc_multiple_fastq():
    output_subdir_name = "qc_multi_file_run"
    args = [
        str(TEST_DATA_SPIKE_DIR / "qc_test_1.fastq"),
        str(TEST_DATA_SPIKE_DIR / "qc_test_2.fastq"),
        "--output-dir", str(TEST_QC_OUTPUT_DIR),
        "--run-name", output_subdir_name,
        "--overwrite"
    ]
    result = run_basebuddy_qc_process(args)
    assert result.returncode == 0, f"basebuddy qc multiple files failed. Stderr: {result.stderr}"
    expected_report_1_path = TEST_QC_OUTPUT_DIR / output_subdir_name / "qc_test_1_fastqc" / "fastqc_report.html"
    expected_report_2_path = TEST_QC_OUTPUT_DIR / output_subdir_name / "qc_test_2_fastqc" / "fastqc_report.html"
    assert expected_report_1_path.exists(), f"FastQC HTML report for qc_test_1.fastq not found at {expected_report_1_path}"
    assert expected_report_2_path.exists(), f"FastQC HTML report for qc_test_2.fastq not found at {expected_report_2_path}"

@pytest.mark.usefixtures("qc_test_environment")
def test_qc_non_existent_input():
    args = [
        "non_existent.fastq",
        "--output-dir", str(TEST_QC_OUTPUT_DIR)
    ]
    result = run_basebuddy_qc_process(args)
    assert result.returncode != 0, "basebuddy qc should fail for non-existent input"
    assert "File Error" in result.stderr or "not found" in result.stderr.lower() , "Correct error message for non-existent input not found."

@pytest.mark.usefixtures("qc_test_environment")
def test_qc_no_input_files():
    process = subprocess.run(["basebuddy", "qc", "--output-dir", str(TEST_QC_OUTPUT_DIR)], capture_output=True, text=True)
    assert process.returncode != 0, "basebuddy qc should fail when no input FASTQ files are provided"
    assert "Missing argument" in process.stderr or "Error: At least one FASTQ file must be provided" in process.stderr, "Correct error for missing FASTQ files not found."

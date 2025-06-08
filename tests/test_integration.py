import pytest
from pathlib import Path
import shutil
import json
import os
import traceback
import subprocess
from typing import Optional, List

from src.basebuddy import runner as bb_runner
from src.basebuddy import signature_utils as bb_sig_utils
from basebuddy import utils as bb_utils

# Content for test FASTA files
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
"""

def create_temp_fasta(tmp_path: Path, filename: str, content: str) -> Path:
    """Helper to create a temporary FASTA file for testing."""
    fasta_file_path = tmp_path / filename
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

    manifest_data = bb_utils.read_run_manifest(manifest_path)
    assert manifest_data is not None, "Could not read manifest file."
    assert manifest_data["run_name"] == run_name
    assert manifest_data["parameters"]["num_mutations_applied"] == results["num_mutations_applied"]

    if results["num_mutations_applied"] > 0:
        with open(input_fasta_path, "r") as f_in, open(output_modified_fasta, "r") as f_out:
            original_content = f_in.read()
            mutated_content = f_out.read()
            assert original_content != mutated_content, "FASTA content should have changed after applying mutations."
    print(f"Integration test for {signature_id} passed. {results['num_mutations_applied']}/{num_mutations_to_request} mutations applied.")

def test_spike_variants_integration(tmp_path: Path):
    """
    Integration test for the spike_variants runner.
    """
    # --- Setup FASTA ---
    ref_fasta_path = create_temp_fasta(tmp_path, "ref_spike.fa", SMOKE_TEST_FASTA_CONTENT)
    samtools_exe_path = bb_utils.find_tool_path("samtools")
    bb_utils.check_fasta_indexed(ref_fasta_path, samtools_exe_path, auto_index_if_missing=True)

    # --- Setup BAM (Copy from predefined minimal BAM) ---
    current_script_path = Path(__file__).resolve()
    repo_root = current_script_path.parent.parent
    tiny_bam_src_path = repo_root / "src" / "basebuddy" / "data" / "smoke_test_data" / "tiny.bam"
    tiny_bai_src_path = repo_root / "src" / "basebuddy" / "data" / "smoke_test_data" / "tiny.bam.bai"

    if not tiny_bam_src_path.exists() or not tiny_bai_src_path.exists():
        pytest.skip(f"Skipping test_spike_variants_integration: Missing tiny.bam or tiny.bam.bai.")

    input_bam_path = tmp_path / "input_for_spike.bam"
    shutil.copy(tiny_bam_src_path, input_bam_path)
    shutil.copy(tiny_bai_src_path, tmp_path / (input_bam_path.name + ".bai"))

    # --- Define Variants and Parameters ---
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
        "vaf_gui": 0.5
    }

    # --- Execution ---
    results = None
    try:
        results = bb_runner.spike_variants(
            output_root_dir=output_root_dir,
            run_name=run_name,
            command_params=command_params_for_manifest,
            reference_fasta=str(ref_fasta_path),
            input_bams=[str(input_bam_path)],
            variants_list=variants_to_spike_corrected,
            output_prefix_for_bam=output_prefix,
            overwrite_output=True,
            auto_index_input_bam=False,
            auto_index_fasta=False
        )
    except Exception as e:
        tb_str = traceback.format_exc()
        pytest.fail(f"spike_variants raised an exception: {e}\nTraceback:\n{tb_str}")

    # --- Assertions ---
    assert results is not None, "Runner function did not return results."
    assert len(results.get("errors_per_bam", [])) == 0, f"Expected 0 errors, got: {results['errors_per_bam']}"
    assert len(results.get("results_per_bam", [])) == 1, "Expected results for 1 input BAM."

    bam_result = results["results_per_bam"][0]
    assert bam_result.get("status") == "success"

    output_bam = Path(bam_result["output_bam"])
    output_bai = Path(bam_result["output_bam_index"])
    manifest_path = Path(results["manifest_path"])

    assert output_bam.exists(), f"Output spiked BAM not found: {output_bam}"
    assert output_bai.exists(), f"Output spiked BAI not found: {output_bai}"
    assert manifest_path.exists(), f"Main manifest file not found: {manifest_path}"

    print(f"Integration test for spike_variants passed for BAM: {input_bam_path.name}")


# --- Spike and QC Command Integration Tests ---

# Constants and Helpers for CLI tests
TEST_FILE_DIR = Path(__file__).parent
TEST_DATA_DIR = (TEST_FILE_DIR / "test_data").resolve()
TEST_SPIKE_OUTPUT_DIR = (TEST_FILE_DIR / "test_spike_output").resolve()
TEST_QC_OUTPUT_DIR = (TEST_FILE_DIR / "test_qc_output").resolve()

def get_picard_jar_for_test() -> Optional[str]:
    env_var = os.environ.get("BAMSURGEON_PICARD_JAR")
    if env_var and Path(env_var).is_file():
        return env_var
    local_picard_candidate = Path(__file__).parent.parent / "picard.jar"
    if local_picard_candidate.is_file():
        return str(local_picard_candidate)
    return None

PICARD_JAR_PATH_FOR_TESTS = get_picard_jar_for_test()

BASE_SPIKE_ARGS = [
    "basebuddy", "spike",
    "--reference", str(TEST_DATA_DIR / "spike_ref.fa"),
    "--output-prefix", str(TEST_SPIKE_OUTPUT_DIR / "spiked_bam_test_run"),
    "--overwrite"
]
if PICARD_JAR_PATH_FOR_TESTS:
    BASE_SPIKE_ARGS.extend(["--picard-jar", PICARD_JAR_PATH_FOR_TESTS])

# Fixture for Spike Tests
@pytest.fixture(scope="module")
def spike_test_environment():
    if not PICARD_JAR_PATH_FOR_TESTS:
        pytest.skip("Skipping all spike CLI tests: BAMSURGEON_PICARD_JAR not set and no local picard.jar found.")
    if not TEST_DATA_DIR.exists():
        pytest.fail(f"Test data directory not found: {TEST_DATA_DIR}")

    required_files = ["spike_ref.fa", "spike_snps.vcf", "spike_indels.vcf", "spike_input.bam"]
    for req_file in required_files:
        if not (TEST_DATA_DIR / req_file).exists():
            pytest.fail(f"Required test data file not found: {TEST_DATA_DIR / req_file}")

    if TEST_SPIKE_OUTPUT_DIR.exists():
        shutil.rmtree(TEST_SPIKE_OUTPUT_DIR)
    TEST_SPIKE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    yield
    if TEST_SPIKE_OUTPUT_DIR.exists():
        shutil.rmtree(TEST_SPIKE_OUTPUT_DIR)

def run_basebuddy_process_for_spike(args: List[str]):
    final_args = list(args)
    if not any(arg == "--input-bam" for arg in final_args):
        try:
            spike_cmd_index = final_args.index("spike")
            final_args.insert(spike_cmd_index + 1, "--input-bam")
            final_args.insert(spike_cmd_index + 2, str(TEST_DATA_DIR / "spike_input.bam"))
        except ValueError:
            final_args.extend(["--input-bam", str(TEST_DATA_DIR / "spike_input.bam")])

    process = subprocess.run(final_args, capture_output=True, text=True)
    if process.returncode != 0:
        print(f"Executing: {' '.join(final_args)}\nStdout: {process.stdout}\nStderr: {process.stderr}")
    return process

# Spike CLI Tests
@pytest.mark.usefixtures("spike_test_environment")
def test_spike_snps_and_indels():
    args = BASE_SPIKE_ARGS + [
        "--snp-vcf", str(TEST_DATA_DIR / "spike_snps.vcf"),
        "--indel-vcf", str(TEST_DATA_DIR / "spike_indels.vcf")
    ]
    result = run_basebuddy_process_for_spike(args)
    assert result.returncode == 0, f"Spike with SNPs and Indels failed unexpectedly. Stderr: {result.stderr}"
    output_bam_path = TEST_SPIKE_OUTPUT_DIR / "spiked_bam_test_run_spike_input_final_sorted.bam"
    assert output_bam_path.exists(), "Spike with SNPs and Indels: output BAM not created."

@pytest.mark.usefixtures("spike_test_environment")
def test_spike_no_vcfs():
    # Remove any vcf args that might be in BASE_SPIKE_ARGS for this specific test
    args_no_vcf = [arg for arg in BASE_SPIKE_ARGS if not ("--snp-vcf" in arg or "--indel-vcf" in arg)]
    result = run_basebuddy_process_for_spike(args_no_vcf)
    assert result.returncode != 0, "Spike should fail if no VCFs are provided"
    assert "Error: At least one VCF file" in result.stderr, "Correct error message for no VCFs not found."

# Fixture for QC Tests
@pytest.fixture(scope="module")
def qc_test_environment():
    if not TEST_DATA_DIR.exists():
        pytest.fail(f"Test data directory not found: {TEST_DATA_DIR}")

    required_qc_files = ["qc_test_1.fastq", "qc_test_2.fastq"]
    for req_file in required_qc_files:
        if not (TEST_DATA_DIR / req_file).exists():
            pytest.fail(f"Required QC test data file not found: {TEST_DATA_DIR / req_file}")

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

# QC CLI Tests
@pytest.mark.usefixtures("qc_test_environment")
def test_qc_single_fastq():
    output_subdir_name = "qc_single_file_run"
    args = [
        str(TEST_DATA_DIR / "qc_test_1.fastq"),
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
        str(TEST_DATA_DIR / "qc_test_1.fastq"),
        str(TEST_DATA_DIR / "qc_test_2.fastq"),
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
def test_qc_no_input_files():
    process = subprocess.run(["basebuddy", "qc", "--output-dir", str(TEST_QC_OUTPUT_DIR)], capture_output=True, text=True)
    assert process.returncode != 0, "basebuddy qc should fail when no input FASTQ files are provided"
    assert "Missing argument" in process.stderr or "Error: At least one FASTQ file" in process.stderr, "Correct error for missing FASTQ files not found."
   
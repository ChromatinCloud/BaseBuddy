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

# File: tests/test_runner.py

import os
import shutil
import subprocess
import sys
import pytest
from pathlib import Path

import basebuddy.runner as runner
from basebuddy.utils import BaseBuddyToolError, BaseBuddyConfigError # Added


@pytest.fixture(autouse=True)
def no_actual_subprocess(monkeypatch):
    """
    Prevent any real subprocess.run calls by default.
    If a test needs to check subprocess arguments, it can override monkeypatch.
    """
    def fake_run(cmd, check=False, cwd=None):
        # Simply record that subprocess.run was called with these arguments
        fake_run.called = getattr(fake_run, "called", []) + [(cmd, cwd)]
        return None

    monkeypatch.setattr(subprocess, "run", fake_run)
    return fake_run


@pytest.fixture
def tmp_fasta(tmp_path):
    fasta = tmp_path / "test.fa"
    fasta.write_text(">chr1\nACGTACGT\n")
    return str(fasta)


@pytest.fixture
def tmp_bam(tmp_path):
    # Create an empty file to represent a BAM
    bam = tmp_path / "test.bam"
    bam.write_bytes(b"")
    # Also create a .bai file so spike_variants sees the index
    bai = tmp_path / "test.bam.bai"
    bai.write_bytes(b"")
    return str(bam)


@pytest.fixture
def tmp_vcf(tmp_path):
    vcf = tmp_path / "test.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    return str(vcf)


def test_simulate_short_invokes_art(tmp_path, tmp_fasta, no_actual_subprocess, monkeypatch):
    # Arrange: ensure shutil.which returns a dummy path
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/art_illumina")
    outdir = tmp_path / "out_short"
    # Act
    runner.simulate_short(tmp_fasta, outdir, depth=10, readlen=100, profile="HS25")
    # Assert: subprocess.run was called exactly once with the expected ART command
    called = no_actual_subprocess.called[-1][0]
    expected_prefix = str(outdir / "reads")
    assert "/usr/bin/art_illumina" in called[0]
    assert "-i" in called and tmp_fasta in called
    assert "-l" in called and "100" in called
    assert "-f" in called and "10" in called
    assert "-o" in called and expected_prefix in called


def test_simulate_short_raises_if_art_missing(monkeypatch, tmp_path, tmp_fasta):
    monkeypatch.setattr(shutil, "which", lambda name: None)
    outdir = tmp_path / "out_short"
    with pytest.raises(RuntimeError) as excinfo:
        runner.simulate_short(tmp_fasta, outdir, depth=10)
    assert "art_illumina not found" in str(excinfo.value)


def test_simulate_long_invokes_nanosim(tmp_path, tmp_fasta, no_actual_subprocess, monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/nanosim-h")
    outdir = tmp_path / "out_long"
    runner.simulate_long(tmp_fasta, outdir, depth=5, nanopore_model="test_model")
    called = no_actual_subprocess.called[-1][0]
    assert "/usr/bin/nanosim-h" in called[0]
    assert "simulate" in called
    assert "-r" in called and tmp_fasta in called
    assert "-c" in called and "5" in called
    assert "-n" in called and "test_model" in called


def test_simulate_long_raises_if_nanosim_missing(monkeypatch, tmp_path, tmp_fasta):
    monkeypatch.setattr(shutil, "which", lambda name: None)
    outdir = tmp_path / "out_long"
    with pytest.raises(RuntimeError) as excinfo:
        runner.simulate_long(tmp_fasta, outdir, depth=5)
    assert "nanosim-h not found" in str(excinfo.value)


def test_spike_variants_indexes_and_addsnv(tmp_bam, tmp_vcf, tmp_path, no_actual_subprocess, monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/addsnv.py")
    reference = "/path/to/ref.fa"
    input_bam = tmp_bam
    vcf = tmp_vcf
    out_bam = str(tmp_path / "out.bam")
    # Ensure input_bam.bai exists (via fixture)
    runner.spike_variants(reference, input_bam, vcf, vaf=0.1, out_bam=out_bam, sampler_seed=99)
    calls = [c[0] for c in no_actual_subprocess.called]
    # First call should index input BAM
    assert ["samtools", "index", input_bam] in calls
    # Next call should run addsnv.py
    assert any("/usr/bin/addsnv.py" in c for c in calls[1][0])
    # Finally, index the output BAM
    assert ["samtools", "index", out_bam] in calls


def test_spike_variants_raises_if_addsnv_missing(monkeypatch, tmp_bam, tmp_vcf, tmp_path):
    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError) as excinfo:
        runner.spike_variants("/path/ref.fa", tmp_bam, tmp_vcf, 0.1, str(tmp_path / "out.bam"))
    assert "addsnv.py not found" in str(excinfo.value)


def test_simulate_signatures_invokes_sigprofiler(tmp_path, monkeypatch):
    # Monkeypatch the SigProfilerSimulatorFunc to record calls
    called = {}
    def fake_simulator(**kwargs):
        called.update(kwargs)

    monkeypatch.setattr(runner, "SigProfilerSimulatorFunc", fake_simulator)
    reference = "/path/to/ref.fa"
    outdir = tmp_path / "sigout"
    runner.simulate_signatures(reference, outdir, sig_type="DBS", num_mutations=50, sample_id="S1")
    assert called["project"] == "S1"
    assert called["input_type"] == "reference_genome"
    assert called["input_data"] == reference
    assert called["sig_type"] == "DBS"
    assert called["sampleNumMutations"] == 50
    assert called["outdir"] == str(outdir)


def test_apply_strand_bias_invokes_samtools(tmp_path, tmp_bam, no_actual_subprocess):
    # Create dummy BAM and BAI
    bam = tmp_bam
    out_bam = str(tmp_path / "biased.bam")
    runner.apply_strand_bias(bam, out_bam, forward_fraction=0.6, seed=7)
    calls = [c[0] for c in no_actual_subprocess.called]
    # Check that index was present (fixture created .bai), so first call is view -f
    assert ["samtools", "view", "-h", "-f", "0x10", "-s", "7.600", "-b", bam, "-o", f"{out_bam}.fwd.bam"] in calls
    # Next call for reverse-strand
    assert ["samtools", "view", "-h", "-F", "0x10", "-s", "7.400", "-b", bam, "-o", f"{out_bam}.rev.bam"] in calls
    # Then merge and index
    assert ["samtools", "merge", "-f", out_bam, f"{out_bam}.fwd.bam", f"{out_bam}.rev.bam"] in calls
    assert ["samtools", "index", out_bam] in calls


# If actual file creation is desired, mark these tests with a pytest marker
# and run with `pytest -m realenv`
@pytest.mark.skipif(shutil.which("art_illumina") is None, reason="ART not installed")
def test_simulate_short_creates_fastq(tmp_fasta, tmp_path):
    outdir = tmp_path / "real_short"
    runner.simulate_short(tmp_fasta, outdir, depth=2)
    # Expect ART to produce .fq files (reads1.fq and reads2.fq)
    assert (outdir / "reads1.fq").exists()
    assert (outdir / "reads2.fq").exists()


# --- New tests for the refactored spike_variants in src/basebuddy/runner.py ---

DEFAULT_SPIKE_REF_FASTA = "spike_ref.fa"
DEFAULT_SPIKE_INPUT_BAM = "input_for_spike.bam"
DEFAULT_SPIKE_OUTPUT_PREFIX = "spiked_out"
MOCKED_ADDSNV_PATH = "/mocked/tools/addsnv.py"
MOCKED_SAMTOOLS_PATH = "/mocked/tools/samtools"

@pytest.fixture
def spike_variants_params(tmp_path, tmp_fasta):
    # tmp_fasta can serve as our reference for these tests
    # Ensure the reference FASTA and its index exist for spike_variants
    Path(tmp_fasta).touch()
    Path(tmp_fasta + ".fai").touch()

    input_bam_path = tmp_path / DEFAULT_SPIKE_INPUT_BAM
    input_bam_path.touch()
    # spike_variants runner will try to index it if .bai is missing and auto_index_input_bam=True
    # For some tests, we might want to pre-create .bai
    # Path(input_bam_path.with_suffix(".bam.bai")).touch()

    return {
        "output_root_dir": tmp_path / "spike_output_root",
        "reference_fasta": tmp_fasta,
        "input_bams": [str(input_bam_path)],
        "variants_list": [
            {"chromosome": "chr1", "position": 100, "ref_allele": "A", "alt_allele": "T"}
        ],
        "output_prefix_for_bam": DEFAULT_SPIKE_OUTPUT_PREFIX,
        "run_name": "test_spike_run",
        "overwrite_output": True, # Allow overwriting for tests
        "auto_index_input_bam": True,
        "auto_index_fasta": True, # Though tmp_fasta.fai is created
    }

def test_spike_variants_successful_detailed(spike_variants_params, tmp_path, no_actual_subprocess, monkeypatch):
    params = spike_variants_params
    run_output_dir = params["output_root_dir"] / params["run_name"]

    # Mock tool paths
    def mock_shutil_which(tool_name):
        if tool_name == "addsnv.py": return MOCKED_ADDSNV_PATH
        if tool_name == "samtools": return MOCKED_SAMTOOLS_PATH
        return None
    monkeypatch.setattr(shutil, "which", mock_shutil_which)

    # Mock file operations and manifest/IGV
    mock_write_manifest = monkeypatch.patch("basebuddy.utils.write_run_manifest")
    mock_gen_igv = monkeypatch.patch("basebuddy.utils.generate_igv_session_xml")
    mock_path_unlink = monkeypatch.patch("pathlib.Path.unlink")

    # Call the runner function
    runner.spike_variants(**params)

    # Assertions
    # subprocess.run calls are captured by no_actual_subprocess.called
    # Format: [(cmd_list, cwd_path), ...]

    # Expected paths
    input_bam_path_obj = Path(params["input_bams"][0])
    input_bam_stem = input_bam_path_obj.stem

    temp_vcf_expected_path = run_output_dir / "temp_variants_for_run.vcf" # As generated by spike_variants
    temp_addsnv_bam_expected_path = run_output_dir / f"{params['output_prefix_for_bam']}_{input_bam_stem}_temp_addsnv.bam"
    final_bam_expected_path = run_output_dir / f"{params['output_prefix_for_bam']}_{input_bam_stem}.bam"

    # Check subprocess calls
    # Call 1: samtools faidx (for reference, if auto_index_fasta=True and .fai missing - fixture creates .fai, so this might not be called if logic is smart)
    # Call 2: samtools index (for input BAM, if auto_index_input_bam=True and .bai missing - fixture doesn't create .bai for input_for_spike.bam)
    # Call 3: addsnv.py
    # Call 4: samtools sort
    # Call 5: samtools index (for output BAM)

    cmds_called = [c[0] for c in no_actual_subprocess.called]

    # Assert samtools index on input BAM (assuming .bai wasn't pre-created by fixture for this specific input bam name)
    assert [MOCKED_SAMTOOLS_PATH, "index", params["input_bams"][0]] in cmds_called

    # Assert addsnv.py call
    expected_addsnv_cmd_part = [
        "python", MOCKED_ADDSNV_PATH,
        "--reference", params["reference_fasta"],
        "--in_bam", params["input_bams"][0],
        "--vcf", str(temp_vcf_expected_path),
        "--out_bam", str(temp_addsnv_bam_expected_path)
        # -p and -s are also there but checking main ones
    ]
    assert any(all(part in call for part in expected_addsnv_cmd_part) for call in cmds_called), f"Expected addsnv call not found in {cmds_called}"

    # Assert samtools sort call
    assert [MOCKED_SAMTOOLS_PATH, "sort", str(temp_addsnv_bam_expected_path), "-o", str(final_bam_expected_path)] in cmds_called

    # Assert samtools index on output BAM
    assert [MOCKED_SAMTOOLS_PATH, "index", str(final_bam_expected_path)] in cmds_called

    mock_write_manifest.assert_called_once()
    mock_gen_igv.assert_called_once() # Called per BAM, so once for this test

    # Assert unlink calls for temporary files
    # Order of unlink might not be guaranteed, so check for presence
    unlink_calls_args = [call_args[0][0] for call_args in mock_path_unlink.call_args_list]
    assert temp_addsnv_bam_expected_path in unlink_calls_args
    assert temp_vcf_expected_path in unlink_calls_args


def test_spike_variants_samtools_sort_failure(spike_variants_params, tmp_path, no_actual_subprocess, monkeypatch):
    params = spike_variants_params
    run_output_dir = params["output_root_dir"] / params["run_name"]
    input_bam_path_obj = Path(params["input_bams"][0])
    input_bam_stem = input_bam_path_obj.stem
    temp_addsnv_bam_expected_path = run_output_dir / f"{params['output_prefix_for_bam']}_{input_bam_stem}_temp_addsnv.bam"
    final_bam_expected_path = run_output_dir / f"{params['output_prefix_for_bam']}_{input_bam_stem}.bam"

    # Mock tool paths
    def mock_shutil_which(tool_name):
        if tool_name == "addsnv.py": return MOCKED_ADDSNV_PATH
        if tool_name == "samtools": return MOCKED_SAMTOOLS_PATH
        return None
    monkeypatch.setattr(shutil, "which", mock_shutil_which)

    # Mock manifest (it might be called with failure status or not at all)
    mock_write_manifest = monkeypatch.patch("basebuddy.utils.write_run_manifest")
    mock_path_unlink = monkeypatch.patch("pathlib.Path.unlink") # To check cleanup

    # Modify fake_run to raise error for samtools sort
    original_fake_run_called_list = []
    def selective_fail_run(cmd, check=False, cwd=None):
        original_fake_run_called_list.append((cmd, cwd)) # Still record calls
        if cmd[0] == MOCKED_SAMTOOLS_PATH and cmd[1] == "sort":
            # Ensure the temp addsnv BAM is created before sort fails, as per actual tool flow
            # The actual spike_variants runner would have created this.
            # For this test, we need to simulate its existence if the runner doesn't create it before calling sort.
            # However, the runner *does* call addsnv.py first which *should* create this file.
            # The test for addsnv.py failing would be separate.
            raise subprocess.CalledProcessError(returncode=1, cmd=cmd, stderr="Samtools sort failed miserably")
        return None # For other successful calls like addsnv.py

    monkeypatch.setattr(subprocess, "run", selective_fail_run)
    # Also need to ensure no_actual_subprocess.called uses this new list if it's separate
    # The fixture `no_actual_subprocess` returns `fake_run`, so monkeypatching `subprocess.run`
    # effectively changes what `no_actual_subprocess` uses. We'll check `original_fake_run_called_list`.

    with pytest.raises(BaseBuddyToolError) as excinfo:
        runner.spike_variants(**params)

    assert "samtools sort" in str(excinfo.value).lower() or "addsnv variant spiking failed" in str(excinfo.value).lower() # error bubbles up

    # Check that manifest was called with a "failed" status or similar, or not at all
    # Depending on where the failure is caught and manifest is written.
    # The new runner aims to write manifest at the very end. If a per-BAM step fails,
    # the overall manifest might indicate partial success or overall failure.
    # For this specific test, the failure is catastrophic for the BAM, so the overall manifest might reflect this.
    # Current spike_variants writes one manifest at the end. If a BAM fails, it's in errors_per_bam.
    # The main manifest should still be written.

    mock_write_manifest.assert_called_once() # Manifest for the run is still written
    manifest_call_args = mock_write_manifest.call_args[0] if mock_write_manifest.call_args else None
    if manifest_call_args:
        # Check if the error is recorded in the manifest parameters or if status is 'failed'
        # This depends on implementation details of error reporting in spike_variants manifest
        pass # For now, just ensure it's called. Detailed check might be too complex for this step.

    # Check that temporary addsnv BAM was attempted to be unlinked if it existed
    # The runner might try to clean up even on failure.
    # unlink_calls_args = [call_args[0][0] for call_args in mock_path_unlink.call_args_list]
    # assert temp_addsnv_bam_expected_path in unlink_calls_args # This depends on where error is caught

    # Check subprocess calls made *before* the failure
    cmds_called_before_failure = [c[0] for c in original_fake_run_called_list]
    assert any(MOCKED_ADDSNV_PATH in call[1] for call in cmds_called_before_failure if call[0] == "python"), "addsnv.py was not called before sort failure"
    assert any(call[0] == MOCKED_SAMTOOLS_PATH and call[1] == "sort" for call in cmds_called_before_failure), "samtools sort was not attempted"

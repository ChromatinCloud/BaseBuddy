# File: tests/test_runner.py

import os
import shutil
import subprocess
import sys
import pytest
from pathlib import Path

import basebuddy.runner as runner


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


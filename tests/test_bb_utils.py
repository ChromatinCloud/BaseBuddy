import pytest
from pathlib import Path
from unittest import mock

# Assuming bb_utils.py is in the same directory or accessible in PYTHONPATH
import bb_utils
import subprocess # For mocker.patch of subprocess.CompletedProcess etc. and type hints
import logging # For mocking loggers

class TestGenerateUniqueRunName:
    def test_generate_unique_run_name_default(self):
        run_name = bb_utils.generate_unique_run_name(command_name="test_command")
        assert isinstance(run_name, str)
        assert "test_command" in run_name
        assert len(run_name) > len("test_command")

    def test_generate_unique_run_name_with_prefix(self):
        command_name = "test_run"
        run_name = bb_utils.generate_unique_run_name(command_name=command_name)
        assert run_name.startswith(command_name + "_")

    def test_generate_unique_run_name_uniqueness(self):
        run_name1 = bb_utils.generate_unique_run_name(command_name="test")
        import time
        time.sleep(1.1)
        run_name2 = bb_utils.generate_unique_run_name(command_name="test")
        assert run_name1 != run_name2

class TestEnsureFileExists:
    def test_ensure_file_exists_when_file_exists(self, tmp_path):
        file_path = tmp_path / "test_file.txt"; file_path.touch()
        bb_utils.ensure_file_exists(str(file_path))
        assert file_path.is_file()

    def test_ensure_file_exists_when_file_does_not_exist(self, tmp_path):
        file_path = tmp_path / "non_existent_file.txt"
        with pytest.raises(bb_utils.BaseBuddyFileError):
            bb_utils.ensure_file_exists(str(file_path))

    def test_ensure_file_exists_when_path_is_directory(self, tmp_path):
        dir_path = tmp_path / "test_dir"; dir_path.mkdir()
        with pytest.raises(bb_utils.BaseBuddyFileError):
            bb_utils.ensure_file_exists(str(dir_path))

class TestEnsureDirectoryExists:
    def test_ensure_directory_exists_when_dir_exists(self, tmp_path):
        dir_path = tmp_path / "existing_dir"; dir_path.mkdir()
        bb_utils.ensure_directory_exists(str(dir_path))
        assert dir_path.is_dir()

    def test_ensure_directory_exists_does_not_exist_create_false(self, tmp_path):
        dir_path = tmp_path / "non_existent_dir"
        with pytest.raises(bb_utils.BaseBuddyFileError):
            bb_utils.ensure_directory_exists(str(dir_path), create=False)

    def test_ensure_directory_exists_does_not_exist_create_true(self, tmp_path):
        dir_path = tmp_path / "new_dir"
        bb_utils.ensure_directory_exists(str(dir_path), create=True)
        assert dir_path.is_dir()

    def test_ensure_directory_exists_when_path_is_file(self, tmp_path):
        file_path = tmp_path / "test_file.txt"; file_path.touch()
        with pytest.raises(bb_utils.BaseBuddyFileError):
            bb_utils.ensure_directory_exists(str(file_path))

    def test_ensure_directory_exists_when_path_is_file_and_create_true(self, tmp_path):
        file_path = tmp_path / "another_test_file.txt"; file_path.touch()
        with pytest.raises(bb_utils.BaseBuddyFileError):
            bb_utils.ensure_directory_exists(str(file_path), create=True)

class TestRunExternalCmd:
    def test_run_external_cmd_success(self, mocker):
        mock_run = mocker.patch('subprocess.run')
        mock_run.return_value = subprocess.CompletedProcess(args=['ls', '-l'], returncode=0, stdout="total 0", stderr="")
        cmd = ['ls', '-l']
        return_code, stdout, stderr = bb_utils.run_external_cmd(cmd)
        assert return_code == 0
        assert stdout == "total 0"
        assert stderr == ""
        mock_run.assert_called_once_with(
            cmd, cwd=None, capture_output=True, text=True,
            encoding='utf-8', errors='replace', timeout=None, env=None, check=False
        )

    def test_run_external_cmd_failure(self, mocker):
        mock_run = mocker.patch('subprocess.run')
        mock_run.return_value = subprocess.CompletedProcess(args=['error_cmd'], returncode=1, stdout="", stderr="Error occurred")
        cmd = ['error_cmd']
        with pytest.raises(bb_utils.BaseBuddyToolError) as excinfo:
            bb_utils.run_external_cmd(cmd)
        assert "External command execution failed" in str(excinfo.value)
        mock_run.assert_called_once()

    def test_run_external_cmd_timeout(self, mocker):
        mock_run = mocker.patch('subprocess.run', side_effect=subprocess.TimeoutExpired(cmd=['sleep', '5'], timeout=5))
        cmd = ['sleep', '5']
        with pytest.raises(bb_utils.BaseBuddyToolError) as excinfo:
            bb_utils.run_external_cmd(cmd, timeout_seconds=5)
        assert "Command execution timed out" in str(excinfo.value)
        mock_run.assert_called_once()

    def test_run_external_cmd_file_not_found(self, mocker):
        mock_run = mocker.patch('subprocess.run', side_effect=FileNotFoundError("Command not found"))
        cmd = ['non_existent_cmd']
        with pytest.raises(bb_utils.BaseBuddyConfigError) as excinfo:
            bb_utils.run_external_cmd(cmd)
        assert "The command 'non_existent_cmd' was not found" in str(excinfo.value)
        mock_run.assert_called_once()

    def test_run_external_cmd_stream_output_success(self, mocker):
        mock_popen_instance = mocker.Mock()
        mock_popen_instance.stdout.readline.side_effect = ["output line 1\n", "output line 2\n", ""]
        mock_popen_instance.stderr.readline.side_effect = [""]
        mock_popen_instance.wait.return_value = None
        mock_popen_instance.returncode = 0
        mock_popen = mocker.patch('subprocess.Popen', return_value=mock_popen_instance)
        mocker.patch('logging.info')
        cmd = ['streaming_cmd']
        return_code, stdout, stderr = bb_utils.run_external_cmd(cmd, stream_output=True, capture_output=False)
        assert return_code == 0
        assert stdout == "output line 1\noutput line 2"
        assert stderr == ""
        mock_popen.assert_called_once_with(
            cmd, cwd=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, encoding='utf-8', errors='replace', env=None, bufsize=1
        )

    def test_run_external_cmd_stream_output_failure(self, mocker):
        mock_popen_instance = mocker.Mock()
        mock_popen_instance.stdout.readline.side_effect = ["output\n", ""]
        mock_popen_instance.stderr.readline.side_effect = ["error line\n", ""]
        mock_popen_instance.wait.return_value = None
        mock_popen_instance.returncode = 1
        mock_popen = mocker.patch('subprocess.Popen', return_value=mock_popen_instance)
        mocker.patch('logging.error')
        cmd = ['failing_streaming_cmd']
        with pytest.raises(bb_utils.BaseBuddyToolError) as excinfo:
            bb_utils.run_external_cmd(cmd, stream_output=True, capture_output=False)
        assert "External command execution failed" in str(excinfo.value)
        mock_popen.assert_called_once()

class TestVerifyFileChecksum:
    def test_verify_checksum_matching_sha256(self, tmp_path):
        file_content = b"Hello BaseBuddy!"
        file_path = tmp_path / "test_file.txt"
        file_path.write_bytes(file_content)
        expected_checksum = "39b730532f70c73ca071cc36b04a1a981361272181e5391a8ee4181f4d7d59b9"
        bb_utils.verify_file_checksum(file_path, expected_checksum)

    def test_verify_checksum_non_matching_sha256(self, tmp_path):
        file_content = b"Hello BaseBuddy!"
        file_path = tmp_path / "test_file.txt"
        file_path.write_bytes(file_content)
        incorrect_checksum = "incorrectchecksum123"
        with pytest.raises(bb_utils.BaseBuddyChecksumError):
            bb_utils.verify_file_checksum(file_path, incorrect_checksum)

    def test_verify_checksum_file_not_found(self, tmp_path):
        non_existent_file = tmp_path / "non_existent.txt"
        with pytest.raises(bb_utils.BaseBuddyFileError):
            bb_utils.verify_file_checksum(non_existent_file, "dummy")

    def test_verify_checksum_default_algorithm_is_sha256(self, tmp_path):
        file_content = b"Test sha256 default"
        file_path = tmp_path / "sha_test.txt"
        file_path.write_bytes(file_content)
        import hashlib
        expected_sha256 = hashlib.sha256(file_content).hexdigest()
        bb_utils.verify_file_checksum(file_path, expected_sha256)
        expected_md5 = hashlib.md5(file_content).hexdigest()
        with pytest.raises(bb_utils.BaseBuddyChecksumError):
            bb_utils.verify_file_checksum(file_path, expected_md5)

class TestManifestOperations:
    def test_write_and_read_run_manifest(self, tmp_path):
        manifest_path = tmp_path / "run_manifest.json"
        run_name, command_name = "test_run_123", "my_command"
        parameters = {"param1": "value1", "input_file": Path("/tmp/input.fa")}
        output_files = [{"name": "output_bam", "path": "results/out.bam", "type": "BAM"}]
        ref_path = "/refs/hg38.fa"
        bb_utils.write_run_manifest(manifest_path, run_name, command_name, parameters, output_files, ref_path)
        read_data = bb_utils.read_run_manifest(manifest_path)
        assert read_data["run_name"] == run_name
        assert read_data["parameters"]["input_file"] == "/tmp/input.fa"

    def test_read_run_manifest_non_existent_file(self, tmp_path, mocker):
        mock_log_warning = mocker.patch('bb_utils.logger.warning')
        assert bb_utils.read_run_manifest(tmp_path / "ghost.json") is None
        mock_log_warning.assert_called_once()

    def test_read_run_manifest_malformed_json(self, tmp_path, mocker):
        mock_log_error = mocker.patch('bb_utils.logger.error')
        malformed_file = tmp_path / "bad.json"; malformed_file.write_text("{")
        assert bb_utils.read_run_manifest(malformed_file) is None
        mock_log_error.assert_called_once()

import xml.etree.ElementTree as ET

class TestCheckFastaIndexed:
    SAMTOOLS_PATH = "samtools"

    def test_index_exists(self, tmp_path):
        fasta_file = tmp_path / "ref.fasta"; fasta_file.touch()
        (tmp_path / "ref.fasta.fai").touch()
        bb_utils.check_fasta_indexed(fasta_file, self.SAMTOOLS_PATH)

    def test_index_missing_no_auto_index(self, tmp_path):
        fasta_file = tmp_path / "ref.fasta"; fasta_file.touch()
        expected_msg = f"Reference FASTA is not indexed \\(\\.fai missing\\): {fasta_file}\\. Please run `{self.SAMTOOLS_PATH} faidx {fasta_file}` first or enable auto-indexing\\."
        with pytest.raises(bb_utils.BaseBuddyFileError, match=expected_msg):
            bb_utils.check_fasta_indexed(fasta_file, self.SAMTOOLS_PATH, auto_index_if_missing=False)

    def test_index_missing_auto_index_success(self, tmp_path, mocker):
        fasta_file = tmp_path / "ref.fasta"; fasta_file.touch()
        fai_path = fasta_file.with_suffix(fasta_file.suffix + ".fai")
        assert not fai_path.exists()

        def mock_samtools_faidx_creates_fai(cmd_list, cwd=None, **kwargs):
            created_fai_path = Path(cmd_list[2] + ".fai")
            created_fai_path.touch()
            return (0, "FASTA index created", "")

        mock_run_cmd = mocker.patch('bb_utils.run_external_cmd', side_effect=mock_samtools_faidx_creates_fai)

        bb_utils.check_fasta_indexed(fasta_file, self.SAMTOOLS_PATH, auto_index_if_missing=True)

        mock_run_cmd.assert_called_once_with([self.SAMTOOLS_PATH, "faidx", str(fasta_file)], cwd=fasta_file.parent)
        assert fai_path.exists()

    def test_index_missing_auto_index_samtools_fails(self, tmp_path, mocker):
        fasta_file = tmp_path / "ref.fasta"; fasta_file.touch()
        fai_path = fasta_file.with_suffix(fasta_file.suffix + ".fai")
        assert not fai_path.exists()

        mock_run_cmd = mocker.patch('bb_utils.run_external_cmd',
                                    side_effect=bb_utils.BaseBuddyToolError("samtools failed", command=["samtools", "faidx"]))

        with pytest.raises(bb_utils.BaseBuddyFileError, match="Failed to auto-index FASTA file"):
            bb_utils.check_fasta_indexed(fasta_file, self.SAMTOOLS_PATH, auto_index_if_missing=True)

        mock_run_cmd.assert_called_once_with([self.SAMTOOLS_PATH, "faidx", str(fasta_file)], cwd=fasta_file.parent)
        assert not fai_path.exists()

    def test_reference_fasta_missing(self, tmp_path, mocker):
        fasta_file_not_existing = tmp_path / "non_existent.fasta"
        assert not fasta_file_not_existing.exists() # Ensure it truly doesn't exist

        mock_run_cmd = mocker.patch('bb_utils.run_external_cmd',
                                    side_effect=bb_utils.BaseBuddyToolError(
                                        message="samtools faidx failed on non-existent file",
                                        command=[self.SAMTOOLS_PATH, "faidx", str(fasta_file_not_existing)]
                                    ))

        with pytest.raises(bb_utils.BaseBuddyFileError, match="Failed to auto-index FASTA file"):
             bb_utils.check_fasta_indexed(fasta_file_not_existing, self.SAMTOOLS_PATH, auto_index_if_missing=True)

        mock_run_cmd.assert_called_once_with([self.SAMTOOLS_PATH, "faidx", str(fasta_file_not_existing)], cwd=fasta_file_not_existing.parent)

class TestCheckBamIndexed:
    SAMTOOLS_PATH = "samtools_bam"

    def test_bam_index_exists_bam_bai(self, tmp_path):
        bam_file = tmp_path / "align.bam"; bam_file.touch()
        (tmp_path / "align.bam.bai").touch()
        bb_utils.check_bam_indexed(bam_file, self.SAMTOOLS_PATH)

    def test_bam_index_exists_bai(self, tmp_path):
        bam_file = tmp_path / "align.bam"; bam_file.touch()
        (tmp_path / "align.bai").touch()
        bb_utils.check_bam_indexed(bam_file, self.SAMTOOLS_PATH)

    def test_bam_index_missing_no_auto_index(self, tmp_path):
        bam_file = tmp_path / "align.bam"; bam_file.touch()
        expected_msg = f"Input BAM is not indexed \\(\\.bai missing\\): {bam_file}\\. Please run `{self.SAMTOOLS_PATH} index {bam_file}` or use an auto-index option\\."
        with pytest.raises(bb_utils.BaseBuddyFileError, match=expected_msg):
            bb_utils.check_bam_indexed(bam_file, self.SAMTOOLS_PATH, auto_index_if_missing=False)

    def test_bam_index_missing_auto_index_success(self, tmp_path, mocker):
        bam_file = tmp_path / "align.bam"; bam_file.touch()
        bai_path = bam_file.with_suffix(bam_file.suffix + ".bai") # .bam.bai
        assert not bai_path.exists()

        def mock_samtools_index_creates_bai(cmd_list, **kwargs):
            created_bai_path = Path(cmd_list[2]).with_suffix(Path(cmd_list[2]).suffix + ".bai")
            created_bai_path.touch()
            return (0, "BAM index created", "")

        mock_run_cmd = mocker.patch('bb_utils.run_external_cmd', side_effect=mock_samtools_index_creates_bai)

        bb_utils.check_bam_indexed(bam_file, self.SAMTOOLS_PATH, auto_index_if_missing=True)

        mock_run_cmd.assert_called_once_with([self.SAMTOOLS_PATH, "index", str(bam_file)])
        assert bai_path.exists()

    def test_bam_index_missing_auto_index_samtools_fails(self, tmp_path, mocker):
        bam_file = tmp_path / "align.bam"; bam_file.touch()
        bai_path = bam_file.with_suffix(bam_file.suffix + ".bai")
        assert not bai_path.exists()

        mock_run_cmd = mocker.patch('bb_utils.run_external_cmd',
                                    side_effect=bb_utils.BaseBuddyToolError("samtools failed", command=["samtools", "index"]))

        with pytest.raises(bb_utils.BaseBuddyFileError, match="Failed to auto-index BAM file"):
            bb_utils.check_bam_indexed(bam_file, self.SAMTOOLS_PATH, auto_index_if_missing=True)

        mock_run_cmd.assert_called_once_with([self.SAMTOOLS_PATH, "index", str(bam_file)])
        assert not bai_path.exists()

class TestGenerateIGVSessionXML:
    def test_generate_igv_session_success(self, tmp_path):
        session_file = tmp_path / "my_session.xml"
        genome_fasta_file = tmp_path / "genome/ref.fasta"
        genome_fasta_file.parent.mkdir(); genome_fasta_file.touch()
        abs_genome_path = str(genome_fasta_file.resolve())
        tracks = [{"name": "My BAM", "path": "alignments/my.bam"}]

        bb_utils.generate_igv_session_xml(session_file, str(genome_fasta_file), tracks, initial_locus="chr1:100-200")

        assert session_file.is_file()
        tree = ET.parse(session_file); root = tree.getroot()
        assert root.tag == "Global"
        assert root.get("genome") == abs_genome_path
        assert root.get("locus") == "chr1:100-200"
        resources_node = root.find("Resources")
        assert resources_node is not None
        xml_tracks = resources_node.findall("Resource")
        assert len(xml_tracks) == len(tracks)
        assert xml_tracks[0].get("name") == tracks[0]["name"]
        assert xml_tracks[0].get("path") == tracks[0]["path"]

    def test_generate_igv_session_genome_path_resolves(self, tmp_path):
        session_file = tmp_path / "session_resolve.xml"
        genome_dir = tmp_path / "genomes"; genome_dir.mkdir()
        relative_genome_path = genome_dir / "human.fasta"; relative_genome_path.touch()
        abs_genome_path = str(relative_genome_path.resolve())

        bb_utils.generate_igv_session_xml(session_file, str(relative_genome_path), [])
        tree = ET.parse(session_file)
        assert tree.getroot().get("genome") == abs_genome_path

    def test_generate_igv_session_io_error_on_write(self, tmp_path, mocker):
        mock_tree_write = mocker.patch.object(ET.ElementTree, 'write', side_effect=IOError("Disk full"))
        with pytest.raises(bb_utils.BaseBuddyFileError, match="Could not write IGV session file"):
            bb_utils.generate_igv_session_xml(tmp_path / "s.xml", "dummy.fa", [])
        mock_tree_write.assert_called_once()

import pytest
from pathlib import Path
from unittest import mock
import copy
import shutil

from src.basebuddy import runner as src_bb_runner # For testing refactored runner
import bb_runners # For testing original runners
import bb_utils

# Default parameters for tests
DEFAULT_ART_PATH = "dummy_art_illumina"
DEFAULT_SAMTOOLS_PATH = "dummy_samtools"
DEFAULT_ADDSNV_PATH = "dummy_addsnv.py"
DEFAULT_CURL_PATH = "dummy_curl"
DEFAULT_REFERENCE_FASTA = "reference.fasta"
DEFAULT_RUN_OUTPUT_DIR_NAME = "test_run"

@pytest.fixture
def default_src_simulate_params(tmp_path): # Renamed fixture
    output_root = tmp_path / "output_data_simulate_src"
    ref_fasta_path_str = str(tmp_path / DEFAULT_REFERENCE_FASTA)
    Path(ref_fasta_path_str).parent.mkdir(parents=True, exist_ok=True)
    Path(ref_fasta_path_str).touch()

    return {
        "output_root_dir": output_root,
        "reference_fasta": ref_fasta_path_str,
        "depth": 50, "read_length": 150, "art_profile": "HS25",
        "run_name": DEFAULT_RUN_OUTPUT_DIR_NAME + "_src_simulate",
        "mean_fragment_length": 400, "std_dev_fragment_length": 50,
        "is_paired_end": True, "overwrite_output": False,
        "art_platform": "illumina", "timeout": 3600.0,
        "auto_index_fasta": True,
        "command_params": {
            "art_exe_path_param": DEFAULT_ART_PATH,
            "samtools_exe_path_param": DEFAULT_SAMTOOLS_PATH,
            "id_prefix": "sim_reads_src", "no_aln_output": False,
            "reference_fasta": ref_fasta_path_str, "depth": 50, "read_length": 150,
            "art_profile": "HS25", "mean_fragment_length": 400, "std_dev_fragment_length": 50,
            "is_paired_end": True, "art_platform": "illumina", "overwrite_output": False,
            "auto_index_fasta": True, # Ensure this is also in command_params if simulate_short copies it
            # Keep run_name and output_root_dir here if they are part of command_params for manifest
            "output_root_dir": str(output_root), "run_name": DEFAULT_RUN_OUTPUT_DIR_NAME + "_src_simulate",
        }
    }

class TestSrcSimulateShortRunner: # Renamed class
    def test_successful_paired_end_simulation(self, tmp_path, mocker, default_src_simulate_params):
        run_params = default_src_simulate_params
        expected_run_output_dir = run_params["output_root_dir"] / run_params["run_name"]
        expected_run_output_dir.mkdir(parents=True, exist_ok=True)
        id_prefix = run_params["command_params"]["id_prefix"]

        mock_find_tool = mocker.patch('bb_utils.find_tool_path', side_effect=lambda tool_name: run_params["command_params"]["art_exe_path_param"] if tool_name.startswith('art_') else run_params["command_params"]["samtools_exe_path_param"])
        mock_ensure_file_exists_util = mocker.patch('bb_utils.ensure_file_exists', return_value=Path(run_params["reference_fasta"]))
        mock_check_fasta_indexed_util = mocker.patch('bb_utils.check_fasta_indexed')
        mock_prepare_dir_util = mocker.patch('bb_utils.prepare_run_output_dir', return_value=expected_run_output_dir)

        def mock_art_cmd_effect(cmd_list, **kwargs):
            output_prefix_path_str = next(cmd_list[i+1] for i, arg in enumerate(cmd_list) if arg == "-o")
            prefix_path = Path(output_prefix_path_str)
            assert prefix_path == expected_run_output_dir / id_prefix
            (prefix_path.parent / (prefix_path.name + "1.fq")).touch()
            if run_params["is_paired_end"]: (prefix_path.parent / (prefix_path.name + "2.fq")).touch()
            if not run_params["command_params"]["no_aln_output"]: (prefix_path.parent / (prefix_path.name + ".aln")).touch()
            return (0, "ART success", "")
        mock_run_external_cmd_util = mocker.patch('bb_utils.run_external_cmd', side_effect=mock_art_cmd_effect)
        mock_write_manifest_util = mocker.patch('bb_utils.write_run_manifest')

        # Pass all params; simulate_short will pick what it needs
        src_bb_runner.simulate_short(**run_params)

        mock_prepare_dir_util.assert_called_once_with(Path(run_params["output_root_dir"]), run_params["run_name"], run_params["overwrite_output"])
        mock_find_tool.assert_any_call(f"art_{run_params['art_platform']}")
        mock_find_tool.assert_any_call("samtools")
        mock_ensure_file_exists_util.assert_any_call(run_params["reference_fasta"], "Reference FASTA")
        mock_check_fasta_indexed_util.assert_called_once_with(Path(run_params["reference_fasta"]), run_params["command_params"]["samtools_exe_path_param"], auto_index_if_missing=run_params["auto_index_fasta"])

        mock_run_external_cmd_util.assert_called_once()
        art_call_args = mock_run_external_cmd_util.call_args[0][0]
        assert str(expected_run_output_dir / id_prefix) in art_call_args

        expected_fq1 = expected_run_output_dir / (id_prefix + "1.fq")
        expected_fq2 = expected_run_output_dir / (id_prefix + "2.fq")
        mock_ensure_file_exists_util.assert_any_call(expected_fq1, "Expected ART output R1 FASTQ")
        mock_ensure_file_exists_util.assert_any_call(expected_fq2, "Expected ART output R2 FASTQ")

        mock_write_manifest_util.assert_called_once()
        manifest_kall = mock_write_manifest_util.call_args
        assert manifest_kall.kwargs['run_name'] == run_params["run_name"]
        assert manifest_kall.kwargs['command_name'] == f"simulate_short_reads_{run_params['art_platform']}"
        assert manifest_kall.kwargs['reference_genome_path'] == run_params["reference_fasta"]
        assert manifest_kall.kwargs['parameters']["id_prefix"] == id_prefix

    def test_successful_single_end_simulation(self, tmp_path, mocker, default_src_simulate_params):
        run_params = copy.deepcopy(default_src_simulate_params)
        run_params["is_paired_end"] = False
        if "command_params" in run_params: run_params["command_params"]["is_paired_end"] = False

        expected_run_output_dir = run_params["output_root_dir"] / run_params["run_name"]
        expected_run_output_dir.mkdir(parents=True, exist_ok=True)
        id_prefix = run_params["command_params"]["id_prefix"]
        expected_fq_path = expected_run_output_dir / (id_prefix + ".fq")

        mocker.patch('bb_utils.find_tool_path', side_effect=lambda tool_name: run_params["command_params"]["art_exe_path_param"] if tool_name.startswith('art_') else run_params["command_params"]["samtools_exe_path_param"])
        mock_ensure_file_exists_util = mocker.patch('bb_utils.ensure_file_exists', return_value=Path(run_params["reference_fasta"]))
        mocker.patch('bb_utils.check_fasta_indexed')
        mocker.patch('bb_utils.prepare_run_output_dir', return_value=expected_run_output_dir)

        def mock_art_cmd_effect_single(cmd_list, **kwargs):
            output_prefix_path_str = next(cmd_list[i+1] for i, arg in enumerate(cmd_list) if arg == "-o")
            prefix_path = Path(output_prefix_path_str)
            (prefix_path.parent / (prefix_path.name + ".fq")).touch() # Single-end
            if not run_params["command_params"]["no_aln_output"]: (prefix_path.parent / (prefix_path.name + ".aln")).touch()
            return (0, "ART SE success", "")
        mock_run_external_cmd_util = mocker.patch('bb_utils.run_external_cmd', side_effect=mock_art_cmd_effect_single)
        mock_write_manifest_util = mocker.patch('bb_utils.write_run_manifest')

        src_bb_runner.simulate_short(**run_params)

        art_call_args_list = mock_run_external_cmd_util.call_args[0][0]
        assert "-p" not in art_call_args_list
        mock_ensure_file_exists_util.assert_any_call(expected_fq_path, "Expected ART output FASTQ")
        manifest_params_written = mock_write_manifest_util.call_args.kwargs['parameters']
        assert manifest_params_written["is_paired_end"] == False

    def test_art_tool_failure(self, tmp_path, mocker, default_src_simulate_params):
        run_params = default_src_simulate_params
        expected_run_output_dir = run_params["output_root_dir"] / run_params["run_name"]
        expected_run_output_dir.mkdir(parents=True, exist_ok=True)

        mocker.patch('bb_utils.find_tool_path', side_effect=lambda tool_name: run_params["command_params"]["art_exe_path_param"] if tool_name.startswith('art_') else run_params["command_params"]["samtools_exe_path_param"])
        mocker.patch('bb_utils.ensure_file_exists', return_value=Path(run_params["reference_fasta"]))
        mocker.patch('bb_utils.check_fasta_indexed')
        mocker.patch('bb_utils.prepare_run_output_dir', return_value=expected_run_output_dir)
        mock_run_external_cmd_util = mocker.patch('bb_utils.run_external_cmd', side_effect=bb_utils.BaseBuddyToolError("ART sim failed", command=["art_exe"]))
        mock_write_manifest_util = mocker.patch('bb_utils.write_run_manifest')

        with pytest.raises(bb_utils.BaseBuddyToolError, match="ART sim failed"):
            src_bb_runner.simulate_short(**run_params)
        mock_run_external_cmd_util.assert_called_once()
        mock_write_manifest_util.assert_not_called()

    def test_missing_art_tool(self, tmp_path, mocker, default_src_simulate_params):
        run_params = default_src_simulate_params
        def find_tool_effect(tool_name):
            if tool_name.startswith('art_'): raise bb_utils.BaseBuddyConfigError("ART not found")
            return run_params["command_params"]["samtools_exe_path_param"]
        mocker.patch('bb_utils.find_tool_path', side_effect=find_tool_effect) # Patched in bb_utils
        mocker.patch('bb_utils.prepare_run_output_dir', return_value=run_params["output_root_dir"] / run_params["run_name"])
        mocker.patch('bb_utils.ensure_file_exists')
        mocker.patch('bb_utils.check_fasta_indexed')
        mocker.patch('bb_utils.run_external_cmd')

        with pytest.raises(bb_utils.BaseBuddyConfigError, match="ART not found"):
            src_bb_runner.simulate_short(**run_params)

    def test_invalid_input_depth(self, tmp_path, default_src_simulate_params):
        run_params = copy.deepcopy(default_src_simulate_params); run_params["depth"] = 0
        if "command_params" in run_params: run_params["command_params"]["depth"] = 0

        with pytest.raises(bb_utils.BaseBuddyInputError, match="Sequencing depth must be a positive integer"):
            src_bb_runner.simulate_short(**run_params)

@pytest.fixture
def default_spike_params(tmp_path): # This fixture remains for the old TestSpikeVariantsRunner
    output_root = tmp_path / "output_data_spike"; ref_fasta_path = tmp_path / "ref_spike.fasta"; ref_fasta_path.touch()
    input_bam_path = tmp_path / "input.bam"; input_bam_path.touch(); (tmp_path / "input.bam.bai").touch()
    vcf_file_path = tmp_path / "variants.vcf"; vcf_file_path.touch()
    return {
        "output_root_dir": output_root, "run_name": DEFAULT_RUN_OUTPUT_DIR_NAME + "_spike",
        "reference_fasta": str(ref_fasta_path), "input_bam": str(input_bam_path), "vcf_file": str(vcf_file_path),
        "output_bam_prefix_rel": "spiked_output", "overwrite_output": False, "auto_index_input_bam": False,
        "timeout": 7200.0,
        "command_params": {"addsnv_path": DEFAULT_ADDSNV_PATH, "samtools_path": DEFAULT_SAMTOOLS_PATH,
                           "reference_fasta": str(ref_fasta_path), "input_bam": str(input_bam_path), "vcf_file": str(vcf_file_path),
                           "output_bam_prefix_rel": "spiked_output", "overwrite_output": False, "auto_index_input_bam": False,
                           "output_root_dir": str(output_root), "run_name": DEFAULT_RUN_OUTPUT_DIR_NAME + "_spike",}}

class TestSpikeVariantsRunner: # This tests the root bb_runners.py version
    def test_successful_spike(self, tmp_path, mocker, default_spike_params):
        run_params = default_spike_params
        run_output_dir = run_params["output_root_dir"] / run_params["run_name"]
        run_output_dir.mkdir(parents=True, exist_ok=True)
        final_bam_path = run_output_dir / (run_params["output_bam_prefix_rel"] + ".bam")
        temp_bam_path = run_output_dir / (run_params["output_bam_prefix_rel"] + "_temp_addsnv.bam")

        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool: DEFAULT_ADDSNV_PATH if tool == "addsnv.py" else DEFAULT_SAMTOOLS_PATH)
        mock_ensure_exists = mocker.patch('bb_runners.ensure_file_exists', side_effect=lambda fp_arg, entity_name_arg="Input file": Path(fp_arg))
        mock_check_fasta_idx = mocker.patch('bb_runners.check_fasta_indexed')
        mock_check_bam_idx = mocker.patch('bb_runners.check_bam_indexed')
        mock_prepare_dir = mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mock_unlink = mocker.patch.object(Path, 'unlink')
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')
        mock_gen_igv = mocker.patch('bb_runners.generate_igv_session_xml')

        def run_cmd_effect(cmd_list, **kwargs):
            if DEFAULT_ADDSNV_PATH in cmd_list: temp_bam_path.touch()
            elif "sort" in cmd_list: final_bam_path.touch()
            return (0, "cmd success", "")
        mock_run_cmd = mocker.patch('bb_runners.run_external_cmd', side_effect=run_cmd_effect)

        runner_args_to_pass = {k: v for k, v in run_params.items() if k in bb_runners.spike_variants_runner.__code__.co_varnames}
        bb_runners.spike_variants_runner(**runner_args_to_pass)

        expected_ensure_calls = [
            mocker.call(run_params["reference_fasta"], "Reference FASTA"),
            mocker.call(run_params["input_bam"], "Input BAM"),
            mocker.call(run_params["vcf_file"], "Input VCF"),
            mocker.call(temp_bam_path, "Temporary output BAM from addsnv.py"),
            mocker.call(final_bam_path, "Final output BAM from variant spiking"),
        ]
        mock_ensure_exists.assert_has_calls(expected_ensure_calls, any_order=True)
        mock_check_fasta_idx.assert_called_once_with(Path(run_params["reference_fasta"]), DEFAULT_SAMTOOLS_PATH)
        mock_check_bam_idx.assert_any_call(Path(run_params["input_bam"]), DEFAULT_SAMTOOLS_PATH, auto_index_if_missing=run_params["auto_index_input_bam"])
        mock_check_bam_idx.assert_any_call(final_bam_path, DEFAULT_SAMTOOLS_PATH, auto_index_if_missing=True)
        assert mock_run_cmd.call_count == 2
        mock_unlink.assert_called_once_with()
        mock_write_manifest.assert_called_once()
        mock_gen_igv.assert_called_once()

    def test_addsnv_tool_failure(self, tmp_path, mocker, default_spike_params):
        run_params = default_spike_params
        run_output_dir = run_params["output_root_dir"] / run_params["run_name"]; run_output_dir.mkdir(parents=True, exist_ok=True)
        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool: DEFAULT_ADDSNV_PATH if tool == "addsnv.py" else DEFAULT_SAMTOOLS_PATH)
        mocker.patch('bb_runners.ensure_file_exists', side_effect=lambda fp_arg, entity_name_arg="Input file": Path(fp_arg))
        mocker.patch('bb_runners.check_fasta_indexed'); mocker.patch('bb_runners.check_bam_indexed')
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        def failing_run_cmd_effect(cmd_list, **kwargs):
            if DEFAULT_ADDSNV_PATH in cmd_list: raise bb_utils.BaseBuddyToolError("addsnv.py failed", command=cmd_list)
            return (0, "", "")
        mock_run_cmd = mocker.patch('bb_runners.run_external_cmd', side_effect=failing_run_cmd_effect)
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')
        runner_args = {k:v for k,v in run_params.items() if k in bb_runners.spike_variants_runner.__code__.co_varnames}
        with pytest.raises(bb_utils.BaseBuddyToolError, match="addsnv.py failed"): bb_runners.spike_variants_runner(**runner_args)
        mock_run_cmd.assert_called_once(); mock_write_manifest.assert_not_called()

    def test_samtools_sort_failure(self, tmp_path, mocker, default_spike_params):
        run_params = default_spike_params
        run_output_dir = run_params["output_root_dir"] / run_params["run_name"]; run_output_dir.mkdir(parents=True, exist_ok=True)
        temp_bam_path = run_output_dir / (run_params["output_bam_prefix_rel"] + "_temp_addsnv.bam")
        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool: DEFAULT_ADDSNV_PATH if tool == "addsnv.py" else DEFAULT_SAMTOOLS_PATH)
        mocker.patch('bb_runners.ensure_file_exists', side_effect=lambda fp_arg, entity_name_arg="Input file": Path(fp_arg))
        mocker.patch('bb_runners.check_fasta_indexed'); mocker.patch('bb_runners.check_bam_indexed')
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mock_unlink = mocker.patch.object(Path, 'unlink')
        def run_cmd_effect(cmd_list, **kwargs):
            if DEFAULT_ADDSNV_PATH in cmd_list: temp_bam_path.touch(); return (0, "addsnv success", "")
            elif "sort" in cmd_list: raise bb_utils.BaseBuddyToolError("samtools sort failed", command=["samtools", "sort"])
            return (0, "", "")
        mock_run_cmd = mocker.patch('bb_runners.run_external_cmd', side_effect=run_cmd_effect)
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')
        runner_args = {k:v for k,v in run_params.items() if k in bb_runners.spike_variants_runner.__code__.co_varnames}
        with pytest.raises(bb_utils.BaseBuddyToolError, match="samtools sort failed"): bb_runners.spike_variants_runner(**runner_args)
        assert mock_run_cmd.call_count == 2;
        mock_unlink.assert_not_called()
        mock_write_manifest.assert_not_called()

    def test_missing_addsnv_tool(self, tmp_path, mocker, default_spike_params):
        run_params = default_spike_params
        def find_tool_effect(tool_name):
            if tool_name == "addsnv.py": raise bb_utils.BaseBuddyConfigError("addsnv.py not found")
            elif tool_name == "samtools": return DEFAULT_SAMTOOLS_PATH
            return None
        mocker.patch('bb_runners.find_tool_path', side_effect=find_tool_effect) # Patched in bb_runners for this test
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_params["output_root_dir"] / run_params["run_name"])
        mocker.patch('bb_runners.ensure_file_exists')
        mocker.patch('bb_runners.check_fasta_indexed')
        mocker.patch('bb_runners.check_bam_indexed')
        mocker.patch('bb_runners.run_external_cmd')

        runner_args = {k:v for k,v in run_params.items() if k in bb_runners.spike_variants_runner.__code__.co_varnames}
        with pytest.raises(bb_utils.BaseBuddyConfigError, match="addsnv.py not found"): bb_runners.spike_variants_runner(**runner_args)

@pytest.fixture
def default_download_params(tmp_path):
    output_root = tmp_path / "output_data_download"
    return {
        "output_root_dir": output_root, "run_name": DEFAULT_RUN_OUTPUT_DIR_NAME + "_download",
        "download_url": "http://example.com/dummy.fasta.gz", "destination_filename": "dummy.fasta.gz",
        "expected_checksum": "abc123checksum", "checksum_algorithm": "sha256",
        "timeout_download": 10800.0, "overwrite_output": False,
        "command_params": {"curl_path": DEFAULT_CURL_PATH, "samtools_path": DEFAULT_SAMTOOLS_PATH,
                           "download_url": "http://example.com/dummy.fasta.gz", "destination_filename": "dummy.fasta.gz",
                           "expected_checksum": "abc123checksum", "checksum_algorithm": "sha256",
                           "output_root_dir": str(output_root), "run_name": DEFAULT_RUN_OUTPUT_DIR_NAME + "_download",
                           "overwrite_output": False,}}

class TestDownloadReferenceRunner: # This tests the root bb_runners.py version
    def test_successful_download_and_index_fasta(self, tmp_path, mocker, default_download_params):
        run_params = copy.deepcopy(default_download_params)
        run_params["destination_filename"] = "dummy.fa"
        run_params["command_params"]["destination_filename"] = "dummy.fa"

        run_output_dir = run_params["output_root_dir"] / run_params["run_name"]
        run_output_dir.mkdir(parents=True, exist_ok=True)
        downloaded_file_path = run_output_dir / run_params["destination_filename"]
        fai_file_path = downloaded_file_path.with_suffix(".fa.fai")

        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool: DEFAULT_CURL_PATH if tool == "curl" else DEFAULT_SAMTOOLS_PATH)
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mock_verify_checksum = mocker.patch('bb_runners.verify_file_checksum')
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')
        mocker.patch('bb_runners.ensure_file_exists', side_effect=lambda fp_arg, entity_name_arg="Input file": Path(fp_arg))
        mock_check_fasta_indexed_runner = mocker.patch('bb_runners.check_fasta_indexed')

        def run_cmd_effect(cmd_list, **kwargs):
            if cmd_list[0] == DEFAULT_CURL_PATH: downloaded_file_path.touch(); return (0, "", "")
            elif cmd_list[0] == DEFAULT_SAMTOOLS_PATH and cmd_list[1] == "faidx": fai_file_path.touch(); return (0, "", "")
            return (1, "unexpected cmd", "error")

        mock_run_cmd_runner = mocker.patch('bb_runners.run_external_cmd', side_effect=run_cmd_effect)

        runner_args = {k:v for k,v in run_params.items() if k in bb_runners.download_reference_runner.__code__.co_varnames}
        bb_runners.download_reference_runner(**runner_args)

        curl_call_expected = mocker.call([DEFAULT_CURL_PATH, "-L", "-o", str(downloaded_file_path), run_params["download_url"]], timeout_seconds=run_params["timeout_download"], stream_output=True, cwd=run_output_dir)
        faidx_call_expected = mocker.call([DEFAULT_SAMTOOLS_PATH, "faidx", str(downloaded_file_path)], cwd=run_output_dir)

        mock_run_cmd_runner.assert_has_calls([curl_call_expected, faidx_call_expected], any_order=False)
        mock_verify_checksum.assert_called_once_with(downloaded_file_path, run_params["expected_checksum"], run_params["checksum_algorithm"])
        mock_check_fasta_indexed_runner.assert_called_once_with(downloaded_file_path, DEFAULT_SAMTOOLS_PATH)
        assert fai_file_path.exists()
        mock_write_manifest.assert_called_once()

    def test_successful_download_non_fasta(self, tmp_path, mocker, default_download_params):
        run_params = copy.deepcopy(default_download_params)
        run_params["destination_filename"] = "dummy_file.txt"
        run_params["command_params"]["destination_filename"] = "dummy_file.txt"

        run_output_dir = run_params["output_root_dir"] / run_params["run_name"]
        run_output_dir.mkdir(parents=True, exist_ok=True)
        downloaded_file_path = run_output_dir / run_params["destination_filename"]

        mocker.patch('bb_runners.find_tool_path', return_value=DEFAULT_CURL_PATH)
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mock_verify_checksum = mocker.patch('bb_runners.verify_file_checksum')
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')
        mocker.patch('bb_runners.ensure_file_exists', side_effect=lambda fp_arg, entity_name_arg="Input file": Path(fp_arg))
        mock_check_fasta_indexed_runner = mocker.patch('bb_runners.check_fasta_indexed')

        def run_cmd_effect(cmd_list, **kwargs):
            if cmd_list[0] == DEFAULT_CURL_PATH: downloaded_file_path.touch(); return (0, "", "")
            return (1, "unexpected cmd", "error")
        mock_run_cmd_runner = mocker.patch('bb_runners.run_external_cmd', side_effect=run_cmd_effect)

        runner_args = {k:v for k,v in run_params.items() if k in bb_runners.download_reference_runner.__code__.co_varnames}
        bb_runners.download_reference_runner(**runner_args)

        mock_run_cmd_runner.assert_called_once_with(
            [DEFAULT_CURL_PATH, "-L", "-o", str(downloaded_file_path), run_params["download_url"]],
            timeout_seconds=run_params["timeout_download"], stream_output=True, cwd=run_output_dir
        )
        mock_verify_checksum.assert_called_once_with(downloaded_file_path, run_params["expected_checksum"], run_params["checksum_algorithm"])
        mock_check_fasta_indexed_runner.assert_not_called()
        mock_write_manifest.assert_called_once()

    def test_download_fails_curl_error(self, tmp_path, mocker, default_download_params):
        run_params = default_download_params
        run_output_dir = run_params["output_root_dir"] / run_params["run_name"]
        if run_output_dir.exists(): shutil.rmtree(run_output_dir) # Clean first
        run_output_dir.mkdir(parents=True, exist_ok=True)
        downloaded_file_path = run_output_dir / run_params["destination_filename"]

        mocker.patch('bb_runners.find_tool_path', return_value=DEFAULT_CURL_PATH)
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)

        curl_command_list = [DEFAULT_CURL_PATH, "-L", "-o", str(downloaded_file_path), run_params["download_url"]]
        mock_run_cmd_runner = mocker.patch('bb_runners.run_external_cmd',
                                           side_effect=bb_utils.BaseBuddyToolError("curl failed", command=curl_command_list))

        mock_path_unlink = mocker.patch.object(Path, 'unlink')

        exists_call_state = {"count": 0}
        def selective_exists_effect(*args):
            nonlocal exists_call_state
            if not args: return True
            path_self_obj = args[0]
            if str(path_self_obj) == str(downloaded_file_path):
                exists_call_state["count"] += 1
                if exists_call_state["count"] == 1: return False # Initial check in runner is False (dir is clean)
                return True  # Second check (in except block) is True (simulating partial file)
            return True # Default for other paths
        mocker.patch.object(Path, 'exists', side_effect=selective_exists_effect)

        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')
        mock_verify_checksum = mocker.patch('bb_runners.verify_file_checksum')
        mock_check_fasta_indexed = mocker.patch('bb_runners.check_fasta_indexed')

        runner_args = {k:v for k,v in run_params.items() if k in bb_runners.download_reference_runner.__code__.co_varnames}
        with pytest.raises(bb_utils.BaseBuddyToolError, match="curl failed"):
            bb_runners.download_reference_runner(**runner_args)

        mock_run_cmd_runner.assert_called_once()
        mock_path_unlink.assert_called_with(downloaded_file_path)
        mock_write_manifest.assert_not_called()
        mock_verify_checksum.assert_not_called()
        mock_check_fasta_indexed.assert_not_called()

    def test_checksum_fails(self, tmp_path, mocker, default_download_params):
        run_params = default_download_params
        run_output_dir = run_params["output_root_dir"] / run_params["run_name"]; run_output_dir.mkdir(parents=True, exist_ok=True)
        downloaded_file_path = run_output_dir / run_params["destination_filename"]

        mocker.patch('bb_runners.find_tool_path', return_value=DEFAULT_CURL_PATH)
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mocker.patch('bb_runners.ensure_file_exists', return_value=downloaded_file_path)
        mocker.patch('bb_runners.run_external_cmd', side_effect=lambda *args, **kwargs: downloaded_file_path.touch() or (0,"",""))

        mocker.patch('bb_runners.verify_file_checksum', side_effect=bb_utils.BaseBuddyChecksumError("checksum mismatch"))
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')

        runner_args = {k:v for k,v in run_params.items() if k in bb_runners.download_reference_runner.__code__.co_varnames}
        with pytest.raises(bb_utils.BaseBuddyChecksumError, match="checksum mismatch"):
            bb_runners.download_reference_runner(**runner_args)
        mock_write_manifest.assert_not_called()

    def test_fasta_indexing_fails_if_fasta(self, tmp_path, mocker, default_download_params):
        run_params = copy.deepcopy(default_download_params)
        run_params["destination_filename"] = "dummy.fa"
        run_params["command_params"]["destination_filename"] = "dummy.fa"
        run_output_dir = run_params["output_root_dir"] / run_params["run_name"]; run_output_dir.mkdir(parents=True, exist_ok=True)
        downloaded_file_path = run_output_dir / run_params["destination_filename"]

        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool: DEFAULT_CURL_PATH if tool == "curl" else DEFAULT_SAMTOOLS_PATH)
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mocker.patch('bb_runners.ensure_file_exists', return_value=downloaded_file_path)
        mocker.patch('bb_runners.verify_file_checksum')
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')

        def run_cmd_effect(cmd_list, **kwargs):
            if cmd_list[0] == DEFAULT_CURL_PATH: downloaded_file_path.touch(); return (0,"","")
            elif cmd_list[0] == DEFAULT_SAMTOOLS_PATH and cmd_list[1] == "faidx":
                raise bb_utils.BaseBuddyToolError("samtools faidx failed", command=cmd_list)
            return (1, "unexpected", "error")
        mocker.patch('bb_runners.run_external_cmd', side_effect=run_cmd_effect)
        mocker.patch('bb_runners.check_fasta_indexed')

        runner_args = {k:v for k,v in run_params.items() if k in bb_runners.download_reference_runner.__code__.co_varnames}
        with pytest.raises(NameError, match="name 'BaseBuddyError' is not defined"):
            bb_runners.download_reference_runner(**runner_args)

        mock_write_manifest.assert_called_once()
        assert mock_write_manifest.call_args[1]['status'] == "failed_indexing"

import xml.etree.ElementTree as ET

import pytest
from pathlib import Path
from unittest import mock
import copy
import shutil

from src.basebuddy import runner as src_bb_runner # For testing refactored runner
from basebuddy import utils as bb_utils

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

        mock_find_tool = mocker.patch('basebuddy.utils.find_tool_path', side_effect=lambda tool_name: run_params["command_params"]["art_exe_path_param"] if tool_name.startswith('art_') else run_params["command_params"]["samtools_exe_path_param"])
        mock_ensure_file_exists_util = mocker.patch('basebuddy.utils.ensure_file_exists', return_value=Path(run_params["reference_fasta"]))
        mock_check_fasta_indexed_util = mocker.patch('basebuddy.utils.check_fasta_indexed')
        mock_prepare_dir_util = mocker.patch('basebuddy.utils.prepare_run_output_dir', return_value=expected_run_output_dir)

        def mock_art_cmd_effect(cmd_list, **kwargs):
            output_prefix_path_str = next(cmd_list[i+1] for i, arg in enumerate(cmd_list) if arg == "-o")
            prefix_path = Path(output_prefix_path_str)
            assert prefix_path == expected_run_output_dir / id_prefix
            (prefix_path.parent / (prefix_path.name + "1.fq")).touch()
            if run_params["is_paired_end"]: (prefix_path.parent / (prefix_path.name + "2.fq")).touch()
            if not run_params["command_params"]["no_aln_output"]: (prefix_path.parent / (prefix_path.name + ".aln")).touch()
            return (0, "ART success", "")
        mock_run_external_cmd_util = mocker.patch('basebuddy.utils.run_external_cmd', side_effect=mock_art_cmd_effect)
        mock_write_manifest_util = mocker.patch('basebuddy.utils.write_run_manifest')

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

        mocker.patch('basebuddy.utils.find_tool_path', side_effect=lambda tool_name: run_params["command_params"]["art_exe_path_param"] if tool_name.startswith('art_') else run_params["command_params"]["samtools_exe_path_param"])
        mock_ensure_file_exists_util = mocker.patch('basebuddy.utils.ensure_file_exists', return_value=Path(run_params["reference_fasta"]))
        mocker.patch('basebuddy.utils.check_fasta_indexed')
        mocker.patch('basebuddy.utils.prepare_run_output_dir', return_value=expected_run_output_dir)

        def mock_art_cmd_effect_single(cmd_list, **kwargs):
            output_prefix_path_str = next(cmd_list[i+1] for i, arg in enumerate(cmd_list) if arg == "-o")
            prefix_path = Path(output_prefix_path_str)
            (prefix_path.parent / (prefix_path.name + ".fq")).touch() # Single-end
            if not run_params["command_params"]["no_aln_output"]: (prefix_path.parent / (prefix_path.name + ".aln")).touch()
            return (0, "ART SE success", "")
        mock_run_external_cmd_util = mocker.patch('basebuddy.utils.run_external_cmd', side_effect=mock_art_cmd_effect_single)
        mock_write_manifest_util = mocker.patch('basebuddy.utils.write_run_manifest')

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

        mocker.patch('basebuddy.utils.find_tool_path', side_effect=lambda tool_name: run_params["command_params"]["art_exe_path_param"] if tool_name.startswith('art_') else run_params["command_params"]["samtools_exe_path_param"])
        mocker.patch('basebuddy.utils.ensure_file_exists', return_value=Path(run_params["reference_fasta"]))
        mocker.patch('basebuddy.utils.check_fasta_indexed')
        mocker.patch('basebuddy.utils.prepare_run_output_dir', return_value=expected_run_output_dir)
        mock_run_external_cmd_util = mocker.patch('basebuddy.utils.run_external_cmd', side_effect=bb_utils.BaseBuddyToolError("ART sim failed", command=["art_exe"]))
        mock_write_manifest_util = mocker.patch('basebuddy.utils.write_run_manifest')

        with pytest.raises(bb_utils.BaseBuddyToolError, match="ART sim failed"):
            src_bb_runner.simulate_short(**run_params)
        mock_run_external_cmd_util.assert_called_once()
        mock_write_manifest_util.assert_not_called()

    def test_missing_art_tool(self, tmp_path, mocker, default_src_simulate_params):
        run_params = default_src_simulate_params
        def find_tool_effect(tool_name):
            if tool_name.startswith('art_'): raise bb_utils.BaseBuddyConfigError("ART not found")
            return run_params["command_params"]["samtools_exe_path_param"]
        mocker.patch('basebuddy.utils.find_tool_path', side_effect=find_tool_effect) # Patched in basebuddy.utils
        mocker.patch('basebuddy.utils.prepare_run_output_dir', return_value=run_params["output_root_dir"] / run_params["run_name"])
        mocker.patch('basebuddy.utils.ensure_file_exists')
        mocker.patch('basebuddy.utils.check_fasta_indexed')
        mocker.patch('basebuddy.utils.run_external_cmd')

        with pytest.raises(bb_utils.BaseBuddyConfigError, match="ART not found"):
            src_bb_runner.simulate_short(**run_params)

    def test_invalid_input_depth(self, tmp_path, default_src_simulate_params):
        run_params = copy.deepcopy(default_src_simulate_params); run_params["depth"] = 0
        if "command_params" in run_params: run_params["command_params"]["depth"] = 0

        with pytest.raises(bb_utils.BaseBuddyInputError, match="Sequencing depth must be a positive integer"):
            src_bb_runner.simulate_short(**run_params)

import xml.etree.ElementTree as ET

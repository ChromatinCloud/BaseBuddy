import pytest
from pathlib import Path
from unittest import mock
import copy
import shutil

# Merged Imports from jules_wip and main (prioritizing jules_wip)
from src.basebuddy import runner as src_bb_runner # For testing refactored runner
import bb_runners # For testing original runners (kept from both)
import bb_utils # Kept from both
import xml.etree.ElementTree as ET

# Merged Default parameters for tests (prioritizing jules_wip where names conflicted)
DEFAULT_ART_PATH = "dummy_art_illumina"
DEFAULT_SAMTOOLS_PATH = "dummy_samtools"
DEFAULT_ADDSNV_PATH = "dummy_addsnv.py"
DEFAULT_CURL_PATH = "dummy_curl"
DEFAULT_REFERENCE_FASTA_FILENAME = "reference.fasta" # Using consistent name, content from both was "reference.fasta"
DEFAULT_RUN_OUTPUT_DIR_NAME = "test_run" # From jules_wip

@pytest.fixture
def default_src_simulate_params(tmp_path): # As per jules_wip intent
    output_root = tmp_path / "output_data_simulate_src"
    # Use DEFAULT_REFERENCE_FASTA_FILENAME for consistency
    ref_fasta_path_str = str(tmp_path / DEFAULT_REFERENCE_FASTA_FILENAME)
    Path(ref_fasta_path_str).parent.mkdir(parents=True, exist_ok=True)
    # Create a dummy reference FASTA file with content
    with open(ref_fasta_path_str, "w") as f:
        f.write(">chr1\nACGTACGTACGTN\n")
        f.write(">chr2\nNNNNNNNNNNNNNN\n")

    current_run_name = DEFAULT_RUN_OUTPUT_DIR_NAME + "_src_simulate"

    return {
        "output_root_dir": str(output_root),
        "reference_fasta": ref_fasta_path_str,
        "depth": 50,
        "read_length": 150,
        "art_profile": "HS25",
        "run_name": current_run_name,
        "mean_fragment_length": 400,
        "std_dev_fragment_length": 50,
        "is_paired_end": True,
        "overwrite_output": False,
        "art_platform": "illumina",
        "timeout": 3600.0,
        "auto_index_fasta": True,
        "command_params": { # Reflecting structure with duplicated params
            "art_exe_path_param": DEFAULT_ART_PATH,
            "samtools_exe_path_param": DEFAULT_SAMTOOLS_PATH,
            "id_prefix": "sim_reads_src",
            "no_aln_output": False,
            "reference_fasta": ref_fasta_path_str, # Duplicated
            "depth": 50, # Duplicated
            "read_length": 150, # Duplicated
            "art_profile": "HS25", # Duplicated
            "mean_fragment_length": 400, # Duplicated
            "std_dev_fragment_length": 50, # Duplicated
            "is_paired_end": True, # Duplicated
            "art_platform": "illumina", # Duplicated
            "overwrite_output": False, # Duplicated
            "auto_index_fasta": True, # Duplicated
            "output_root_dir": str(output_root), # Duplicated
            "run_name": current_run_name, # Duplicated
        }
    }

class TestSrcSimulateShortRunner: # From jules_wip
    def test_successful_paired_end_simulation(self, tmp_path, mocker, default_src_simulate_params):
        run_params = default_src_simulate_params
        # Convert output_root_dir to Path for internal test logic, though fixture provides string
        run_params_output_root_dir_path = Path(run_params["output_root_dir"])
        expected_run_output_dir = run_params_output_root_dir_path / run_params["run_name"]
        
        id_prefix = run_params["command_params"]["id_prefix"]

        mock_find_tool = mocker.patch('bb_utils.find_tool_path', side_effect=lambda tool_name: run_params["command_params"]["art_exe_path_param"] if tool_name.startswith('art_') else run_params["command_params"]["samtools_exe_path_param"])
        mock_ensure_file_exists_util = mocker.patch('bb_utils.ensure_file_exists', return_value=Path(run_params["reference_fasta"]))
        mock_check_fasta_indexed_util = mocker.patch('bb_utils.check_fasta_indexed')
        # mock_prepare_dir_util will be called by the runner to get the output directory path
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

        src_bb_runner.simulate_short(**run_params)

        mock_prepare_dir_util.assert_called_once_with(run_params_output_root_dir_path, run_params["run_name"], run_params["overwrite_output"])
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
        manifest_call_args = mock_write_manifest_util.call_args # Corrected variable name from manifest_kall
        assert manifest_call_args.kwargs['run_name'] == run_params["run_name"]
        assert manifest_call_args.kwargs['command_name'] == f"simulate_short_reads_{run_params['art_platform']}"
        assert manifest_call_args.kwargs['reference_genome_path'] == run_params["reference_fasta"]
        assert manifest_call_args.kwargs['parameters']["id_prefix"] == id_prefix

    def test_successful_single_end_simulation(self, tmp_path, mocker, default_src_simulate_params):
        run_params = copy.deepcopy(default_src_simulate_params)
        run_params["is_paired_end"] = False
        if "command_params" in run_params: # Keep command_params consistent if used by runner
             run_params["command_params"]["is_paired_end"] = False

        run_params_output_root_dir_path = Path(run_params["output_root_dir"])
        expected_run_output_dir = run_params_output_root_dir_path / run_params["run_name"]
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
        run_params_output_root_dir_path = Path(run_params["output_root_dir"])
        expected_run_output_dir = run_params_output_root_dir_path / run_params["run_name"]

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
        mocker.patch('bb_utils.find_tool_path', side_effect=find_tool_effect)
        mocker.patch('bb_utils.prepare_run_output_dir', return_value=Path(run_params["output_root_dir"]) / run_params["run_name"])
        mocker.patch('bb_utils.ensure_file_exists')
        mocker.patch('bb_utils.check_fasta_indexed')
        # run_external_cmd mock is not needed if find_tool_path fails before it

        with pytest.raises(bb_utils.BaseBuddyConfigError, match="ART not found"):
            src_bb_runner.simulate_short(**run_params)

    def test_invalid_input_depth(self, tmp_path, default_src_simulate_params):
        run_params = copy.deepcopy(default_src_simulate_params)
        run_params["depth"] = 0
        if "command_params" in run_params and "depth" in run_params["command_params"]: # if depth is also in command_params
             run_params["command_params"]["depth"] = 0

        with pytest.raises(bb_utils.BaseBuddyInputError, match="Sequencing depth must be a positive integer"):
            src_bb_runner.simulate_short(**run_params)

@pytest.fixture
def default_spike_params(tmp_path): # This fixture remains for the old TestSpikeVariantsRunner (from jules_wip)
    output_root = tmp_path / "output_data_spike"
    ref_fasta_path = tmp_path / "ref_spike.fasta"; ref_fasta_path.touch()
    input_bam_path = tmp_path / "input.bam"; input_bam_path.touch(); (tmp_path / "input.bam.bai").touch()
    vcf_file_path = tmp_path / "variants.vcf"; vcf_file_path.touch()
    run_name = DEFAULT_RUN_OUTPUT_DIR_NAME + "_spike"
    return {
        "output_root_dir": str(output_root),
        "run_name": run_name,
        "reference_fasta": str(ref_fasta_path),
        "input_bam": str(input_bam_path),
        "vcf_file": str(vcf_file_path),
        "output_bam_prefix_rel": "spiked_output",
        "overwrite_output": False,
        "auto_index_input_bam": False,
        "timeout": 7200.0, # Top-level timeout
        "command_params": { # Specific paths for tools
            "addsnv_path": DEFAULT_ADDSNV_PATH,
            "samtools_path": DEFAULT_SAMTOOLS_PATH,
        }
    }

class TestSpikeVariantsRunner: # This tests the root bb_runners.py version (from jules_wip)
    def test_successful_spike(self, tmp_path, mocker, default_spike_params):
        run_params = default_spike_params
        run_params_output_root_dir_path = Path(run_params["output_root_dir"]) # Convert to Path for consistency
        run_output_dir = run_params_output_root_dir_path / run_params["run_name"]
        final_bam_path = run_output_dir / (run_params["output_bam_prefix_rel"] + ".bam")
        temp_bam_path = run_output_dir / (run_params["output_bam_prefix_rel"] + "_temp_addsnv.bam")

        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool: default_spike_params["command_params"]["addsnv_path"] if tool == "addsnv.py" else default_spike_params["command_params"]["samtools_path"])
        mock_ensure_exists = mocker.patch('bb_runners.ensure_file_exists', side_effect=lambda fp_arg, entity_name_arg="Input file": Path(fp_arg))
        mock_check_fasta_idx = mocker.patch('bb_runners.check_fasta_indexed')
        mock_check_bam_idx = mocker.patch('bb_runners.check_bam_indexed')
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mock_unlink = mocker.patch.object(Path, 'unlink', return_value=None)
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')
        mock_gen_igv = mocker.patch('bb_runners.generate_igv_session_xml')

        def run_cmd_effect(cmd_list, **kwargs):
            if default_spike_params["command_params"]["addsnv_path"] in cmd_list: temp_bam_path.touch()
            elif "sort" in cmd_list: final_bam_path.touch()
            return (0, "cmd success", "")
        mock_run_cmd = mocker.patch('bb_runners.run_external_cmd', side_effect=run_cmd_effect)

        runner_arg_names = bb_runners.spike_variants_runner.__code__.co_varnames
        runner_args_to_pass = {k: v for k, v in run_params.items() if k in runner_arg_names}
        if "command_params" in runner_arg_names and "command_params" not in runner_args_to_pass :
             runner_args_to_pass["command_params"] = run_params.get("command_params", {})

        bb_runners.spike_variants_runner(**runner_args_to_pass)

        expected_ensure_calls = [
            mocker.call(run_params["reference_fasta"], "Reference FASTA"),
            mocker.call(run_params["input_bam"], "Input BAM"),
            mocker.call(run_params["vcf_file"], "Input VCF"),
            mocker.call(str(temp_bam_path), "Temporary output BAM from addsnv.py"),
            mocker.call(str(final_bam_path), "Final output BAM from variant spiking"),
        ]
        mock_ensure_exists.assert_has_calls(expected_ensure_calls, any_order=True) # any_order for flexibility
        mock_check_fasta_idx.assert_called_once_with(Path(run_params["reference_fasta"]), default_spike_params["command_params"]["samtools_path"])
        mock_check_bam_idx.assert_any_call(Path(run_params["input_bam"]), default_spike_params["command_params"]["samtools_path"], auto_index_if_missing=run_params["auto_index_input_bam"])
        mock_check_bam_idx.assert_any_call(final_bam_path, default_spike_params["command_params"]["samtools_path"], auto_index_if_missing=True)
        assert mock_run_cmd.call_count == 2
        mock_unlink.assert_called_once_with()
        mock_write_manifest.assert_called_once()
        mock_gen_igv.assert_called_once()

    def test_addsnv_tool_failure(self, tmp_path, mocker, default_spike_params):
        run_params = default_spike_params
        run_params_output_root_dir_path = Path(run_params["output_root_dir"])
        run_output_dir = run_params_output_root_dir_path / run_params["run_name"]
        
        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool: default_spike_params["command_params"]["addsnv_path"] if tool == "addsnv.py" else default_spike_params["command_params"]["samtools_path"])
        mocker.patch('bb_runners.ensure_file_exists', side_effect=lambda fp_arg, entity_name_arg="Input file": Path(fp_arg))
        mocker.patch('bb_runners.check_fasta_indexed'); mocker.patch('bb_runners.check_bam_indexed')
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        def failing_run_cmd_effect(cmd_list, **kwargs):
            if default_spike_params["command_params"]["addsnv_path"] in cmd_list: raise bb_utils.BaseBuddyToolError("addsnv.py failed", command=cmd_list)
            return (0, "", "") # Should not be reached for addsnv call
        mock_run_cmd = mocker.patch('bb_runners.run_external_cmd', side_effect=failing_run_cmd_effect)
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')
        
        runner_arg_names = bb_runners.spike_variants_runner.__code__.co_varnames
        runner_args_to_pass = {k: v for k, v in run_params.items() if k in runner_arg_names}
        if "command_params" in runner_arg_names and "command_params" not in runner_args_to_pass :
             runner_args_to_pass["command_params"] = run_params.get("command_params", {})

        with pytest.raises(bb_utils.BaseBuddyToolError, match="addsnv.py failed"):
            bb_runners.spike_variants_runner(**runner_args_to_pass)
        mock_run_cmd.assert_called_once()
        mock_write_manifest.assert_not_called()

    def test_samtools_sort_failure(self, tmp_path, mocker, default_spike_params):
        run_params = default_spike_params
        run_params_output_root_dir_path = Path(run_params["output_root_dir"])
        run_output_dir = run_params_output_root_dir_path / run_params["run_name"]
        temp_bam_path = run_output_dir / (run_params["output_bam_prefix_rel"] + "_temp_addsnv.bam")

        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool: default_spike_params["command_params"]["addsnv_path"] if tool == "addsnv.py" else default_spike_params["command_params"]["samtools_path"])
        mocker.patch('bb_runners.ensure_file_exists', side_effect=lambda fp_arg, entity_name_arg="Input file": Path(fp_arg))
        mocker.patch('bb_runners.check_fasta_indexed'); mocker.patch('bb_runners.check_bam_indexed')
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mock_unlink = mocker.patch.object(Path, 'unlink', return_value=None)
        def run_cmd_effect(cmd_list, **kwargs):
            if default_spike_params["command_params"]["addsnv_path"] in cmd_list: temp_bam_path.touch(); return (0, "addsnv success", "")
            elif "sort" in cmd_list: raise bb_utils.BaseBuddyToolError("samtools sort failed", command=["samtools", "sort"])
            return (1, "unexpected cmd", "error") # Fallback
        mock_run_cmd = mocker.patch('bb_runners.run_external_cmd', side_effect=run_cmd_effect)
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')

        runner_arg_names = bb_runners.spike_variants_runner.__code__.co_varnames
        runner_args_to_pass = {k: v for k, v in run_params.items() if k in runner_arg_names}
        if "command_params" in runner_arg_names and "command_params" not in runner_args_to_pass :
             runner_args_to_pass["command_params"] = run_params.get("command_params", {})

        with pytest.raises(bb_utils.BaseBuddyToolError, match="samtools sort failed"):
            bb_runners.spike_variants_runner(**runner_args_to_pass)
        assert mock_run_cmd.call_count == 2
        mock_unlink.assert_not_called() # Temp file should not be unlinked if sort fails
        mock_write_manifest.assert_not_called()

    def test_missing_addsnv_tool(self, tmp_path, mocker, default_spike_params):
        run_params = default_spike_params
        def find_tool_effect(tool_name):
            if tool_name == "addsnv.py": raise bb_utils.BaseBuddyConfigError("addsnv.py not found")
            elif tool_name == "samtools": return default_spike_params["command_params"]["samtools_path"]
            return None 
        mocker.patch('bb_runners.find_tool_path', side_effect=find_tool_effect)
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=Path(run_params["output_root_dir"]) / run_params["run_name"])
        mocker.patch('bb_runners.ensure_file_exists')
        mocker.patch('bb_runners.check_fasta_indexed')
        mocker.patch('bb_runners.check_bam_indexed')
        # run_external_cmd mock not needed if find_tool_path fails early

        runner_arg_names = bb_runners.spike_variants_runner.__code__.co_varnames
        runner_args_to_pass = {k: v for k, v in run_params.items() if k in runner_arg_names}
        if "command_params" in runner_arg_names and "command_params" not in runner_args_to_pass :
             runner_args_to_pass["command_params"] = run_params.get("command_params", {})
             
        with pytest.raises(bb_utils.BaseBuddyConfigError, match="addsnv.py not found"):
            bb_runners.spike_variants_runner(**runner_args_to_pass)

@pytest.fixture
def default_download_params(tmp_path): # From jules_wip
    output_root = tmp_path / "output_data_download"
    run_name = DEFAULT_RUN_OUTPUT_DIR_NAME + "_download"
    return {
        "output_root_dir": str(output_root),
        "run_name": run_name,
        "download_url": "http://example.com/dummy.fasta.gz",
        "destination_filename": "dummy.fasta.gz",
        "expected_checksum": "abc123checksum",
        "checksum_algorithm": "sha256",
        "timeout_download": 10800.0,
        "overwrite_output": False,
        "command_params": {
            "curl_path": DEFAULT_CURL_PATH,
            "samtools_path": DEFAULT_SAMTOOLS_PATH,
        }
    }

class TestDownloadReferenceRunner: # This tests the root bb_runners.py version (from jules_wip)
    def test_successful_download_and_index_fasta(self, tmp_path, mocker, default_download_params):
        run_params = copy.deepcopy(default_download_params)
        run_params["destination_filename"] = "dummy.fa" # Make it a .fa for indexing

        run_params_output_root_dir_path = Path(run_params["output_root_dir"])
        run_output_dir = run_params_output_root_dir_path / run_params["run_name"]
        downloaded_file_path = run_output_dir / run_params["destination_filename"]
        fai_file_path = downloaded_file_path.with_suffix(downloaded_file_path.suffix + ".fai")

        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool: run_params["command_params"]["curl_path"] if tool == "curl" else run_params["command_params"]["samtools_path"])
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mock_verify_checksum = mocker.patch('bb_runners.verify_file_checksum')
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')
        mocker.patch('bb_runners.ensure_file_exists', side_effect=lambda fp_arg, entity_name_arg="Input file": Path(fp_arg))
        mock_check_fasta_indexed_in_bb_runners = mocker.patch('bb_runners.check_fasta_indexed')

        def run_cmd_effect(cmd_list, **kwargs):
            if cmd_list[0] == run_params["command_params"]["curl_path"]: downloaded_file_path.touch(); return (0, "", "")
            elif cmd_list[0] == run_params["command_params"]["samtools_path"] and cmd_list[1] == "faidx": fai_file_path.touch(); return (0, "", "")
            return (1, "unexpected cmd", "error")
        mock_run_cmd_runner = mocker.patch('bb_runners.run_external_cmd', side_effect=run_cmd_effect)

        runner_arg_names = bb_runners.download_reference_runner.__code__.co_varnames
        runner_args_to_pass = {k: v for k, v in run_params.items() if k in runner_arg_names}
        if "command_params" in runner_arg_names and "command_params" not in runner_args_to_pass :
             runner_args_to_pass["command_params"] = run_params.get("command_params", {})

        bb_runners.download_reference_runner(**runner_args_to_pass)

        curl_call_expected = mocker.call([run_params["command_params"]["curl_path"], "-L", "-o", str(downloaded_file_path), run_params["download_url"]], timeout_seconds=run_params["timeout_download"], stream_output=True, cwd=run_output_dir)
        faidx_call_expected = mocker.call([run_params["command_params"]["samtools_path"], "faidx", str(downloaded_file_path)], cwd=run_output_dir)

        mock_run_cmd_runner.assert_has_calls([curl_call_expected, faidx_call_expected], any_order=False)
        mock_verify_checksum.assert_called_once_with(downloaded_file_path, run_params["expected_checksum"], run_params["checksum_algorithm"])
        mock_check_fasta_indexed_in_bb_runners.assert_called_once_with(downloaded_file_path, run_params["command_params"]["samtools_path"])
        assert fai_file_path.exists()
        mock_write_manifest.assert_called_once()

    def test_successful_download_non_fasta(self, tmp_path, mocker, default_download_params):
        run_params = copy.deepcopy(default_download_params)
        run_params["destination_filename"] = "dummy_file.txt"

        run_params_output_root_dir_path = Path(run_params["output_root_dir"])
        run_output_dir = run_params_output_root_dir_path / run_params["run_name"]
        downloaded_file_path = run_output_dir / run_params["destination_filename"]

        mocker.patch('bb_runners.find_tool_path', return_value=run_params["command_params"]["curl_path"])
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mock_verify_checksum = mocker.patch('bb_runners.verify_file_checksum')
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')
        mocker.patch('bb_runners.ensure_file_exists', side_effect=lambda fp_arg, entity_name_arg="Input file": Path(fp_arg))
        mock_check_fasta_indexed_in_bb_runners = mocker.patch('bb_runners.check_fasta_indexed')

        def run_cmd_effect(cmd_list, **kwargs):
            if cmd_list[0] == run_params["command_params"]["curl_path"]: downloaded_file_path.touch(); return (0, "", "")
            return (1, "unexpected cmd", "error")
        mock_run_cmd_runner = mocker.patch('bb_runners.run_external_cmd', side_effect=run_cmd_effect)

        runner_arg_names = bb_runners.download_reference_runner.__code__.co_varnames
        runner_args_to_pass = {k: v for k, v in run_params.items() if k in runner_arg_names}
        if "command_params" in runner_arg_names and "command_params" not in runner_args_to_pass :
             runner_args_to_pass["command_params"] = run_params.get("command_params", {})

        bb_runners.download_reference_runner(**runner_args_to_pass)

        mock_run_cmd_runner.assert_called_once_with(
            [run_params["command_params"]["curl_path"], "-L", "-o", str(downloaded_file_path), run_params["download_url"]],
            timeout_seconds=run_params["timeout_download"], stream_output=True, cwd=run_output_dir
        )
        mock_verify_checksum.assert_called_once_with(downloaded_file_path, run_params["expected_checksum"], run_params["checksum_algorithm"])
        mock_check_fasta_indexed_in_bb_runners.assert_not_called()
        mock_write_manifest.assert_called_once()

    def test_download_fails_curl_error(self, tmp_path, mocker, default_download_params):
        run_params = default_download_params
        run_params_output_root_dir_path = Path(run_params["output_root_dir"])
        run_output_dir = run_params_output_root_dir_path / run_params["run_name"]
        downloaded_file_path = run_output_dir / run_params["destination_filename"]

        mocker.patch('bb_runners.find_tool_path', return_value=run_params["command_params"]["curl_path"])
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)

        curl_command_list = [run_params["command_params"]["curl_path"], "-L", "-o", str(downloaded_file_path), run_params["download_url"]]
        mock_run_cmd_runner = mocker.patch('bb_runners.run_external_cmd',
                                           side_effect=bb_utils.BaseBuddyToolError("curl failed", command=curl_command_list))

        mock_path_unlink = mocker.patch.object(Path, 'unlink', return_value=None)
        # Simulate file was partially created before error, so unlink is attempted
        downloaded_file_path.touch() 

        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')
        mock_verify_checksum = mocker.patch('bb_runners.verify_file_checksum')
        mock_check_fasta_indexed = mocker.patch('bb_runners.check_fasta_indexed')

        runner_arg_names = bb_runners.download_reference_runner.__code__.co_varnames
        runner_args_to_pass = {k: v for k, v in run_params.items() if k in runner_arg_names}
        if "command_params" in runner_arg_names and "command_params" not in runner_args_to_pass :
             runner_args_to_pass["command_params"] = run_params.get("command_params", {})

        with pytest.raises(bb_utils.BaseBuddyToolError, match="curl failed"):
            bb_runners.download_reference_runner(**runner_args_to_pass)

        mock_run_cmd_runner.assert_called_once()
        if downloaded_file_path.exists(): # Runner should attempt to unlink if file exists after error
            mock_path_unlink.assert_called_with(downloaded_file_path)
        else: # If runner successfully unlinks it, this also fine
            mock_path_unlink.assert_called_once()


        mock_write_manifest.assert_not_called()
        mock_verify_checksum.assert_not_called()
        mock_check_fasta_indexed.assert_not_called()

    def test_checksum_fails(self, tmp_path, mocker, default_download_params):
        run_params = default_download_params
        run_params_output_root_dir_path = Path(run_params["output_root_dir"])
        run_output_dir = run_params_output_root_dir_path / run_params["run_name"]
        downloaded_file_path = run_output_dir / run_params["destination_filename"]

        mocker.patch('bb_runners.find_tool_path', return_value=run_params["command_params"]["curl_path"])
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mocker.patch('bb_runners.ensure_file_exists', return_value=downloaded_file_path)
        mocker.patch('bb_runners.run_external_cmd', side_effect=lambda *args, **kwargs: downloaded_file_path.touch() or (0,"",""))

        mocker.patch('bb_runners.verify_file_checksum', side_effect=bb_utils.BaseBuddyChecksumError("checksum mismatch"))
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')

        runner_arg_names = bb_runners.download_reference_runner.__code__.co_varnames
        runner_args_to_pass = {k: v for k, v in run_params.items() if k in runner_arg_names}
        if "command_params" in runner_arg_names and "command_params" not in runner_args_to_pass :
             runner_args_to_pass["command_params"] = run_params.get("command_params", {})

        with pytest.raises(bb_utils.BaseBuddyChecksumError, match="checksum mismatch"):
            bb_runners.download_reference_runner(**runner_args_to_pass)
        mock_write_manifest.assert_not_called()

    def test_fasta_indexing_fails_if_fasta(self, tmp_path, mocker, default_download_params):
        run_params = copy.deepcopy(default_download_params)
        run_params["destination_filename"] = "dummy.fa"

        run_params_output_root_dir_path = Path(run_params["output_root_dir"])
        run_output_dir = run_params_output_root_dir_path / run_params["run_name"]
        downloaded_file_path = run_output_dir / run_params["destination_filename"]

        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool: run_params["command_params"]["curl_path"] if tool == "curl" else run_params["command_params"]["samtools_path"])
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=run_output_dir)
        mocker.patch('bb_runners.ensure_file_exists', return_value=downloaded_file_path)
        mocker.patch('bb_runners.verify_file_checksum') # Assume checksum passes
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')

        # Mock bb_runners.check_fasta_indexed to raise an error that the runner should handle
        # by writing manifest with status="failed_indexing"
        mock_check_fasta_idx = mocker.patch('bb_runners.check_fasta_indexed', 
                                            side_effect=bb_utils.BaseBuddyToolError("Simulated faidx error for testing manifest status"))

        # Simulate curl command succeeding
        mocker.patch('bb_runners.run_external_cmd', side_effect=lambda cmd_list, **kwargs: downloaded_file_path.touch() or (0,"","") if cmd_list[0] == run_params["command_params"]["curl_path"] else (1, "unexpected", "error"))


        runner_arg_names = bb_runners.download_reference_runner.__code__.co_varnames
        runner_args_to_pass = {k: v for k, v in run_params.items() if k in runner_arg_names}
        if "command_params" in runner_arg_names and "command_params" not in runner_args_to_pass :
             runner_args_to_pass["command_params"] = run_params.get("command_params", {})
        
        # Assuming download_reference_runner catches BaseBuddyToolError from check_fasta_indexed
        # and proceeds to write the manifest with a 'failed_indexing' status.
        bb_runners.download_reference_runner(**runner_args_to_pass)
        
        mock_write_manifest.assert_called_once()
        # Ensure the status keyword argument was passed correctly to write_run_manifest
        assert 'status' in mock_write_manifest.call_args.kwargs
        assert mock_write_manifest.call_args.kwargs['status'] == "failed_indexing"


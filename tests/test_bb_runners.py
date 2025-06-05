import pytest
from pathlib import Path
from unittest import mock

# Assuming bb_runners.py and bb_utils.py are accessible in PYTHONPATH
# PYTHONPATH will be set to /app in the execution environment.
import bb_runners
import bb_utils # Used for exception types and potentially other utils if not mocked

# Default parameters for tests
DEFAULT_ART_PATH = "dummy_art_illumina"
DEFAULT_SAMTOOLS_PATH = "dummy_samtools"
DEFAULT_REFERENCE_FASTA = "reference.fasta"
DEFAULT_RUN_OUTPUT_DIR_NAME = "test_run_paired"

@pytest.fixture
def default_runner_params(tmp_path):
    """Parameters that are direct arguments to simulate_short_reads_runner."""
    output_root = tmp_path / "output_data"
    ref_fasta = tmp_path / DEFAULT_REFERENCE_FASTA
    return {
        "output_root_dir": output_root, # Path object
        "run_name": DEFAULT_RUN_OUTPUT_DIR_NAME,
        "reference_fasta": str(ref_fasta), # String path
        "depth": 50,
        "read_length": 150,
        "art_profile": "HS25", # Example, ensure it's a valid ART seq system
        "mean_fragment_length": 400,
        "std_dev_fragment_length": 50,
        "is_paired_end": True,
        "overwrite_output": False,
        "art_platform": "illumina",
        "timeout": 3600.0,
        # command_params will hold other settings for manifest and ART construction
        "command_params": {
            # These might be used by the runner to construct ART cmd or for manifest
            "art_illumina_path": DEFAULT_ART_PATH, # For find_tool_path mock
            "samtools_path": DEFAULT_SAMTOOLS_PATH, # For find_tool_path mock
            # ART specific params not directly in runner signature but needed for cmd construction
            # (The runner currently derives them or uses defaults, but good to have for flexibility)
            "quality_shift": None,
            "quality_shift2": None,
            "random_seed": 12345,
            "profile_type": "empirical",
            "custom_profile_path1": None,
            "custom_profile_path2": None,
            "id_prefix": "sim_reads", # Used by runner for output naming
            "no_aln_output": False,
            "log_level": "INFO",
            # Include original params for manifest that were direct args to runner
            "reference_fasta": str(ref_fasta),
            "depth": 50,
            "read_length": 150,
            "art_profile": "HS25",
            "mean_fragment_length": 400,
            "std_dev_fragment_length": 50,
            "is_paired_end": True,
            "art_platform": "illumina",
            "output_root_dir": str(output_root), # Store as string for manifest
            "run_name": DEFAULT_RUN_OUTPUT_DIR_NAME,
            "overwrite_output": False,
        }
    }

class TestSimulateShortReadsRunner:

    def test_successful_paired_end_simulation(self, tmp_path, mocker, default_runner_params):
        # 1. Setup paths and expected outputs
        # Use resolved paths from fixture where appropriate
        output_run_dir = default_runner_params["output_root_dir"] / default_runner_params["run_name"]
        id_prefix = default_runner_params["command_params"]["id_prefix"]
        expected_fq1_path = output_run_dir / (id_prefix + "1.fq")
        expected_fq2_path = output_run_dir / (id_prefix + "2.fq")
        expected_aln_path = output_run_dir / (id_prefix + ".aln")

        # Ensure the output directory that _prepare_run_output_directory is supposed to create, actually exists for the mock_art_cmd_effect
        output_run_dir.mkdir(parents=True, exist_ok=True)

        # 2. Mock external dependencies
        # Patch find_tool_path where it's looked up (in bb_runners module)
        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool_name:
                     default_runner_params["command_params"]["art_illumina_path"] if tool_name.startswith('art_') else
                     default_runner_params["command_params"]["samtools_path"] if tool_name == 'samtools' else None)

        mock_ensure_file_exists = mocker.patch('bb_runners.ensure_file_exists', return_value=Path(default_runner_params["reference_fasta"])) # Target bb_runners
        mock_check_fasta_indexed = mocker.patch('bb_runners.check_fasta_indexed') # Target bb_runners

        # Mock _prepare_run_output_directory
        mock_prepare_dir = mocker.patch('bb_runners._prepare_run_output_directory', return_value=output_run_dir)

        # Mock run_external_cmd for ART
        def mock_art_cmd_effect(cmd_list, **kwargs):
            # Based on ART's output naming: {id_prefix}[1/2].fq and {id_prefix}.aln
            # Extract the output prefix from the -o argument
            output_prefix_path_str = None
            # The runner calculates output_prefix as: run_output_dir / f"simulated_{art_platform}_reads"
            # It then passes this full path string to ART's -o argument.
            # So, cmd_list will contain this full path.
            for i, arg in enumerate(cmd_list):
                if arg == "-o" and i + 1 < len(cmd_list):
                    output_prefix_path_str = cmd_list[i+1] # This is output_run_dir / id_prefix
                    break

            if output_prefix_path_str:
                prefix_path = Path(output_prefix_path_str) # e.g., .../output_data/test_run_paired/sim_reads
                # Create dummy FASTQ files based on the prefix name
                (prefix_path.parent / (prefix_path.name + "1.fq")).touch()
                if default_runner_params["is_paired_end"]: # Check original params
                     (prefix_path.parent / (prefix_path.name + "2.fq")).touch()

                # Create dummy ALN file if not no_aln_output (from command_params)
                if not default_runner_params["command_params"]["no_aln_output"]:
                    (prefix_path.parent / (prefix_path.name + ".aln")).touch()
                return (0, "ART simulation successful", "")
            raise ValueError(f"ART command mock did not find -o argument correctly in {cmd_list}")

        mock_run_art_cmd = mocker.patch('bb_runners.run_external_cmd', side_effect=mock_art_cmd_effect) # Target bb_runners
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest') # Target bb_runners

        # Create dummy reference fasta
        Path(default_runner_params["reference_fasta"]).parent.mkdir(parents=True, exist_ok=True)
        Path(default_runner_params["reference_fasta"]).touch()

        # 3. Run the function - pass only known direct arguments
        runner_args_to_pass = {
            k: v for k, v in default_runner_params.items()
            if k in bb_runners.simulate_short_reads_runner.__code__.co_varnames
        }
        bb_runners.simulate_short_reads_runner(**runner_args_to_pass)


        # 4. Assertions
        # find_tool_path is already mocked via bb_runners.find_tool_path
        # No direct assertion on bb_utils.find_tool_path needed if already covered by the mock effect.
        # We can assert on the mock object if needed, e.g. mock_find_tool_path.assert_any_call(...)
        # For now, the side_effect implicitly tests it's called.

        # Check calls to ensure_file_exists
        # The runner calls ensure_file_exists(reference_fasta, "Reference FASTA")
        # where "Reference FASTA" is a positional argument for entity_name.
        mock_ensure_file_exists.assert_any_call(default_runner_params["reference_fasta"], "Reference FASTA")

        art_output_filename_prefix = f"simulated_{default_runner_params['art_platform']}_reads"
        # For these calls, entity_name is also passed positionally by the runner.
        mock_ensure_file_exists.assert_any_call(output_run_dir / (art_output_filename_prefix + "1.fq"), "Expected ART output R1 FASTQ")
        if default_runner_params["is_paired_end"]:
            mock_ensure_file_exists.assert_any_call(output_run_dir / (art_output_filename_prefix + "2.fq"), "Expected ART output R2 FASTQ")

        mock_check_fasta_indexed.assert_called_once_with(
            Path(default_runner_params["reference_fasta"]), # First positional argument
            default_runner_params["command_params"]["samtools_path"]  # Second positional argument
            # auto_index_if_missing uses its default value (False) in the runner's call
        )
        mock_prepare_dir.assert_called_once_with(
            default_runner_params["output_root_dir"],
            default_runner_params["run_name"],
            default_runner_params["overwrite_output"] # Passed positionally
        )

        assert mock_run_art_cmd.call_count == 1
        art_call_args_list = mock_run_art_cmd.call_args[0][0]
        assert default_runner_params["command_params"]["art_illumina_path"] in art_call_args_list
        assert "-l" in art_call_args_list and str(default_runner_params["read_length"]) in art_call_args_list
        assert "-f" in art_call_args_list and str(default_runner_params["depth"]) in art_call_args_list # ART uses -f for depth/coverage
        assert "-m" in art_call_args_list and str(default_runner_params["mean_fragment_length"]) in art_call_args_list
        assert "-s" in art_call_args_list and str(default_runner_params["std_dev_fragment_length"]) in art_call_args_list
        if default_runner_params["is_paired_end"]:
            assert "-p" in art_call_args_list

        # The output prefix for ART is `run_output_dir / f"simulated_{art_platform}_reads"`
        # The runner uses id_prefix from command_params for the manifest, but a fixed "simulated_{art_platform}_reads" for ART output files
        art_output_filename_prefix = f"simulated_{default_runner_params['art_platform']}_reads"
        assert "-o" in art_call_args_list and str(output_run_dir / art_output_filename_prefix) in art_call_args_list

        # ensure_file_exists calls are checked above with assert_has_calls

        mock_write_manifest.assert_called_once()
        manifest_call_args = mock_write_manifest.call_args

        assert manifest_call_args[0][0] == output_run_dir / "manifest.json" # manifest_path
        assert manifest_call_args[0][1] == default_runner_params["run_name"] # run_name
        assert manifest_call_args[0][2] == "short" # command_name in manifest (hardcoded in runner)

        manifest_params_written = manifest_call_args[0][3] # parameters
        assert manifest_params_written["reference_fasta"] == default_runner_params["reference_fasta"]
        assert manifest_params_written["read_length"] == default_runner_params["read_length"]
        assert manifest_params_written["is_paired_end"] == default_runner_params["is_paired_end"]

        manifest_outputs_written = manifest_call_args[0][4] # output_files

        # Expected paths for manifest are relative to output_run_dir, not output_run_dir.parent
        # And filenames are based on "simulated_{art_platform}_reads"
        expected_manifest_fq1 = {"name": "Simulated Reads (R1)", "path": art_output_filename_prefix + "1.fq", "type": "FASTQ"}
        assert any(o == expected_manifest_fq1 for o in manifest_outputs_written)
        if default_runner_params["is_paired_end"]:
            expected_manifest_fq2 = {"name": "Simulated Reads (R2)", "path": art_output_filename_prefix + "2.fq", "type": "FASTQ"}
            assert any(o == expected_manifest_fq2 for o in manifest_outputs_written)
        if not default_runner_params["command_params"]["no_aln_output"]:
            expected_manifest_aln = {"name": "ART Alignment ALN", "path": art_output_filename_prefix + ".aln", "type": "ALN"}
            assert any(o == expected_manifest_aln for o in manifest_outputs_written)
        else:
            # Ensure ALN is not in manifest if no_aln_output is True
            assert not any(o["type"] == "ALN" for o in manifest_outputs_written)

        # reference_genome_path is passed as a keyword argument
        assert manifest_call_args[1]['reference_genome_path'] == str(Path(default_runner_params["reference_fasta"]))
        # status is not explicitly passed by the runner, so it takes its default value in write_run_manifest.
        # Thus, it won't appear in call_args[1] (kwargs). We don't assert it here.

    # Placeholder for other tests
    def test_successful_single_end_simulation(self, tmp_path, mocker, default_runner_params):
        single_end_params = default_runner_params.copy()
        single_end_params["is_paired_end"] = False
        # Update command_params as well if it's used to derive is_paired_end for manifest or other logic
        single_end_params["command_params"] = default_runner_params["command_params"].copy()
        single_end_params["command_params"]["is_paired_end"] = False

        output_run_dir = single_end_params["output_root_dir"] / single_end_params["run_name"]
        # id_prefix = single_end_params["command_params"]["id_prefix"] # Not used for ART output file names
        art_output_filename_prefix = f"simulated_{single_end_params['art_platform']}_reads"
        expected_fq_path = output_run_dir / (art_output_filename_prefix + ".fq") # Single end uses .fq

        output_run_dir.mkdir(parents=True, exist_ok=True)

        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool_name:
                     single_end_params["command_params"]["art_illumina_path"] if tool_name.startswith('art_') else
                     single_end_params["command_params"]["samtools_path"] if tool_name == 'samtools' else None)

        mock_ensure_file_exists = mocker.patch('bb_runners.ensure_file_exists', return_value=Path(single_end_params["reference_fasta"]))
        mocker.patch('bb_runners.check_fasta_indexed')
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=output_run_dir)

        def mock_art_cmd_effect_single(cmd_list, **kwargs):
            output_prefix_path_str = None
            for i, arg in enumerate(cmd_list):
                if arg == "-o" and i + 1 < len(cmd_list):
                    output_prefix_path_str = cmd_list[i+1]
                    break
            if output_prefix_path_str:
                prefix_path = Path(output_prefix_path_str)
                (prefix_path.parent / (prefix_path.name + ".fq")).touch() # Single-end .fq
                if not single_end_params["command_params"]["no_aln_output"]:
                    (prefix_path.parent / (prefix_path.name + ".aln")).touch()
                return (0, "ART SE simulation successful", "")
            raise ValueError("ART SE command mock did not find -o argument correctly")

        mock_run_art_cmd = mocker.patch('bb_runners.run_external_cmd', side_effect=mock_art_cmd_effect_single)
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')

        Path(single_end_params["reference_fasta"]).touch() # Ensure dummy fasta exists

        runner_args_to_pass = {
            k: v for k, v in single_end_params.items()
            if k in bb_runners.simulate_short_reads_runner.__code__.co_varnames
        }
        bb_runners.simulate_short_reads_runner(**runner_args_to_pass)

        art_call_args_list = mock_run_art_cmd.call_args[0][0]
        assert "-p" not in art_call_args_list # Should not have paired-end flag
        mock_ensure_file_exists.assert_any_call(expected_fq_path, "Expected ART output FASTQ")

        manifest_call_args = mock_write_manifest.call_args
        manifest_params_written = manifest_call_args[0][3]
        assert manifest_params_written["is_paired_end"] == False

        manifest_outputs_written = manifest_call_args[0][4]
        expected_manifest_fq = {"name": "Simulated Reads", "path": art_output_filename_prefix + ".fq", "type": "FASTQ"}
        assert any(o == expected_manifest_fq for o in manifest_outputs_written)
        assert not any("R1" in o["name"] or "R2" in o["name"] for o in manifest_outputs_written if o["type"] == "FASTQ")


    def test_art_tool_failure(self, tmp_path, mocker, default_runner_params):
        output_run_dir = default_runner_params["output_root_dir"] / default_runner_params["run_name"]
        output_run_dir.mkdir(parents=True, exist_ok=True)

        mocker.patch('bb_runners.find_tool_path', side_effect=lambda tool_name:
                     default_runner_params["command_params"]["art_illumina_path"] if tool_name.startswith('art_') else
                     default_runner_params["command_params"]["samtools_path"] if tool_name == 'samtools' else None)

        mocker.patch('bb_runners.ensure_file_exists', return_value=Path(default_runner_params["reference_fasta"]))
        mocker.patch('bb_runners.check_fasta_indexed')
        mocker.patch('bb_runners._prepare_run_output_directory', return_value=output_run_dir)

        # Mock run_external_cmd for ART to raise BaseBuddyToolError
        art_command_str_part = default_runner_params["command_params"]["art_illumina_path"] # ART exe path
        mock_run_art_cmd = mocker.patch('bb_runners.run_external_cmd',
                                        side_effect=bb_utils.BaseBuddyToolError(
                                            message="ART simulation failed",
                                            command=[art_command_str_part, "-ss", "..."] # Dummy command list
                                        ))

        # write_run_manifest should not be called if ART fails before manifest generation
        mock_write_manifest = mocker.patch('bb_runners.write_run_manifest')

        Path(default_runner_params["reference_fasta"]).touch()

        runner_args_to_pass = {
            k: v for k, v in default_runner_params.items()
            if k in bb_runners.simulate_short_reads_runner.__code__.co_varnames
        }

        with pytest.raises(bb_utils.BaseBuddyToolError) as excinfo:
            bb_runners.simulate_short_reads_runner(**runner_args_to_pass)

        assert "ART simulation failed" in str(excinfo.value)
        mock_run_art_cmd.assert_called_once() # Ensure ART was attempted
        mock_write_manifest.assert_not_called() # Manifest should not be written on tool failure


    def test_missing_art_tool(self, tmp_path, mocker, default_runner_params):
        # Mock find_tool_path to raise BaseBuddyConfigError for art_illumina
        def find_tool_side_effect(tool_name):
            if tool_name.startswith('art_'):
                raise bb_utils.BaseBuddyConfigError(f"Tool {tool_name} not found")
            if tool_name == 'samtools':
                return default_runner_params["command_params"]["samtools_path"]
            return None

        mocker.patch('bb_runners.find_tool_path', side_effect=find_tool_side_effect)

        # Other essential mocks that might be called before find_tool_path for ART
        output_run_dir = default_runner_params["output_root_dir"] / default_runner_params["run_name"]
        # output_run_dir might not be created if find_tool_path fails early, so no need to mkdir it here.
        # mock_prepare_dir = mocker.patch('bb_runners._prepare_run_output_directory', return_value=output_run_dir)
        # mock_ensure_file_exists = mocker.patch('bb_runners.ensure_file_exists', return_value=Path(default_runner_params["reference_fasta"]))
        # mock_check_fasta_indexed = mocker.patch('bb_runners.check_fasta_indexed')
        # Path(default_runner_params["reference_fasta"]).touch() # Dummy fasta might not be needed if it fails before this stage

        runner_args_to_pass = {
            k: v for k, v in default_runner_params.items()
            if k in bb_runners.simulate_short_reads_runner.__code__.co_varnames
        }

        with pytest.raises(bb_utils.BaseBuddyConfigError, match="Tool art_illumina not found"):
            bb_runners.simulate_short_reads_runner(**runner_args_to_pass)


    def test_invalid_input_depth(self, tmp_path, default_runner_params):
        pytest.skip("Not yet implemented")

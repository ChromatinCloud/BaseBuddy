#!/usr/bin/env python3
"""
Unit tests for newly improved features in BaseBuddy
"""

import unittest
from unittest.mock import patch, MagicMock, call
from pathlib import Path
import sys
import tempfile
import json

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from basebuddy import runner, utils as bb_utils
from basebuddy.signature_utils import parse_signature_matrix_tsv


class TestLongReadSimulation(unittest.TestCase):
    """Test improvements to long read simulation."""
    
    def test_num_reads_parameter_validation(self):
        """Test that num_reads parameter is properly validated."""
        with self.assertRaises(bb_utils.BaseBuddyInputError) as context:
            runner.simulate_long(
                output_root_dir=Path("test_output"),
                reference_fasta="dummy.fa",
                depth=-5,  # Invalid negative depth
                model="nanopore_R9.4.1",
                command_params={"num_reads": None}
            )
        self.assertIn("Either depth", str(context.exception))
        self.assertIn("num_reads must be positive", str(context.exception))
    
    def test_num_reads_overrides_depth(self):
        """Test that num_reads properly overrides depth when specified."""
        with patch('basebuddy.utils.Command') as mock_command:
            with patch('basebuddy.utils.run_external_cmd'):
                with patch('basebuddy.utils.prepare_run_output_dir') as mock_prep:
                    with patch('basebuddy.utils.ensure_file_exists') as mock_ensure:
                        with patch('basebuddy.utils.check_fasta_indexed'):
                            with patch('pathlib.Path.exists', return_value=True):
                                mock_prep.return_value = Path("test_output/run")
                                mock_ensure.return_value = Path("ref.fa")
                                
                                # Create mock command builder
                                mock_cmd_instance = MagicMock()
                                mock_command.return_value = mock_cmd_instance
                                mock_cmd_instance.get_command_parts.return_value = ["nanosim-h", "simulate"]
                                
                                runner.simulate_long(
                                    output_root_dir=Path("test_output"),
                                    reference_fasta="ref.fa",
                                    depth=30,
                                    model="nanopore_R9.4.1",
                                    command_params={"num_reads": 1000}
                                )
                                
                                # Check that -N flag was used instead of -c
                                mock_cmd_instance.add_option.assert_any_call("-N", "1000")
                                # Ensure -c was not called with depth
                                calls = mock_cmd_instance.add_option.call_args_list
                                depth_calls = [c for c in calls if c[0][0] == "-c"]
                                self.assertEqual(len(depth_calls), 0)
    
    def test_output_detection_multiple_patterns(self):
        """Test that various NanoSim output patterns are detected."""
        with patch('basebuddy.utils.Command'):
            with patch('basebuddy.utils.run_external_cmd'):
                with patch('basebuddy.utils.prepare_run_output_dir') as mock_prep:
                    with patch('basebuddy.utils.ensure_file_exists') as mock_ensure:
                        with patch('basebuddy.utils.check_fasta_indexed'):
                            with patch('pathlib.Path.glob') as mock_glob:
                                with patch('pathlib.Path.exists') as mock_exists:
                                    mock_prep.return_value = Path("test_output/run")
                                    mock_ensure.return_value = Path("ref.fa")
                                    
                                    # Mock different output file patterns
                                    mock_exists.side_effect = lambda self: str(self).endswith("_unaligned_reads.fastq")
                                    mock_glob.return_value = [Path("test_output/run/nanosim_reads_unaligned_reads.fastq")]
                                    
                                    result = runner.simulate_long(
                                        output_root_dir=Path("test_output"),
                                        reference_fasta="ref.fa",
                                        depth=30,
                                        model="nanopore_R9.4.1"
                                    )
                                    
                                    # Check that unaligned reads were found
                                    self.assertTrue(any("unaligned" in f.get("name", "").lower() 
                                                      for f in result.get("output_files", [])))


class TestMutationalSignatures(unittest.TestCase):
    """Test improvements to mutational signatures functionality."""
    
    def test_signature_type_validation(self):
        """Test that signature types are properly validated."""
        with self.assertRaises(bb_utils.BaseBuddyInputError) as context:
            runner.simulate_signatures(
                output_root_dir=Path("test_output"),
                reference_fasta=None,
                sig_type="INVALID",  # Invalid signature type
                num_mutations=1000,
                sample_id="test"
            )
        self.assertIn("Invalid signature type", str(context.exception))
        self.assertIn("SBS", str(context.exception.details))
        self.assertIn("DBS", str(context.exception.details))
        self.assertIn("ID", str(context.exception.details))
    
    def test_default_reference_handling(self):
        """Test that default GRCh38 is used when no reference is provided."""
        with patch('basebuddy.utils.prepare_run_output_dir') as mock_prep:
            with patch('SigProfilerSimulator.SigProfilerSimulator') as mock_sps:
                mock_prep.return_value = Path("test_output/run")
                
                # Mock SigProfilerSimulator to avoid actual execution
                mock_sps_instance = MagicMock()
                mock_sps.return_value = mock_sps_instance
                
                runner.simulate_signatures(
                    output_root_dir=Path("test_output"),
                    reference_fasta=None,  # No reference provided
                    sig_type="SBS",
                    num_mutations=1000,
                    sample_id="test"
                )
                
                # Check that GRCh38 was used as genome build
                mock_sps.assert_called_once()
                call_args = mock_sps.call_args[1]
                self.assertEqual(call_args['genome_build'], 'GRCh38')
                self.assertIsNone(call_args['custom_genome'])
    
    def test_signature_data_parsing(self):
        """Test parsing of signature TSV files."""
        # Create a temporary TSV file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("MutationType\tSBS1\tSBS2\n")
            f.write("A[C>A]A\t0.011\t0.00018\n")
            f.write("A[C>A]C\t0.009\t0.00021\n")
            temp_path = Path(f.name)
        
        try:
            signatures = parse_signature_matrix_tsv(temp_path)
            
            # Verify structure
            self.assertIn("SBS1", signatures)
            self.assertIn("SBS2", signatures)
            self.assertIn("A[C>A]A", signatures["SBS1"])
            self.assertIn("A[C>A]C", signatures["SBS1"])
            
            # Verify values
            self.assertAlmostEqual(signatures["SBS1"]["A[C>A]A"], 0.011)
            self.assertAlmostEqual(signatures["SBS2"]["A[C>A]C"], 0.00021)
        finally:
            temp_path.unlink()


class TestImprovedErrorMessages(unittest.TestCase):
    """Test improved error messages across features."""
    
    def test_short_read_error_messages(self):
        """Test improved error messages for short read simulation."""
        # Test fragment length validation
        with self.assertRaises(bb_utils.BaseBuddyInputError) as context:
            runner.simulate_short(
                output_root_dir=Path("test_output"),
                reference_fasta="dummy.fa",
                depth=30,
                read_length=150,
                art_profile="HS25",
                mean_fragment_length=100,  # Too short for 150bp reads
                is_paired_end=True
            )
        self.assertIn("Fragment length too short", str(context.exception))
        self.assertIn("150bp", str(context.exception.details))
        self.assertIn("100bp", str(context.exception.details))
    
    def test_invalid_platform_error(self):
        """Test error message for invalid ART platform."""
        with self.assertRaises(bb_utils.BaseBuddyInputError) as context:
            runner.simulate_short(
                output_root_dir=Path("test_output"),
                reference_fasta="dummy.fa",
                depth=30,
                read_length=150,
                art_profile="HS25",
                art_platform="pacbio"  # Invalid platform
            )
        self.assertIn("Invalid ART platform", str(context.exception))
        self.assertIn("illumina", str(context.exception.details))
        self.assertIn("454", str(context.exception.details))
        self.assertIn("solid", str(context.exception.details))


class TestCLIIntegration(unittest.TestCase):
    """Test CLI integration for new features."""
    
    @patch('basebuddy.runner.simulate_signatures')
    def test_signature_cli_command(self, mock_sim_sig):
        """Test that signature CLI command properly calls runner."""
        from basebuddy.cli import app
        from typer.testing import CliRunner
        
        runner_cli = CliRunner()
        mock_sim_sig.return_value = {
            "run_name": "test_run",
            "output_directory": "/tmp/test",
            "manifest_path": "/tmp/test/manifest.json",
            "output_files": []
        }
        
        result = runner_cli.invoke(app, [
            "signature",
            "--sig-type", "SBS",
            "--num-mutations", "5000",
            "--sample-id", "test_sample"
        ])
        
        self.assertEqual(result.exit_code, 0)
        mock_sim_sig.assert_called_once()
        
        # Verify parameters passed
        call_args = mock_sim_sig.call_args[1]
        self.assertEqual(call_args['sig_type'], 'SBS')
        self.assertEqual(call_args['num_mutations'], 5000)
        self.assertEqual(call_args['sample_id'], 'test_sample')


if __name__ == '__main__':
    unittest.main()
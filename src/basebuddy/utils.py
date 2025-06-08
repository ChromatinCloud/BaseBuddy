# File 1: bb_utils.py (Updated)

import subprocess
import hashlib
import logging
import shutil
import sys
import time
import json
import datetime
from pathlib import Path
from typing import List, Union, Optional, Tuple, Dict, Any
from xml.etree.ElementTree import Element, SubElement, ElementTree, indent # For IGV XML
import pysam # For FASTA manipulation
import typer # Added

# --- Logger Setup ---
logger = logging.getLogger(__name__)

# --- Custom Exceptions ( 그대로 유지 ) ---
class BaseBuddyError(Exception):
    """Base class for all custom exceptions in BaseBuddy."""
    def __init__(self, message: str, details: Optional[str] = None):
        super().__init__(message)
        self.details = details

    def __str__(self):
        if self.details:
            return f"{super().__str__()} \nDetails: {self.details}"
        return super().__str__()

class BaseBuddyConfigError(BaseBuddyError):
    """For errors related to configuration or setup (e.g., missing tools)."""
    pass

class BaseBuddyInputError(BaseBuddyError, ValueError):
    """For errors related to invalid user inputs or parameters."""
    pass

class BaseBuddyFileError(BaseBuddyError, IOError):
    """For errors related to file operations (missing, permissions, format)."""
    pass

class BaseBuddyToolError(BaseBuddyError):
    """For errors arising from external tool execution."""
    def __init__(self, message: str, command: List[str], return_code: Optional[int] = None,
                 stdout: Optional[str] = None, stderr: Optional[str] = None):
        details = f"Command: {' '.join(command)}\n"
        if return_code is not None:
            details += f"Return Code: {return_code}\n"
        if stderr:
            details += f"Stderr: {stderr.strip()}\n"
        if stdout and not stderr and return_code !=0 :
            details += f"Stdout: {stdout.strip()}\n"
        super().__init__(message, details.strip())
        self.command = command
        self.return_code = return_code
        self.stdout = stdout
        self.stderr = stderr

class BaseBuddyChecksumError(BaseBuddyFileError):
    """For checksum verification failures."""
    pass

# --- External Tool Path Checking ( 그대로 유지 ) ---
def find_tool_path(tool_name: str) -> str:
    tool_path = shutil.which(tool_name)
    if not tool_path:
        logger.error(f"External tool '{tool_name}' not found in PATH.")
        raise BaseBuddyConfigError(
            f"External tool '{tool_name}' not found in PATH. "
            "Please ensure it is installed and your system's PATH environment variable is configured correctly."
        )
    logger.debug(f"Found tool '{tool_name}' at '{tool_path}'.")
    return tool_path

# --- Subprocess Runner ( 그대로 유지 ) ---
def run_external_cmd(
    cmd: List[str],
    cwd: Optional[Union[str, Path]] = None,
    capture_output: bool = True,
    stream_output: bool = False,
    timeout_seconds: Optional[float] = None,
    encoding: str = 'utf-8',
    env: Optional[dict] = None
) -> Tuple[int, str, str]:
    cmd_str = " ".join(map(str, cmd))
    logger.info(f"Running command: {cmd_str}")
    if cwd:
        logger.info(f"Working directory: {cwd}")

    stdout_lines = []
    stderr_lines = []

    try:
        if stream_output and not capture_output:
            process = subprocess.Popen(
                cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, encoding=encoding, errors='replace', env=env, bufsize=1
            )
            if process.stdout:
                for line in iter(process.stdout.readline, ''):
                    line_strip = line.strip()
                    logger.info(f"[TOOL_STDOUT] {line_strip}")
                    stdout_lines.append(line_strip)
                process.stdout.close()
            if process.stderr:
                for line in iter(process.stderr.readline, ''):
                    line_strip = line.strip()
                    logger.error(f"[TOOL_STDERR] {line_strip}")
                    stderr_lines.append(line_strip)
                process.stderr.close()
            
            process.wait(timeout=timeout_seconds)
            return_code = process.returncode
            stdout_str = "\n".join(stdout_lines)
            stderr_str = "\n".join(stderr_lines)
        else:
            completed_process = subprocess.run(
                cmd, cwd=cwd, capture_output=True, text=True,
                encoding=encoding, errors='replace', timeout=timeout_seconds, env=env, check=False
            )
            return_code = completed_process.returncode
            stdout_str = completed_process.stdout.strip() if completed_process.stdout else ""
            stderr_str = completed_process.stderr.strip() if completed_process.stderr else ""
            if stdout_str: logger.debug(f"Command stdout: {stdout_str}")
            if stderr_str: logger.debug(f"Command stderr: {stderr_str}")

        if return_code != 0:
            logger.error(f"Command failed with exit code {return_code}: {cmd_str}")
            raise BaseBuddyToolError(
                f"External command execution failed.",
                command=cmd, return_code=return_code, stdout=stdout_str, stderr=stderr_str
            )
        
        logger.info(f"Command executed successfully: {cmd_str}")
        return return_code, stdout_str, stderr_str

    except subprocess.TimeoutExpired as e:
        logger.error(f"Command timed out after {timeout_seconds}s: {cmd_str}")
        raise BaseBuddyToolError(
            f"Command execution timed out after {timeout_seconds} seconds.",
            command=cmd, stderr=str(e)
        )
    except FileNotFoundError:
        logger.error(f"Command not found: {cmd[0]}. Is it installed and in PATH?")
        raise BaseBuddyConfigError(f"The command '{cmd[0]}' was not found. Ensure it is installed and in your PATH.")
    except PermissionError:
        logger.error(f"Permission denied for command: {cmd[0]}. Check execute permissions.")
        raise BaseBuddyConfigError(f"Permission denied when trying to execute '{cmd[0]}'.")
    except Exception as e:
        logger.exception(f"An unexpected error occurred while running command '{cmd_str}'.")
        raise BaseBuddyToolError(f"An unexpected error occurred during command execution.", command=cmd, stderr=str(e))

# --- File System Utilities ( 그대로 유지, 약간 수정된 ensure_directory_exists ) ---
def ensure_file_exists(file_path: Union[str, Path], entity_name: str = "Input file") -> Path:
    path = Path(file_path)
    if not path.exists():
        raise BaseBuddyFileError(f"{entity_name} not found at specified path: {path}")
    if not path.is_file():
        raise BaseBuddyFileError(f"Path exists but is not a file for {entity_name}: {path}")
    try:
        with open(path, 'rb') as f:
            f.read(1)
    except IOError as e:
        raise BaseBuddyFileError(f"Cannot read {entity_name} at {path}. Check permissions. Error: {e}")
    logger.debug(f"{entity_name} verified at {path}")
    return path

def ensure_directory_exists(dir_path: Union[str, Path], entity_name: str = "Directory", create: bool = False, check_write_perm: bool = True) -> Path:
    path = Path(dir_path)
    if path.exists():
        if not path.is_dir():
            raise BaseBuddyFileError(f"Path exists for {entity_name} but is not a directory: {path}")
    else:
        if create:
            logger.info(f"{entity_name} not found at {path}. Creating it.")
            try:
                path.mkdir(parents=True, exist_ok=True)
            except OSError as e:
                raise BaseBuddyFileError(f"Could not create {entity_name} at {path}. Error: {e}")
        else:
            raise BaseBuddyFileError(f"{entity_name} not found: {path}")
    
    if check_write_perm:
        try:
            temp_file = path / f".__write_test_{time.time_ns()}"
            temp_file.touch()
            temp_file.unlink(missing_ok=True) # missing_ok for robustness if test is interrupted
        except OSError as e:
            raise BaseBuddyFileError(f"Cannot write to {entity_name} at {path}. Check permissions. Error: {e}")
    logger.debug(f"{entity_name} verified at {path}")
    return path

def prepare_run_output_dir(
    output_root_dir: Path,
    run_name: str,
    overwrite_output: bool
) -> Path:
    """Creates and returns the specific output directory for the current run."""
    run_output_dir = output_root_dir / run_name
    if run_output_dir.exists() and not overwrite_output and any(run_output_dir.iterdir()):
        # Check if directory is not empty. any(run_output_dir.iterdir()) is a simple way.
        raise BaseBuddyFileError(
            f"Run output directory '{run_output_dir}' is not empty and overwrite is not permitted. "
            "Use --overwrite, --run-name with a new name, or clear the directory."
        )
    # ensure_directory_exists will handle creation and also check if it's a directory if it exists.
    # It also handles parent creation.
    ensure_directory_exists(run_output_dir, "Run output directory", create=True, check_write_perm=True) # check_write_perm is True by default
    logger.info(f"Using run output directory: {run_output_dir.resolve()}")
    return run_output_dir

def check_fasta_indexed(
    reference_path: Path,
    samtools_path: str,
    auto_index_if_missing: bool = False # New parameter, default to False to not change existing explicit calls without it
) -> None:
    """
    Checks if a FASTA file is indexed ('.fai' file exists).
    Optionally attempts to index it if `auto_index_if_missing` is True.
    """
    fai_path = reference_path.with_suffix(reference_path.suffix + ".fai")
    if not fai_path.exists():
        logger.warning(f"FASTA index (.fai) not found for {reference_path}.")
        if auto_index_if_missing:
            logger.info(f"Attempting to auto-index FASTA file: {reference_path} using {samtools_path}")
            try:
                # Ensure the directory containing the reference_path is writable for the .fai file
                # This is usually the case if the user provided the FASTA path.
                run_external_cmd([samtools_path, "faidx", str(reference_path)], cwd=reference_path.parent)
                if not fai_path.exists(): # Re-check if index was created
                     # This specific error message for tool failure is better than generic BaseBuddyFileError here
                     raise BaseBuddyToolError(
                         f"FASTA index (.fai) still not found after attempting to create it with 'samtools faidx'. "
                         f"The command might have failed silently or the output was unexpected.",
                         command=[samtools_path, "faidx", str(reference_path)]
                     )
                logger.info(f"Successfully auto-indexed FASTA file: {reference_path}")
            except BaseBuddyToolError as e: # Catch specific tool errors from run_external_cmd
                # Re-raise as a BaseBuddyFileError or a more specific indexing error if preferred
                raise BaseBuddyFileError(
                    f"Failed to auto-index FASTA file {reference_path}. Please index it manually. Original error: {e.args[0]}",
                    details=e.details
                )
            except Exception as e: # Catch other unexpected errors during indexing attempt
                 raise BaseBuddyFileError(
                    f"An unexpected error occurred while attempting to auto-index FASTA file {reference_path}. Please index it manually. Error: {e}"
                )
        else: # auto_index_if_missing is False and index is missing
            raise BaseBuddyFileError(
                f"Reference FASTA is not indexed (.fai missing): {reference_path}. "
                f"Please run `{samtools_path} faidx {reference_path}` first or enable auto-indexing."
            )
    logger.debug(f"FASTA index verified for {reference_path}")


def check_bam_indexed(bam_path: Path, samtools_path: str, auto_index_if_missing: bool = False) -> None:
    bai_path_1 = bam_path.with_suffix(bam_path.suffix + ".bai")
    bai_path_2 = bam_path.with_suffix(".bai")

    if not (bai_path_1.exists() or bai_path_2.exists()):
        logger.warning(f"BAM index (.bai) not found for {bam_path}.")
        if auto_index_if_missing:
            logger.info(f"Attempting to auto-index BAM file: {bam_path}")
            try:
                run_external_cmd([samtools_path, "index", str(bam_path)])
                if not (bai_path_1.exists() or bai_path_2.exists()):
                     raise BaseBuddyToolError(f"BAM index (.bai) still not found after attempting to create it for {bam_path}.",
                                              command=[samtools_path, "index", str(bam_path)])
                logger.info(f"Successfully auto-indexed BAM file: {bam_path}")
            except BaseBuddyToolError as e:
                raise BaseBuddyFileError(
                    f"Failed to auto-index BAM file {bam_path}. Please index it manually using `samtools index`. Original error: {e}"
                )
        else:
            raise BaseBuddyFileError(
                f"Input BAM is not indexed (.bai missing): {bam_path}. "
                f"Please run `{samtools_path} index {bam_path}` or use an auto-index option."
            )
    logger.debug(f"BAM index verified for {bam_path}")

# --- Checksum Verification ( 그대로 유지 ) ---
def verify_file_checksum(file_path: Path, expected_checksum: str, algorithm: str = "sha256"):
    logger.info(f"Verifying {algorithm} checksum for {file_path}...")
    ensure_file_exists(file_path, f"File for checksum '{file_path.name}'")
    
    hasher = hashlib.new(algorithm)
    try:
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                hasher.update(chunk)
    except IOError as e:
        raise BaseBuddyFileError(f"Could not read file {file_path} for checksum calculation. Error: {e}")

    calculated_checksum = hasher.hexdigest()
    if calculated_checksum.lower() != expected_checksum.lower():
        logger.error(f"Checksum mismatch for {file_path}. Expected: {expected_checksum}, Got: {calculated_checksum}")
        raise BaseBuddyChecksumError(
            f"{algorithm.upper()} checksum mismatch for {file_path}.",
            details=f"Expected: {expected_checksum}\nCalculated: {calculated_checksum}"
        )
    logger.info(f"Checksum verified successfully for {file_path}.")

# --- New: Run Naming and Manifest Utilities ---
def generate_unique_run_name(command_name: str) -> str:
    """Generates a unique run name using command name and timestamp."""
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{command_name}_{timestamp}"

def write_run_manifest(
    manifest_path: Path,
    run_name: str,
    command_name: str,
    parameters: Dict[str, Any],
    output_files: List[Dict[str, str]], # Each dict: {"name": "desc", "path": "rel_path", "type": "BAM|VCF|FASTQ|IGV_SESSION|LOG"}
    reference_genome_path: Optional[str] = None,
    status: str = "completed"
) -> None:
    """Writes a manifest.json file for a run."""
    manifest_data = {
        "run_name": run_name,
        "command": command_name,
        "timestamp_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "status": status,
        "parameters": {k: str(v) if isinstance(v, Path) else v for k, v in parameters.items()}, # Ensure paths are strings
        "reference_genome_path": str(reference_genome_path) if reference_genome_path else None,
        "outputs": output_files
    }
    try:
        with open(manifest_path, 'w') as f:
            json.dump(manifest_data, f, indent=4)
        logger.info(f"Run manifest written to {manifest_path}")
    except IOError as e:
        logger.error(f"Failed to write run manifest to {manifest_path}: {e}")
        # Depending on policy, this could raise an error or just be a warning

def read_run_manifest(manifest_path: Path) -> Optional[Dict[str, Any]]:
    """Reads a manifest.json file."""
    if not manifest_path.is_file():
        logger.warning(f"Manifest file not found: {manifest_path}")
        return None
    try:
        with open(manifest_path, 'r') as f:
            return json.load(f)
    except (IOError, json.JSONDecodeError) as e:
        logger.error(f"Failed to read or parse manifest file {manifest_path}: {e}")
        return None

# --- New: IGV Session Generation ---
def generate_igv_session_xml(
    session_file_path: Path,
    genome_fasta_path: str, # Absolute path to FASTA
    tracks: List[Dict[str, str]], # Each dict: {"name": "Track Name", "path": "relative/path/to/file.bam", "format": "bam", "type": "alignment/variant"}
    # IGV uses 'locus' for initial view, e.g., "chr1:1000-2000" or "All"
    initial_locus: str = "All"
) -> None:
    """
    Generates an IGV session XML file.
    Paths for tracks in the XML should be relative to the session file if possible, or absolute.
    Here, we assume track paths are already relative to the session file's directory.
    """
    logger.info(f"Generating IGV session file at: {session_file_path}")
    
    # Ensure genome_fasta_path is absolute for IGV stability
    abs_genome_fasta_path = str(Path(genome_fasta_path).resolve())

    root = Element("Global")
    root.set("genome", abs_genome_fasta_path)
    root.set("locus", initial_locus) # IGV 2.8+ uses locus, older versions might use other tags.
    root.set("version", "4") # Common IGV session version

    resources_node = SubElement(root, "Resources")
    for track in tracks:
        track_path = track.get("path")
        track_name = track.get("name", Path(track_path).name) # Default to filename if name not provided
        
        resource_node = SubElement(resources_node, "Resource")
        resource_node.set("name", track_name)
        resource_node.set("path", track_path) # This path should be findable by IGV
        # IGV often infers type from extension, but explicit hints can be added if needed:
        # if "format" in track: resource_node.set("format", track["format"])
        # if "type" in track: resource_node.set("type", track["type"])

    # Pretty print XML
    try:
        # For Python 3.9+
        indent(root) # type: ignore
    except AttributeError: # For older Python versions that don't have xml.etree.ElementTree.indent
        pass # XML will not be pretty-printed but still valid

    tree = ElementTree(root)
    try:
        tree.write(session_file_path, encoding="UTF-8", xml_declaration=True)
        logger.info(f"IGV session file generated: {session_file_path}")
    except IOError as e:
        logger.error(f"Failed to write IGV session file {session_file_path}: {e}")
        raise BaseBuddyFileError(f"Could not write IGV session file to {session_file_path}", details=str(e))


# --- FASTA Manipulation ---
def apply_variants_to_fasta(
    original_fasta_path: Path,
    variants_list: List[Dict[str, Any]],
    output_mutated_fasta_path: Path
) -> None:
    """
    Applies a list of SNV variants to a FASTA file and writes a new FASTA file.
    Assumes variants are 1-based. Converts to 0-based for sequence manipulation.
    Currently handles only SNVs (single nucleotide variants where ref and alt are single bases).
    """
    ensure_file_exists(original_fasta_path, "Original FASTA for mutation")
    logger.info(f"Applying {len(variants_list)} variants from list to {original_fasta_path}, outputting to {output_mutated_fasta_path}")

    variants_by_chrom: Dict[str, List[Dict[str, Any]]] = {}
    for var in variants_list:
        chrom = var.get("chromosome")
        if not chrom:
            logger.warning(f"Skipping variant with missing chromosome: {var}")
            continue
        if chrom not in variants_by_chrom:
            variants_by_chrom[chrom] = []
        variants_by_chrom[chrom].append(var)

    # Sort variants by position for each chromosome (ascending)
    for chrom in variants_by_chrom:
        variants_by_chrom[chrom].sort(key=lambda v: v.get("position", 0))

    try:
        with pysam.FastaFile(str(original_fasta_path)) as fasta_in, \
             open(output_mutated_fasta_path, 'w') as fasta_out:

            processed_chroms_with_variants = set()

            for chrom_name in fasta_in.references:
                original_sequence = fasta_in.fetch(chrom_name)

                if chrom_name not in variants_by_chrom:
                    # No variants for this chromosome, write original sequence
                    fasta_out.write(f">{chrom_name}\n{original_sequence}\n")
                    continue

                logger.debug(f"Processing chromosome {chrom_name} for mutations.")
                seq_list = list(original_sequence)
                num_applied_on_chrom = 0

                for variant in variants_by_chrom[chrom_name]:
                    pos_1_based = variant.get("position")
                    ref_allele = variant.get("ref_allele")
                    alt_allele = variant.get("alt_allele")

                    if not all([pos_1_based, ref_allele, alt_allele]):
                        logger.warning(f"Skipping incomplete variant on {chrom_name}: {variant}")
                        continue

                    # Basic SNV check (single base ref and alt)
                    if not (len(ref_allele) == 1 and len(alt_allele) == 1):
                        logger.warning(f"Skipping non-SNV variant on {chrom_name} at pos {pos_1_based}: Ref '{ref_allele}', Alt '{alt_allele}'. Only SNVs are currently supported by apply_variants_to_fasta.")
                        continue

                    pos_0_based = int(pos_1_based) - 1

                    if 0 <= pos_0_based < len(seq_list):
                        actual_ref_at_pos = seq_list[pos_0_based]
                        if actual_ref_at_pos.upper() != ref_allele.upper():
                            logger.warning(
                                f"Reference mismatch on {chrom_name} at 1-based position {pos_1_based}. "
                                f"Expected FASTA ref: '{ref_allele}', Actual FASTA ref: '{actual_ref_at_pos}'. "
                                f"Proceeding with substitution using alt allele '{alt_allele}'."
                            )

                        seq_list[pos_0_based] = alt_allele
                        num_applied_on_chrom +=1
                    else:
                        logger.warning(
                            f"Variant position {pos_1_based} on chromosome {chrom_name} is out of bounds "
                            f"(sequence length {len(seq_list)}). Skipping variant."
                        )

                fasta_out.write(f">{chrom_name}\n{''.join(seq_list)}\n")
                processed_chroms_with_variants.add(chrom_name)
                if num_applied_on_chrom > 0:
                    logger.info(f"Applied {num_applied_on_chrom} SNV(s) to chromosome {chrom_name}.")

            # Check if any variants were for chromosomes not in the FASTA
            for chrom_name in variants_by_chrom:
                if chrom_name not in fasta_in.references and chrom_name not in processed_chroms_with_variants : # fasta_in.references might be empty if file is bad
                     logger.warning(f"Chromosome '{chrom_name}' specified in variants list was not found in the reference FASTA file '{original_fasta_path}'.")


    except FileNotFoundError:
        raise BaseBuddyFileError(f"Original FASTA file not found at {original_fasta_path}")
    except (IOError, ValueError) as e: # ValueError for pysam issues like bad fasta
        logger.error(f"Error during FASTA processing or writing mutated FASTA: {e}")
        raise BaseBuddyFileError(f"Error processing FASTA file {original_fasta_path} or writing to {output_mutated_fasta_path}. Details: {e}")
    except Exception as e: # Catch-all for unexpected pysam errors or other issues
        logger.exception(f"An unexpected error occurred in apply_variants_to_fasta: {e}")
        raise BaseBuddyError(f"An unexpected error occurred while applying variants to FASTA. Details: {e}")

    logger.info(f"Finished applying variants. Mutated FASTA written to {output_mutated_fasta_path}")


def extract_ranges_to_fasta(
    source_fasta_path: Path,
    ranges_list: List[str], # List of strings like "chr1:1000-2000"
    output_subset_fasta_path: Path,
    samtools_path: Optional[str] = None # Optional: if not provided, will try to find it
) -> None:
    """
    Extracts specified genomic ranges from a source FASTA file and writes them
    to a new FASTA file using `samtools faidx`.
    """
    ensure_file_exists(source_fasta_path, "Source FASTA for range extraction")
    if not ranges_list:
        raise BaseBuddyInputError("Genomic ranges list cannot be empty for extraction.")

    if not samtools_path:
        samtools_path = find_tool_path("samtools")

    # Ensure the source FASTA is indexed. faidx requires it.
    # check_fasta_indexed will raise an error if it can't index or find the index.
    check_fasta_indexed(source_fasta_path, samtools_path, auto_index_if_missing=True)

    cmd = [samtools_path, "faidx", str(source_fasta_path)] + ranges_list
    cmd_str = " ".join(cmd)
    logger.info(f"Extracting ranges using command: {cmd_str} > {output_subset_fasta_path}")

    try:
        with open(output_subset_fasta_path, 'w') as outfile_handle:
            # Using subprocess.run directly to redirect stdout to a file.
            # `run_external_cmd` is more for capturing output as strings or streaming to logs.
            process = subprocess.run(
                cmd,
                stdout=outfile_handle,
                stderr=subprocess.PIPE, # Capture stderr to check for errors
                text=True,
                check=False # We will check returncode manually
            )
            if process.returncode != 0:
                stderr_output = process.stderr.strip()
                logger.error(f"'samtools faidx' command failed with exit code {process.returncode}. Stderr: {stderr_output}")
                # Attempt to remove potentially incomplete output file
                try:
                    output_subset_fasta_path.unlink(missing_ok=True)
                except OSError:
                    logger.warning(f"Could not remove incomplete subset FASTA file: {output_subset_fasta_path}")
                raise BaseBuddyToolError(
                    message=f"'samtools faidx' for range extraction failed.",
                    command=cmd,
                    return_code=process.returncode,
                    stderr=stderr_output
                )

        # Verify that the output file was created and is not empty
        if not output_subset_fasta_path.exists() or output_subset_fasta_path.stat().st_size == 0:
            # This might happen if samtools faidx ran but produced no output for valid (but empty) regions,
            # or if regions were invalid but samtools didn't error in a way caught above.
            stderr_msg_if_any = process.stderr.strip() if process.returncode !=0 else "No specific error from samtools, but output is empty."
            logger.warning(f"Subset FASTA file {output_subset_fasta_path} is empty or was not created after 'samtools faidx'. Samtools stderr (if any): '{stderr_msg_if_any}'")
            # Decide if this is an error or just a warning. For now, a warning, but could be an error.
            # If it should be an error:
            # raise BaseBuddyFileError(f"Output subset FASTA {output_subset_fasta_path} is empty or missing after successful samtools faidx command. Ranges might have been invalid or produced no sequence.")

        logger.info(f"Subset FASTA with specified ranges written to {output_subset_fasta_path}")

    except FileNotFoundError: # For samtools_path itself, though find_tool_path should catch this.
        logger.error(f"Command not found: {samtools_path}. Is samtools installed and in PATH?")
        raise BaseBuddyConfigError(f"The command '{samtools_path}' was not found. Ensure it is installed and in your PATH.")
    except subprocess.TimeoutExpired: # Should not happen with default run() unless timeout is added
        logger.error(f"samtools faidx command timed out: {cmd_str}")
        raise BaseBuddyToolError(f"samtools faidx command timed out.", command=cmd)
    except Exception as e:
        logger.exception(f"An unexpected error occurred while extracting ranges with samtools faidx: {e}")
        raise BaseBuddyError(f"An unexpected error occurred during FASTA range extraction. Details: {e}")


# --- New Command Class and Typer Runner ---

class Command:
    """A helper class to build command-line commands programmatically."""
    def __init__(self, base_command: str):
        self.base_command_name = base_command # Store for error messages
        # find_tool_path will raise BaseBuddyConfigError if not found.
        # This error should be caught by the calling CLI layer and handled with Typer.
        self.executable_path = find_tool_path(base_command)
        self.parts = [self.executable_path]

    def add_option(self, option: str, value: Any):
        """Adds an option with a value, e.g., `-o output.txt`."""
        if value is not None:
            self.parts.extend([option, str(value)])
        return self # Allow chaining

    def add_flag(self, flag: str, condition: bool = True):
        """Adds a flag if a condition is met, e.g., `--verbose`."""
        if condition:
            self.parts.append(flag)
        return self # Allow chaining

    def get_command_parts(self) -> List[str]:
        """Returns the complete command as a list of strings."""
        return self.parts

    def get_command_str(self) -> str:
        """Returns the complete command as a single string (for shell=True or display)."""
        return " ".join(self.parts)


def run_external_command_typer(cmd_parts: List[str], use_shell: bool = False, cwd: Optional[Union[str, Path]] = None) -> subprocess.CompletedProcess:
    """
    Executes an external command using subprocess.run, tailored for Typer CLI output.
    Uses a list of command parts to avoid issues with shell=True and string formatting if not careful.
    If use_shell is True, cmd_parts will be joined to a string.
    """
    if not cmd_parts:
        typer.secho("Error: No command provided to run_external_command_typer.", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    if use_shell:
        cmd_to_run: Union[str, List[str]] = " ".join(cmd_parts)
    else:
        cmd_to_run = cmd_parts

    cmd_display_str = " ".join(cmd_parts) # For display purposes

    typer.echo(f"Executing: {cmd_display_str}")
    if cwd:
        typer.echo(f"Working directory: {cwd}")

    try:
        result = subprocess.run(
            cmd_to_run,
            shell=use_shell, # Be cautious with shell=True if cmd_parts contains user input not properly sanitized.
            capture_output=True,
            text=True,
            check=False, # We check the returncode manually
            cwd=cwd
        )

        if result.returncode != 0:
            typer.secho(f"Error running command: {cmd_display_str}", fg=typer.colors.RED, err=True)
            # Log full details for debugging, but typer output might be more concise
            logger.error(f"Command failed: {cmd_display_str}\nReturn Code: {result.returncode}\nStdout: {result.stdout}\nStderr: {result.stderr}")
            typer.secho(f"Stderr:\n{result.stderr}", fg=typer.colors.RED, err=True) # Show stderr
            raise typer.Exit(code=1)

        if result.stdout:
            logger.info(f"Command stdout: {result.stdout}")
        if result.stderr: # Some tools print to stderr even on success
            logger.info(f"Command stderr (non-error): {result.stderr}")

        return result
    except FileNotFoundError: # Should be caught by Command class, but as a safeguard if called directly
        typer.secho(f"Error: Command not found '{cmd_parts[0]}'. Is it installed and in PATH?", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except PermissionError:
        typer.secho(f"Error: Permission denied for command '{cmd_parts[0]}'. Check execute permissions.", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    except Exception as e: # Catch other potential subprocess errors
        typer.secho(f"An unexpected error occurred while trying to run command: {cmd_display_str}", fg=typer.colors.RED, err=True)
        typer.secho(str(e), fg=typer.colors.RED, err=True)
        logger.exception(f"Unexpected error in run_external_command_typer for command: {cmd_display_str}")
        raise typer.Exit(code=1)

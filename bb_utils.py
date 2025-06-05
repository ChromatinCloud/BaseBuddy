# File 1: bb_utils.py

import subprocess
import hashlib
import logging
import shutil
import sys
import time
from pathlib import Path
from typing import List, Union, Optional, Tuple

# --- Logger Setup ---
# Configure logger for this module. It will inherit the root logger's configuration
# if set up in the main CLI script.
logger = logging.getLogger(__name__)

# --- Custom Exceptions ---
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
        if stdout and not stderr and return_code !=0 : # Sometimes errors go to stdout
            details += f"Stdout: {stdout.strip()}\n"
        super().__init__(message, details.strip())
        self.command = command
        self.return_code = return_code
        self.stdout = stdout
        self.stderr = stderr

class BaseBuddyChecksumError(BaseBuddyFileError):
    """For checksum verification failures."""
    pass

# --- External Tool Path Checking ---
def find_tool_path(tool_name: str) -> str:
    """
    Checks if an external tool exists in PATH and returns its full path.
    Raises BaseBuddyConfigError if not found.
    """
    tool_path = shutil.which(tool_name)
    if not tool_path:
        logger.error(f"External tool '{tool_name}' not found in PATH.")
        raise BaseBuddyConfigError(
            f"External tool '{tool_name}' not found in PATH. "
            "Please ensure it is installed and your system's PATH environment variable is configured correctly."
        )
    logger.debug(f"Found tool '{tool_name}' at '{tool_path}'.")
    return tool_path

# --- Subprocess Runner ---
def run_external_cmd(
    cmd: List[str],
    cwd: Optional[Union[str, Path]] = None,
    capture_output: bool = True,
    stream_output: bool = False,
    timeout_seconds: Optional[float] = None,
    encoding: str = 'utf-8',
    env: Optional[dict] = None
) -> Tuple[int, str, str]:
    """
    Runs an external command with enhanced error handling, logging, and options.

    Args:
        cmd: Command and arguments as a list of strings.
        cwd: Current working directory.
        capture_output: If True, stdout and stderr are captured. Mutually exclusive with stream_output for capture.
        stream_output: If True, streams stdout/stderr to BaseBuddy's logger in real-time.
        timeout_seconds: Optional timeout for the command.
        encoding: Encoding for decoding stdout/stderr.
        env: Optional environment variables for the subprocess.

    Returns:
        A tuple (return_code, stdout_str, stderr_str).
        stdout_str and stderr_str are captured output if capture_output=True,
        otherwise they are minimal messages indicating streaming.

    Raises:
        BaseBuddyToolError: If the command fails (non-zero exit or timeout) or cannot be started.
    """
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
                text=True, encoding=encoding, errors='replace', env=env, bufsize=1 # Line buffered
            )
            # Stream stdout
            if process.stdout:
                for line in iter(process.stdout.readline, ''):
                    line_strip = line.strip()
                    logger.info(f"[TOOL_STDOUT] {line_strip}")
                    stdout_lines.append(line_strip)
                process.stdout.close()
            # Stream stderr
            if process.stderr:
                for line in iter(process.stderr.readline, ''):
                    line_strip = line.strip()
                    logger.error(f"[TOOL_STDERR] {line_strip}") # Log tool's stderr as error
                    stderr_lines.append(line_strip)
                process.stderr.close()
            
            process.wait(timeout=timeout_seconds)
            return_code = process.returncode
            stdout_str = "\n".join(stdout_lines)
            stderr_str = "\n".join(stderr_lines)

        else: # Default: capture output
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
    except Exception as e: # Catch-all for other subprocess or unexpected errors
        logger.exception(f"An unexpected error occurred while running command '{cmd_str}'.") # Includes stack trace
        raise BaseBuddyToolError(f"An unexpected error occurred during command execution.", command=cmd, stderr=str(e))

# --- File System Utilities ---
def ensure_file_exists(file_path: Union[str, Path], entity_name: str = "Input file") -> Path:
    path = Path(file_path)
    if not path.exists():
        raise BaseBuddyFileError(f"{entity_name} not found at specified path: {path}")
    if not path.is_file():
        raise BaseBuddyFileError(f"Path exists but is not a file for {entity_name}: {path}")
    # Basic read permission check (more robust checks are OS-dependent)
    try:
        with open(path, 'rb') as f:
            f.read(1) # Try to read one byte
    except IOError as e:
        raise BaseBuddyFileError(f"Cannot read {entity_name} at {path}. Check permissions. Error: {e}")
    logger.debug(f"{entity_name} verified at {path}")
    return path

def ensure_directory_exists(dir_path: Union[str, Path], entity_name: str = "Directory", create: bool = False) -> Path:
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
    # Basic write permission check
    try:
        temp_file = path / f".__write_test_{time.time_ns()}"
        temp_file.touch()
        temp_file.unlink()
    except OSError as e:
        raise BaseBuddyFileError(f"Cannot write to {entity_name} at {path}. Check permissions. Error: {e}")
    logger.debug(f"{entity_name} verified at {path}")
    return path


def check_fasta_indexed(reference_path: Path, samtools_path: str) -> None:
    fai_path = reference_path.with_suffix(reference_path.suffix + ".fai")
    if not fai_path.exists():
        logger.warning(f"FASTA index (.fai) not found for {reference_path}.")
        raise BaseBuddyFileError(
            f"Reference FASTA is not indexed (.fai missing): {reference_path}. "
            f"Please run `{samtools_path} faidx {reference_path}` first."
        )
    logger.debug(f"FASTA index verified for {reference_path}")

def check_bam_indexed(bam_path: Path, samtools_path: str, auto_index_if_missing: bool = False) -> None:
    bai_path_1 = bam_path.with_suffix(bam_path.suffix + ".bai")
    bai_path_2 = bam_path.with_suffix(".bai") # For cases like .cram.bai -> .crai, .bam -> .bai

    if not (bai_path_1.exists() or bai_path_2.exists()):
        logger.warning(f"BAM index (.bai) not found for {bam_path}.")
        if auto_index_if_missing:
            logger.info(f"Attempting to auto-index BAM file: {bam_path}")
            try:
                run_external_cmd([samtools_path, "index", str(bam_path)])
                if not (bai_path_1.exists() or bai_path_2.exists()): # Re-check
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


# --- Checksum Verification ---
def verify_file_checksum(file_path: Path, expected_checksum: str, algorithm: str = "sha256"):
    logger.info(f"Verifying {algorithm} checksum for {file_path}...")
    ensure_file_exists(file_path, f"File for checksum '{file_path.name}'")
    
    hasher = hashlib.new(algorithm)
    try:
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""): # Read in chunks
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

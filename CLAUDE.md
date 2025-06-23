# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Development Setup
```bash
# Install in development mode
pip install -e .

# Or use conda environment
conda env create -f environment.yml
conda activate basebuddy
```

### Testing
```bash
# Run smoke test
bash tests/smoke_test.sh

# Run integration tests
python tests/test_integration.py

# Run specific test
python -m pytest tests/test_integration.py::test_specific_function
```

### Docker
```bash
# Build Docker image
DOCKER_BUILDKIT=1 docker build --progress=plain -t basebuddy:latest .

# Run Docker container
docker run -v /path/to/data:/data basebuddy:latest <command>
```

## Architecture

BaseBuddy is a bioinformatics CLI/GUI tool for simulating and manipulating sequencing data. It wraps several external tools into a unified interface.

### Core Components

1. **CLI Layer** (`src/basebuddy/cli.py`): Typer-based command interface defining all commands (short, long, spike, signature, strand-bias, qc)

2. **Runner Layer** (`src/basebuddy/runner.py`, `src/basebuddy/bb_runners.py`): Core business logic that orchestrates external tools:
   - `Runner` base class with common functionality
   - Specialized runners for each command type
   - Handles tool execution, file management, and error handling

3. **External Tool Integration**: Wraps bioinformatics tools:
   - ART → short read simulation
   - NanoSim-h → long read simulation  
   - BAMSurgeon → variant spiking
   - SigProfilerSimulator → mutational signatures
   - SAMtools → BAM/FASTA manipulation
   - FastQC → quality control

4. **GUI Layer** (`src/basebuddy/gui/main_app.py`): Optional customtkinter GUI

### Key Design Patterns

- **Exception Hierarchy**: Custom exceptions in `utils.py` (BaseBuddyConfigError, BaseBuddyToolError, BaseBuddyFileError)
- **Automatic File Indexing**: Creates .fai/.bai indices when needed
- **Output Management**: Creates timestamped run directories with manifest files
- **Tool Availability Checking**: Validates external tools before execution

### Important Notes

- The codebase has an issue with missing `Optional` import in `cli.py` that needs fixing
- External tools must be available in PATH (handled automatically in Docker)
- Reference genomes can be cached or provided via FASTA files
- All runners inherit from base `Runner` class for consistent behavior
# External Tool Dependencies

BaseBuddy relies on several external tools that need to be installed and available on your system's `$PATH`. Here's a list of these dependencies and how to install them:

*   **ART (Illumina read simulator):**
    *   Installation: `conda install -c bioconda art_illumina`
    *   Alternatively, download from [ART's website](https://www.niehs.nih.gov/research/resources/software/art/index.cfm) and add to `$PATH`.
*   **NanoSim-h (Nanopore read simulator):**
    *   Installation: `conda install -c bioconda nanosim-h`
    *   Alternatively, download from [NanoSim-h GitHub](https://github.com/bcgsc/NanoSim) and add to `$PATH`.
*   **SAMtools:**
    *   Installation: `conda install -c bioconda samtools`
*   **BAMSurgeon (for spiking in variants):**
    *   Installation: `conda install -c bioconda bamsurgeon`
    *   This provides `addsnv.py` among other tools.
*   **SigProfilerSimulator (for mutational signature simulation):**
    *   Installation: `pip install SigProfilerSimulator`

**Note:** If you are using the provided Docker image, these dependencies are already included.

*   **FastQC (for read quality control):**
    *   Installation: `conda install -c bioconda fastqc`
    *   Website: [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
*   **Perl (required for SimuG):**
    *   Usually pre-installed on Linux/macOS.
    *   Installation on Conda: `conda install -c conda-forge perl`
    *   Windows users might need to install Strawberry Perl or use WSL.
*   **SimuG (for modifying FASTA with germline variants):**
    *   SimuG is a Perl script and needs Perl to run.
    *   Installation:
        ```bash
        git clone https://github.com/yjx1217/simuG.git
        # Add the directory containing simuG.pl to your PATH, or call it with its full path.
        # Example: perl /path/to/simuG/simuG.pl -h
        ```
    *   GitHub: [https://github.com/yjx1217/simuG](https://github.com/yjx1217/simuG)

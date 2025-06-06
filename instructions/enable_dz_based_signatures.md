Comprehensive Feature Plan: Disease-Based Signature Simulation
1. Objective

To make BaseBuddy more intuitive for its target audience—clinical professionals—who think in terms of diseases (e.g., "Lung Cancer") rather than abstract mutational signatures (e.g., "SBS4"). This feature will allow a user to apply a realistic, multi-class (SBS, DBS, ID) mutational landscape to a FASTA file by simply selecting a cancer type.
2. Data Asset Summary

The implementation relies on two distinct categories of bundled data files:

    Canonical Signature Matrices (/data/Signatures/COSMIC...): These files define the fundamental signatures (e.g., SBS1, DBS2). Each file is a matrix where rows are mutation types and columns are signature IDs, containing the probability of each mutation type for that signature. These are the "building blocks."
    Patient Signature Profiles (/data/SP_Signatures_in_Samples/...): These files contain real-world data linking diseases to signature compositions. The columns represent individual patient samples, grouped by Cancer Types. The values are the counts of each canonical signature found in that patient. This data provides the plausible "recipe" for a given disease.

3. Proposed User Workflow

    The user selects a "Simulate by Disease" option in the GUI.
    A dropdown menu appears, populated exclusively with the cancer types that have complete profiles across all three (SBS, DBS, ID) ..._in_samples.csv files.
    The user selects a cancer type (e.g., "ColoRect-AdenoCA").
    The tool applies a plausible spectrum of SBS, DBS, and ID mutations simultaneously to the user's input FASTA file.

4. Technical Implementation Plan

This feature will be implemented by creating a new backend runner that leverages the SigProfilerSimulator Python library.

    Step 1: GUI Enhancement & Dynamic Dropdown
        Task: Update the GUI with the "Simulate by Disease" workflow. The cancer type dropdown menu must be dynamically populated at startup by parsing the available cancer types from the ..._in_samples.csv files.
        Tidbit: To ensure a robust simulation is always possible, the list should be filtered to include only those cancer types present in all three patient profile files (SBS, DBS, and ID). This list generation should be a backend utility function.

    Step 2: Backend Logic for Profile Selection
        Task: When a user selects a cancer type, the backend will randomly select one patient sample (i.e., one column) from the patient profile tables that matches the chosen disease. It will then extract the full signature profile for that patient—a list of signature IDs and their corresponding mutation counts (e.g., {'SBS1': 1496, 'SBS5': 1825, 'DBS2': 91, ...}).

    Step 3: Profile Transformation & Simulation
        Task: The patient profile from Step 2 must be converted into the precise input format required by SigProfilerSimulator.
        Key Tidbit - The Transformation Logic: SigProfilerSimulator is designed to handle this complex step. It takes the canonical signature matrices (the "what") and a weights/counts file (the "how much") as input. Our task is simplified to formatting the randomly selected patient profile into the specific file format the simulator expects for its signatures parameter.

    Step 4: Invoking the Simulator
        Task: Invoke the SigProfilerSimulator Python library, providing it with:
            The user's input FASTA file.
            The paths to our local Canonical Signature Matrix files.
            The Patient Signature Profile (from Step 2), formatted as a text file or pandas DataFrame specifying the number of mutations to generate from each signature.

5. Key Prerequisite: Prototyping the Simulator Input

The research you provided clarifies the path forward. The new prerequisite is to write a small proof-of-concept script that directly uses SigProfilerSimulator with the required inputs:

    A sample FASTA file.
    The paths to our local COSMIC... matrix files.
    A manually created patient profile file (e.g., a simple text file with the format SBS1,1496\nSBS5,1825\n...) that mimics the data from Step 2.

The goal is to confirm the exact function calls and data formats, validating the approach before full integration.
6. Dependencies and Documentation

    Primary Tool: SigProfilerSimulator (Python library).
    Key Documentation:
        OSF Wiki (Comprehensive): https://osf.io/usxjz/wiki/home/
        Source Code (GitHub): https://github.com/AlexandrovLab/SigProfilerSimulator

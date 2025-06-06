import csv
import logging
import random
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

logger = logging.getLogger(__name__)

class SignatureError(Exception):
    """Base class for errors in signature processing."""
    pass

class SignatureFormatError(SignatureError):
    """Error for malformed signature files."""
    pass

# Dictionary mapping signature type prefixes to their master filenames
BUNDLED_SIGNATURE_MASTER_FILES: Dict[str, str] = {
    "sbs": "sbs_grch37_cosmic_v3.3.tsv", # Example, actual filename might differ
    "dbs": "dbs_grch37_cosmic_v3.3.tsv", # Example
    "id": "id_grch37_cosmic_v3.3.tsv",   # Example
}

def parse_signature_matrix_tsv(file_path: Path) -> Dict[str, Dict[str, float]]:
    """
    Parses a tab-separated (TSV) file containing a generic signature matrix.

    The expected format is:
    Header: MutationType SignatureName1 SignatureName2 ...
            (e.g., MutationType SBS1 SBS5a DBS2 ID8 ...)
    Rows:   <MutationTypeID> <Prob_Sig1> <Prob_Sig2> ...
            (e.g., A[C>A]A      0.01         0.005        ...)
            (e.g., AC>TT        0.003        0.001        ...)

    Args:
        file_path: Path to the signature matrix TSV file.

    Returns:
        A dictionary where outer keys are signature names (e.g., "SBS1"),
        and inner keys are mutation types (e.g., "A[C>A]A"), with values
        being the probabilities (float).
        Example: {"SBS1": {"A[C>A]A": 0.01, ...}, "DBS2": {"AC>TT": 0.003, ...}}

    Raises:
        FileNotFoundError: If the file_path does not exist.
        SignatureFormatError: If the file is malformed.
    """
    if not file_path.exists():
        logger.error(f"Signature matrix file not found: {file_path}")
        raise FileNotFoundError(f"Signature matrix file not found: {file_path}")

    logger.info(f"Parsing signature matrix from: {file_path}")
    signatures: Dict[str, Dict[str, float]] = {}

    try:
        with open(file_path, 'r', newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')

            header = next(reader, None)
            if not header or len(header) < 2 or header[0].lower().strip() != "mutationtype":
                raise SignatureFormatError(
                    f"Invalid header in signature file {file_path}. Expected 'MutationType' "
                    f"as the first column (case-insensitive), followed by signature names. Got: {header}"
                )

            signature_names = header[1:]
            for sig_name in signature_names:
                signatures[sig_name.strip()] = {} # Ensure signature names are stripped

            for i, row in enumerate(reader):
                if not row: # Skip empty rows
                    continue
                if len(row) != len(header):
                    raise SignatureFormatError(
                        f"Row {i+2} in {file_path} has an inconsistent number of columns. "
                        f"Expected {len(header)}, got {len(row)}. Row content: {row}"
                    )

                mutation_type = row[0].strip()
                if not mutation_type: # Skip rows with empty mutation type
                    logger.warning(f"Skipping row {i+2} in {file_path} due to empty mutation type.")
                    continue

                for j, sig_name in enumerate(signature_names):
                    try:
                        probability_str = row[j+1].strip()
                        if not probability_str: # Handle empty probability strings if necessary
                            logger.warning(f"Empty probability for {sig_name} at {mutation_type}, row {i+2}. Assuming 0.0 or error.")
                            probability = 0.0 # Or raise error/skip
                        else:
                            probability = float(probability_str)

                        signatures[sig_name.strip()][mutation_type] = probability
                    except ValueError:
                        raise SignatureFormatError(
                            f"Non-numeric probability found in {file_path} for signature '{sig_name.strip()}', "
                            f"mutation type '{mutation_type}', value '{row[j+1]}' at row {i+2}, column {j+2}."
                        )

            if not signatures:
                 raise SignatureFormatError(f"No signature data found in {file_path} beyond the header.")
            for sig_name, mut_probs in signatures.items():
                if not mut_probs:
                     raise SignatureFormatError(f"Signature '{sig_name}' in {file_path} contains no mutation probabilities.")

    except StopIteration:
        raise SignatureFormatError(f"Signature file {file_path} is empty or contains only a header.")
    except csv.Error as e:
        raise SignatureFormatError(f"CSV parsing error in {file_path}: {e}")
    except FileNotFoundError:
        logger.error(f"Attempted to parse non-existent file: {file_path}")
        raise
    except Exception as e:
        logger.exception(f"An unexpected error occurred while parsing signature file {file_path}: {e}")
        raise SignatureError(f"An unexpected error occurred parsing {file_path}: {str(e)}")

    logger.info(f"Successfully parsed {len(signatures)} signature(s) from {file_path}: {list(signatures.keys())}")
    return signatures

def get_bundled_signature_master_file_path(sig_type_prefix: str, bundled_signatures_dir: Path) -> Optional[Path]:
    """
    Constructs the path to a bundled master signature TSV file based on type prefix.

    Args:
        sig_type_prefix: The prefix of the signature type (e.g., "sbs", "dbs", "id"). Case-insensitive.
        bundled_signatures_dir: The directory where bundled signature master files are stored.

    Returns:
        A Path object to the master signature file if a corresponding entry exists,
        otherwise None.
    """
    if not sig_type_prefix:
        return None

    master_filename = BUNDLED_SIGNATURE_MASTER_FILES.get(sig_type_prefix.lower())
    if not master_filename:
        logger.warning(f"No bundled master file defined for signature type prefix: '{sig_type_prefix.lower()}'")
        return None

    file_path = bundled_signatures_dir / master_filename
    # This function only constructs the path; existence should be checked by the caller.
    return file_path


# --- SBS-Specific Mutation Application Functions ---

def get_sbs_mutation_contexts(sequence: str) -> List[Tuple[int, str, str]]:
    """
    Finds all possible 3-base contexts for SBS mutations in a DNA sequence.
    (Implementation as provided previously)
    """
    contexts: List[Tuple[int, str, str]] = []
    seq_upper = sequence.upper()
    valid_bases = "ACGT"

    if len(seq_upper) < 3:
        return contexts

    for i in range(len(seq_upper) - 2):
        trinucleotide = seq_upper[i : i + 3]
        if all(base in valid_bases for base in trinucleotide):
            central_base_pos = i + 1
            original_central_base = trinucleotide[1]
            contexts.append((central_base_pos, trinucleotide, original_central_base))
    return contexts

SBS_MUTATION_TYPE_PATTERN = re.compile(r"([ACGTN])\[([ACGTN])>([ACGTN])\]([ACGTN])") # Remains SBS specific

def _parse_mutation_type_string(mutation_type_str: str) -> Tuple[str, str, str]: # Remains SBS specific for now
    """
    Parses a standard SBS mutation type string like "A[C>T]G".
    (Implementation as provided previously)
    """
    match = SBS_MUTATION_TYPE_PATTERN.fullmatch(mutation_type_str)
    if not match:
        raise SignatureFormatError(f"Invalid SBS mutation type string format: '{mutation_type_str}'. Expected format like 'N[N>N]N'.")

    prefix, original_base, new_base, suffix = match.groups()
    trinucleotide_context = f"{prefix}{original_base}{suffix}"
    return trinucleotide_context, original_base, new_base

def choose_sbs_mutation(signature_profile: Dict[str, float]) -> str: # Name implies SBS, but logic is generic for weighted choice
    """
    Stochastically selects an SBS mutation type based on its probability in the signature.
    (Implementation as provided previously)
    """
    if not signature_profile:
        logger.warning("Cannot choose mutation from empty signature profile.")
        return ""

    mutation_types = list(signature_profile.keys())
    probabilities = list(signature_profile.values())

    if not mutation_types or not probabilities or len(mutation_types) != len(probabilities):
        logger.warning("Invalid signature profile for choosing mutation (empty or mismatched keys/values).")
        return ""

    if any(p < 0 for p in probabilities) or sum(probabilities) <= 0:
        logger.warning(f"Probabilities in signature profile are invalid (negative or sum to zero/negative): {probabilities}")
        if sum(probabilities) <= 0:
            logger.warning("Sum of probabilities is <= 0, cannot make a weighted choice.")
            # If all probabilities are zero, could select uniformly:
            # return random.choice(mutation_types) if mutation_types else ""
            return "" # Or raise error

    try:
        chosen_list = random.choices(mutation_types, weights=probabilities, k=1)
        return chosen_list[0]
    except ValueError as e:
        logger.error(f"Error in random.choices, possibly due to invalid weights (e.g., all zero): {e}. Profile: {signature_profile}")
        return ""

def apply_sbs_mutations_to_sequence(
    sequence: str,
    num_mutations: int,
    signature_profile: Dict[str, float],
    available_contexts: List[Tuple[int, str, str]]
) -> Tuple[str, int]:
    """
    Applies a specified number of SBS mutations to a sequence according to a signature.
    (Implementation as provided previously, still SBS specific due to _parse_mutation_type_string and context matching logic)
    """
    if not sequence or num_mutations == 0 or not signature_profile :
        return sequence, 0

    if not available_contexts:
        logger.warning("No available contexts to apply mutations.")
        return sequence, 0

    sequence_list = list(sequence)
    applied_count = 0

    for _ in range(num_mutations):
        if not available_contexts:
            logger.warning("Ran out of available contexts during mutation application.")
            break

        chosen_mutation_type_str = choose_sbs_mutation(signature_profile)
        if not chosen_mutation_type_str:
            logger.warning("choose_sbs_mutation returned no choice, skipping this mutation attempt.")
            continue

        try:
            required_trinucleotide_context, original_central_base, new_central_base = \
                _parse_mutation_type_string(chosen_mutation_type_str) # SBS specific parsing
        except SignatureFormatError as e:
            logger.warning(f"Could not parse chosen mutation type '{chosen_mutation_type_str}': {e}. Skipping mutation.")
            continue

        candidate_sites = [
            site for site in available_contexts
            if site[1] == required_trinucleotide_context and \
               site[2] == original_central_base and \
               sequence_list[site[0]] == original_central_base
        ]

        if not candidate_sites:
            continue

        chosen_site_context_tuple = random.choice(candidate_sites)
        position_to_mutate_0_based = chosen_site_context_tuple[0]

        sequence_list[position_to_mutate_0_based] = new_central_base
        applied_count += 1

    return "".join(sequence_list), applied_count

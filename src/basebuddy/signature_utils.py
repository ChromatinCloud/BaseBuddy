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

def parse_sbs_signature_matrix(file_path: Path) -> Dict[str, Dict[str, float]]:
    """
    Parses a tab-separated (TSV) file containing an SBS signature matrix.
    (Implementation as provided previously)
    """
    if not file_path.exists():
        logger.error(f"Signature matrix file not found: {file_path}")
        raise FileNotFoundError(f"Signature matrix file not found: {file_path}")

    logger.info(f"Parsing SBS signature matrix from: {file_path}")
    signatures: Dict[str, Dict[str, float]] = {}

    try:
        with open(file_path, 'r', newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')

            header = next(reader, None)
            if not header or len(header) < 2 or header[0].lower() != "mutationtype":
                raise SignatureFormatError(
                    f"Invalid header in signature file {file_path}. Expected 'MutationType' "
                    f"as the first column, followed by signature names. Got: {header}"
                )

            signature_names = header[1:]
            for sig_name in signature_names:
                signatures[sig_name] = {}

            for i, row in enumerate(reader):
                if len(row) != len(header):
                    raise SignatureFormatError(
                        f"Row {i+2} in {file_path} has an inconsistent number of columns. "
                        f"Expected {len(header)}, got {len(row)}. Row content: {row}"
                    )

                mutation_type = row[0]
                for j, sig_name in enumerate(signature_names):
                    try:
                        probability = float(row[j+1])
                        signatures[sig_name][mutation_type] = probability
                    except ValueError:
                        raise SignatureFormatError(
                            f"Non-numeric probability found in {file_path} for signature '{sig_name}', "
                            f"mutation type '{mutation_type}', value '{row[j+1]}' at row {i+2}, column {j+2}."
                        )

            if not signatures: # File only had a header
                 raise SignatureFormatError(f"No signature data found in {file_path} beyond the header.")
            for sig_name, mut_probs in signatures.items():
                if not mut_probs: # Signature column had no data rows
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

def get_bundled_sbs_signature_path(signature_id: str, bundled_signatures_dir: Path) -> Optional[Path]:
    """
    Constructs the path to a bundled SBS signature file.
    (Implementation as provided previously)
    """
    if not signature_id:
        return None
    filename = f"{signature_id.lower()}.tsv"
    file_path = bundled_signatures_dir / filename
    return file_path

# --- New SBS Mutation Application Functions ---

def get_sbs_mutation_contexts(sequence: str) -> List[Tuple[int, str, str]]:
    """
    Finds all possible 3-base contexts for SBS mutations in a DNA sequence.

    Args:
        sequence: A single DNA sequence (e.g., for one chromosome).

    Returns:
        A list of tuples: (position_0_based_py, trinucleotide, original_central_base).
        - position_0_based_py: 0-based index of the central base in the input sequence.
        - trinucleotide: The 3-base string (e.g., "ACA").
        - original_central_base: The central base of the trinucleotide (e.g., 'C').
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

SBS_MUTATION_TYPE_PATTERN = re.compile(r"([ACGTN])\[([ACGTN])>([ACGTN])\]([ACGTN])")

def _parse_mutation_type_string(mutation_type_str: str) -> Tuple[str, str, str]:
    """
    Parses a standard SBS mutation type string like "A[C>T]G".

    Args:
        mutation_type_str: The string to parse.

    Returns:
        A tuple: (trinucleotide_context, original_central_base, new_central_base).
        Example: For "A[C>T]G", returns ("ACG", "C", "T").

    Raises:
        SignatureFormatError: If the string format is invalid.
    """
    match = SBS_MUTATION_TYPE_PATTERN.fullmatch(mutation_type_str)
    if not match:
        raise SignatureFormatError(f"Invalid SBS mutation type string format: '{mutation_type_str}'. Expected format like 'N[N>N]N'.")

    prefix, original_base, new_base, suffix = match.groups()
    trinucleotide_context = f"{prefix}{original_base}{suffix}"
    return trinucleotide_context, original_base, new_base

def choose_sbs_mutation(signature_profile: Dict[str, float]) -> str:
    """
    Stochastically selects an SBS mutation type based on its probability in the signature.

    Args:
        signature_profile: A dictionary where keys are mutation types (e.g., "A[C>A]A")
                           and values are their probabilities.

    Returns:
        The chosen mutation type string (e.g., "A[C>A]A").
        Returns an empty string if the profile is empty or probabilities are invalid.
    """
    if not signature_profile:
        logger.warning("Cannot choose mutation from empty signature profile.")
        return "" # Or raise error

    mutation_types = list(signature_profile.keys())
    probabilities = list(signature_profile.values())

    if not mutation_types or not probabilities or len(mutation_types) != len(probabilities):
        logger.warning("Invalid signature profile for choosing mutation (empty or mismatched keys/values).")
        return "" # Or raise error

    # Ensure probabilities are non-negative and sum to a positive value
    if any(p < 0 for p in probabilities) or sum(probabilities) <= 0:
        # Attempt to normalize if sum is not ~1.0 but all are positive, or handle as error
        # For now, if sum is zero or negative, it's problematic.
        logger.warning(f"Probabilities in signature profile are invalid (negative or sum to zero/negative): {probabilities}")
        # Fallback: could pick uniformly if all weights are zero, or raise error.
        # For simplicity, if probabilities are messed up, we might not be able to make a weighted choice.
        # random.choices will raise ValueError if sum of weights is 0.
        # Let it raise, or handle more gracefully:
        if sum(probabilities) <= 0:
            # If all probabilities are zero, could select uniformly or return error/empty.
            # For now, let random.choices handle it (it will likely error).
            # A simple uniform choice if all are zero: return random.choice(mutation_types) if mutation_types else ""
            logger.warning("Sum of probabilities is <= 0, cannot make a weighted choice.")
            return ""


    try:
        chosen_list = random.choices(mutation_types, weights=probabilities, k=1)
        return chosen_list[0]
    except ValueError as e: # Handles issues like sum of weights being zero
        logger.error(f"Error in random.choices, possibly due to invalid weights (e.g., all zero): {e}. Profile: {signature_profile}")
        return "" # Fallback or re-raise

def apply_sbs_mutations_to_sequence(
    sequence: str,
    num_mutations: int,
    signature_profile: Dict[str, float],
    available_contexts: List[Tuple[int, str, str]]
) -> Tuple[str, int]:
    """
    Applies a specified number of SBS mutations to a sequence according to a signature.

    Args:
        sequence: The original DNA sequence.
        num_mutations: The total number of mutations to attempt.
        signature_profile: Probabilities for each mutation type for the desired signature.
        available_contexts: Pre-calculated list of mutable sites from
                            get_sbs_mutation_contexts(original_sequence).

    Returns:
        A tuple: (mutated_sequence_str, applied_mutations_count).
        The available_contexts list is NOT dynamically updated after each mutation in this version.
    """
    if not sequence or num_mutations == 0 or not signature_profile :
        return sequence, 0

    if not available_contexts: # No valid sites on the original sequence
        logger.warning("No available contexts to apply mutations.")
        return sequence, 0

    sequence_list = list(sequence)
    applied_count = 0

    for _ in range(num_mutations):
        if not available_contexts: # Should ideally not be reached if initial check is done, but as safeguard
            logger.warning("Ran out of available contexts during mutation application.")
            break

        chosen_mutation_type_str = choose_sbs_mutation(signature_profile)
        if not chosen_mutation_type_str:
            logger.warning("choose_sbs_mutation returned no choice, skipping this mutation attempt.")
            continue # Skip this attempt if no mutation type could be chosen

        try:
            # Example: chosen_mutation_type_str = "A[C>T]G"
            # required_trinucleotide_context = "ACG"
            # original_central_base = "C"
            # new_central_base = "T"
            required_trinucleotide_context, original_central_base, new_central_base = \
                _parse_mutation_type_string(chosen_mutation_type_str)
        except SignatureFormatError as e:
            logger.warning(f"Could not parse chosen mutation type '{chosen_mutation_type_str}': {e}. Skipping mutation.")
            continue

        # Find sites in the *original* sequence that match this chosen mutation type's context
        # Note: available_contexts stores (0_based_pos, trinuc_str, original_center_base_in_trinuc)
        candidate_sites = [
            site for site in available_contexts
            if site[1] == required_trinucleotide_context and \
               site[2] == original_central_base and \
               sequence_list[site[0]] == original_central_base # Ensure current base at site is still the original one
        ]

        if not candidate_sites:
            # logger.debug(f"No available sites match the chosen mutation type context: {chosen_mutation_type_str}")
            continue # Try choosing another mutation type in the next iteration

        chosen_site_context_tuple = random.choice(candidate_sites)
        position_to_mutate_0_based = chosen_site_context_tuple[0]

        # Apply mutation
        # logger.debug(f"Applying mutation: {chosen_mutation_type_str} at pos {position_to_mutate_0_based} (original: {sequence_list[position_to_mutate_0_based]}, new: {new_central_base})")
        sequence_list[position_to_mutate_0_based] = new_central_base
        applied_count += 1

        # As per subtask, available_contexts is not dynamically updated here.
        # This means a site could theoretically be chosen again if the new base still allows for
        # a valid context according to the *original* available_contexts list.
        # However, the check `sequence_list[site[0]] == original_central_base` helps prevent
        # re-mutating an already mutated site *if the original base it expects has changed*.

    return "".join(sequence_list), applied_count

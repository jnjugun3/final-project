#!/usr/bin/env python3

import sys
from typing import Dict

'''
Step 1. Reads and validates a nucleotide FASTA file.

- Ensures file is not empty
- Ensures headers begin with ">"
- Ensures sequences contain only valid bases (A, T, C, G)
- Ensures no duplicate sequence IDs
- Returns a dictionary {seq_id : sequence}
'''

# Allowed nucleotide bases 
VALID_BASES = set("ATCGNRYKMSWBDHV") 


def read_fasta(filepath: str) -> Dict[str, str]:
    """
    Read and validate a FASTA file.

    Parameters
    ----------
    filepath : str
        Path to FASTA file.

    Returns
    -------
    Dict[str, str]
        Dictionary {seq_id : validated nucleotide sequence}.
    """

    # Read all lines from FASTA file
    with open(filepath, "r") as fh:
        lines = fh.readlines()
        
    # Raise error if FASTA file contains no lines
    if not lines:
        raise ValueError("FASTA file is empty.")

    # Storage for sequences and tracking 
    seq_dict = {}
    current_header = None
    current_seq = []
    header_found = False

    # Reads line by line, skips if empty, and removes whitespace
    for line in lines:
        line = line.strip()
        if not line:
            continue

        # Looks for header lines
        if line.startswith(">"):
            header_found = True

            # Save previous sequence if one is being built
            if current_header and current_seq:
                seq_dict[current_header] = "".join(current_seq).upper()

            # Extract header 
            header = line[1:].strip()

            # Checks for duplicate IDs
            if header in seq_dict:
                raise ValueError(f"Duplicate FASTA header found: '{header}'")

            # Starts recording a new sequence 
            current_header = header
            current_seq = []
            continue

        # If we find a sequence before any header, raise an error
        if current_header is None:
            raise ValueError("Found sequence before FASTA header ('>').")

        # Remove spaces and convert sequence to uppercase
        seq_line = line.replace(" ", "").upper()

        # Validate each nucleotide base
        fixed_line = []  
        for base in seq_line:

           
            if base not in VALID_BASES:
                print(f"[WARNING] Nonstandard nucleotide '{base}' in '{current_header}'. Converting to 'N'.")
                base = "N"

            fixed_line.append(base)

        # Append validated base line into a list 
        current_seq.append("".join(fixed_line))

    # Save last sequence after the loop ends
    if current_header and current_seq:
        seq_dict[current_header] = "".join(current_seq).upper()

    # Make sure that at least 1 header is present
    if not header_found:
        raise ValueError("No FASTA headers ('>') found.")

    # Make sure at least one valid sequence exists
    if not seq_dict:
        raise ValueError("No valid sequences found in FASTA.")

    return seq_dict

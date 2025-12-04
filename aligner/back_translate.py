#!/usr/bin/env python3

"""
 Codon-Aware Back-Translation

Converts an amino-acid multiple sequence alignment (MSA) produced in Step 5
back into a codon-aware nucleotide alignment.

Inputs
------
aligned_aa : dict[str, str]
    { seq_id : aligned amino-acid string }
    Example: { "s1": "M-AK--T*", "s2": "MTAK--T-" }

longest_orfs : dict[str, str]
    { seq_id : nucleotide ORF string (ungapped) }
    Example: { "s1": "ATGGCTAAA..." }

Rules
-----
- Each non-gap amino acid corresponds to the next codon (3 nt) in the ORF.
- Alignment gaps '-' become codon-sized gaps: '---'.
- Optional stop '*' at the end of AA sequence:
      - consumes one codon
      - you can keep or drop the stop codon depending on downstream needs.

Output
------
codon_aligned : dict[str, str]
    { seq_id : codon-aware nucleotide alignment string }
    Example: { "s1": "ATG---GCT...", "s2": "ATGAAA---..." }

This nucleotide MSA can then be used in Step 7 for codon-position statistics.
"""

from typing import Dict, List


def back_translate(
    aligned_aa,            # mapping: seq_id → aligned amino-acid string
    longest_orfs,          # mapping: seq_id → original ungapped ORF nucleotides
    drop_terminal_stop=True
):  
    # type: (...) -> Dict[str, str]

    codon_aligned = {}     # dictionary to hold final codon-aware nucleotide alignments

    # Iterate through each aligned amino-acid sequence
    for seq_id, aa_aln in aligned_aa.items():

        # Ensure ORF exists for this sequence
        if seq_id not in longest_orfs:
            raise KeyError(
                "Sequence '%s' is present in aligned amino acids but missing from longest_orfs."
                % seq_id
            )

        # Fetch original nucleotide ORF for this sequence
        nt_seq = longest_orfs[seq_id]    
        nt_len = len(nt_seq)             # should always be divisible by 3

        # Validate ORF triplet structure
        if nt_len % 3 != 0:
            raise ValueError(
                "ORF for '%s' has length %d, which is not a multiple of 3."
                % (seq_id, nt_len)
            )

        pos = 0                # pointer to next codon in nt_seq
        nt_out = []            # accumulator list for codon-aligned NT output

        # NOTE:
        # We ALWAYS iterate over the FULL aligned AA string.
        # This fixes the length-mismatch bug when '*' was sliced off.

        for aa in aa_aln:

            if aa == "-":
                nt_out.append("---")   # gap in AA alignment → codon-sized gap

            elif aa == "*":
                # Stop codon position is ALWAYS represented in the alignment
                if drop_terminal_stop:
                    # Keep column but represent '*' as gap codon
                    nt_out.append("---")
                else:
                    # Use actual stop codon from ORF (if available)
                    if pos + 3 <= nt_len:
                        stop_codon = nt_seq[pos : pos + 3]
                        nt_out.append(stop_codon)
                        pos += 3
                    else:
                        nt_out.append("---")

            else:
                # Ensure enough nucleotides remain for a full codon
                if pos + 3 <= nt_len:
                    codon = nt_seq[pos : pos + 3]  # extract next codon
                    nt_out.append(codon)
                    pos += 3                       # advance ORF pointer
                else:
                    nt_out.append("---")            # AA alignment longer than ORF → codon gap

        # Join list of codons into a single aligned NT string
        codon_aligned[seq_id] = "".join(nt_out)

    return codon_aligned

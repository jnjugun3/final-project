#!/usr/bin/env python3

'''
Step 3. Translate nucleotide ORFs into amino-acid sequences.

- Input:  longest_orfs  {seq_id : nucleotide_ORF_string}
- Uses Biopython Seq.translate()
- User may choose preferred translation table (default table = 1)
- Keeps '*' stop codon if present (back-translation can handle it)
- Returns amino-acid dictionary for k-mer similarity step

Common translation tables used in practice:
    1  = Standard Genetic Code
    2  = Vertebrate Mitochondrial Code
    4  = Mold / Protozoan Mitochondrial Code
    11 = Bacterial / Archaeal / Plant Plastid Code
'''

from Bio.Seq import Seq
from typing import Dict


def translate_orfs(longest_orfs: Dict[str, str], table: int = 1) -> Dict[str, str]:
    """
    Translate each nucleotide ORF into an amino-acid sequence
    using the specified NCBI translation table.

    Returns:
        aa_dict = { seq_id : amino_acid_string }
    """

    aa_dict = {}   # storage for translated amino-acid sequences

    # Loop through each sequence ID and its nucleotide ORF
    for seq_id, nt_orf in longest_orfs.items():

        # Convert nucleotide string into a Biopython Seq object
        seq = Seq(nt_orf)

        # Translate the ORF using the user-selected translation table
        aa = seq.translate(
            table=table,   # NCBI translation table ID
            to_stop=False  # do not stop at stop codon; include '*' if present
        )

        # Keep '*' if present — handled later during back-translation
        aa_clean = str(aa)

        # Store the translated sequence in the dictionary
        aa_dict[seq_id] = aa_clean

    # Return all amino-acid sequences for downstream k-mer similarity step
    return aa_dict

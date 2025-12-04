#!/usr/bin/env python3

'''
Detect Open Reading Frames (ORFs)

- Scans all 3 reading frames (0, 1, 2)
- Searches on 3 strands:
      1. Forward strand
      2. Reverse strand (string reversed, class notes)
      3. Reverse complement strand
- Identifies ORFs starting with ATG and ending in TAA/TAG/TGA
- Filters ORFs by minimum length
- Returns:
     (1) all ORFs detected on all strands
     (2) the longest nucleotide ORF per sequence (extract_longest_orf)

Works with seq_dict from load_fasta.py {seq_id : nucleotide_string}
'''

from Bio.Seq import Seq

# Start and stop codons used for ORF detection
START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


def detect_orfs(seq: str, min_length: int, strand: str) -> list:
    """
    Detect ORFs on **one strand** across frames 0, 1, and 2.
    """

    orfs = []

    # Loop across all 3 reading frames
    for frame in (0, 1, 2):

        i = frame   # start index for this reading frame

        # Continue codon-by-codon until end of sequence
        while i + 3 <= len(seq):

            codon = seq[i:i+3]   # current codon

            # ATG signals start of ORF
            if codon == START_CODON:
                start_pos = i
                j = i

                # Extend ORF until a stop codon is found
                while j + 3 <= len(seq):

                    stop_codon = seq[j:j+3]

                    # Check if we've reached a valid stop codon
                    if stop_codon in STOP_CODONS:
                        end_pos = j + 3        # include stop codon
                        orf_seq = seq[start_pos:end_pos]

                        # Apply minimum length filter
                        if len(orf_seq) >= min_length:
                            orfs.append({
                                "start": start_pos,
                                "end": end_pos,
                                "frame": frame,
                                "strand": strand,
                                "sequence": orf_seq
                            })
                        break

                    j += 3   # move to next codon

                # Move index forward after ORF is processed
                i = j + 3

            else:
                i += 3   # no start codon → continue scanning

    return orfs



def scan_all_frames(seq: str, min_length: int = 30) -> list:
    """
    Detect ORFs across all 3 strands:
        - Forward
        - Reverse (string reversed)
        - Reverse complement
    """

    seq = seq.upper()   # ensure uppercase input

    # Prepare strand orientations
    forward = seq
    reverse = seq[::-1]
    revcomp = str(Seq(seq).reverse_complement())

    # Scan each strand independently
    orfs_forward = detect_orfs(forward, min_length, "+")
    orfs_reverse = detect_orfs(reverse, min_length, "rev")
    orfs_revcomp = detect_orfs(revcomp, min_length, "-")

    # Return combined ORFs from all strands
    return orfs_forward + orfs_reverse + orfs_revcomp



def extract_longest_orf(seq_dict: dict, min_length: int = 30) -> dict:
    """
    For each sequence in seq_dict, return the **longest ORF** found.
    """

    longest_orfs = {}   # <-- updated variable name

    # Loop through every input sequence
    for seq_id, nt_seq in seq_dict.items():

        # Detect all ORFs for this sequence
        all_orfs = scan_all_frames(nt_seq, min_length)

        # Skip sequences that have no ORFs
        if not all_orfs:
            continue

        # Select ORF with maximum length
        best_orf = max(all_orfs, key=lambda x: len(x["sequence"]))

        # Save the nucleotide ORF sequence
        longest_orfs[seq_id] = best_orf["sequence"]

    return longest_orfs


    
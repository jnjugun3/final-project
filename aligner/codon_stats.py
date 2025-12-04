#!/usr/bin/env python3
"""
Step 7. Codon Statistics 
-------------------------

Computes codon-level statistics from a codon-aware, gap-preserving
nucleotide alignment.

This script provides:

    • Percent identity at nucleotide positions 1, 2, 3
    • Codon-level percent identity (avg of the three)
    • Global percent identity across ALL aligned columns
    • GC fraction at positions 1, 2, 3
    • Alignment-wide summary (clean, sys.stdout-based)
    • CSV output for plotting in Step 8

This file is FAIR-compliant, PEP-8 clean, and fully modular.
"""

import csv
import sys
from typing import Dict, List




def compute_percent_identity(chars: List[str]) -> float:
    """
    Compute percent identity for a list of nucleotides in a single
    position across sequences.
    Only A, T, C, G are considered valid.
    """
    filtered = [c.upper() for c in chars if c.upper() in "ATCG"]
    if not filtered:
        return 0.0
    most_common = max(set(filtered), key=filtered.count)
    return filtered.count(most_common) / len(filtered)


def compute_codon_identity_total(pos1: List[str], pos2: List[str], pos3: List[str]) -> float:
    """Average percent identity across codon positions 1, 2, and 3."""
    return (
        compute_percent_identity(pos1) +
        compute_percent_identity(pos2) +
        compute_percent_identity(pos3)
    ) / 3.0


def compute_global_identity(sequences: List[str]) -> float:
    """
    Compute alignment-wide percent identity across ALL nucleotide columns.
    """
    aln_len = len(sequences[0])
    total = 0.0

    for col in range(aln_len):
        chars = [seq[col].upper() for seq in sequences if seq[col].upper() in "ATCG"]
        if not chars:
            continue
        most_common = max(set(chars), key=chars.count)
        total += chars.count(most_common) / len(chars)

    return total / aln_len



#    Codon Statistics Function

def step7_codon_statistics(
    aligned_nt_dict: Dict[str, str],
    output_csv: str = "codon_stats.csv"
) -> None:
    """
    Compute codon-level percent identity + GC content for every codon
    in a codon-aware nucleotide alignment.

    Outputs:
        • CSV file with full stats
        • Clean, robust terminal summary using sys.stdout.write()
    """

    seq_ids = list(aligned_nt_dict.keys())
    sequences = list(aligned_nt_dict.values())

    # Validate equal lengths
    aln_len = len(sequences[0])
    if any(len(seq) != aln_len for seq in sequences):
        raise ValueError("Aligned sequences must all have identical length.")

    # Validate codon alignment
    if aln_len % 3 != 0:
        raise ValueError("Aligned sequences must be codon-aligned (length % 3 == 0).")

    num_codons = aln_len // 3
    rows = []

    # Global percent identity across ALL nucleotide columns
    global_identity = compute_global_identity(sequences)

    # Running sums for alignment-wide position averages
    sum_id1 = sum_id2 = sum_id3 = 0.0
    sum_id_total = 0.0
    sum_gc1 = sum_gc2 = sum_gc3 = 0.0

    # Per-codon stats
    for codon_index in range(num_codons):

        i = codon_index * 3
        pos1 = [seq[i] for seq in sequences]
        pos2 = [seq[i + 1] for seq in sequences]
        pos3 = [seq[i + 2] for seq in sequences]

        id1 = compute_percent_identity(pos1)
        id2 = compute_percent_identity(pos2)
        id3 = compute_percent_identity(pos3)
        id_total = compute_codon_identity_total(pos1, pos2, pos3)

        gc1 = (pos1.count("G") + pos1.count("C")) / len(pos1)
        gc2 = (pos2.count("G") + pos2.count("C")) / len(pos2)
        gc3 = (pos3.count("G") + pos3.count("C")) / len(pos3)

        # accumulate for summary
        sum_id1 += id1
        sum_id2 += id2
        sum_id3 += id3
        sum_id_total += id_total
        sum_gc1 += gc1
        sum_gc2 += gc2
        sum_gc3 += gc3

        rows.append([
            codon_index + 1,
            id1, id2, id3,
            id_total,
            gc1, gc2, gc3,
            global_identity,
        ])

    
    #   Save CSV
    
    header = [
        "codon_index",
        "identity_pos1", "identity_pos2", "identity_pos3",
        "identity_codon_total",
        "gc_pos1", "gc_pos2", "gc_pos3",
        "identity_global"
    ]

    with open(output_csv, "w", newline="") as fh:
        csv.writer(fh).writerows([header] + rows)

    
    #   Terminal Summary  
     

    sys.stdout.write("\n[Codon Statistics Summary]\n")

    sys.stdout.write(f"Total codons analyzed: {num_codons}\n")
    sys.stdout.write(f"Global percent identity (all nucleotide columns): {global_identity:.4f}\n\n")

    sys.stdout.write("Average Percent Identity (across ALL codons):\n")
    sys.stdout.write(f"  • Position 1 identity:      {sum_id1 / num_codons:.4f}\n")
    sys.stdout.write(f"  • Position 2 identity:      {sum_id2 / num_codons:.4f}\n")
    sys.stdout.write(f"  • Position 3 identity:      {sum_id3 / num_codons:.4f}\n")
    sys.stdout.write(f"  • Combined codon identity:  {sum_id_total / num_codons:.4f}\n\n")

    sys.stdout.write("Average GC Fraction (across ALL codons):\n")
    sys.stdout.write(f"  • GC at position 1: {sum_gc1 / num_codons:.4f}\n")
    sys.stdout.write(f"  • GC at position 2: {sum_gc2 / num_codons:.4f}\n")
    sys.stdout.write(f"  • GC at position 3: {sum_gc3 / num_codons:.4f}\n\n")

    sys.stdout.write(f"[Step 7] Codon statistics saved → {output_csv}\n")
    sys.stdout.write("==========================================================\n\n")

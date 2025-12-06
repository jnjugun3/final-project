#!/usr/bin/env python3
"""
Step 5. Progressive Needleman–Wunsch Alignment

This module performs progressive amino-acid multiple sequence alignment (MSA)
using a right-branching guide tree.

For each internal node:
    - Align its two children with Needleman–Wunsch.
    - Store the aligned children.
    - Build a consensus node.
    - Insert any new gaps into all previously-aligned nodes.

The result is a set of aligned amino-acid sequences (strings) ready for
codon-aware back-translation.
"""

from typing import Dict, List, Any
import numpy as np
from newu3.num import init_mat, fill_matrix, trace_matrix



def array_to_string(arr: np.ndarray) -> str:
    """Convert array-of-sets (e.g., [{'A'}, {'-'}]) into a simple string."""
    return "".join(next(iter(x)) for x in arr)


def build_consensus(seq1: np.ndarray, seq2: np.ndarray) -> np.ndarray:
    """
    Build a consensus sequence for an internal tree node.
    Each column is the union of residues from both children.
    """
    return np.array([a.union(b) for a, b in zip(seq1, seq2)], dtype=object)


 
# Gap patching


def patch_gaps(old_seq: np.ndarray, aligned_seq: np.ndarray) -> np.ndarray:
    """
    Insert gaps into a previously aligned sequence so that its length matches
    the newest alignment.

    Steps:
        1. Find where the new alignment has gaps.
        2. Insert gaps in the old sequence at the same positions.
        3. Pad or trim so final lengths match exactly.

    Parameters
    ----------
    old_seq : np.ndarray
        Previously aligned sequence.
    aligned_seq : np.ndarray
        Newly aligned sequence (defines final column structure).

    Returns
    -------
    np.ndarray
        Updated old_seq, equal in length to aligned_seq.
    """
    # Gap positions from the new alignment
    gap_positions = sorted(
        [i for i, col in enumerate(aligned_seq) if col == {'-'}],
        reverse=True
    )

    seq = old_seq.copy()

    # Insert new gaps right→left to avoid shifting
    for pos in gap_positions:
        if pos >= len(seq):
            seq = np.append(seq, {'-'})
        else:
            seq = np.insert(seq, pos, {'-'})

    # Ensure exact length match
    target_len = len(aligned_seq)

    if len(seq) < target_len:
        seq = np.append(seq, [{'-'}] * (target_len - len(seq)))
    elif len(seq) > target_len:
        seq = seq[:target_len]

    return seq


 
# Progressive NW alignment
 

def progressive_nw(
    sequences: Dict[Any, np.ndarray],
    terminals_dic: Dict[int, str],
    internal_nodes: Dict[int, List[Any]],
    *,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -1,
) -> Dict[str, str]:
    """
    Progressive Needleman–Wunsch alignment.

    Parameters
    ----------
    sequences : dict
        Node ID → sequence (array-of-sets).
    terminals_dic : dict
        Terminal node index → sequence name.
    internal_nodes : dict
        Internal node → [childA, childB].
    match, mismatch, gap : int
        NW scoring parameters.

    Returns
    -------
    dict[str, str]
        Final aligned amino-acid sequences (strings) for terminal IDs.
    """

    alignment_history: List[Any] = []

    for node_id in sorted(internal_nodes.keys()):
        childA, childB = internal_nodes[node_id]

        # Retrieve sequences
        seqA = sequences[childA]
        seqB = sequences[childB]

        # NW DP matrix
        mat = init_mat(seqA, seqB)
        mat = fill_matrix(
            matrix=mat,
            seq1_array=seqA,
            seq2_array=seqB,
            match=match,
            mismatch=mismatch,
            indel=gap
        )

        # Traceback
        alignedA, alignedB = trace_matrix(
            matrix=mat,
            seq1_array=seqA,
            seq2_array=seqB,
            match=match,
            mismatch=mismatch,
            indel=gap
        )

        # Save aligned children
        sequences[childA] = alignedA
        sequences[childB] = alignedB

        # Make a consensus internal node
        sequences[node_id] = build_consensus(alignedA, alignedB)

        # Patch gaps into all previous sequences
        for prev_node in alignment_history:
            sequences[prev_node] = patch_gaps(
                old_seq=sequences[prev_node],
                aligned_seq=alignedA
            )

        # Add nodes to alignment history
        for x in (childA, childB, node_id):
            if x not in alignment_history:
                alignment_history.append(x)

        # Length check
        final_len = len(alignedA)
        for nid in alignment_history:
            if len(sequences[nid]) != final_len:
                raise RuntimeError(
                    f"Inconsistent alignment lengths after node {node_id}: "
                    f"{nid} length={len(sequences[nid])}, expected={final_len}"
                )

    # Convert terminal sequences to strings
    output: Dict[str, str] = {}
    for term_index, seq_name in terminals_dic.items():
        output[seq_name] = array_to_string(sequences[seq_name])

    return output

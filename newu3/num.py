#!/usr/bin/env python3
"""
num.py

Matrix and sequence utilities for the newu3 module.

This file implements the core Needleman–Wunsch dynamic programming engine:
    • DP matrix initialization
    • DP matrix scoring (with customizable match/mismatch/gap)
    • Full scoring-aware traceback (NEW — supports arbitrary scoring)
    • FASTA parsing into array-of-set format

This engine is used by your translation-based progressive MSA pipeline.
"""

import sys
import numpy as np
from typing import Tuple, Dict, Any


def parse_fa(fasta_str: str) -> Dict[str, np.ndarray]:
    """
    Parse a formatted FASTA-like string into headers and sequences stored
    as numpy arrays of sets.

    Parameters
    ----------
    fasta_str : str
        FASTA-like string containing multiple sequences.

    Returns
    -------
    dict[str, np.ndarray]
        Mapping from header → sequence-as-array-of-sets.
    """
    dic = {}

    lines = fasta_str.strip().splitlines()
    if len(lines) < 4:
        raise ValueError("Formatted FASTA string must contain at least two sequences.")

    lines = [line.strip() for line in lines if line.strip()]

    headers = []
    sequences = []

    for line in lines:
        if line.startswith(">"):
            headers.append(line[1:])
        else:
            sequences.append([set(i) for i in list(line)])

    for i in range(len(headers)):
        dic[headers[i]] = np.array(sequences[i], dtype=object)

    return dic


def init_mat(seq1: np.ndarray, seq2: np.ndarray) -> np.ndarray:
    """
    Initialize a DP scoring matrix for Needleman–Wunsch alignment.

    Parameters
    ----------
    seq1 : np.ndarray
        Sequence 1 as array-of-sets (columns).
    seq2 : np.ndarray
        Sequence 2 as array-of-sets (rows).

    Returns
    -------
    np.ndarray
        Zero-initialized matrix of shape (len(seq2)+1, len(seq1)+1).
    """
    rows = len(seq2) + 1
    cols = len(seq1) + 1
    return np.zeros((rows, cols), dtype=int)


def fill_matrix(
    matrix: np.ndarray,
    seq1_array: np.ndarray,
    seq2_array: np.ndarray,
    match: int = 1,
    mismatch: int = -1,
    indel: int = -1,
) -> np.ndarray:
    """
    Fill a DP matrix using the Needleman–Wunsch scoring scheme.

    Parameters
    ----------
    matrix : np.ndarray
        Zero-initialized DP matrix.
    seq1_array : np.ndarray
        Sequence 1 (columns).
    seq2_array : np.ndarray
        Sequence 2 (rows).
    match : int
        Reward for matching residues.
    mismatch : int
        Penalty for mismatched residues.
    indel : int
        Penalty for insertion/deletion (gap).

    Returns
    -------
    np.ndarray
        Filled DP scoring matrix.
    """

    rows, cols = matrix.shape

    # Initialize top row and left column with indel penalties
    for i in range(1, rows):
        matrix[i, 0] = matrix[i - 1, 0] + indel
    for j in range(1, cols):
        matrix[0, j] = matrix[0, j - 1] + indel

    # Fill DP matrix
    for i in range(1, rows):
        for j in range(1, cols):

            # Determine match vs mismatch scoring
            is_match = bool(seq1_array[j - 1] & seq2_array[i - 1])
            diag_score = match if is_match else mismatch

            score_diag = matrix[i - 1, j - 1] + diag_score
            score_left = matrix[i, j - 1] + indel
            score_up = matrix[i - 1, j] + indel

            matrix[i, j] = max(score_diag, score_left, score_up)

    return matrix


def trace_matrix(
    matrix: np.ndarray,
    seq1_array: np.ndarray,
    seq2_array: np.ndarray,
    match: int = 1,
    mismatch: int = -1,
    indel: int = -1,
) -> np.ndarray:
    """
    Trace back through the DP matrix to reconstruct the optimal global alignment.

    This is a FULL scoring-aware traceback implementation. It correctly respects
    custom match, mismatch, and gap penalties.

    Parameters
    ----------
    matrix : np.ndarray
        DP scoring matrix.
    seq1_array : np.ndarray
        Sequence 1 (columns).
    seq2_array : np.ndarray
        Sequence 2 (rows).
    match : int
        Match reward.
    mismatch : int
        Mismatch penalty.
    indel : int
        Gap penalty.

    Returns
    -------
    np.ndarray
        2-row array: [aligned_seq1, aligned_seq2].
    """

    i = matrix.shape[0] - 1
    j = matrix.shape[1] - 1

    aligned_seq1 = []
    aligned_seq2 = []

    while i > 0 or j > 0:

        # Boundary cases
        if i == 0:
            aligned_seq1.append(seq1_array[j - 1])
            aligned_seq2.append({'-'})
            j -= 1
            continue

        if j == 0:
            aligned_seq1.append({'-'})
            aligned_seq2.append(seq2_array[i - 1])
            i -= 1
            continue

        current_score = matrix[i, j]

        # Determine expected diagonal score
        is_match = bool(seq1_array[j - 1] & seq2_array[i - 1])
        diag_score = match if is_match else mismatch

        score_diag = matrix[i - 1, j - 1] + diag_score
        score_left = matrix[i, j - 1] + indel
        score_up = matrix[i - 1, j] + indel

        # Priority 1: Diagonal (match or mismatch)
        if current_score == score_diag:
            aligned_seq1.append(seq1_array[j - 1])
            aligned_seq2.append(seq2_array[i - 1])
            i -= 1
            j -= 1
            continue

        # Priority 2: Left (gap in seq2)
        if current_score == score_left:
            aligned_seq1.append(seq1_array[j - 1])
            aligned_seq2.append({'-'})
            j -= 1
            continue

        # Priority 3: Up (gap in seq1)
        if current_score == score_up:
            aligned_seq1.append({'-'})
            aligned_seq2.append(seq2_array[i - 1])
            i -= 1
            continue

        raise RuntimeError("Traceback error: no valid move found.")

    aligned_seq1.reverse()
    aligned_seq2.reverse()

    return np.array([aligned_seq1, aligned_seq2], dtype=object)

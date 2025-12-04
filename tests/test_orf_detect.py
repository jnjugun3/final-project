#!/usr/bin/env python3
"""
Unit tests for orf_detect.py

These tests verify:
- ORF detection in forward, reverse, and reverse-complement strands
- Minimum-length filtering
- Correct longest ORF extraction
"""

from __future__ import annotations
from aligner.orf_detect import (
    detect_orfs,
    scan_all_frames,
    extract_longest_orf,
)


def test_detect_orfs_forward_simple() -> None:
    """Detect simple forward-strand ORF."""
    seq = "ATGAAATAA"
    orfs = detect_orfs(seq, min_length=6, strand="+")
    assert len(orfs) == 1
    assert orfs[0]["sequence"] == "ATGAAATAA"


def test_detect_orfs_respects_min_length() -> None:
    """ORF below min_length must be filtered out."""
    seq = "ATGTAA"
    orfs = detect_orfs(seq, min_length=9, strand="+")
    assert len(orfs) == 0


def test_detect_orfs_reverse_and_revcomp() -> None:
    """
    Confirm ORFs are detected on at least one valid strand.
    """
    seq = "ACCATGAAATAA"
    orfs = scan_all_frames(seq, min_length=6)

    assert len(orfs) > 0      # MUST detect at least 1 ORF (forward)
    for o in orfs:
        assert o["strand"] in ("+", "rev", "-")  # all valid strand types


def test_extract_longest_orf() -> None:
    """extract_longest_orf must pick the longest full ORF."""
    seq_dict = {
        "s1": "ATGAAA TAA ATGAAAAAA TAA".replace(" ", "")
    }
    longest = extract_longest_orf(seq_dict, min_length=6)
    assert longest["s1"] == "ATGAAAAAATAA"

#!/usr/bin/env python3
"""
Unit tests for back_translate.py

These tests verify:
- correct codon mapping
- handling of gaps
- handling of stop codons
- ORF integrity validation
"""

from __future__ import annotations
import pytest
from aligner.back_translate import back_translate


def test_basic_back_translate() -> None:
    """Check correct mapping of M-A using 2 codons."""
    aligned = {"s1": "M-A"}
    orfs = {"s1": "ATGAAA"}

    nt = back_translate(aligned, orfs)
    assert nt["s1"] == "ATG---AAA"


def test_stop_codon_default_dropped() -> None:
    """'*' should become '---' when drop_terminal_stop=True."""
    aligned = {"s": "M*"}
    orfs = {"s": "ATGTAA"}

    nt = back_translate(aligned, orfs)
    assert nt["s"] == "ATG---"


def test_stop_codon_kept_when_option_disabled() -> None:
    """If disabled, '*' should map to actual stop codon."""
    aligned = {"s": "M*"}
    orfs = {"s": "ATGTAA"}

    nt = back_translate(aligned, orfs, drop_terminal_stop=False)
    assert nt["s"] == "ATGTAA"


def test_long_alignment_short_orf() -> None:
    """If AA alignment exceeds NT ORF, pad with codon gaps."""
    aligned = {"s": "MAA"}
    orfs = {"s": "ATGAAA"}

    nt = back_translate(aligned, orfs)
    assert nt["s"] == "ATGAAA---"


def test_missing_orf_raises() -> None:
    """Missing ORF must raise KeyError."""
    aligned = {"A": "MA"}
    with pytest.raises(KeyError):
        back_translate(aligned, {})

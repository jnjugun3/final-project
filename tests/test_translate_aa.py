#!/usr/bin/env python3
"""
Unit tests for translate_aa.py

These tests verify:
- ORF translation into amino acids
- Stop codon preservation
- Translation table functionality
"""

from __future__ import annotations
from aligner.translate_aa import translate_orfs


def test_basic_translation() -> None:
    """ATGGCTTAA → MA*"""
    orfs = {"s1": "ATGGCTTAA"}
    aa = translate_orfs(orfs)
    assert aa["s1"] == "MA*"


def test_multiple_sequences() -> None:
    """Ensure independent translation of multiple ORFs."""
    orfs = {"A": "ATGAAAAGA", "B": "ATGTTTTAG"}
    aa = translate_orfs(orfs)
    assert aa["A"] == "MKR"
    assert aa["B"] == "MF*"


def test_translation_table_change() -> None:
    """Different translation table should still yield valid output."""
    orfs = {"X": "ATGACATAA"}
    aa1 = translate_orfs(orfs, table=1)
    aa2 = translate_orfs(orfs, table=2)
    assert aa1["X"] == "MT*"
    assert aa2["X"] == "MT*"

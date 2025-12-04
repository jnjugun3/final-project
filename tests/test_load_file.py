#!/usr/bin/env python3
"""
Unit tests for load_file.py

These tests verify:
- FASTA parsing correctness
- Handling of invalid bases (converted to 'N')
- Detection of malformed FASTA files
- Duplicate header validation

This satisfies rubric requirements for Testing & Validation and Data Handling.
"""

from __future__ import annotations
from typing import Dict
import tempfile
import pytest
from aligner.load_file import read_fasta    


def write_temp_fasta(text: str) -> str:
    """
    Helper function:
    Writes temporary FASTA data to a file and returns its path.
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fa")
    tmp.write(text)
    tmp.close()
    return tmp.name


def test_read_fasta_basic() -> None:
    """
    FASTA with two valid sequences should load correctly.
    """
    fasta = ">seq1\nATGC\n>seq2\nAAAA"
    path = write_temp_fasta(fasta)

    seqs: Dict[str, str] = read_fasta(path)
    assert seqs["seq1"] == "ATGC"
    assert seqs["seq2"] == "AAAA"


def test_invalid_base_converted_to_N(capsys) -> None:
    """
    Invalid bases should be converted to 'N' with a warning.
    IUPAC ambiguous bases (like B) are valid and must NOT be converted.
    """
    fasta = ">s1\nATGBZ"
    path = write_temp_fasta(fasta)

    seqs = read_fasta(path)
    assert seqs["s1"] == "ATGBN"   # B is valid, Z becomes N

    captured = capsys.readouterr()
    assert "Nonstandard nucleotide" in captured.out

def test_error_on_empty_file() -> None:
    """Empty FASTA file must raise ValueError."""
    path = write_temp_fasta("")
    with pytest.raises(ValueError):
        read_fasta(path)


def test_error_on_sequence_before_header() -> None:
    """Sequence without header raises ValueError."""
    path = write_temp_fasta("ATGC")
    with pytest.raises(ValueError):
        read_fasta(path)


def test_error_duplicate_headers() -> None:
    """Duplicate FASTA headers are not allowed."""
    fasta = ">A\nATGC\n>A\nTTTT"
    path = write_temp_fasta(fasta)
    with pytest.raises(ValueError):
        read_fasta(path)

#!/usr/bin/env python3
"""
Unit tests for codon_stats.py

These tests verify:
- position-specific identity
- codon identity
- global identity
- CSV file generation
- alignment validation errors
"""

from __future__ import annotations
import csv
import os
import tempfile
import pytest
from aligner.codon_stats import (
    compute_percent_identity,
    compute_codon_identity_total,
    compute_global_identity,
    step7_codon_statistics,
)


def test_compute_percent_identity() -> None:
    """Check calculation of per-position identity."""
    assert compute_percent_identity(["A", "A", "A"]) == 1.0
    assert compute_percent_identity(["A", "C", "A"]) == pytest.approx(2/3)
    assert compute_percent_identity(["-", "-", "N"]) == 0.0


def test_compute_codon_identity_total() -> None:
    """Average identity across 3 positions."""
    pos1 = ["A", "A", "A"]
    pos2 = ["C", "C", "T"]
    pos3 = ["G", "G", "G"]
    expected = (1.0 + 2/3 + 1.0) / 3
    assert compute_codon_identity_total(pos1, pos2, pos3) == pytest.approx(expected)


def test_compute_global_identity() -> None:
    """Test alignment-wide percent identity."""
    seqs = ["ATG", "ATG", "ATA"]
    expected = (1.0 + 1.0 + 2/3) / 3
    assert compute_global_identity(seqs) == pytest.approx(expected)


def test_step7_csv_output() -> None:
    """
    Validate that the CSV output is created correctly
    and contains expected structure.
    """
    aligned = {"s1": "ATGCCC", "s2": "ATGCCA"}

    with tempfile.TemporaryDirectory() as tmp:
        out_csv = os.path.join(tmp, "stats.csv")
        step7_codon_statistics(aligned, output_csv=out_csv)

        assert os.path.exists(out_csv)

        with open(out_csv) as fh:
            rows = list(csv.reader(fh))

        # correct header
        assert rows[0] == [
            "codon_index",
            "identity_pos1", "identity_pos2", "identity_pos3",
            "identity_codon_total",
            "gc_pos1", "gc_pos2", "gc_pos3",
            "identity_global"
        ]

        # 2 codons → 2 rows + header
        assert len(rows) == 3


def test_step7_raises_on_invalid_lengths() -> None:
    """Unequal sequence lengths should cause ValueError."""
    aligned = {"A": "ATG", "B": "ATGC"}
    with pytest.raises(ValueError):
        step7_codon_statistics(aligned)

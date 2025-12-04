#!/usr/bin/env python3
"""
Unit tests for kmer_sim.py

These tests verify:
- k-mer extraction (without stop codons)
- Z-distance computation
- Greedy sequence ordering
- Guide tree structure
"""

from __future__ import annotations
from aligner.kmer_sim import (
    sequence_to_kmers,
    compute_z,
    compute_pair_distances,
    order_sequences,
    build_tree_dict_from_kmers,
    kmer_guide_tree,
)


def test_sequence_to_kmers_basic() -> None:
    """Stop codons (*) must be removed before k-mer extraction."""
    kmers = sequence_to_kmers("MA*KTL", k=2)
    assert kmers == ["MA", "AK", "KT", "TL"]


def test_compute_z_distance() -> None:
    """Z-distance equals sum of absolute frequency differences."""
    assert compute_z(["AAA", "BBB"], ["AAA"]) == 1


def test_pairwise_distance_computation() -> None:
    kmers = {"A": ["AA", "AA"], "B": ["AA"]}
    dist = compute_pair_distances(kmers)
    assert dist[("A", "B")] == 1


def test_order_sequences_simple_case() -> None:
    """Closest pair should appear first."""
    dist = {("A", "B"): 0, ("A", "C"): 5, ("B", "C"): 5}
    ordered = order_sequences(dist, ["A", "B", "C"])
    assert ordered[-1] == "C"


def test_full_kmer_guide_tree() -> None:
    """Guide tree must produce correct number of terminals & internal nodes."""
    aa = {"A": "AAA", "B": "AAA", "C": "CCC"}
    ordered, terminals, internal_nodes = kmer_guide_tree(aa, k=2)

    assert len(ordered) == 3
    assert len(terminals) == 3
    assert len(internal_nodes) == 2

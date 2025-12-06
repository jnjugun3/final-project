#!/usr/bin/env python3

"""
Step 4. K-mer Similarity & Guide Tree Construction (String Version)

This module performs:

1. Extracts k-mers from amino-acid sequences (default k=3)
2. Computes pairwise distances using the Z-difference metric:
       Z = Σ |freq_seq1(k) - freq_seq2(k)|
3. Orders sequences so that the closest sequences appear first
4. Builds a right-branching guide tree for progressive NW alignment

Inputs:
    aa_dict = { seq_id : "MAKTLLA..." }

Outputs:
    ordered        : List[str]             # order used for progressive alignment
    terminals_dic  : Dict[int, str]        # terminal node → sequence ID
    internal_nodes : Dict[int, List[Any]]  # internal node → [child1, child2]
"""

from typing import Dict, List, Tuple, Any


def sequence_to_kmers(seq, k=3):
    """
    Convert an amino-acid sequence into overlapping k-mers.

    Notes:
    - Removes '*' stop codons (common in translated ORFs)
    - Uses sliding-window extraction of k characters
    """
    seq = seq.replace("*", "")        # remove stop codons to avoid false k-mers
    return [seq[i:i+k]                # create substring of length k
            for i in range(len(seq) - k + 1)]   # sliding window over sequence


def compute_z(k1, k2):
    """
    Compute Z-distance between two k-mer lists.
    Z = Σ |freq1(k) - freq2(k)|
    """

    freq1, freq2 = {}, {}             # frequency tables

    # Count frequency of each k-mer in first sequence
    for k in k1:
        freq1[k] = freq1.get(k, 0) + 1

    # Count frequency of each k-mer in second sequence
    for k in k2:
        freq2[k] = freq2.get(k, 0) + 1

    # All unique k-mers from both sequences
    all_k = set(freq1) | set(freq2)

    # Sum absolute differences → Z-distance
    return sum(abs(freq1.get(k, 0) - freq2.get(k, 0)) for k in all_k)


def compute_pair_distances(kmer_dict):
    """
    Compute Z-distance for each unordered pair of sequences.
    """
    ids = list(kmer_dict.keys())      # list of sequence IDs
    distances = {}                    # dictionary to store pairwise distances

    # Build pairwise combinations (A,B)
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            a, b = ids[i], ids[j]
            distances[(a, b)] = compute_z(kmer_dict[a], kmer_dict[b])  # compute Z-distance

    return distances


def order_sequences(distance_dict, original_ids):
    """
    Order sequences using a greedy nearest-neighbor strategy.
    """

    if len(original_ids) <= 1:
        return original_ids[:]        # return input if only 0–1 sequences

    seqs = list(original_ids)         # ensure stable ordering

    # If there are pairwise distances, find closest pair
    if distance_dict:
        best_pair = min(distance_dict, key=distance_dict.get)
        ordered = [best_pair[0], best_pair[1]]     # initial cluster
    else:
        return seqs[:]                # fallback if all sequences identical

    remaining = set(seqs) - set(ordered)

    # Iteratively attach closest remaining sequence
    while remaining:
        best_next = None
        best_dist = float("inf")

        for seq in remaining:
            # distance from this seq to current cluster
            d = min(
                distance_dict.get((min(seq, o), max(seq, o)), float("inf"))
                for o in ordered
            )

            if d < best_dist:               # keep best match
                best_dist = d
                best_next = seq

        ordered.append(best_next)           # attach next closest sequence
        remaining.remove(best_next)

    return ordered


def build_tree_dict_from_kmers(kmer_dict):
    """
    Build a right-branching binary guide tree using ordered sequence IDs.
    """

    ids = list(kmer_dict.keys())               # list of sequence names
    dist = compute_pair_distances(kmer_dict)   # compute all pairwise distances
    ordered_ids = order_sequences(dist, ids)   # compute progressive order

    terminals_dic = {                          # map terminal node numbers → names
        i: ordered_ids[i - 1]
        for i in range(1, len(ordered_ids) + 1)
    }

    n = len(ordered_ids)
    internal_nodes = {}                        # dictionary for tree structure

    if n <= 1:
        return terminals_dic, internal_nodes   # simple tree

    # First internal node joins first two sequences
    node_id = n + 1
    internal_nodes[node_id] = [terminals_dic[1], terminals_dic[2]]

    # Right-branch: attach each new terminal to previous internal node
    for idx in range(3, n + 1):
        prev = node_id          # previous internal node
        node_id += 1            # new internal node ID
        internal_nodes[node_id] = [prev, terminals_dic[idx]]

    return terminals_dic, internal_nodes


def kmer_guide_tree(aa_dict, k=3):
    """
    Run the entire Step-4 pipeline:
        amino acids → k-mers → pairwise distances → ordering → guide tree
    """

    # Convert each AA sequence to list of k-mers
    kmer_dict = {sid: sequence_to_kmers(seq, k)
                 for sid, seq in aa_dict.items()}

    # Build the guide tree
    terminals_dic, internal_nodes = build_tree_dict_from_kmers(kmer_dict)

    # Extract progressive order directly from terminals dictionary (1..N)
    ordered_ids = [terminals_dic[i] for i in range(1, len(terminals_dic) + 1)]

    return ordered_ids, terminals_dic, internal_nodes

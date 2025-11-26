#!/usr/bin/env python3

from typing import Dict, List, Tuple
from collections import Counter


def sequence_to_kmers(seq: str, k: int = 3) -> List[str]:
    """
    Convert an amino-acid sequence into a list of k-mers.

    Parameters
    ----------
    seq : str
        Amino-acid sequence (string of residues).
    k : int, default=3
        Length of each k-mer.

    Returns
    -------
    List[str]
        A list of k-mer strings extracted using a sliding window.
    """
    if seq is None:
        return []
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]


def unique_pairs(ids: List[str]) -> List[Tuple[str, str]]:
    """
    Generate all unique unordered pairs of sequence IDs.

    Parameters
    ----------
    ids : List[str]
        Sequence ID list.

    Returns
    -------
    List[Tuple[str, str]]
        List of unique (A, B) pairs where A < B.
    """
    return [(ids[i], ids[j]) for i in range(len(ids)) for j in range(i + 1, len(ids))]


def compute_z_value(kmers1: Counter, kmers2: Counter) -> int:
    """
    Compute z-distance between two sequences based on k-mer counts:

        z = Σ | countA(k) - countB(k) |

    Parameters
    ----------
    kmers1 : Counter
        Counter of k-mer frequencies for the first sequence.
    kmers2 : Counter
        Counter of k-mer frequencies for the second sequence.

    Returns
    -------
    int
        The z-distance between the two sequences.
    """
    all_kmers = set(kmers1.keys()) | set(kmers2.keys())
    return sum(abs(kmers1.get(k, 0) - kmers2.get(k, 0)) for k in all_kmers)


# -------------------------------------------------------
# Compute all pairwise z-distances
# -------------------------------------------------------

def compute_all_z_distances(kmer_dict: Dict[str, List[str]]) -> Dict[Tuple[str, str], int]:
    """
    Compute z-distance for every pair of sequences using extracted k-mers.

    Parameters
    ----------
    kmer_dict : Dict[str, List[str]]
        Maps seq_id -> list of k-mers.

    Returns
    -------
    Dict[Tuple[str, str], int]
        Pairwise z-distances: {(id1, id2): z_value}
    """
    kmer_counts = {sid: Counter(klist) for sid, klist in kmer_dict.items()}

    ids = list(kmer_dict.keys())
    distances: Dict[Tuple[str, str], int] = {}

    for a, b in unique_pairs(ids):
        distances[(a, b)] = compute_z_value(kmer_counts[a], kmer_counts[b])

    return distances


def order_by_kmer_similarity(kmer_dict: Dict[str, List[str]]) -> List[str]:
    """
    Order sequence IDs from most similar to least similar using z-distance.
    Lower z-distance = higher similarity.

    Parameters
    ----------
    kmer_dict : Dict[str, List[str]]
        seq_id -> list of k-mers.

    Returns
    -------
    List[str]
        Ordered list of sequence IDs for guide-tree construction.
    """
    distances = compute_all_z_distances(kmer_dict)
    ids = list(kmer_dict.keys())

    if len(ids) <= 1:
        return ids

    # Step 1 — find closest pair (min z-distance)
    (best_a, best_b), best_z = min(distances.items(), key=lambda x: x[1])
    ordered = [best_a, best_b]
    remaining = set(ids) - set(ordered)

    # Step 2 — add remaining sequences by smallest distance to cluster
    while remaining:
        best_cand = None
        best_score = float("inf")

        for cand in remaining:
            score = min(
                distances.get((cand, existing),
                              distances.get((existing, cand), float("inf")))
                for existing in ordered
            )
            if score < best_score:
                best_score = score
                best_cand = cand

        ordered.append(best_cand)
        remaining.remove(best_cand)

    return ordered



def build_tree_dict(
    ordered_ids: List[str]
) -> Tuple[Dict[int, str], Dict[int, List]]:
    """
    Build a right-branching guide tree from ordered sequence IDs.

    Terminal nodes: IDs 1..N  
    Internal nodes: IDs N+1 .. (2N-1)

    Parameters
    ----------
    ordered_ids : List[str]
        Sequence IDs ordered by similarity (most similar first).

    Returns
    -------
    terminals : Dict[int, str]
        Maps terminal node numbers to sequence IDs.
    internal_nodes : Dict[int, List]
        Maps internal node numbers to their children.
        Children may be raw sequence names or internal node IDs.
    """
    terminals = {i + 1: ordered_ids[i] for i in range(len(ordered_ids))}
    internal_nodes: Dict[int, List] = {}

    n = len(ordered_ids)
    if n <= 1:
        return terminals, internal_nodes

    # First internal node: join terminal 1 and terminal 2
    next_node = n + 1
    internal_nodes[next_node] = [terminals[1], terminals[2]]
    last_node = next_node

    # Add remaining terminals
    for i in range(3, n + 1):
        next_node += 1
        internal_nodes[next_node] = [last_node, terminals[i]]
        last_node = next_node

    return terminals, internal_nodes

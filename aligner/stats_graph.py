#!/usr/bin/env python3
"""
Step 8 – Codon Statistics Plotting


This module generates biologically meaningful visualizations from the
codon_stats.csv file produced in Step 7.

    1. Percent identity at codon positions 1, 2, and 3
    2. A skyline plot of codon-level identity (identity_codon_total)
    3. GC content at codon positions 1, 2, and 3
    4. A pairwise percent identity heatmap (global similarity)

All plots are saved as PNG images into the specified output directory.


"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict



#  Pairwise Identity Heatmap (Plot 4)


def plot_pairwise_identity_heatmap(aligned_nt_dict: Dict[str, str], outdir: str) -> None:
    """
    Generate a heatmap showing pairwise global percent identity between
    all sequences in the codon-aware nucleotide alignment.

    This provides a simple, visual overview of which sequences are more
    similar or more divergent based on the full MSA.

    Saves:
        - pairwise_identity_heatmap.png
        - heatmap_labels.txt (mapping short labels → real IDs)
    """

    seq_ids = list(aligned_nt_dict.keys())
    seqs = list(aligned_nt_dict.values())
    n = len(seq_ids)

    # Use short readable labels for plotting
    short_ids = [f"Seq{i+1}" for i in range(n)]

    # Save label mapping for reference
    map_path = os.path.join(outdir, "heatmap_labels.txt")
    with open(map_path, "w") as fp:
        for sid, short in zip(seq_ids, short_ids):
            fp.write(f"{short} → {sid}\n")

    # Compute pairwise percent identity matrix
    identity_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):

            s1, s2 = seqs[i], seqs[j]
            matches = 0
            valid = 0

            for a, b in zip(s1, s2):
                if a in "ATCG" and b in "ATCG":
                    valid += 1
                    if a == b:
                        matches += 1

            identity_matrix[i, j] = matches / valid if valid > 0 else 0.0

    # Convert to DataFrame for clean labeling
    df = pd.DataFrame(identity_matrix, index=short_ids, columns=short_ids)

    # Make heatmap
    plt.figure(figsize=(8, 6))
    plt.imshow(df, cmap="viridis", vmin=0, vmax=1)
    plt.colorbar(label="Percent Identity")

    plt.title("Pairwise Global Percent Identity Heatmap")
    plt.xticks(np.arange(n), short_ids, rotation=45, ha="right")
    plt.yticks(np.arange(n), short_ids)

    plt.tight_layout()
    out_path = os.path.join(outdir, "pairwise_identity_heatmap.png")
    plt.savefig(out_path)
    plt.close()

    print(f"[PLOT] → {out_path}")
    print(f"[INFO] Heatmap labels → {map_path}")



#  Step 8 Plotting Function  


def step8_generate_plots(csv_file: str, aligned_nt_dict: Dict[str, str], outdir: str) -> None:
    """
    Generate all Step 8 plots:
        1) Positional identity (pos1/pos2/pos3)
        2) Codon identity skyline
        3) GC content at codon positions
        4) Pairwise global identity heatmap

    Parameters
    ----------
    csv_file : str
        Path to codon_stats.csv created in Step 7.
    aligned_nt_dict : dict
        Final nucleotide alignment dictionary (for heatmap).
    outdir : str
        Directory where plots will be saved.
    """

    # Ensure directory exists
    os.makedirs(outdir, exist_ok=True)

    # Load codon statistics
    df = pd.read_csv(csv_file)
    codons = df["codon_index"]

    
    # 1. Positional Identity (identity_pos1, pos2, pos3)
    
    plt.figure(figsize=(12, 5))
    plt.plot(codons, df["identity_pos1"], label="Position 1", linewidth=2)
    plt.plot(codons, df["identity_pos2"], label="Position 2", linewidth=2)
    plt.plot(codons, df["identity_pos3"], label="Position 3", linewidth=2)

    plt.xlabel("Codon Index")
    plt.ylabel("Percent Identity")
    plt.title("Positional Percent Identity per Codon")
    plt.legend()
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "identity_positions.png"))
    plt.close()

    
    # 2. Codon Identity Skyline (identity_codon_total)
   
    plt.figure(figsize=(12, 4))
    plt.plot(
        codons,
        df["identity_codon_total"],
        color="purple",
        linewidth=2
    )
    plt.xlabel("Codon Index")
    plt.ylabel("Codon Identity")
    plt.title("Codon Identity Skyline Plot")
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "identity_codon_total.png"))
    plt.close()

     
    # 3. GC Content Plot (gc1, gc2, gc3)
  
    plt.figure(figsize=(12, 5))
    plt.plot(codons, df["gc_pos1"], label="GC1", linewidth=2)
    plt.plot(codons, df["gc_pos2"], label="GC2", linewidth=2)
    plt.plot(codons, df["gc_pos3"], label="GC3", linewidth=2)

    plt.xlabel("Codon Index")
    plt.ylabel("GC Fraction")
    plt.title("GC Content at Codon Positions")
    plt.legend()
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "gc_positions.png"))
    plt.close()

    
    # 4. Pairwise Identity Heatmap
   
    plot_pairwise_identity_heatmap(aligned_nt_dict, outdir)

    print(f"[Step 8] Plots saved → {outdir}")

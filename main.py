#!/usr/bin/env python3
"""
main.py — Translation-Based Multiple Sequence Alignment Pipeline
----------------------------------------------------------------

This script coordinates all 9 stages of the translation-based
multiple sequence alignment (MSA) pipeline. It is intentionally modular
and each stage lives in its own module.

Pipeline Steps
--------------
1. Load nucleotide FASTA
2. Detect longest ORF in each sequence
3. Translate ORFs to amino acids
4. Extract k-mers (AA)
5. Compute k-mer similarity matrix + build guide tree
6. Progressive Needleman–Wunsch AA alignment
7. Codon-aware back-translation to nucleotide space
8. Codon position statistics (percent identity, GC%)
9. Plot generation

Now includes:
- CPU core count  
- Total runtime  
- CLI usage summary printed before pipeline runs
"""

import argparse
import os
import shutil
import sys
import time
import multiprocessing
from typing import Dict, List, Any
import numpy as np

# Import modules for each step
from aligner.load_file import read_fasta
from aligner.orf_detect import extract_longest_orf
from aligner.translate_aa import translate_orfs
from aligner.kmer_sim import kmer_guide_tree
from aligner.needleman_wunch import progressive_nw
from aligner.back_translate import back_translate
from aligner.codon_stats import step7_codon_statistics
from aligner.stats_graph import step8_generate_plots


# Translation tables
TRANSLATION_TABLES: Dict[int, str] = {
    1: "Standard (Universal) Code",
    2: "Vertebrate Mitochondrial",
    4: "Mold / Protozoan Mitochondrial",
    11: "Bacterial / Archaeal / Plastid"
}


# -------------------------------------------------------------
# CLEAN CLI parser  
# -------------------------------------------------------------
class CleanParser(argparse.ArgumentParser):
    """Custom ArgumentParser that prints cleaner error messages."""

    def error(self, message: str) -> None:
        sys.stderr.write("\nERROR: Incorrect CLI input → " + message + "\n")
        sys.stderr.write("Use: python3 main.py -h   for full help.\n\n")
        self.print_usage(sys.stderr)
        sys.exit(2)


# -------------------------------------------------------------
# Pipeline Controller Class
# -------------------------------------------------------------
class PipelineRunner:
    """Runs the full translation-based MSA workflow."""

    def __init__(self, args: argparse.Namespace) -> None:
        """Initialize state and prepare output folders."""
        self.args = args

        fasta_name = os.path.basename(args.input)
        self.run_name = args.run_name if args.run_name else os.path.splitext(fasta_name)[0]

        # Create result directories
        self.dirs = self._make_results_dirs(self.run_name)
        sys.stdout.write(f"[INFO] Created results directory: results/{self.run_name}/\n")

        # Pipeline storage dictionaries
        self.seq_dict: Dict[str, str] = {}
        self.longest_orfs: Dict[str, str] = {}
        self.aa_dict: Dict[str, str] = {}
        self.terminals_dic: Dict[int, str] = {}
        self.internal_nodes: Dict[int, List[Any]] = {}
        self.aligned_aa: Dict[str, str] = {}
        self.nt_alignment: Dict[str, str] = {}
        self.stats_csv: str = ""

    # ----------------------------
    # Directory setup
    # ----------------------------
    def _make_results_dirs(self, run_name: str) -> Dict[str, str]:
        """Create results/<run_name>/ and all subdirectories."""
        root = os.path.join("results", run_name)

        if os.path.exists(root):
            shutil.rmtree(root)

        subdirs = {
            "root": root,
            "input": os.path.join(root, "input"),
            "orfs": os.path.join(root, "orfs"),
            "aa": os.path.join(root, "aa"),
            "nt_backtranslated": os.path.join(root, "nt_backtranslated"),
            "stats": os.path.join(root, "stats"),
            "plots": os.path.join(root, "plots"),
            "logs": os.path.join(root, "logs"),
        }

        for d in subdirs.values():
            os.makedirs(d, exist_ok=True)

        return subdirs

    # ----------------------------
    # CLI Summary
    # ----------------------------
    def print_cli_summary(self) -> None:
        """Print available command-line options and translation tables."""
        sys.stdout.write("\n========== COMMAND-LINE OPTIONS ==========\n")
        sys.stdout.write("  -i, --input <file>         Input FASTA file\n")
        sys.stdout.write("      --run-name <name>      Optional output folder name\n")
        sys.stdout.write("      --min-orf-length <n>   Minimum ORF length\n")
        sys.stdout.write("  -t, --table <id>           Translation table ID\n")
        sys.stdout.write("  -k, --kmer <n>             K-mer size\n")
        sys.stdout.write("      --match <n>            NW match score\n")
        sys.stdout.write("      --mismatch <n>         NW mismatch penalty\n")
        sys.stdout.write("      --gap <n>              NW gap penalty\n")

        sys.stdout.write("\nAvailable Translation Tables:\n")
        for code, desc in TRANSLATION_TABLES.items():
            sys.stdout.write(f"  {code:<3} → {desc}\n")

        sys.stdout.write("==========================================\n\n")

    # ----------------------------
    # Step 1 — Load FASTA
    # ----------------------------
    def run_load_fasta(self) -> None:
        """Step 1 — Load FASTA."""
        sys.stdout.write("[LOAD] Reading FASTA...\n")
        self.seq_dict = read_fasta(self.args.input)
        shutil.copy(self.args.input, os.path.join(self.dirs["input"], "input_sequences.fasta"))

    # ----------------------------
    # Step 2 — ORF detection
    # ----------------------------
    def run_extract_orfs(self) -> None:
        """Step 2 — Detect longest ORF."""
        sys.stdout.write("[ORF] Extracting longest ORFs...\n")
        self.longest_orfs = extract_longest_orf(self.seq_dict, self.args.min_orf_length)

        out = os.path.join(self.dirs["orfs"], "longest_orfs_nt.fasta")
        with open(out, "w") as f:
            for sid, seq in self.longest_orfs.items():
                f.write(f">{sid}\n{seq}\n")

    # ----------------------------
    # Step 3 — Translation
    # ----------------------------
    def run_translate_orfs(self) -> None:
        """Step 3 — Translate ORFs to amino acids."""
        sys.stdout.write("[AA] Translating ORFs...\n")
        self.aa_dict = translate_orfs(self.longest_orfs, self.args.table)

        out = os.path.join(self.dirs["aa"], "translated_orfs_aa.fasta")
        with open(out, "w") as f:
            for sid, aa in self.aa_dict.items():
                f.write(f">{sid}\n{aa}\n")

    # ----------------------------
    # Step 4–5 — K-mer similarity tree
    # ----------------------------
    def run_kmer_tree(self) -> None:
        """Step 4–5 — Build k-mer similarity guide tree."""
        sys.stdout.write("[KMER] Computing guide tree...\n")
        (_, self.terminals_dic, self.internal_nodes) = kmer_guide_tree(
            self.aa_dict, self.args.kmer
        )

    # ----------------------------
    # Step 6 — NW alignment
    # ----------------------------
    def run_needleman_wunsch(self) -> None:
        """Step 6 — Progressive NW alignment."""
        sys.stdout.write("[ALIGN] Running NW alignment...\n")

        seq_arrays = {
            sid: np.array([set(res) for res in aa], dtype=object)
            for sid, aa in self.aa_dict.items()
        }

        self.aligned_aa = progressive_nw(
            sequences=seq_arrays,
            terminals_dic=self.terminals_dic,
            internal_nodes=self.internal_nodes,
            match=self.args.match,
            mismatch=self.args.mismatch,
            gap=self.args.gap
        )

        out = os.path.join(self.dirs["aa"], "progressive_alignment_aa.fasta")
        with open(out, "w") as f:
            for sid, aln in self.aligned_aa.items():
                f.write(f">{sid}\n{aln}\n")

    # ----------------------------
    # Step 7 — Back-translation
    # ----------------------------
    def run_back_translate(self) -> None:
        """Step 7 — Convert AA alignment to codon alignment."""
        sys.stdout.write("[BACK] Back-translating to nucleotides...\n")

        self.nt_alignment = back_translate(
            aligned_aa=self.aligned_aa,
            longest_orfs=self.longest_orfs,
            drop_terminal_stop=False
        )

        out = os.path.join(self.dirs["nt_backtranslated"], "codon_alignment_nt.fasta")
        with open(out, "w") as f:
            for sid, aln in self.nt_alignment.items():
                f.write(f">{sid}\n{aln}\n")

    # ----------------------------
    # Step 8 — Codon statistics
    # ----------------------------
    def run_codon_statistics(self) -> None:
        """Step 8 — Compute codon statistics."""
        sys.stdout.write("[STATS] Computing codon statistics...\n")
        csv_path = os.path.join(self.dirs["stats"], "codon_stats.csv")
        step7_codon_statistics(self.nt_alignment, csv_path)
        self.stats_csv = csv_path

    # ----------------------------
    # Step 9 — Graphing
    # ----------------------------
    def run_generate_plots(self) -> None:
        """Step 9 — Generate identity, GC%, and heatmap plots."""
        sys.stdout.write("[PLOT] Generating plots...\n")
        step8_generate_plots(
            csv_file=self.stats_csv,
            aligned_nt_dict=self.nt_alignment,
            outdir=self.dirs["plots"]
        )

    # ----------------------------
    # Execute entire pipeline
    # ----------------------------
    def run(self) -> None:
        """Execute all pipeline stages + runtime + CPU report."""
        self.print_cli_summary()

        start = time.time()

        self.run_load_fasta()
        self.run_extract_orfs()
        self.run_translate_orfs()
        self.run_kmer_tree()
        self.run_needleman_wunsch()
        self.run_back_translate()
        self.run_codon_statistics()
        self.run_generate_plots()

        runtime = time.time() - start
        cpu_cores = multiprocessing.cpu_count()

        sys.stdout.write(f"\n[INFO] CPU cores: {cpu_cores}\n")
        sys.stdout.write(f"[INFO] Total runtime: {runtime:.2f} seconds\n")
        sys.stdout.write(
            f"\n[DONE] Pipeline complete.\n"
            f"Results saved to: results/{self.run_name}/\n"
        )

        # Save run metadata
        log_file = os.path.join(self.dirs["logs"], "run_parameters.txt")
        with open(log_file, "w") as f:
            f.write("Translation-Based MSA Pipeline — Run Log\n\n")
            f.write(f"Run Name: {self.run_name}\n")
            f.write(f"Input FASTA: {self.args.input}\n\n")
            f.write("Parameters:\n")
            f.write(f"  ORF Min Length: {self.args.min_orf_length}\n")
            f.write(f"  Translation Table: {self.args.table}\n")
            f.write(f"  K-mer Size: {self.args.kmer}\n")
            f.write(f"  Match Score: {self.args.match}\n")
            f.write(f"  Mismatch Penalty: {self.args.mismatch}\n")
            f.write(f"  Gap Penalty: {self.args.gap}\n\n")
            f.write("System Info:\n")
            f.write(f"  CPU Cores: {cpu_cores}\n\n")
            f.write("Performance:\n")
            f.write(f"  Total Runtime: {runtime:.2f} seconds\n")


# -------------------------------------------------------------
# CLI Parser
# -------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Define command-line arguments."""
    parser = CleanParser(
        description="Translation-Based MSA Pipeline",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-i", "--input", required=True, help="Input FASTA file.")
    parser.add_argument("--run-name", default=None, help="Optional name for results folder.")
    parser.add_argument("--min-orf-length", type=int, default=30, help="Minimum ORF length.")
    parser.add_argument("-t", "--table", type=int, default=1, help="NCBI translation table ID.")
    parser.add_argument("-k", "--kmer", type=int, default=3, help="K-mer size.")
    parser.add_argument("--match", type=int, default=1, help="Match score for NW.")
    parser.add_argument("--mismatch", type=int, default=-1, help="Mismatch penalty.")
    parser.add_argument("--gap", type=int, default=-1, help="Gap penalty.")

    return parser.parse_args()


# -------------------------------------------------------------
# Program Entry Point
# -------------------------------------------------------------
def main() -> None:
    """Start the pipeline."""
    args = parse_args()
    runner = PipelineRunner(args)
    runner.run()


if __name__ == "__main__":
    main()

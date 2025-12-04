# 🧬 **Translation-Based Multiple Sequence Alignment Pipeline**

### *Final Project Report — Jakes Njuguna*

------

### 📝 **Abstract**

This project implements a complete, modular, and FAIR-compliant translation-based multiple sequence alignment (MSA) pipeline.
The workflow processes raw nucleotide FASTA sequences, detects open reading frames (ORFs), translates them into amino acids, aligns them using a progressive Needleman–Wunsch algorithm guided by k-mer similarity, back-translates the aligned amino-acid sequences into codon-respecting nucleotide alignments, and computes codon-position statistics with visualization outputs.

The pipeline is organized into standalone Python modules inside `aligner/`, orchestrated through `main.py`, and fully validated using a comprehensive unit-test suite located in `tests`

### 📂 **Project Directory Layout**

Final_Project_Prog2/
│
├── aligner/
│   ├── load_file.py            # Step 1: FASTA reader + validation
│   ├── orf_detect.py           # Step 2: ORF scanning
│   ├── translate_aa.py         # Step 3: Translate ORFs → AA
│   ├── kmer_sim.py             # Step 4: k-mer similarity + guide tree
│   ├── needleman_wunch.py      # Step 5: Progressive NW alignment
│   ├── back_translate.py       # Step 6: Codon-aware back-translation
│   ├── codon_stats.py          # Step 7: Codon identity & GC stats
│   └── stats_graph.py          # Step 8: Plotting + heatmap
│
├── tests/                      # Complete pytest suite
│   ├── test_load_file.py
│   ├── test_orf_detect.py
│   ├── test_translate_aa.py
│   ├── test_kmer_sim.py
│   ├── test_needleman_wunch.py
│   ├── test_back_translate.py
│   ├── test_codon_stats.py
│   └── __init__.py
│
├── results/                    # Auto-generated outputs
│   ├── aa/
│   ├── nt_backtranslated/
│   ├── stats/
│   ├── plots/
│   ├── logs/
│   └── input/
│
├── test.fasta                  # Example input FASTA
├── main.py                     # Pipeline driver
└── README.md

### 🧰 **Environment Setup**

#### Recommended: Conda environment

conda create -n bioinfo_env python3
conda activate bioinfo_env
pip install biopython numpy pandas matplotlib pytest

### ▶️ **Running the Pipeline**

Run : 

python3 main.py -i test.fasta -o results/

### All CLI arguments:

| Flag         | Meaning             | Default    |
| ------------ | ------------------- | ---------- |
| `-i`         | Input FASTA file    | required   |
| `-o`         | Output directory    | `results/` |
| `--k`        | k-mer size          | `3`        |
| `--match`    | NW match score      | `1`        |
| `--mismatch` | NW mismatch penalty | `-1`       |
| `--gap`      | NW gap penalty      | `-1`       |
| `--no-plots` | Disable plotting    | False      |

#### See help:

python3 main.py --help

## 🧵 **Pipeline Overview**

FASTA → ORF Detection → AA Translation → k-mer Similarity
     → Guide Tree → Progressive NW → Codon Back-Translation
     → Codon Stats → Plots & CSVs

### 🧵  Full Pipeline Description**

##### **Step 1 — Loading FASTA (`load_file.py`)**

This module:

- Reads nucleotide FASTA files

- Ensures every header starts with `>`

- Ensures only valid IUPAC nucleotide characters

- Converts invalid characters → `"N"`

- Prevents duplicate sequence IDs

- Print warnings when needed

  Input FASTA:

  seq1
  ATGCZ

  Output:

  {"s1": "ATGCN"}

#####  **Step 2 — ORF Detection (`orf_detect.py`)**

For each input sequence:

- Scan **3 reading frames**
- Scan **forward**, **reverse**, **reverse complement**
- Identify ORFs that start with `"ATG"` and end with `TAA/TAG/TGA`
- Extract **longest ORF per sequence**

Input :

ACCATGAAATAA

Output:

{"sequence": "ATGAAATAA", "strand": "-"}

#### **Step 3 — Amino-Acid Translation (`translate_aa.py`)**

uses Biopython translation:

Seq(nt).translate(table=1, to_stop=False)

Input :

"ATGGCTTAA"

Output:

"MA*"

#### Step 4 — K-mer Similarity & Guide Tree (`kmer_sim.py`)**

- Extract overlapping k-mers from AA sequences
- Compute **Z-distance**
- Greedy ordering to cluster similar sequences
- Build **right-branching guide tree**
- Remove stop codons (`*`)

Sequence:

MAKTL

2-mers expected:

["MA", "AK", "KT", "TL"]

#### **Step 5 — Progressive NW Alignment (`needleman_wunch.py`)**

- Build DP matrix
- Score with match/mismatch/gap
- Traceback to aligned arrays
- Insert matching gaps into earlier sequences
- Build consensus nodes

Output: 

String-named terminals remain unaligned 

"A" → "ACG"
"B" → "AG"

Numeric node alignment (used internally)

[{'A'}, {'C'}, {'G'}]
[{'A'}, {'-'}, {'G'}]

#### Step 6 — Codon Back-Translation (`back_translate.py`)

AA alignment:

M-A

ORF:

ATGAAA

Output:

ATG---AAA

## Step 7 — Codon Position Statistics (`codon_stats.py`)

Computes:

- Identity at codon positions 1, 2, 3
- GC content per position
- Global identity
- Per-codon identity
- Writes `codon_stats.csv`
  - containing identity + GC content
  - codon_index, identity_pos1, identity_pos2, identity_pos3,
    identity_codon_total, gc_pos1, gc_pos2, gc_pos3, identity_global

#### Step 8 — Plotting (`stats_graph.py`) 

Generates:

- Identity skyline plot
- Positional identity plot
- GC position plot
- Pairwise identity heatmap

#### 🧪 **Testing & Validation**

The pipeline includes unit  tests all located in:

Final_Project_Prog2/tests/

To run all tests:

pytest -v

Testing covers:

- FASTA loading
- ORF detection
- Translation
- K-mer similarity
- Guide tree
- Progressive NW
- Gap patching
- Back-translation
- Codon statistics

#### 📤 **Outputs Generated**

Inside the `results/` folder:

results/
├── aa/                     # aligned amino acids
├── nt_backtranslated/      # aligned codon-aware nucleotides
├── stats/                  # codon_stats.csv
├── plots/                  # PNG figures
├── logs/                   # pipeline logs
└── input/                  # copy of user input

📚 **References**

I used AI to assist with code structuring, pseudocode generation, debugging, ORF detection logic, k-mer similarity design, Needleman–Wunsch alignment structure, codon back-translation design, and unit test creation.

#### 🎓 **Author**

**Jakes Njuguna**
Master of Science in Bioinformatics
University of North Carolina at Charlotte
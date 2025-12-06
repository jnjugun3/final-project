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

conda create -n msa_env python=3.10
conda activate msa_env

conda install biopython numpy pandas matplotlib pytest


### ▶️ **Running the Pipeline**

Run : 

python3 main.py -i test.fasta  

This will automatically create:

results/<run_name>/

## Command-Line Arguments

| Flag / Argument      | Description                                   | Default        |
|----------------------|-----------------------------------------------|----------------|
| `-i`, `--input`      | Input FASTA file (required)                   | —              |
|----------------------|-----------------------------------------------|----------------|
| `--run-name`         | Optional name for results folder              | `<input_name>` |
|----------------------|-----------------------------------------------|----------------|
| `--min-orf-length`   | Minimum ORF length to keep                    | 30             |
|----------------------|-----------------------------------------------|----------------|
| `-t`, `--table`      | NCBI translation table ID                     | 1              |
|----------------------|-----------------------------------------------|----------------|
| `-k`, `--kmer`       | k-mer size for amino-acid similarity          | 3              |
|----------------------|-----------------------------------------------|----------------|
| `--match`            | Needleman–Wunsch match score                  | 1              |
|----------------------|-----------------------------------------------|----------------|
| `--mismatch`         | Needleman–Wunsch mismatch penalty             | -1             |
|----------------------|-----------------------------------------------|----------------|
| `--gap`              | Needleman–Wunsch gap penalty                  | -1             |


#### See help:

python3 main.py --help

## 🧵 **Pipeline Overview**

FASTA → ORF Detection → AA Translation → k-mer Similarity
     → Guide Tree → Progressive NW → Codon Back-Translation
     → Codon Stats → Plots & CSVs

### 🧵  Full Pipeline Description**

##### **Step 1 — Loading FASTA (`load_file.py`)**

- Reads nucleotide FASTA files
- Ensures every header starts with `>`
- Ensures only valid IUPAC nucleotide characters
- Converts invalid characters → `"N"`
- Prevents duplicate sequence IDs
- Print warnings when needed

Example 
 
Input FASTA:

>seq1
ATGCZ

Output:

{"s1": "ATGCN"}

#####  **Step 2 — ORF Detection (`orf_detect.py`)**

For each input sequence:

- Scan all 3 reading frames 
- Scan **forward**, **reverse**, **reverse complement**
- Identify ORFs that start with `"ATG"` and end with `TAA/TAG/TGA`
- Extract **longest ORF per sequence**

Example

Input:

ACCATGAAATAA

Output:

{"sequence": "ATGAAATAA", "strand"}

#### **Step 3 — Amino-Acid Translation (`translate_aa.py`)**

- Translate each ORF using Biopython’s Seq.translate()
- Use user-selected NCBI translation table
- Keep "*" stop codons (handled later in back-translation)
- Return dictionary {seq_id : AA_string}

Example

Input:

"ATGGCTTAA"

Output:

"MA*"

#### **Step 4 — K-mer Similarity & Guide Tree (`kmer_sim.py`)**

- Remove terminal "*" from amino-acid sequences
- Extract overlapping k-mers (ex: MAKTL → MA, AK, KT, TL)
- Compute Z-distance between sequences
- Perform greedy ordering to cluster similar sequences
- Build a right-branching guide tree used for progressive alignment

Example 

Sequence:

MAKTL

3-mers extracted::

["MAK", "AKT", "KTL"]

#### **Step 5 — Progressive NW Alignment (`needleman_wunch.py`)**

- Align pairs of amino-acid sequences using the Needleman–Wunsch algorithm
- Use match / mismatch / gap scores to fill the DP matrix
- Trace back to get two aligned sequences (with gaps)
- Build a consensus sequence from the aligned pair
- Add the same gaps to all previously aligned sequences so every sequence stays the same length
- Repeat until all sequences are aligned

Example:

A → "ACG"
B → "AG"


[{'A'}, {'C'}, {'G'}]
[{'A'}, {'-'}, {'G'}]

#### **Step 6 — Codon Back-Translation (`back_translate.py`)**

- Map each aligned amino acid to the next codon in its original ORF
- Convert alignment gaps (-) into codon-sized gaps (---)
- Preserve the exact alignment structure produced in Step 5
- Advance to the next ORF codon only when the AA is not a gap

Example:

AA alignment:
M-A

Original ORF:
ATGAAA

Back-translated output:
ATG---AAA


## **Step 7 — Codon Position Statistics (`codon_stats.py`)**

- Take the codon-aligned nucleotide sequences from Step 6
- Split the alignment into codons (groups of 3 nt)
- Compute percent identity at positions 1, 2, and 3
- Compute average codon identity (mean of pos1 + pos2 + pos3)
- Compute GC content at positions 1, 2, and 3
- Compute global identity across the entire nucleotide alignment
- Save all results into codon_stats.csv
The CSV file contains:
codon_index, identity_pos1, identity_pos2, identity_pos3, 
identity_codon_total, gc_pos1, gc_pos2, gc_pos3, identity_global


#### **Step 8 — Plotting (`stats_graph.py`)** 

- Read the values from codon_stats.csv
- Plot identity at codon positions 1, 2, and 3
- Plot the combined codon identity across the whole alignment
- Plot GC content at positions 1, 2, and 3
- Create a pairwise percent-identity heatmap using the nucleotide alignment
- Save all plots as PNG images in the plots/ folder

#### 🧪 **Testing & Validation**

The pipeline includes unit  tests all located in:

final-project/tests/

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

Inside the results/<input_name>/ folder:

results/
├── aa/                 # Amino-acid alignments
├── input/              # Copy of input FASTA
├── logs/               # Run logs and parameters
├── nt_backtranslated/  # Codon-aware nucleotide alignments
├── plots/              # All generated plots
└── stats/              # codon_stats.csv and statistics


📚 **References**

I used AI to assist with refining Python code. All debugging and modifications, were done by me.

#### 🎓 **Author**

**Jakes Njuguna**
Master of Science in Bioinformatics
University of North Carolina at Charlotte
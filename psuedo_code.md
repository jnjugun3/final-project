## Final Project Goal

To design and implement a modular Python pipeline for **translation-based multiple sequence alignment (MSA)**.

The pipeline must:

- Detect **open reading frames (ORFs)** from nucleotide sequences
- Translate ORFs into **amino-acid sequences**
- Align amino acids using **k-mer–guided progressive alignment** (via the Needleman–Wunsch algorithm)
- **Back-translate** the amino-acid alignment into a **codon-aware nucleotide MSA**
- Compute **codon-position statistics** and generate **visualizations**
- Be fully runnable from the **terminal** via clear command-line interfaces

The code is designed to:

- Follow **FAIR principles** (Findable, Accessible, Interoperable, Reusable)
- Include a **comprehensive unit-test suite**
- Handle **ambiguous or imperfect input** gracefully, with informative warnings and errors

#### **Step 1 — Load FASTA File (`load_file.py`)**

- Open the FASTA file.
- Read every line.
- If the file is empty → raise an error.
- For each line:
  - If it starts with `>` → it is a header.
    - Store the header name (without `>`).
    - If this header already exists → raise an error.
  - Otherwise → it is a sequence line.
    - Remove whitespace.
    - Convert all characters to uppercase.
    - For each nucleotide:
      - If it is a valid IUPAC base → keep it.
      - If it is invalid → print a warning and convert it to `'N'`.
- After reading all lines:
  - Save the final sequence for the last header.
- Return a dictionary:
  **{ sequence_id : cleaned_nucleotide_string }**

#### **Step 2 — Detect ORFs (`orf_detect.py`)**

#### For each sequence:

- Convert sequence to uppercase.
- Generate:
  - forward strand
  - reverse (string reversed)
  - reverse-complement (Biopython)
- For **each strand**:
  - For each frame (0, 1, 2):
    - Start scanning codons 3 nt at a time.
    - When you see `"ATG"` → mark an ORF start.
    - Continue scanning until a stop codon (`TAA`, `TAG`, `TGA`) is found.
    - Extract that ORF.
    - If its length ≥ minimum length:
      - Record the ORF with:
        - start
        - end
        - frame
        - strand
        - sequence
- Combine ORFs from all strands.
- Return full list of ORFs.

#### Extract longest ORF per sequence:

- For each sequence ID:
  - Run ORF detection.
  - If no ORFs found → skip.
  - Choose the ORF with the longest nucleotide length.
  - Store it.

Return: **{ seq_id : longest_orf_sequence }**

------

#### **Step 3 — Translate ORFs to Amino Acids (`translate_aa.py`)**

For each sequence ID in longest_orfs:

- Convert nucleotide string → Biopython `Seq` object.
- Translate using `Seq.translate(table=<user_input>)`.
- Do **not** stop at the first stop codon.
- Keep `"*"` characters for stops.
- Store AA sequence in dictionary.

Return: **{ seq_id : amino_acid_sequence }**

#### **Step 4 — K-mer Similarity & Guide Tree (`kmer_sim.py`)**

#### Convert AA sequences → k-mers

- Remove `"*"` characters.
- Slide a window of size `k` across the sequence.
- Extract k-mers (overlapping).
- Store: **{ seq_id : list_of_kmers }**

#### Compute pairwise Z-distance

For every unordered pair:

- Count frequency of each k-mer in sequence A.
- Count frequency in sequence B.
- Z = sum of absolute differences in frequencies.
- Store distance in dictionary.

#### Order sequences by similarity

- Find the closest pair = smallest Z-distance.
- Start ordered list with these two.
- For remaining sequences:
  - For each remaining sequence:
    - Compute its distance to all sequences already in the ordered list.
    - Pick the one with lowest distance.
  - Append it to the ordered list.

#### Build a right-branching guide tree

- Terminals:
  `{1: seq1, 2: seq2, 3: seq3, ...}`
- Internal nodes:
  - Node `n+1` joins terminal 1 and 2.
  - Node `n+2` joins previous node with terminal 3.
  - Repeat.

Return:
`ordered_ids, terminals_dic, internal_nodes`

#### **Step 5 — Progressive Needleman–Wunsch Alignment (`needleman_wunch.py`)**

#### Initialize starting sequences

- Convert each AA sequence into **array-of-sets**:
  `"A"` → `{"A"}`, gap → `{"-"}`.

#### For each internal node in the tree:

1. Retrieve children A and B (numeric node IDs).
2. Build a DP matrix using `init_mat`.
3. Fill matrix using match / mismatch / gap scoring.
4. Trace back to produce:
   - alignedA (array-of-sets)
   - alignedB (array-of-sets)
5. Save aligned versions:
   - `sequences[childA] = alignedA`
   - `sequences[childB] = alignedB`
6. Build a consensus sequence:
   - For each column → union of sets in A and B
7. Patch gaps into all previously aligned sequences:
   - Look for gap positions in alignedA
   - Insert matching gaps into earlier sequences
8. Add nodes to alignment history.
9. Check that all aligned sequences now have same length.

#### After finishing tree:

- Convert terminal sequences ONLY for string keys (e.g., `"A"`, `"B"``).
- Numeric nodes remain internal and unused.

Return:
**{ seq_name : aligned_aa_string }**

#### **Step 6 — Back-Translate AA Alignment (`back_translate.py`)**

For each aligned AA sequence:

- Get its original ORF nucleotides.
- Set pointer at codon start.
- For each AA character:
  - If AA is:
    - `"-"` → append `"---"`
    - `"*"`:
      - If drop_terminal_stop=True → `"---"`
      - Else → append actual stop codon
    - any other AA → append next codon (`nt[pos:pos+3]`) and advance pointer
- Join codons together.

Return:
**{ seq_id : codon-aware_aligned_nt_sequence }**

#### **Step 7 — Codon Statistics (`codon_stats.py`)**

Given aligned nucleotide sequences:

- Ensure all have same length.
- Length must be divisible by 3.
- For each codon index:
  - Extract position 1 bases from all sequences.
  - Extract position 2 bases.
  - Extract position 3 bases.
  - Compute:
    - percent identity at each position
    - combined codon identity (mean of 3 positions)
    - GC fraction (per position)
- Compute global identity across entire alignment.
- Write all values to `codon_stats.csv`.

#### **Step 8 — Plotting (`stats_graph.py`)**

Based on CSV + aligned nucleotides:

- Plot identity at positions 1, 2, 3.
- Plot codon total identity (skyline).
- Plot GC content at positions 1, 2, 3.
- Compute pairwise identity matrix.
- Create heatmap.
- Save PNG files.

#### **Step 9 — Main Execution (`main.py`)**

Main program:

1. Parse command-line arguments.
2. Create results directory structure.
3. Step 1: Load FASTA.
4. Step 2: Detect longest ORFs.
5. Step 3: Translate ORFs.
6. Step 4: Compute k-mer similarity & build guide tree.
7. Prepare sequences for NW alignment.
8. Step 5: Run progressive NW.
9. Step 6: Back-translate AA alignment.
10. Step 7: Compute codon statistics.
11. Save outputs.
12. If plotting enabled → Step 8 plots.
13. Print completion message.
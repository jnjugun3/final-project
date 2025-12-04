# Needleman–Wunsch Alignment and Mutation Rates: A Demonstration

This project implements a **Needleman–Wunsch global sequence alignment** in Python 3 using a modular structure that includes separate components for reading FASTA files, initializing and filling scoring matrices, performing traceback, writing aligned sequences, and analyzing mutation/indel rates using **pandas**.

The code is designed to run inside a **Conda environment** named **`bioenv`** (defined in `environment.yml`) that ensures consistent dependency versions.

## Features

- Reads nucleotide sequences from FASTA input.
- Performs pairwise progressive alignment using the **Needleman–Wunsch** algorithm.
- Generates a final multiple sequence alignment in FASTA format.
- Computes **per-site and overall mismatch/indel rates** using `pandas`.
- Uses **NumPy** for efficient numerical operations and **Biopython** for FASTA handling.

## Requirements

- **Python 3.9+**
- **Conda** package manager (from [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html))

Dependencies are handled automatically through the provided `environment.yml` file and include:

- `numpy`
- `pandas`
- `biopython`

## Environment Setup

To reproduce the exact environment used for development:

1. **Clone or copy** this project folder to your local system.

2. **Create the Conda environment** from the provided `environment.yml` file:

   ```bash
   conda env create -f environment.yml
   ```

3. **Activate the environment**:

   ```bash
   conda activate bioenv
   ```

4. **Verify installation** (optional):

   ```bash
   python --version
   conda list
   ```

------

## Running the Program

The main driver script is `main.py`.

Once your environment is active, execute the program directly from the terminal:

```bash
python3 main.py
```

Alternatively, you can specify input and output files explicitly:

```bash
python3 main.py -i sequences.fasta -o alignment.fasta -r s1
```

### Command-line options

| Option            | Description                                              | Default           |
| ----------------- | -------------------------------------------------------- | ----------------- |
| `-i, --input`     | Path to input FASTA file with unaligned sequences        | `sequences.fasta` |
| `-o, --output`    | Path to output FASTA alignment file                      | `alignment.fasta` |
| `-r, --reference` | Sequence ID used as reference for mutation-rate analysis | `s1`              |
| `--match`         | Match score                                              | `1`               |
| `--mismatch`      | Mismatch penalty                                         | `-1`              |
| `--indel`         | Insertion/deletion penalty                               | `-1`              |

To view all options:

```bash
python3 main.py --help
```

------

## Example Execution

Example input FASTA file (`sequences.fasta`):

```
>s1
GATACAGATACAGAGATACA
>s2
GATACAGATTAGAGATACA
>s3
GATACACTCTAGAGATACA
>s4
GATACACCCTACAGATACA
>s5
GATACACACTACAGATACA
```

Example command:

```bash
conda activate bioenv
python3 main.py
```

### Example Output (truncated)

```
Alignment:
>s1
GATACAGATACAG-AGATACA
>s2
--GATACAGATTA-AGATACA
>s3
--GATACACTCTA-AGATACA
>s4
--GATACACCCTACAGATACA
>s5
--GATACACACTACAGATACA
Alignment saved into alignment.fasta
//
Mutation rates per position:
 Position  MismatchRate  IndelRate
        1          0.00        1.0
        2          0.00        1.0
        3          1.00        0.0
        ...
//
Average mismatch rate per position: 0.3214
Average indel rate per position: 0.1190
```

## Output Files

After running, the script produces:

- **`alignment.fasta`** — the aligned sequences in FASTA format.
- **Console summary** — alignment preview and per-position mutation rates.

## Notes

- The workflow is modular, allowing reuse of alignment, parsing, and analysis modules independently.
- The script expects valid FASTA-formatted nucleotide sequences of similar length (aligned or nearly aligned).
- Results are printed to standard output and can be redirected to a log file if needed.

## 📚 References

- Chapman, B., & Chang, J. (2000). Biopython: Python tools for computational biology. *ACM Sigbio Newsletter*, *20*(2), 15-19.
- Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & De Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, *25*(11), 1422.
- Gupta, P., & Bagchi, A. (2024). Introduction to pandas. In *Essentials of python for artificial intelligence and machine learning* (pp. 161-196). Cham: Springer Nature Switzerland.
- Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of molecular biology*, *48*(3), 443-453.
- Van Der Walt, S., Colbert, S. C., & Varoquaux, G. (2011). The NumPy array: a structure for efficient numerical computation. *Computing in science & engineering*, *13*(2), 22-30.

## Author

Dr. Denis Jacob Machado is an Assistant Professor at the Department of Bioinformatics and Genomics of the University of North Carolina at Charlotte.
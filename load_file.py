#!/usr/bin/env python3

import sys
from Bio.Seq import Seq

def load_fasta():
    fasta_file = input("Enter FASTA file path: ").strip()
    
    try: 
        with open(fasta_file, "r") as fh:
            lines = fh.readlines()
    except FileNotFoundError:
        print("Error: The file does not exist.")
        sys.exit(1)
        
    if len(lines) == 0:
        print("The input file is empty. Please input a valid FASTA file.")
        sys.exit(1)
        
    seq_dict = {}                 # store sequences
    valid_bases = "ATCG"          # allowed bases
    
    current_header = None         # track current header
    current_seq = []              # accumulate sequence lines
    header_found = False          # check if at least one header exists

    # parse fasta file
    for line in lines:
        line = line.strip()

        # check for header
        if line.startswith(">"):
            header_found = True

            # save previous sequence
            if current_header and current_seq:
                seq = "".join(current_seq)
                seq_dict[current_header] = seq

            header = line[1:].strip()

            # Check for duplicate IDs
            if header in seq_dict:
                print(f"Duplicate header found ({header}). Skipping this sequence.")
                current_header = None
                current_seq = []
                continue

            # reset sequences
            current_header = header
            current_seq = []
            continue   # move on to next line

        # if sequence appears before a header
        if current_header is None:
            print("Error: Found sequence data without a FASTA header '>'.")
            sys.exit(1)

        # clean and uppercase sequence
        clean_seq = ""
        for char in line:
            if char != " ":
                clean_seq += char.upper()

        # validate sequences
        skip_sequence = False
        for base in clean_seq:
            if base not in valid_bases:
                print(f"Invalid characters '{base}' found in sequence for {current_header}. Skipping this sequence.")
                skip_sequence = True
                break

        # skip invalid sequences   
        if skip_sequence:
            current_header = None
            current_seq = []
            continue
        
        # append cleaned up sequences to the list
        current_seq.append(clean_seq)

    # save last sequence if exists
    if current_header and current_seq:
        seq = "".join(current_seq)
        seq_dict[current_header] = seq

    if not header_found:
        print("No FASTA header was found. Please check file format.")
        sys.exit(1)
        
    if len(seq_dict) == 0:
        print("No valid sequences found. Exiting pipeline.")
        sys.exit(1)

    print("FASTA file correctly uploaded and all sequences are validated.")
    return seq_dict


# call step 1 in pipeline
if __name__ == "__main__":
    seq_dict = load_fasta()
    print(seq_dict)

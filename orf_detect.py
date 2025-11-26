#!/usr/bin/env python3

import sys
from Bio.Seq import Seq

def detect_orfs(seq_dict):
    # ORF detection parameters
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    min_length = 30
    
     # dictionary to store all ORFs
    orf_dict = {}
    
	longest_orfs = {}
   
	forward_frames = [0, 1, 2]
	reverse_frames = [0, 1, 2]

	for seq_id, nuc_seq in seq_dict.items():
		forward = Seq(nt_seq)
 		reverse = forward.reverse_complement()
 		
		all_orfs = []
		
		strands = {
            "+": seq,
            "-": rev_seq
            }
	for strand_symbol, seq_obj in strands.items():
            seq_len = len(seq_obj)

            # Reading frames 0, 1, 2
            for frame in range(3):
                i = frame  # starting index of frame

                while i < seq_len - 2:
                    codon = seq_obj[i:i+3]

                    # START codon found
                    if codon == start_codon:
                        start_pos = i
                        j = i + 3

                        # Extend ORF until STOP codon
                        while j < seq_len - 2:
                            next_codon = seq_obj[j:j+3]

                            if next_codon in stop_codons:
                                stop_pos = j + 3
                                orf_len = stop_pos - start_pos

                                # Check minimum length requirement
                                if orf_len >= min_length:
                                    orf_entry = {
                                        "strand": strand_symbol,
                                        "frame": frame,
                                        "start": start_pos,
                                        "stop": stop_pos,
                                        "nt_sequence": str(seq_obj[start_pos:stop_pos])
                                    }
                                    all_orfs.append(orf_entry)

                                break  # Stop scanning this ORF when stop codon found

                            j += 3

                        i = j  # Jump to after the ORF to continue scanning
                    else:
                        i += 3

        # Store results for this sequence
        orf_dict[seq_id] = all_orfs

        # Find the longest ORF in this sequence
        if all_orfs:
            longest = max(all_orfs, key=lambda x: len(x["nt_sequence"]))
            longest_orf_dict[seq_id] = longest
        else:
            longest_orf_dict[seq_id] = None

    print("ORF detection complete.")
    return orf_dict, longest_orf_dict
    
    
    
    
    
    
    
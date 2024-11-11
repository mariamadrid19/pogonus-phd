#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def filter_fasta(fasta_file, ids_file, output_file):
    # Read IDs to be filtered from the text file into a set
    with open(ids_file, "r") as f:
        ids_to_filter = set(line.strip() for line in f if line.strip())

    # Parse the input FASTA file
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    # Filter sequences by excluding those with IDs in ids_to_filter
    filtered_sequences = [seq for seq in sequences if seq.id not in ids_to_filter]

    # Write the filtered sequences to a new FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(filtered_sequences, output_handle, "fasta")

    # Check for any IDs that were in the list but not in the FASTA headers
    found_ids = {seq.id for seq in sequences}
    missing_ids = ids_to_filter - found_ids

    if missing_ids:
        print("Warning: The following IDs were not found in the FASTA file:")
        for missing_id in missing_ids:
            print(missing_id)
    else:
        print("All specified IDs were found and filtered.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter sequences from a FASTA file based on a list of IDs.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file")
    parser.add_argument("ids_file", help="Path to the text file with IDs to remove")
    parser.add_argument("output_file", help="Path to the output filtered FASTA file")

    args = parser.parse_args()
    filter_fasta(args.fasta_file, args.ids_file, args.output_file)

#!/usr/bin/env python3
#
# This script generates a nucleotide fasta sequence of a requested length.  
import argparse
import random


__author__ = "Ekaterina Osipova, 2024."


def generate_random_fasta(sequence_length):
    # Define the possible nucleotides
    nucleotides = ['A', 'T', 'C', 'G']
    fasta = ""
    
    # Generate a random sequence
    sequence = ''.join(random.choice(nucleotides) for _ in range(sequence_length))
    
    # Wrap sequence lines to a maximum of 80 characters per line (FASTA format standard)
    for i in range(0, len(sequence), 80):
        fasta += sequence[i:i+80] + "\n"
    
    return fasta


def main():
    parser = argparse.ArgumentParser(description="Simulate nucleotide fasta")
    parser.add_argument("--length", "-l", type=int, default=100, help="Sequence length to generate; default: 100")
    args = parser.parse_args()

    fasta = generate_random_fasta(args.length)
    print('>random_seq')
    print(fasta)


if __name__ == "__main__":
    main()
#!/usr/bin/env python
# coding: utf-8

# import required modules and packages
from Bio import SeqIO

# read the protein sequence
#seq = [str(record.seq) for record in SeqIO.parse("sequence.fasta", "fasta")][0]
sequence = input("Enter the protein sequence filename: ")
seq = [str(record.seq) for record in SeqIO.parse(sequence, "fasta")][0]
# length of the protein sequence
N = len(seq)
print(f"Length of Sequence: {N}\n{seq}")

# Count Lysine (K)
k = seq.count("K")
# calculate percentage of K
k_pcnt = round((k/N)*100, 2)

# Count Aspartic Acid (D)
d = seq.count("D")
# calculate percentage of D
d_pcnt = round((d/N)*100, 2)

print(f"Percentage of Lysine (K): {k_pcnt}%")
print(f"Percentage of Aspartic Acid (D): {d_pcnt}%")

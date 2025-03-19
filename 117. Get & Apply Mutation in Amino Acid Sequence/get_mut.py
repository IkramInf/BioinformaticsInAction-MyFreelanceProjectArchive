#!/usr/bin/env python
# coding: utf-8

# import required libraries
from Bio import pairwise2  # perform pairwise alignment
from Bio.Align import substitution_matrices  # contains blosum62 matrix

# take inputs from user
sequence1 = input("Enter the first amino acid sequence: ")
sequence2 = input("Enter the second amino acid sequence: ")

# for testing purpose
#sequence1 = "ASTVRHIL"
#sequence2 = "ASTVREPEHVL"

# get blosum62 matrix
blosum62 = substitution_matrices.load("BLOSUM62")
# perform pairwise alignment between two input sequences
alignment = pairwise2.align.globalds(sequence1, sequence2, blosum62, -10, -10)[0]
#print(f"{alignment.seqA}\n{alignment.seqB}")

mutations = []  # this will contain all mutations
# iterate over two aligned sequence together
i = 0
for x, y in zip(alignment.seqA, alignment.seqB):
    # x: first sequence base
    # y: second sequence base
    # if the two bases are different
    if x != y:
        # if first base is deleted...
        if x == "-":
            mutations.append(f"ins{i}{y}")
        # if insertion in second base
        elif y == "-":
            mutations.append(f"del{i+1}")
        # point mutation
        else:
            mutations.append(f"{x}{i+1}{y}")
       
    # increment the index value
    if x != "-":
        i += 1
            
# display the list of mutations
print(" ".join(mutations))
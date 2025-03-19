#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
from random import randint

class FastAreader:
    def __init__(self, fname=''):
        """contructor: saves attribute fname"""
        self.fname = fname
        self.fileH = None
        
    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
    
    def readFasta(self):
        records = []
        with self.doOpen() as self.fileH:
            header = ''
            sequence = ''
            # skip to first fasta header
            line = self.fileH.readline()
            while not line.startswith('>'):
                line = self.fileH.readline()
            header = line[1:].rstrip()
            
            for line in self.fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
        records.append((header, sequence))
        return records

def Count(Motifs):
    transpose = [*zip(*Motifs)] # transpose the motif matrix
    count = {nucleotide: [index.count(nucleotide) for index in transpose] for nucleotide in "ATCG"}
    return count

def Motifs(Profile, dna):
    return [ProfileMostProbableKmer(text,len(Profile["A"]),Profile) for text in dna]

def Pr(sequence, profile):
    """gives the probability that sequence fits the profile"""
    probability = 1.0
    for position, base in enumerate(sequence):
        probability *= profile[base][position]
    return probability

def Consensus(Motifs):
    matrix = Count(Motifs)
    consensus = ""
    for i in range(len(Motifs[0])):
        consensus += max(matrix, key=lambda e: matrix[e][i])
    return consensus

def HammingDistance(dna_str_1, dna_str_2):
    return sum([1 for x, y in zip(dna_str_1, dna_str_2) if x != y])

def ProfileMostProbableKmer(text, k, profile):
    MostProbable = text[0:k]
    p = Pr(text[0:k], profile)
    for i in range(len(text)-k+1):
        if p < Pr(text[i:i+k], profile):
            p = Pr(text[i:i+k], profile)
            MostProbable = text[i:i+k]
    return MostProbable

def Score(motifs):
    consensusString = Consensus(motifs)
    return sum([HammingDistance(consensusString, line) for line in motifs])

def CountWithPseudocounts(Motif):
    t = len(Motif)
    k = len(Motif[0])
    pseudocount = Count(Motif)
    for symbol in "ACGT":
        for j in range(k): #For each symbol it loops over the length of the motif
            pseudocount[symbol][j] += 1
    return(pseudocount)

def ProfileWithPseudocounts(Motifs):
    profile = CountWithPseudocounts(Motifs)
    t = len(Motifs)
    k = len(Motifs[0])
    return {k: [i/(t+4) for i in v] for (k, v) in profile.items()}

def RandomMotifs(dna, k, t):
    motifs = []
    for i in range(t):
        ind = randint(0, len(dna[0]) - k)
        motif = dna[i][ind:ind + k]
        motifs.append(motif)
    return motifs

def RandomizedMotifSearch(dna, k, t):

    M = RandomMotifs(dna, k, t)
    BestMotifs = M

    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Randomized Search")

    parser.add_argument('-i', '--iteration', type=int, default=1000, help='Number of iterations')
    parser.add_argument('-k', '--ks', type=int, default=13, help='kmer size')
    parser.add_argument('-p', '--ps', type=int, default=1, help='profile size')
    parser.add_argument('input', help="input filename")
    parser.add_argument('output', help="output filename")
    args = parser.parse_args()
    
    fr = FastAreader(args.input)
    records = list(fr.readFasta())
    records = {header: seq for header, seq in records}

    results = "\n".join(RandomizedMotifSearch(list(records.values()), args.ks, args.ps))
    with open(args.output, "w") as f:
        f.write(results)
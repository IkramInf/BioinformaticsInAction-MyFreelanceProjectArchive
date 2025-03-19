# import required libraries
from itertools import combinations
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator

# extract all sequences from 'original_sequences_exercise.fasta' into a dictionary
seqs = {record.id:record.seq for record in SeqIO.parse("original_sequences_exercise.fasta", "fasta")}

# calculate pairwise alignment scores
aligner = PairwiseAligner(mode = 'global')
scores = {}
sequences = []
for pair in combinations(seqs.keys(), 2):
    alignments = aligner.align(seqs[pair[0]], seqs[pair[1]])
    optimal = next(alignments)
    scores.update({pair : optimal.score})

# calculate distance matrices
dm = []
for k1, v1 in scores.items():
    temp = []
    for k2, v2 in scores.items():
        if k1 == k2:
            temp.append(0)
        else:
            temp.append(abs(v1-v2))
    dm.append(temp)
    
print(list(scores.values()))
print(dm)
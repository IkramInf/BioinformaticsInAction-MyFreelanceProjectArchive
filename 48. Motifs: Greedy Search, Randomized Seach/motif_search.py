# import randint from random module
from random import randint

# input for first exercise
dna = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
k = 5
profile = [[0.2, 0.2, 0.3, 0.2, 0.3],
            [0.4, 0.3, 0.1, 0.5, 0.1],
            [0.3, 0.3, 0.5, 0.2, 0.4],
            [0.1, 0.2, 0.1, 0.1, 0.2]]


# calculate probability with kmer and profile matrix
def calc_probalility(kmer, profile):
    # if profile is passed as list, convert it into dict
    if isinstance(profile, list):
        profile = {base : count for base, count in zip('ACGT', profile)}
    # let initial probability as 1    
    probability = 1.0
    # iterate over each base of kmer and 
    # multiply the corresponding bases's position-wise value with initial probability
    for position, base in enumerate(kmer):
        probability *= profile[base][position]
    return probability

# find out most probable kmer in a dna sequence
def mostProbableKmer(dna, k, profile):
    kmers = {}
    for i in range(len(dna)-k+1):
        kmer = dna[i:i+k]
        kmers.update({kmer : calc_probalility(kmer, profile)})
    kmers = dict(sorted(kmers.items(), key=lambda x:x[0].lower()))
    return sorted(kmers.items(), key=lambda x:x[1], reverse=True)[0][0]


print("Profile-most probable k-mer is ", mostProbableKmer(dna, k, profile))


# inputs for second exercise
k, t = 3, 5
dna = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']


# calculate hamming distance between two sequences
def hammingDistance(dna1, dna2):
    return sum([1 for x, y in zip(dna1, dna2) if x != y])


# calculate score of motifs
def score(motifs):
    # count bases as position-wise for motifs
    counts = {base: [bases.count(base) for bases in [*zip(*motifs)]] for base in "ACGT"}
    # create a variable to store consensus sequence
    consensus = ""
    # find out consensus sequence among all motifs
    for i in range(len(motifs[0])):
        consensus += max(counts, key=lambda x: counts[x][i])
    # iterate over all motifs, calculate hamming distance with consensus sequence and return sum of them
    return sum([hammingDistance(consensus, motif) for motif in motifs])


# calculate profile from motifs
def Profile(motifs):
    L = len(motifs)
    return {base: [bases.count(base)/L for bases in [*zip(*dna)]] for base in "ATCG"}


# Find out motif with greedy motif search
def GreedyMotifSearch(dna, k, t):
    # create a list to store best motif
    best_motif = []
    # store k-length bases from each motif into best_motif variable
    for i in range(0, t):
        best_motif.append(dna[i][0:k])
    # iterate over dna and find out most probable kmers
    for i in range(len(dna[0])-k+1):
        motifs = []
        motifs.append(dna[0][i:i+k])
        for j in range(1, t):
            profile = Profile(motifs[0:j])
            motifs.append(mostProbableKmer(dna[j], k, profile))
        # compare between score     
        if score(motifs) < score(best_motif):
            best_motif = motifs
    return best_motif


print("Best motifs by Greedy Motif Search are : ", GreedyMotifSearch(dna, k, t))


# Exercise 3

# count motifs with pseudocount
def profileWithPseudocounts(motifs):
    t, k = len(motifs), len(motifs[0])
    # count bases
    pseudocount = {base: [bases.count(base) for bases in [*zip(*motifs)]] for base in "ACGT"}
    for symbol in "ACGT":
        for j in range(k):
            # add 1 at position j
            pseudocount[symbol][j] += 1
    # return pseudocount by dividing each element by t
    return {base: [count/t for count in counts] for base, counts in pseudocount.items()}


# generate motifs randomly from dna list
def randomMotifs(dna, k, t):
    motifs = []
    for i in range(t):
        # randomly generate a index r
        r = randint(0,len(dna[0])-k)
        # append motif from r to r+k from dna
        motifs.append(dna[i][r:r+k])
    return motifs


# search randomly to find out motifs
def RandomizedMotifSearch(dna, k, t):
    # generate random motifs
    motifs = randomMotifs(dna, k, t)
    # keep the motifs into best_motif varibale
    best_motif = motifs

    while True:
        # create profile with pseudo count
        profile = profileWithPseudocounts(motifs)
        # find out most probable motifs
        motifs = [mostProbableKmer(kmer,len(profile["A"]),profile) for kmer in dna]
        # compare score
        if score(motifs) < score(best_motif):
            best_motif = motifs
        else:
            return best_motif



print("Best motifs by Randomized Search are :", RandomizedMotifSearch(dna, k, t))





#! /usr/bin/python3

from Bio import Seq, SeqIO
from Bio.Seq import MutableSeq
import sys


# Global dictionary for our kmers
kmers = {}

# Given a position determine if any high coverage kmers can be found
# which include that base.  If even one 'real' kmer using that base
# can be identified then we assume the base call is correct.
#
# Input
#   pos - the base position in the read to check
#   size - kmer size
#   threshold - the kmer count threshold to determine noise or not
#   s - the Seq object we are checking for errors
#
# Return - False is returned if a kmer count above the threshold is found. This
#          means the base is NOT an error.  True is returned if no kmers involving
#          the base have counts >= the threshold (the base IS ERRONEOUS)
#
def isErroneous(pos, size, threshold, s):
    # setup a sliding window with pos as 3' as possible
    # HINT: start = pos-size+1.  If start is less than 0, then set it to zero
    # make s as uppercase
    start = pos-size+1
    if start < 0:
        start = 0
    
    # find out kmers that contains the base, s[pos]
    base_involved_kmers = [s[i : i+size] for i in range(start, pos+1)]
    # read all kmers from test_kmers.fasta
    read_kmers = [str(record.seq) for record in SeqIO.parse("test_kmers.fasta", "fasta")]
    
    count = []
    # check for each base involved kmer with all read kmers from file and count their differences in base
    for kmer in base_involved_kmers:
        for checker in read_kmers:
            # if the base is in file read kmer, compare with that
            if s[pos] in checker:
                # count base difference
                count.append(sum([1 for x, y in zip(kmer, checker) if x!=y]))
            else:
                continue  
    # compare the count with threshold value            
    if max(count) > threshold:
        return True
    else:
        return False
        
# Given the position of a likely erroneous base, test if substituting
# the base with the other three results in a high coverage kmer
# Input
#   pos - the read position to interrogate
#   size - kmer size
#   threshold - If there are less than this many kmers supporting the base
#               then the base is an error
#
# Return: the highest rated alternative base
#
def getAlternativeBase(pos,size,threshold,s):
    bases = ['A','G','C','T']
    # There is no method to deepcopy a biopython sequence so I copied the sequence
    # by converting it to a list and later putting it back together as a string for
    # the call to isErroneous
    newseq = list(s)
    
    counts = {}
    # iterate over all 4 bases A,C,G,T and don't consider the erroneous base
    for base in bases:
        # to check if the base is erroneous
        if base != s[pos]:
            # replace erroneous base position with other base
            s[pos] = base
            # add count to dictionary as base:count
            counts[base] = isErroneous(pos,size,threshold,s)
        else:
            continue
    # return the base that has smallest base difference    
    return min(counts, key=counts.get)

# Set up the kmer dictionary from the Jellyfish dump command.  It is provided
# with this assignment.
# Dictionary - kmers
#       key - the kmer sequence
#       value - the count for the kmer
# HINT: Make sure that the kmer is cast to a string and the value is cast to an int
def setUpKmers():
    
    for record in SeqIO.parse("test_data.fastq", "fastq"):
        seq = str(record.seq)
        for i in range(len(seq)):
            kmer = seq[i : i+17]
            if len(kmer) == 17:
                kmers[kmer] = kmers.get(kmer, 0) + 1
    
    # Unique k-mers are those that appear only once
    unique = len([1 for k, v in kmers.items() if v==1])
    # distinct k-mers are counted only once, even if they appear more times
    distinct = len(kmers)
    # add all count value of kmers dictionary
    total = sum([v for k, v in kmers.items()])
    # maximum count value of kmer
    max_count = kmers[max(kmers,key=kmers.get)]
    print(f"Unique: {unique}\nDistinct: {distinct}\nTotal: {total}\nMax_count: {max_count}")


def main():
    setUpKmers()
    #kmer_size = int(sys.argv[3])
    #threshold = int(sys.argv[4])
    
    kmer_size = 17
    threshold = 0
    
    outfile = open("outfile.fastq","w")
    #with open(sys.argv[2],'r') as f:
    
    with open("test_data.fastq",'r') as f:
        seqIO = SeqIO.parse(f,'fastq')
        for sr in seqIO:
            s = MutableSeq(sr.seq)
            s = s.upper()
            for pos in range(0,len(s)):
                if isErroneous(pos, kmer_size, threshold, s):
                    newbase = getAlternativeBase(pos,kmer_size,threshold,s)
                    s[pos] = newbase
            sr.seq = s
            SeqIO.write(sr,outfile,'fastq')
    outfile.close()
    
if __name__ == "__main__":
    main()


#!/usr/bin/python3
# ---------------------------------------------------
# MINLENGTH = 100
# ---------------------------------------------------
MINLENGTH = 100
thresholds = [0, 1, 2, 3, 4, 5, 6]

# ---------------------------------------------------
# e.g.,  findpattern('ATG', sequence, 0, 100, 3, 0)
# ---------------------------------------------------
def findpattern(pattern, seq, start, end, increment, threshold=0) :
    '''
    Assumes that pattern is a string and seq is a string,
    start, end, increment, and threshold are integers.
    Finds the first match of the string pattern in seq
    when looking at positions in range(start, end, increment),
    where a match means mismatches <= threshold.
    Returns the index of that first match, or -1 if no match.
    '''
    n = len(pattern)
    for i in range(start, end, increment) :
        # try to match the pattern to seq[i:i + npattern]
        if matches(pattern, seq[i:i+n], threshold) :
            return i
    return -1
# ---------------------------------------------------
# e.g.,  findstop(sequence, 0, 100, 3)
# ---------------------------------------------------
def findstop(seq, start, end, increment) :
    '''
    Assumes seq is a string, start, end, increment are integers.
    Finds the first exact match to a stop codon in seq
    when looking at positions in range(start, end, increment).
    Returns the index of that first match, or -1 if no match.
    '''
    stops = ('TAG', 'TGA', 'TAA')
    for i in range(start, end, increment) :
        triplet = seq[i:i+3]
        if triplet in stops :
            return i
    return -1

def matches(s1, s2, threshold=0) :
    '''
    Assumes that s1 and s2 are equal length sequences.
    Returns the number of mismatches when aligning them directly.
    '''
    n1 = len(s1)
    n2 = len(s2)
    if n1 != n2 :
        print( 'Error: the input sequences are not the same length:', s1, s2)
        exit(0)
    count = 0
    for i in range(n1) :
        if s1[i] != s2[i] :
            count += 1
    return count <= threshold

def reversecomplement(seq) :
    '''
    Assumes that seq is a DNA sequence, 5' to 3'
    Returns the reverse complement, 
    which is the complementary strand but 5' to 3'
    '''
    n = len(seq)
    result = ''
    for i in range(n-1, -1, -1) :
        if seq[i] == 'A' :
            result += 'T'
        elif seq[i] == 'T' :
            result += 'A'
        elif seq[i] == 'G' :
            result += 'C'
        elif seq[i] == 'C' :
            result += 'G'
        else :
            print('We encountered a letter that is not A,T,G, or C in sequence')
            exit()
    return result

# ---------------------------------------------------
# Look for all of the ORFs in a single frame
# but only report an ORF if the length is >= minlength
# ---------------------------------------------------
def findorfs(seq, frame, minlength, threshold) :
    rejects = []
    nseq = len(seq)
    # start at position frame, which is 0, 1 or 2
    i = frame
    countORF = 0
    
    while i < nseq :
        # locate the next ATG in this frame
        atg = findpattern('ATG', seq, i, nseq - 3, 3, threshold=0)
        # a return value of -1 means none found
        if atg == -1 :
            break # since there are no more start codons in this frame
        else :
          # locate the nearest stop codon to the ATG
            stop = findstop(seq, atg + 3, nseq - 3, 3)
            # a return value of -1 means there are no stop codons
            if stop == -1 :
                break # since there are no more stop codons in this frame
            else :  # found an ORF
                # only count it if the length is at least minlength
                length = stop + 3 - atg # if we include the stop codon in the length
                if length >= minlength :
                    #for th in thresholds:
                    aggagg = findpattern('AGGAGG', seq, atg - 13, atg - 9, 1, threshold=threshold)
                    #print('aggagg = ', aggagg)
                    # if the Shine-Dalgarno pattern is found, then we'll accept this ORF
                    # as a gene
                    if aggagg > -1 :
                        AGGAGG = seq[aggagg : aggagg + 6]
                        spacer = seq[aggagg + 6 : atg]
                        countORF += 1
                        orf = seq[atg : stop + 3] # to include the stop codon
                        #print( 'ORF' + str(countORF), 'at position', atg, 'in frame', frame, 'with length', length)
                        #print('Shine-Dalgarno', AGGAGG, 'spacer', spacer)
                        #print('----------------------------------------------------------GENE:')
                        #print(orf)
                        #print('----------------------------------------------------------')

                        i = stop + 3 # continue looking after this ORF for the next one
                    else :
                        rejects.append(length) # rejected no AGGAGG
                        i = atg + 3 # look after the ATG since another ATG might have an AGGAGG
                else : # the length was too short
                    rejects.append(length) # rejected too short
                    i = stop + 3 # continue looking after the stop codon
    #print('Lengths of ORFs rejected as too short or without AGGAGG:', rejects)
    return countORF

# ---------------------------------------------------
# Step 3. 
# Look for all of the ORFs in all of the frames
# ---------------------------------------------------
def findallorfs(seq) :
#  print( seq)
    n = len(seq)
    #MINLENGTH = 90   # Setting this variable at the top as global
    gene_counts = []
    for th in thresholds:
        counts = []
        for frame in range(3) :
            #print('----------------------------------- Frame:', frame) 
            counts.append(findorfs(seq, frame, MINLENGTH, th))
        #print( '------------------------------------------')
        #print( 'Doing the reverse complement sequence')
        #print( '------------------------------------------REVERSE COMPLEMENT:')
        reverse = reversecomplement(seq)
        #print( reverse[0:10], '...', reverse[-10 : -1])
        #print('-------------------------------')
        for i in range(3) :
            #print('-------------------------------- Frame: -' + str(i + 1))
            counts.append(findorfs(reverse, i, MINLENGTH, th))
        
        gene_counts.append(sum(counts))
    #print( '------------------------------------------')
    return gene_counts
# ---------------------------------------------------------------
# Try it out
# -----------------------------------------------------------------
'''
To find genes in a given sequence 
by trying all frames forward and backward.
1. Look for ORFs
2. Restrict length to >= MINLENGTH
3. Require AGGAGG with 3 to 7 bases between AGGAGG and ATG
'''
# --------------------------------------------------------
# Data
# ---------------------------------------------------------
f = open('plasmid.txt')
line = f.readline()
#print(line[:-1])
seq = f.read()
seq = seq.replace('\n', '')
seq = seq.upper()
n = len(seq)    
#print('Sequence length:', n)
#print(seq)
#print('-------------------------------------')

print("Threshold\tNumber of genes")
gene_counts = findallorfs(seq)

for t, c in zip(thresholds, gene_counts):
    print(f"{t}\t{c}")
    
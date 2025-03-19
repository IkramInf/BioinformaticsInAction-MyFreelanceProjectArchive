#!/usr/bin/python3
import math
# --------------------------------------------------------
# Construct the substitution matrix from the training alignments
# ---------------------------------------------------------
f = open('alignments.txt')
# --------------------------------------------------------
# Construct the substitution matrix from the training alignments
# ---------------------------------------------------------
def remove(s, i) :
    '''Returns the string s after removing position i'''
    return s[0:i] + s[i+1:]

def read_training_alignments(f):
    # read the files as a list of strings
    lines = f.readlines()
    # store all sequences
    sequences = []
    for line in lines:
        # to separate the lines that contains sequence
        if not line.startswith((' ', '>', '.', ':', '|', '\n')):
            # split the line into 4 parts
            line = line.strip().split()
            # check the line has 4 parts
            if len(line) == 4:
                # as each alignment line hase 4 parts, sequence in 2 as python index
                sequences.append(line[2])
            else:
                # if not 4 parts, skip that line
                continue
            
    # create two string object to store sequence        
    seq1, seq2 = "", ""
    
    # iterate over sequences
    for i in range(len(sequences)):
        # if even lines as python index i.e. 0, 2, 4...
        if i%2 == 0:
            seq1 += sequences[i]
        else:
            seq2 += sequences[i]
    
    # remove all the gapes from seq1 and seq2
    seq1 = seq1.replace("-", "")
    seq2 = seq2.replace("-", "")
    
    # return gapless seq1 and seq2
    return seq1, seq2

# calling the read_training_alignments function with file object f
seq1, seq2 = read_training_alignments(f)   
print(seq1, len(seq1))
print(seq2, len(seq2))
#exit(0)
# ----------------------------------------------
# A function to count alignments between letter1 and letter2
# For simplicity we assume that seq1 and seq2 are global variables
# ----------------------------------------------
def count_align(letter1, letter2) :
    '''Assume we include the pseudocount of 1'''
    count = 0
    for a, b in zip(seq1, seq2) :
        if a == letter1 and b == letter2 or a == letter2 and b == letter1 :
            count += 1
    return count + 1
print('AC = ', count_align('A', 'C'))
print('GA = ', count_align('G', 'A'))
#exit(0)
# ----------------------------------------------
# Function to count occurrences of given letter in both sequences
# ----------------------------------------------
def count_letter(letter) :
    '''Assumes we include the pseudocount'''
    return (seq1+seq2).count(letter) + 1
    
print('A', count_letter('A'), 'expect 18 + 1 = 19')
print('C', count_letter('C'), 'expect 3 + 1 = 4')
#exit(0)
# ----------------------------------------------
# Probability for a letter
# ----------------------------------------------
def p(letter) :
    '''Assumes it is count/2n'''
    return count_letter(letter)/(2 * len(seq1))
    
print('pA', round(p('A'),3), 'expect 19/102 = 0.186')
print('pC', round(p('C'),3), 'expect 4/102 = 0.039')
#exit(0)
# ----------------------------------------------
# Calculate S
# ----------------------------------------------
def calculate_S(letter1, letter2) :
    '''
    S = log_2(M/E)
    M = pair alignments / n
    E = 2 * p1 * p2 or p1 ** 2
    p1 = count / (2n)
    '''
    n = len(seq1)
    M = count_align(letter1, letter2) / n
    E = p(letter1) * p(letter2)
    if letter1 != letter2 :
        E = 2 * E
    R = M / E
    return math.log(R, 2)
    

# ----------------------
# Check a couple of values
# ----------------------
S = calculate_S('A', 'A')
print('S_AA = ', round(S, 1), 'Expect 2.0, see table')
S = calculate_S('A', 'C')
print('S_AC = ', round(S, 1), 'Expect 2.4, see table')
#exit(0)

# ----------------------
# Calculate the Row of the table for A
# ----------------------
AAs = 'ACDFGNSWY'
row = []
for letter in AAs :
    S = calculate_S('A', letter)
    row.append(round(S,1))
print('Row for A = ', row)
print('Expecting: ', [2.0, 2.4, 0.4, -1.8, 2.4, -2.4, -1.3, -1.6, 0.4])
#exit(0)
# ----------------------
# Store the whole matrix as a 2D dictionary
# ----------------------
S = {}
for a in AAs :
    # ------------------------------
    # Make a dictionary for the letter, i.e., the row
    # ------------------------------
    S[a] = {}
    # ------------------------------
    # Another loop to calculate S[a][b]
    # ------------------------------
    for b in AAs :
        S[a][b] = round(calculate_S(a, b), 1)
# ----------------------
# Check a few values
# ----------------------
print(f"S[A][A] = {S['A']['A']}, Expecting 2.0")
print(f"S[A][D] = {S['A']['D']}, Expecting 0.4")
print(f"S[F][A] = {S['F']['A']}, Expecting -1.8")

print(S)
# --------------------------------------------------------
# The End
# ---------------------------------------------------------

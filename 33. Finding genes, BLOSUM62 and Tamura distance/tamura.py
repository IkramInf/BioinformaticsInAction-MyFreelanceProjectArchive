#!/usr/bin/python3
import math
# -----------------------------------------------------
# Implementing Tamura's formula from Figure 7 of images.
# -----------------------------------------------------
def tamura(s1, s2, S, V) :
    '''
    Assume that s1, s2 are two sequences of DNA without gaps,
    which are of equal length and aligned directly starting 
    at the first letter,
    and S is the number of transition mismatches in the alignment,
    and V is the number of transversion mismatches in the alignment,

    Return the value of K using the Tamura's formula.
    
    The formula is on the given image sheet.
    Put your code in the indicated region below
    '''
    # ---------------------------------------------------
    # Step 1. calculate the GC content of each string,
    # where GC content is (number of G's + number of C's) / length
    # ---------------------------------------------------
    GCs1 = float(s1.count('G')+s1.count('C'))/len(s1)*100
    GCs2 = float(s2.count('G')+s2.count('C'))/len(s2)*100
    
    # ---------------------------------------------------
    # Step 2. calculate C using the formula
    # ---------------------------------------------------
    C = (GCs1 + GCs2) - (2 * GCs1 * GCs2)
    
    # ---------------------------------------------------
    # Step 3. calculate K and return it.
    # ---------------------------------------------------
    K = -C*math.log(1-(S/C)-V) - 0.5*(1-C) * math.log(1-2*V)
    
    return K
    
# -----------------------------------------------------
# Try it out
# -----------------------------------------------------
s1 = 'TAATTCTCGTTGGATCCCACCAGCCG'
s2 = 'TAGTTCTCGATGTACGACATCGGCCG'
S = 0.2
V = 0.1
K = tamura(s1, s2, S, V)
print('S = ', S, 'V =', V, 'K = ', round(K, 2))
# -----------------------------------------------------
# The End
# -----------------------------------------------------

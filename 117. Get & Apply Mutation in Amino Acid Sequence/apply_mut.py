#!/usr/bin/env python
# coding: utf-8

# take inputs from user
sequence = list(input("Enter the amino acid sequence: "))  # e.g., AGSTV
mutations = input("Enter list of mutations separating by space: ")  # e.g., A1V T4R

# for testing 
#sequence=list("ASTVR")
#mutations="del1 del2"

# iterate over each mutation and apply them into original sequence
# keep track deleted and inserted
k = 0
for mutation in mutations.strip().split():
    # apply mutation for deletion
    if "del" in mutation.lower():  # deln
        index = int(mutation[3:])
        del sequence[index-k-1]
        k +=1
        #print(sequence)
    
    # apply mutation for insertion    
    elif "ins" in mutation.lower():  # insnaa
        e = [i for i, s in enumerate(mutation) if s.isdigit()][-1]
        index, ins_aa = int(mutation[3:e+1]), mutation[e+1:]
        for aa in ins_aa:
            sequence.insert(index-k, aa)
            #print(sequence)
            k -= 1
            index += 1
    
    # apply other mutation
    else:
        original, mut = mutation[0], mutation[-1]
        index = mutation [1:-1]
        sequence[int(index)-k-1] =  mut
        #print(sequence)
        
# display the mutated sequence
sequence = "".join(sequence)
print(sequence)

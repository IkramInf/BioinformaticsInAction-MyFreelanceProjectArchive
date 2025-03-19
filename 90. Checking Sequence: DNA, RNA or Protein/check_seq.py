# read the sequence and remove the newline characters and make it uppercase
with open("seq2.txt", "r") as f:
    sequence = f.read().replace("\n", "").upper()
# find the unique bases in sequence
bases = set(sequence)

# dna, rna and protein bases for comparison
dna = {'A', 'C', 'G', 'T'}
rna = {'A', 'C', 'G', 'U'}
protein = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

# check if the sequence is dna, rna or protein
#<TEST CASE 1> Nucleotide sequence
if bases.issubset(dna):
    print('The Given Sequence is a Nucleotide Sequence')
elif bases.issubset(rna):
    print('The Given Sequence is a Nucleotide Sequence')
    
#<TEST CASE 2> Protein sequences
elif bases.issubset(protein):
    print('The Given Sequence is a Protein Sequence')
    
#<TEST CASE 3> Neither aminoacids nor nucleotides
else:
    print('Please enter valid Input')
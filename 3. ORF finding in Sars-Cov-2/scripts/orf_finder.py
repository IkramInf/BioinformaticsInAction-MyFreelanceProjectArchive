# Finding ORFs in Sars cov2 
f = open('nc_045512.txt')
description = f.readline()
#print(description)

seq = f.read()
seq = seq.replace('\n', '')
#print(seq)
#print('\nLength of SARS-CoV-2 viral sequence is', len(seq), 'nucleotides')

def firstORF(sequence, start_position) :
    '''
    Returns the (start, stop) index of ORF in the given sequence.
    '''
    n = len(sequence)
    
    for i in range(start_position, n, 3) :
        triplet = sequence[i:i + 3]
        # Check for a start codon
        if triplet == 'ATG' :
            # Look at triplets up to a stop codon
            for j in range(i, n, 3) :
                triplet = sequence[j:j + 3]
                if triplet in ('TAG', 'TGA', 'TAA') :
                    return (i, j+2) #j+2 is for including last 2 bases of stop codon as j points only T of TGA
            
            # Reaching this point, no stop codon was found, so return 0
            return 0
    # Reaching this point, no start codon was found, so return 0
    return 0

# it will contain all the start and stop index of orfs
positions = []

for frame in range(1, 4):
    start_position = frame - 1
    firstorfposition = firstORF(seq, start_position)

    if (firstorfposition[1] - firstorfposition[0]) > 100:
        # NCBI index starts from 1. So, add 1 to orf position
        positions.append((firstorfposition[0]+1, firstorfposition[1]+1))
    start_position = firstorfposition[1]

    while start_position:
        firstorfposition = firstORF(seq, start_position)

        if firstorfposition == 0:
            break
        else:
            if (firstorfposition[1] - firstorfposition[0]) > 100:
                # NCBI index starts from 1. So, add 1 to orf position
                positions.append((firstorfposition[0]+1, firstorfposition[1]+1))
            start_position = firstorfposition[1]


# removes the duplicate ranges from the positions variable            
positions = list(set(positions))

with open("output.txt", "w") as f:
    f.write("start..stop\n")
    f.write("-----------\n")
    for position in sorted(positions):
        f.write(f"{position[0]}..{position[1]}\n")
        
    f.write("\nMy output orf matches with NCBI CDS range are:\n25393..26220\n26523..27191")


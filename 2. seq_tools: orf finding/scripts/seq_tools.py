import re

def load_fasta(fasta_filename):
    fasta = []
    test = []
    with open(fasta_filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                seq = line[1:]
                if seq not in fasta:
                    test.append(''.join(fasta))
                    fasta = []
                continue
            sequence = line
            fasta.append(sequence)

    if fasta:
        test.append(''.join(fasta))

    fasta = test[1:]
    return fasta

# creating a class named 'Sequence' to store a sequence
class Sequence:
    
    def __init__(self, seq): # seq is an attribute of Sequence class
        self.seq = seq
        
    # get_orf is a method to return largest orf as a protein sequence
    def get_orf(self):
        
        seq = self.seq
        
        # complement the sequence and then reverse complement it
        complement_pair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        complement_seq = "".join(complement_pair.get(base) for base in seq)
        reverse_complement_seq = complement_seq[::-1]
    
        # this pattern will search for orfs with start and stop codon using regular expression
        pattern = "ATG.*?(?:TGA|TAA|TAG)"
        extracted1 = re.findall(pattern, seq) # for 5'-->3' sequence
        extracted2 = re.findall(pattern, reverse_complement_seq) # for 3'-->5' sequence
        extracted_seq = list(set(extracted1).union(set(extracted2))) # concatenating 2 extracted orfs
        
        # "codon : amino acid" pair table for translation
        table = {
                'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
                'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
                }
        
        if extracted_seq: # if extracted_seq is not empty
            maximum_orf = max(extracted_seq, key = len) # finding out largest orfs

            # prepare sequence for translation
            index = 3 * (len(maximum_orf) // 3)
            sequence = maximum_orf[0:index]

            # translate the sequence into protein
            protein = "".join([table[sequence[i:i + 3]] for i in range(0, len(sequence), 3)])
            print(protein)
        else:
            print("No orf is found in the given sequence!!!")
        
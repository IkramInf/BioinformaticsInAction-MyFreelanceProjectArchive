# you may run the ORF_finder.py as: python ORF_finder.py genome.fasta -m 100 -o output.txt

# Import the library
import argparse
import re

# Create the parser
parser = argparse.ArgumentParser()
# Add arguments
parser.add_argument('filenames', metavar='input', type=str)
parser.add_argument('-m', '--minlength', type=int, default=100)
parser.add_argument('-o', '--output', type=str, default='output.txt')

# Parse the argument
args = parser.parse_args()

class ORFfinder():
    def __init__(self, input_filename, min_length, output_filename):
        # set all the variables
        self.inp_filename = input_filename
        self.min_length = min_length
        self.out_filename = output_filename
          
    def readSequence(self):
        # open file in read mode
        with open(self.inp_filename, "r") as f:
            description = f.readline()
            sequence = f.read()
        # replace all newline characters
        description = description.replace('\n', '')
        sequence = sequence.replace('\n', '')
        # create set object with a, c, g, t
        valid_base = {'a', 'c', 'g', 't'}
        # find out bases differ than a, c, g, t
        invalid_base = set(sequence).difference(valid_base)
        # remove all invalid bases
        for invalid in invalid_base:
            sequence = sequence.replace(invalid, '')
        # return the sequence
        return description, sequence.upper()

    def getORF(self, sequence):
        # to find out orf between M and *
        pattern = "(M.*?)(?:\*)"
        orf = re.finditer(pattern, sequence)
        return orf
    
    def getReverseStrand(self, forward_strand):
        # create dictionary of base:complement pair
        complement_pair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        # complement the sequence
        complement = "".join(complement_pair[base] for base in forward_strand)
        # reverse the complementary sequence
        reverse_strand = complement[::-1]
        return reverse_strand
    
    def orfExtract(self, sequence):
        # dictionary of codon : amino acid pair
        codonTable = {
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
                    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }

        # create a dictionary to store proteins
        ORFs = []

        # for iterating over forward sequence
        for frame in range(1, 7):
            # make the sequence length as multiple of 3
            # if sequence length = 25 and frame = 0 then, 25-0=25, 25//3=8 and finally length=3*8
            L = 3 * ((len(sequence)-frame) // 3)
            # take the part of sequence to perform translation
            seq = sequence[frame:frame+L]
            protein = "".join([codonTable[seq[i:i + 3]] for i in range(0, len(seq), 3)])
            orf = self.getORF(protein)
            for m in orf:
                s, e = m.start(), m.end()
                if (e-s) > self.min_length:
                    ORFs.append((protein[s:e-1], frame, e-s, 3*s+frame))

        # reverse the sequence
        reverse_strand = self.getReverseStrand(sequence)

        # for iterating over reverse sequence
        for frame in range(1, 7):
            L = 3 * ((len(reverse_strand)-frame) // 3) # Multiple of three
            rev = reverse_strand[frame:frame+L]
            protein = "".join([codonTable[rev[j:j + 3]] for j in range(0, len(rev), 3)])
            orf = self.getORF(protein)
            for n in orf:
                s, e = n.start(), n.end()
                if (e-s) > self.min_length:
                    ORFs.append((protein[s:e-1], frame, e-s, 3*s+frame))

        return ORFs       
    
    def write(self):
        
        description, genome = self.readSequence()
        orfs = self.orfExtract(genome)
        
        with open(self.out_filename, "w") as f:
            for value in orfs:
                f.write(f"{description} {value[1]} {value[2]} {value[3]}\n")
                L = len(value[0])
                # format sequence as 60 bases per line
                for i in range(0, L, 60):
                    if (i+60) < L:
                        f.write(f"{value[0][i:i+60]}\n")
                    else:
                        f.write(f"{value[0][i:]}\n")
                        
                        
if __name__ == '__main__':
    orf_finder = ORFfinder(args.filenames, args.minlength, args.output)
    orf_finder.write()
    
    
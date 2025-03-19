# you may run the ORF_finder.py as: python ORF_finder.py genome.fasta -m 100 -o output.txt

# Import the library
import argparse
import re

# Create the parser
parser = argparse.ArgumentParser()
# Add arguments
# filenames as sys.argv[1]
parser.add_argument('filenames', metavar='input', type=str)
# -m as sys.argv[2]
parser.add_argument('-m', '--minlength', type=int, default=100)
# -o as sys.argv[3]
parser.add_argument('-o', '--output', type=str, default='output.txt')

# Parse the argument
args = parser.parse_args()

# find out orfs and write to a fasta file
class ORFfinder():
    # __init__ is used to initialize the class
    def __init__(self, input_filename, min_length, output_filename):
        # set all the variables
        self.inp_filename = input_filename
        self.min_length = min_length
        self.out_filename = output_filename
    
    # read a sequence from a fasta file and return description and sequence
    def readSequence(self):
        # open file in read mode
        with open(self.inp_filename, "r") as f:
            # read the first line as description
            description = f.readline()
            # read the other lines as sequence
            sequence = f.read()
        # replace all newline characters
        description = description.replace('\n', '')
        # join each word of description with _ (underscore) and remove newline character
        description = "_".join(description.split())
        sequence = sequence.replace('\n', '')
        # create set object with a, c, g, t
        valid_base = {'a', 'c', 'g', 't'}
        # find out unique bases in the sequence
        unique_base = set(sequence)
        # find out bases differ than a, c, g, t
        invalid_base = unique_base.difference(valid_base)
        # remove all invalid bases
        for invalid in invalid_base:
            sequence = sequence.replace(invalid, '')
        # return description and sequence in uppercase
        return description, sequence.upper()
    
    # find out iterator yielding orf instances 
    def getORF(self, sequence):
        # to find out orf between M and *
        pattern = "(M.*?)(?:\*)"
        # use finditer from re module
        orf = re.finditer(pattern, sequence)
        return orf
    
    # complement a sequence first, then reverse it
    def getReverseStrand(self, forward_strand):
        # create dictionary of base:complement pair
        complement_pair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        # complement the sequence
        # create a list to store complement base
        complement = []
        for base in forward_strand:
            # append the complemented base to the complement list
            complement.append(complement_pair[base])
        # make complement list a string sequence
        complement = "".join(complement)
        # reverse the complementary sequence
        reverse_strand = complement[::-1]
        return reverse_strand
    
    # extract all orfs from forward and reverse strand
    def orfExtract(self, sequence):
        # dictionary of codon : amino acid pair
        # represent stop codon with *
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

        # create a list to store all orfs of forward and reverse strand
        ORFs = []

        # for iterating over forward sequence
        for frame in range(1, 7):
            # make the sequence length as multiple of 3
            # if sequence length = 25 and frame = 0 then, 25-0=25, 25//3=8 and finally length=3*8
            L = 3 * ((len(sequence)-frame) // 3)
            # take the part of sequence to perform translation
            seq = sequence[frame:frame+L]
            # create a list to store forward strand amino acids
            protein = []
            # iterate over all forward sequences
            for i in range(0, len(seq), 3):
                # take codon as 3 bases
                codon = seq[i:i + 3]
                # take corresponding amino acid of codon and append to protein list
                protein.append(codonTable[codon])
            # make the protein list to protein sequence    
            protein = "".join(protein)
            # get all orfs of the protein
            orf = self.getORF(protein)
            # iterate over iterable orf instances and extract information
            for i1, m in enumerate(orf):
                # get orf start and end index
                s, e = m.start(), m.end()
                # check the length of orf is greater than min_length or not
                if (e-s) > self.min_length:
                    # append orf sequence, frame, length, start position and index to ORFs list
                    ORFs.append((protein[s:e-1], frame, e-s, 3*s+frame, i1+1))

        # reverse the forward strand
        reverse_strand = self.getReverseStrand(sequence)

        # for iterating over reverse sequence
        for frame in range(1, 7):
            # make the sequence length as multiple of 3
            L = 3 * ((len(reverse_strand)-frame) // 3)
            # take the part of sequence to perform translation
            rev = reverse_strand[frame:frame+L]
                    
            # create a list to store reverse strand amino acids
            protein = []        
            # iterate over all reverse strand sequences
            for j in range(0, len(rev), 3):
                # take codon as 3 bases
                codon = rev[j:j + 3]
                # take corresponding amino acid of codon and append to protein list
                protein.append(codonTable[codon])
            # make the protein list to protein sequence    
            protein = "".join(protein)
            # get all orfs of the protein
            orf = self.getORF(protein)
            # iterate over iterable orf instances and extract information
            for i2, n in enumerate(orf):
                # get orf start and end index
                s, e = n.start(), n.end()
                # check the length of orf is greater than min_length or not
                if (e-s) > self.min_length:
                    # append orf sequence, frame, length, start position and index to ORFs list
                    ORFs.append((protein[s:e-1], frame, e-s, 3*s+frame, i2+1))
        # return all orfs as a list
        return ORFs       
    
    # write the orfs into a file
    def write(self):
        # calling the readSequence and orfExtract function
        description, genome = self.readSequence()
        orfs = self.orfExtract(genome)
        # write the orfs into a file
        with open(self.out_filename, "w") as f:
            # iterate over all orfs
            for value in orfs:
                # write the description line in a fasta file followed by > sign
                f.write(f"{description}_{value[4]} {value[1]} {value[2]} {value[3]}\n")
                L = len(value[0])
                # format sequence as 60 bases per line
                for i in range(0, L, 60):
                    # if the length of orf greater than 60, write 60 bases and proceed for next 60
                    if (i+60) < L:
                        f.write(f"{value[0][i:i+60]}\n")
                    else:
                        # if orf length less than 60, write the entire orf in a line
                        f.write(f"{value[0][i:]}\n")
                        
# to ensure run this section only when run the file directly                        
if __name__ == '__main__':
    # create object of ORFfinder class with arguments
    orf_finder = ORFfinder(args.filenames, args.minlength, args.output)
    # run the class method 'write' with orf_finder object to extract and write orfs into fasta file
    orf_finder.write()
    
    
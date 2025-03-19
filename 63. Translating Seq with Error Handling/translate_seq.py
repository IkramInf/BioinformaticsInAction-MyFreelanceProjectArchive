import os
import sys

def file_handler(filename):
    """
    Handle File opening errors
    
    Parameters:
        filename : input filename
    Returns:
        returns file handler object
    """
    if os.path.exists(filename):
        if not os.stat(filename).st_size == 0:
            input_file = open(filename, "r")
            return input_file            
        else:
            sys.exit("The file is empty!!!")
    else:
        sys.exit("The file doesn't exist!!! Please enter a valid filename.")
        
def extract_codon_table(filename):    
    """
    Extract codon table from text file
    
    Parameters:
        filename : input filename
    Returns:
        returns codon table as dictionary
    """
    Codon_Table = file_handler(filename)
    codon_table = {}
    for line in Codon_Table.readlines():
        codon, aa = line.split()
        codon_table.update({codon.upper():aa})

    Codon_Table.close()        
    return codon_table

def read_fasta(filename):
    """
    Read fasta file
    
    Parameters:
        filename : input filename
    Returns:
        returns dna sequence from fasta file
    """
    fasta_reader = file_handler(filename)  
    header = fasta_reader.readline()
    seq = fasta_reader.read().replace("\n", "")
    fasta_reader.close()
    return seq.upper().replace("U", "T")

def complement(seq):
    """
    Complement a dna sequence
    
    Parameters:
        seq : a dna sequence
    Returns:
        returns complement of dna sequence
    """
    complement_pair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement_pair[base] for base in seq)

def triplets(seq):
    """
    Make list of codons from a dna sequence
    
    Parameters:
        seq : a dna sequence
    Returns:
        returns list of codons from a dna sequence
    """
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]

def grouping(seq):
    """
    Make list of sixty characters substring from a dna sequence
    
    Parameters:
        seq : a dna sequence
    Returns:
        returns list of sixty characters substring from a dna sequence
    """   
    return [seq[i:i + 60] for i in range(0, len(seq), 60)]

def check_length(proteins, sequence):
    """
    check the length of each protein and return all of them as a list
    
    Parameters:
        proteins : list of proteins
        sequence : a dna sequence
    Returns:
        returns a list of proteins and sequence together
    """
    N = len(sequence)
    prots = []
    for i, prot in enumerate(proteins.values()):
        if i > 2:
            prot = "".join(triplets(prot)[::-1])
        if len(prot) != N:
            prot = prot + "   "
        prots.append(prot)
    return prots[:3] + [sequence] + [prots[4]] + [prots[5]] + [prots[3]]

def count_orf(proteins, min_size=1):
    """
    Count orf in a protein sequence
    
    Parameters:
        proteins : a protein sequence
        min_size : minimum length of orf
    Returns:
        returns the list of number of orfs in the sequence
    """
    return ([len([orf for orf in v.split("***") if len(orf)>min_size]) for v in proteins.values()], min_size)


def translate(sequence):
    """
    Translate dna sequence into protein in six frames
    
    Parameters:
        sequence : a dna sequence
    Returns:
        returns a dictionary of proteins as {frame:protein}
    """
    # dictionary of codon : amino acid pair
    codon_table = extract_codon_table("Codon_Table.txt")
    
    # create a dictionary to store proteins
    proteins = {}
    
    # for iterating over forward sequence, frame: f1, f2, f3
    for frame in range(3):
        # make the sequence length as multiple of 3
        # if sequence length = 25 and frame = 0 then, 25-0=25, 25//3=8 and finally length=3*8
        length = 3 * ((len(sequence)-frame) // 3)
        # take the part of sequence to perform translation
        seq = sequence[frame:frame+length]
        proteins["F"+str(frame+1)] = "".join([codon_table[seq[i:i + 3]] for i in range(0, len(seq), 3)])
        
    # reverse the sequence
    rev_seq = complement(sequence[::-1])
    
    # for iterating over reverse sequence
    for frame in range(3):
        #length = 3 * ((len(rev_seq)-frame) // 3) # Multiple of three
        length = len(rev_seq)
        if frame%2==1:
            rev = rev_seq[frame:length] + rev_seq[0:frame]
        else:
            rev = rev_seq[frame:length]
        proteins["F"+str(6-frame)] = "".join([codon_table.get(rev[j:j + 3],'XXX') for j in range(0, len(rev), 3)])
        
    orfs, ms = count_orf(proteins)
    results = check_length(proteins, sequence)
    group_seq = [grouping(seq) for seq in results]
    group_seq_iter = list(zip(*group_seq))
    format_output(group_seq_iter, orfs, ms)

def format_output(group_seq_iter, orfs, min_size):   
    """
    Format proteins and sequence to print on screen
    
    Parameters:
        group_seq_iter : output of grouping function
        orfs : a list of number of orfs
        min_size : minimum length of orf
    Returns:
        None
    """
    L = len(group_seq_iter)
    p = 1
    three_to_one_letter = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                           'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 
                           'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 
                           'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', '***': '*'}
    
    for n, seqs in enumerate(group_seq_iter):
        if n+1 != L:
            for i, ss in enumerate(zip([triplets(seq) for seq in seqs])):
                print("   ", end="")
                if i in [0, 1, 2]:
                    print("\t  " + " "*i, end="")
                    print("  ".join([three_to_one_letter.get(s, 'X') for s in ss[0]]) +" "*(2-i)+("F"+str(i+1)).rjust(8, " "))
                elif i == 3:
                    print(str(p).rjust(6, " ") + " " + "".join(ss[0]).upper() + " " + str(p+59).ljust(7, " "))
                    print("\t  " + "----:----|"*6)
                    print("   " + str(p).rjust(6, " ") + " " +complement("".join(ss[0]).upper()) + " " + str(p+59).ljust(7, " "))
                elif i == 4:
                    print("\t  " + " " + "  ".join([three_to_one_letter.get(s, 'X') for s in ss[0]]) +" "+ ("F"+str(i+2)).rjust(8, " "))
                elif i == 5:
                    print("\t  " + "  ".join([three_to_one_letter.get(s, 'X') for s in ss[0]]) +"  "+ ("F"+str(i)).rjust(8, " "))
                else:
                    print("\t  " + "  " + "  ".join([three_to_one_letter.get(s, 'X') for s in ss[0]]) +""+ ("F"+str(i-2)).rjust(8, " "))
        else:
            for i, ss in enumerate(zip([triplets(seq) for seq in seqs])):
                print("   ", end="")
                if i in [0, 1, 2]:
                    print("\t  " + " "*i, end="")
                    print("  ".join([three_to_one_letter.get(s, 'X') for s in ss[0]])+" "*(2-i)+("F"+str(i+1)).rjust(8, " "))
                elif i == 3:
                    SEQ = "".join(ss[0]).upper()
                    print(str(p).rjust(6, " ") + " " + SEQ + " " + str(p+len(SEQ)-1).ljust(7, " "))
                    print("\t  " + ("----:----|"*6)[:len(SEQ)])
                    print("   " + str(p).rjust(6, " ") + " " + complement("".join(ss[0]).upper()) + " " + str(p+len(SEQ)-1).ljust(7, " "))
                elif i == 4:
                    print("\t  " + " " + "  ".join([three_to_one_letter.get(s, 'X') for s in ss[0]]) +"  "+ ("F"+str(i+2)).rjust(7, " "))
                elif i == 5:
                    print("\t  " + "  ".join([three_to_one_letter.get(s, 'X') for s in ss[0]]) +"  "+ ("F"+str(i)).rjust(8, " "))
                else:
                    print("\t  " + "  " + "  ".join([three_to_one_letter.get(s, 'X') for s in ss[0]]) +""+ ("F"+str(i-2)).rjust(8, " "))
        print("\n\n")
        p += 60
        
    print(f"##############################\nMinimum size of ORFs : {min_size}\n\n")
    orfs = orfs[:3] + [orfs[3]] + [orfs[5]] + [orfs[4]]
    for i, orf in enumerate(orfs):
        print(f"Total ORFs in frame {i+1} :    {orf}")
    print(f"\nTotal ORFs :    {sum(orfs)}\n##############################\n")


if __name__ == "__main__":
    sequence = read_fasta("myoglobin.txt")
    translate(sequence.upper())
    
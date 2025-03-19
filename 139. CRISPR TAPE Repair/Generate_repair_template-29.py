import sys
import random
import argparse
import collections
import numpy as np
import itertools
import pandas as pd

codon_dict = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S1', 'AGT': 'S1', 'AGA': 'R1', 'AGG': 'R1',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R2', 'CGC': 'R2', 'CGG': 'R2', 'CGT': 'R2',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S2', 'TCC': 'S2', 'TCG': 'S2', 'TCT': 'S2',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TGC': 'C', 'TGT': 'C',
        'TGG': 'W', 'TAA': '*', 'TAG': '*', 'TGA': '*'}

aa_to_codon = {
        'I':['ATA', 'ATC', 'ATT'],
        'M': ['ATG'],
        'T':['ACA', 'ACC', 'ACG', 'ACT'],
        'N':['AAC', 'AAT'], 
        'K':['AAA', 'AAG'],
        'S1':['AGC', 'AGT'],
        'R1':['AGA', 'AGG'],
        'L':['CTA', 'CTC', 'CTG', 'CTT'],
        'P':['CCA', 'CCC', 'CCG', 'CCT'],
        'H':['CAC', 'CAT'],
        'Q':['CAA', 'CAG'],
        'R2':['CGA', 'CGC', 'CGG', 'CGT'],
        'V':['GTA', 'GTC', 'GTG', 'GTT'],
        'A':['GCA', 'GCC', 'GCG', 'GCT'],
        'D':['GAC', 'GAT'],
        'E':['GAA', 'GAG'],
        'G':['GGA', 'GGC', 'GGG', 'GGT'],
        'S2':['TCA', 'TCC','TCG', 'TCT'],
        'F':['TTC', 'TTT'],
        'L':['TTA', 'TTG'],
        'Y':['TAC', 'TAT'],
        'C':['TGC', 'TGT'],
        'W': ["TGG"]}
    
human_frequency_dict = {
        'ATA': 0.16, 'ATC': 0.48, 'ATT': 0.36, 'ATG': 1.00,
        'ACA': 0.28, 'ACC': 0.36, 'ACG': 0.12, 'ACT': 0.24,
        'AAC': 0.54, 'AAT': 0.46, 'AAA': 0.42, 'AAG': 0.58,
        'AGC': 0.24, 'AGT': 0.15, 'AGA': 0.20, 'AGG': 0.20,
        'CTA': 0.07, 'CTC': 0.20, 'CTG': 0.41, 'CTT': 0.13,
        'CCA': 0.27, 'CCC': 0.33, 'CCG': 0.11, 'CCT': 0.28,
        'CAC': 0.59, 'CAT': 0.41, 'CAA': 0.25, 'CAG': 0.75,
        'CGA': 0.11, 'CGC': 0.19, 'CGG': 0.21, 'CGT': 0.08,
        'GTA': 0.11, 'GTC': 0.24, 'GTG': 0.47, 'GTT': 0.18,
        'GCA': 0.23, 'GCC': 0.40, 'GCG': 0.11, 'GCT': 0.26,
        'GAC': 0.54, 'GAT': 0.46, 'GAA': 0.42, 'GAG': 0.58,
        'GGA': 0.25, 'GGC': 0.34, 'GGG': 0.25, 'GGT': 0.16,
        'TCA': 0.15, 'TCC': 0.22, 'TCG': 0.06, 'TCT': 0.18,
        'TTC': 0.55, 'TTT': 0.45, 'TTA': 0.07, 'TTG': 0.13,
        'TAC': 0.57, 'TAT': 0.43, 'TGC': 0.55, 'TGT': 0.45,
        'TGG': 1.00, 'TAA': 0, 'TAG': 0, 'TGA': 0}

def get_exons(seq):
    """
    splicing introns to get cds
    """
    return "".join([base for base in seq if base.isupper()])

# 生成反向互补序列
def reverse_complement(dna_sequence):
    """
    Reverse complement a dna sequence
    """
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    reverse_sequence = dna_sequence[::-1]
    complement_sequence = ''.join([complement_dict[base] for base in reverse_sequence])
    return complement_sequence

def read_seq(filename):
    """
    Read Fasta File (Only a sequence)
    """
    with open(filename, "r") as reader:
        header = reader.readline()
        sequence = reader.read().strip().replace("\n", "")
    return sequence

# PAM变异
def mutate_pam(pam, gene_sequence):
    """
    Perform mutation in pam sequence
    """
    pam = pam.upper()
    index = gene_sequence.upper().find(pam)
    # True: coding sequence
    if gene_sequence[index : index + 3].isupper():
        amino_acid = codon_dict[pam]
        synonymous_codons = aa_to_codon[amino_acid]
        synonymous_codons = [codon for codon in synonymous_codons if (not codon.endswith("AG")) and (codon != pam)]
        if synonymous_codons:
            mutated_pam = random.choice(synonymous_codons)
        else:
            mutated_pam = ""
    else:
        while True:
            second_base = random.choice(['C', 'G', 'T'])
            third_base = random.choice(['A', 'C', 'G', 'T'])
            if not pam.endswith(f"{second_base}{third_base}"):
                mutated_pam = f"{pam[0]}{second_base}{third_base}"
                break
    return mutated_pam

def recodonize_sequence(original_dna_sequence):
    """
    Recodonise dna sequence to most abundant codons but not changing protein sequence
    """
    dna_sequence_exon = get_exons(original_dna_sequence)

    # Iterate through codons
    codon_list = [dna_sequence_exon[i:i+3] for i in range(0, len(dna_sequence_exon) - 2, 3)]

    # Initialize variables
    current_list_check = []
    recodnized_sequence = ""

    # Iterate through codons and find possible amino acids
    for current_codon in codon_list:
        current_aa = codon_dict[current_codon]
        
        if current_aa in ['M', 'W', '*']:
            recodnized_sequence += current_codon
        else:
            freq_dict = {}
            for codon in aa_to_codon[current_aa]:
                freq_dict[codon] = human_frequency_dict[codon]
                
            max_freq_codon = max(freq_dict, key=lambda x: freq_dict[x])
            #print(freq_dict, max_freq_codon)
            recodnized_sequence += max_freq_codon

    #print("Recodonized Sequence:", recodnized_sequence)
    return recodnized_sequence

# 修复片段基因提取
def get_donor_template(foreign_DNA: str, forward_sgRNA: str, reverse_sgRNA: str):
    """
    Get the donor template for 300 bp sequence
    """
    forward_pam = mutate_pam(forward_sgRNA[:3], foreign_DNA)
    forward_sgRNA = forward_sgRNA[3:]
    reverse_pam = mutate_pam(reverse_sgRNA[-3:], foreign_DNA)
    reverse_sgRNA = reverse_sgRNA[:-3]
    
    forward_sgRNA_len = len(forward_sgRNA)
    reverse_sgRNA_len = len(reverse_sgRNA)

    up_index = foreign_DNA.find(forward_sgRNA)
    down_index = foreign_DNA.find(reverse_complement(reverse_sgRNA))

    # Adjust indices if reverse sgRNA is found before forward sgRNA
    if down_index < up_index:
        gene_start_index = down_index + reverse_sgRNA_len
        gene_end_index = up_index
        forward_seq = f"{reverse_pam}{reverse_sgRNA}"
        reverse_seq = f"{forward_sgRNA}{forward_pam}"
    else:
        gene_start_index = up_index + forward_sgRNA_len
        gene_end_index = down_index
        forward_seq = f"{forward_pam}{forward_sgRNA}"
        reverse_seq = f"{reverse_sgRNA}{reverse_pam}"

    internal_seq = foreign_DNA[gene_start_index:gene_end_index]
    internal_seq = recodonize_sequence(internal_seq)
    
    return f"{forward_seq}{internal_seq}{reverse_seq}"

def get_pri_hom(foreign_DNA: str, forward_sgRNA: str, reverse_sgRNA: str, gene_len: int):
    """
    Get up and down primer, up and down homology sequence
    """
    forward_sgRNA = forward_sgRNA[3:]
    reverse_sgRNA = reverse_sgRNA[:-3]
    forward_sgRNA_len = len(forward_sgRNA)
    reverse_sgRNA_len = len(reverse_sgRNA)

    forward_sgRNA_index = foreign_DNA.find(forward_sgRNA)
    
    # Check if forward_sgRNA is found
    if forward_sgRNA_index == -1:
        #raise ValueError("Forward sgRNA not found in foreign DNA.")
        return ('', '', '', '')

    # Find the reverse complement of reverse_sgRNA
    reverse_complement_reverse_sgRNA = reverse_complement(reverse_sgRNA)
    
    reverse_sgRNA_index = foreign_DNA.find(reverse_complement(reverse_sgRNA))
    
    # Check if reverse_sgRNA reverse complement is found
    if reverse_sgRNA_index == -1:
        #raise ValueError("Reverse complement of reverse sgRNA not found in foreign DNA.")
        return ('', '', '', '')

    # Adjust indices if reverse sgRNA is found before forward sgRNA
    if reverse_sgRNA_index < forward_sgRNA_index:
        gene_start_index = reverse_sgRNA_index
        gene_end_index = forward_sgRNA_index + forward_sgRNA_len
    else:
        gene_start_index = forward_sgRNA_index
        gene_end_index = reverse_sgRNA_index + reverse_sgRNA_len

    # calculation of homology length
    hom_len = int((260 - gene_len) / 2)

    # primer sequences
    up_primer = foreign_DNA[gene_start_index - hom_len - 20: gene_start_index - hom_len]
    down_primer = foreign_DNA[gene_end_index + hom_len: gene_end_index + hom_len + 20]

    # homology sequences
    up_homolog = foreign_DNA[gene_start_index - hom_len: gene_start_index]
    down_homolog = foreign_DNA[gene_end_index: gene_end_index + hom_len]

    return (up_primer, down_primer, up_homolog, down_homolog)


# 检查sgRNA
def check_dna(tag,temp_DNA):
    if tag in temp_DNA:
        return tag
    else:
        rv_tag = reverse_complement(tag)
        if rv_tag in temp_DNA:
            return rv_tag
        else:
            sys.exit('Please enter the correct sgRNA!')

def check_grna(grna_df, foreign_dna):
    """
    Check gRNA present in foreign DNA sequence 
    """
    #foreign_dna = foreign_dna.upper()
    forward_grna, reverse_grna = [], []
    for i in range(grna_df.shape[0]):
        row = grna_df.loc[i]
        distance = int(row['Distance from Amino Acid (bp)'])
        grna_seq = row['gRNA Sequence']
        pam_seq = row['PAM']
        strand = row['Strand']
        
        if strand == "forward":
            if grna_seq in foreign_dna.upper():
                forward_grna.append(pam_seq + grna_seq)
        else:
            if reverse_complement(grna_seq) in foreign_dna:
                reverse_grna.append(grna_seq + pam_seq)
    return (forward_grna, reverse_grna)

def get_sgRNA(gRNA_path, cut_site:int, foreign_DNA):
    """
    Get each possible pair of forward and reverse sgRNA
    """
    cut_site = cut_site*3
    grna_df = pd.read_csv(gRNA_path, index_col=0, encoding='utf-8')
    # Drop rows with all NaN values
    grna_df.dropna(how='all', inplace=True)

    forward_grna, reverse_grna = check_grna(grna_df, foreign_DNA)
    
    if forward_grna and reverse_grna:
        grna_list_pairs = []
        for fgrna in forward_grna:
            for rgrna in reverse_grna:
                grna_list_pairs.append((fgrna, rgrna))
        return grna_list_pairs
    else:
        sys.exit('The gRNA and PAM in the table do not match the entered DNA template. Please check!')

def get_args():
    parser = argparse.ArgumentParser(description="Generate repair donor.py")
    parser.add_argument("-f", "--foreign_DNA", help="foreign DNA", required=True)
    parser.add_argument("-g", "--gRNA_list", help="gRNA list", required=True)
    parser.add_argument("-c", "--cut_site", help="cut site", type=int, required=True)
    parser.add_argument("-t", "--thres", help="Threshold for white list", type=int, required=True)
    parser.add_argument("-o", "--outfile", help="outfile path", required=True)
    return parser.parse_args()

if __name__ == '__main__':

    args = get_args()
    foreign_DNA = read_seq(args.foreign_DNA)
    gRNA_list = args.gRNA_list
    cut_site = args.cut_site
    thres = args.thres
    outfile = args.outfile
    
    # check sgrna whether is complemet with dna
    # up_sgRNA = check_dna(up_sgRNA,foreign_DNA)
    # down_sgRNA = check_dna(down_sgRNA,foreign_DNA)

    out_obj = open(outfile,'w')
    for forward_sgRNA, reverse_sgRNA in get_sgRNA(gRNA_list, cut_site, foreign_DNA):
        # extract dna
        repir_gene = get_donor_template(foreign_DNA, forward_sgRNA, reverse_sgRNA)
        # extraxt primer and homo
        # DNA length
        DNA_length = len(repir_gene)
        if DNA_length <= 94:
            up_primer, down_primer, up_homolog, down_homolog = get_pri_hom(foreign_DNA,forward_sgRNA,reverse_sgRNA,DNA_length)
        else:
            up_primer, down_primer, up_homolog, down_homolog = '','','',''
            print('The target gene must be less than 94bp, please re-enter the new sgRNA!')

        # mutate pam
        #mutation_forward_sgRNA = mutate_pam(forward_sgRNA[-3:], repir_gene)
        #mutation_reverse_sgRNA = mutate_pam(reverse_sgRNA[-3:], repir_gene)

        # output 
        out_obj.write(f'Donor Template: {repir_gene}\n')
        out_obj.write(f'Upstream primer: {up_primer}\n')
        out_obj.write(f'Downstream primer: {down_primer}\n')
        out_obj.write(f'Upstream homologous: {up_homolog}\n')
        out_obj.write(f'Downstream homologous: {down_homolog}\n')
        out_obj.write(f'Upstream gRNA: {forward_sgRNA[3:]}\n')
        out_obj.write(f'Downstream gRNA: {reverse_sgRNA[:-3]}\n')
        #out_obj.write(f'Upstream PAM: {mutation_forward_sgRNA}\n')
        #out_obj.write(f'Downstream PAM: {mutation_reverse_sgRNA}\n')
        out_obj.write(f"Sequence: {up_primer}{up_homolog}{repir_gene}{down_homolog}{down_primer}\n\n")
        
        #whitelist = mutation_main(repir_gene, thres)
        #out_obj.write(f'DNA white:\n')
        #for d in whitelist:
        #    out_obj.write(f'{d}\n')
    out_obj.close()




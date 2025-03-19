import re
import sys
import random
import argparse
import pandas as pd
from Bio.Data import CodonTable
from Bio.SeqUtils import MeltingTemp
from Bio.Seq import complement, reverse_complement, translate
from CRISPR_TAPE.shared_functions import get_codon_index
from CRISPR_TAPE.specific_function import specific_function

codon_dict = CodonTable.standard_dna_table.forward_table
aa_to_codon = {}
for codon, aa in codon_dict.items():
    aa_to_codon.setdefault(aa, []).append(codon)
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

def read_seq(filename):
    """
    Read Fasta File (Only a sequence)
    """
    with open(filename, "r") as reader:
        header = reader.readline()
        sequence = reader.read().strip().replace("\n", "")
    return sequence

def mutate_pam(pam):
    return f"{pam[0]}C{pam[1]}"

def get_target_position(reference_genome, coding_sequence, aa_position):
    #aa = translate(coding_sequence)[aa_position - 1]
    i = (aa_position - 1) * 3
    target_codon = coding_sequence[i : i + 3]
    target_pos = reference_genome[i+3:].index(codon) + i + 3
    return (target_codon, target_pos)

def mutate_target_codon(target_codon):
    amino_acid = translate(target_codon)

    abundant_codons = {}
    for codon in aa_to_codon[amino_acid]:
        abundant_codons[codon] = human_frequency_dict[codon]

    most_abundant_codon = max(abundant_codons, key=lambda x: abundant_codons[x])
    return most_abundant_codon

def recodonize_sequence(dna_seq):
    """
    Recodonise dna sequence to most abundant codons but not changing protein sequence
    """
    dna_sequence = ""
    for base in dna_seq:
        if base == "N":
            alternative_base = random.choice(['A', 'C', 'G', 'T'])
            dna_sequence += alternative_base
        else:
            dna_sequence += base

    # Iterate through codons
    codon_list = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence) - 2, 3)]

    # Initialize variables
    current_list_check = []
    recodnized_sequence = ""

    # Iterate through codons and find possible amino acids
    for current_codon in codon_list:
        if current_codon in ['ATG', 'TGG', 'TAA', 'TAG', 'TGA']:
            recodnized_sequence += current_codon
        else:
            current_aa = codon_dict[current_codon]
            freq_dict = {}
            for codon in aa_to_codon[current_aa]:
                freq_dict[codon] = human_frequency_dict[codon]

            max_freq_codon = max(freq_dict, key=lambda x: freq_dict[x])
            #print(freq_dict, max_freq_codon)
            recodnized_sequence += max_freq_codon

    #print("Recodonized Sequence:", recodnized_sequence)
    return recodnized_sequence

def get_excised_sequences(unmodified_gene_sequence, gRNA_df):

    sgrna_pair_list = get_sgrna_pairs(gRNA_df)

    unmodified_gene_sequence = unmodified_gene_sequence.upper()
    excised_sequences = []

    for fgrna, rgrna in sgrna_pair_list:
        fgrna, fpam = fgrna
        rgrna, rpam = rgrna

        fsequence = fgrna + fpam
        rsequence = reverse_complement(rgrna + rpam)

        try:
            up_index = unmodified_gene_sequence.index(fsequence) + len(fsequence)
            down_index = unmodified_gene_sequence.index(rsequence) - 3  # -3 for pam
        except ValueError:
            continue

        length = len(unmodified_gene_sequence) - (up_index + down_index)
        total_length = len(fgrna) + length + len(rgrna)
        if (total_length >= 70) and (total_length <= 94):
            seq = unmodified_gene_sequence[up_index : up_index + length]
            excised_sequences.append((fpam, fgrna, up_index, length, rgrna, rpam, seq))

    return excised_sequences

def get_donor_template(excised_sequence, mutated_target_codon):
    """
    Get the donor template for 300 bp sequence
    """

    fpam, fgrna, up_index, length, rgrna, rpam, internal_seq = excised_sequence
    fpam, rpam = reverse_complement(mutate_pam(fpam)), mutate_pam(rpam)

    internal_seq = recodonize_sequence(internal_seq)
    target_pos = len(internal_seq) // 2
    internal_seq = internal_seq[:target_pos] + mutated_target_codon + internal_seq[target_pos+3:]
    return (f"{fpam} {reverse_complement(fgrna)} {internal_seq} {rgrna} {rpam}", fpam, rpam)

def get_homology_regions(excised_sequence, reference_genome, homology_arm_length=80):
    fpam, fgrna, up_index, length, rgrna, rpam, seq = excised_sequence
    up_end = up_index - len(fgrna) - 3
    up_start = up_end - homology_arm_length
    if up_start < 0:
        up_start = 0
    up_homolog = reference_genome[up_start : up_end]

    down_start = length + len(rgrna) + 3
    down_end = down_start + homology_arm_length
    down_homolog = reference_genome[down_start : down_end]

    return (up_homolog, down_homolog)

def get_primer(reference_genome, coding_sequence, aa_position=10, primer_length=20):
    up_end = get_target_position(reference_genome, coding_sequence, aa_position)[1]
    up_start = up_end - primer_length
    if up_start < 0:
        up_start = 0
    up_primer = reference_genome[up_start : up_end]

    down_start = up_end + 3
    down_end = down_start + primer_length
    down_primer = reference_genome[down_start : down_end]

    return (up_primer, down_primer)

def filter_grna(grna_df):
    # Convert 'G/C Content (%)' to numeric values
    grna_df['G/C Content (%)'] = pd.to_numeric(grna_df['G/C Content (%)'], errors='coerce')
    
    # Calculate the melting temperature
    # tm = MeltingTemp.Tm_staluc(dna_sequence)
    
    condition = (
        (grna_df['gRNA Sequence'].str.startswith('G')) &
        (grna_df['Off Target Count'] == 0) &
        (grna_df['G/C Content (%)'].between(60, 70))
    )
    
    filtered_grna_df = grna_df[condition]
    return filtered_grna_df

def get_sgrna_df(amino_acid_position, maximum_gRNA_distance, PAM,
                 genomic_loci, coding_sequence, reference_genome):
    # Run the cripr-tape specific function
    gRNA_df = specific_function(amino_acid_position, maximum_gRNA_distance, PAM,
                                   genomic_loci, coding_sequence, "", "", reference_genome)
    #print(gRNA_df)
    #gRNA_df = pd.read_csv(gRNA_filename, index_col=0, encoding='utf-8')
    # Drop rows with all NaN values
    gRNA_df.dropna(how='all', inplace=True)
    #gRNA_df = filter_grna(gRNA_df)  # filter gRNA_df
    return gRNA_df

def check_grna(genomic_loci, coding_sequence, grna_df, cut_site):
    """
    Check grna in genomic loci
    """
    protein_dict = get_codon_index(genomic_loci, coding_sequence)
    target_codon_start = protein_dict[cut_site - 1]['base_1'][0]
    target_codon_end = protein_dict[cut_site - 1]['base_3'][0]
    dna_edited = genomic_loci.upper()
    
    grnas = []  # (sgrna, cut_site, strand)
    
    for i in range(grna_df.shape[0]):
        row = grna_df.loc[i]
        distance = row['Distance from Amino Acid (bp)']
        strand = row['Strand']
        grna = row['gRNA Sequence']
        pam = row['PAM']

        if strand == "forward":
            if distance > 0:
                position = int(target_codon_end + distance)
                start, end = position - 16, position + 7
            else:
                position = int(target_codon_start + distance)
                start, end = position - 17, position + 6
            cut_grna = dna_edited[start : end]
            if cut_grna == grna + pam:
                grnas.append((cut_grna, start, 'forward'))
        else:  # strand == "reverse"
            if distance > 0:
                position = int(target_codon_end + distance)
                start, end = position - 5, position + 18
            else:
                position = int(target_codon_start + distance)
                start, end = position - 6, position + 17
            cut_grna = reverse_complement(dna_edited[start : end])
            if cut_grna == grna + pam:
                grnas.append((cut_grna, start, 'reverse'))
            
    return grnas

def get_sgrna_pairs(gRNA_df):

    forward_grna = gRNA_df[gRNA_df['Strand'] == 'forward'][['gRNA Sequence', 'PAM']].values
    reverse_grna = gRNA_df[gRNA_df['Strand'] == 'reverse'][['gRNA Sequence', 'PAM']].values

    grna_list_pairs = []
    if len(forward_grna)>=1 and len(reverse_grna)>=1:
        grna_list_pairs = []
        for fgrna in forward_grna:
            for rgrna in reverse_grna:
                grna_list_pairs.append((fgrna, rgrna))

    return grna_list_pairs

def get_args():
    parser = argparse.ArgumentParser(description="Generate repair donor.py")
    parser.add_argument("-f", "--foreign_DNA", help="foreign DNA", required=True)
    parser.add_argument("-r", "--reference_genome", help="reference genome", required=True)
    #parser.add_argument("-g", "--gRNA_list", help="gRNA list", required=True)
    parser.add_argument("-c", "--cut_site", help="cut site", type=int, default=10)
    #parser.add_argument("-t", "--thres", help="Threshold for white list", type=int, default=10)
    parser.add_argument("-o", "--outfile", help="outfile path", default="output.txt")
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    genomic_loci = read_seq(args.foreign_DNA)
    coding_sequence = get_exons(genomic_loci)
    reference_genome = read_seq(args.reference_genome)
    #df = get_sgrna_df(args.gRNA_list)
    amino_acid_position, maximum_gRNA_distance, PAM = args.cut_site, 1000, "NGG"
    df = get_sgrna_df(amino_acid_position, maximum_gRNA_distance, PAM,
                                   genomic_loci, coding_sequence, reference_genome)
    #check_grna(genomic_loci, coding_sequence, df, args.cut_site)
    grna_list_pairs = get_sgrna_pairs(df)
    excised_sequences = get_excised_sequences(genomic_loci, df)
    target_codon, target_pos = get_target_position(reference_genome, coding_sequence, aa_position=args.cut_site)
    mutated_target_codon = mutate_target_codon(target_codon)
    up_primer, down_primer = get_primer(reference_genome, coding_sequence)

    with open(args.outfile, "w") as out_obj:
        for excised_sequence in excised_sequences:
            donor_template, fpam, rpam = get_donor_template(excised_sequence, mutated_target_codon)
            up_homolog, down_homolog = get_homology_regions(excised_sequence, reference_genome)
            core_mutational_template = f"{up_primer} {up_homolog} {donor_template} {down_homolog} {up_primer}"

            # output 
            out_obj.write(f'Donor Template: {donor_template}\n')
            out_obj.write(f'Upstream primer: {up_primer}\n')
            out_obj.write(f'Downstream primer: {down_primer}\n')
            out_obj.write(f'Upstream homologous: {up_homolog}\n')
            out_obj.write(f'Downstream homologous: {down_homolog}\n')
            out_obj.write(f'Upstream gRNA: {reverse_complement(excised_sequence[1])}\n')
            out_obj.write(f'Downstream gRNA: {excised_sequence[4]}\n')
            out_obj.write(f'Upstream PAM: {fpam}\n')
            out_obj.write(f'Downstream PAM: {rpam}\n')
            out_obj.write(f"Core Mutational Template: {core_mutational_template}\n\n")

if __name__ == '__main__':
    # execute the whole script
    main()

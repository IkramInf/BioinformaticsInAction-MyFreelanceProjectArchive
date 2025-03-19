import os
import re
import sys
import random
import argparse
from itertools import product
import pandas as pd
from Bio.Data import CodonTable
from Bio.SeqUtils import MeltingTemp, gc_fraction
from Bio.Seq import complement, reverse_complement, translate
from CRISPR_TAPE.shared_functions import get_codon_index, PAMposition
from CRISPR_TAPE.specific_function import specific_function

codon_dict = CodonTable.standard_dna_table.forward_table
#print(codon_dict)
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
    """
    Perform mutation in pam sequence
    """
    # True: coding sequence
    if pam.isupper():
        amino_acid = codon_dict[pam]
        synonymous_codons = aa_to_codon[amino_acid]
        synonymous_codons = [codon for codon in synonymous_codons if (not codon.endswith("AG")) and (codon != pam)]
        if synonymous_codons:
            mutated_pam = random.choice(synonymous_codons)
        else:
            mutated_pam = f"{pam[0]}GA"
    else:
        while True:
            second_base = random.choice(['C', 'T'])
            third_base = random.choice(['A', 'C', 'T'])
            if not pam.endswith(f"{second_base}{third_base}"):
                mutated_pam = f"{pam[0]}{second_base}{third_base}"
                break
            else:
                continue
    return mutated_pam.upper()

def get_target_position(genomic_loci, coding_sequence, aa_position):
    protein_dict = get_codon_index(genomic_loci, coding_sequence)
    target_dict = protein_dict[aa_position - 1]
    target_codon = f"{target_dict['base_1'][1]}{target_dict['base_2'][1]}{target_dict['base_3'][1]}"
    target_pos = (target_dict['base_1'][0], target_dict['base_2'][0], target_dict['base_3'][0])
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

    L = len(dna_sequence)
    # Iterate through codons
    codon_list = [dna_sequence[i:i+3] for i in range(0, L - 2, 3)]

    # Initialize variables
    current_list_check = []
    recodnized_sequence = ""

    # Iterate through codons and find possible amino acids
    for current_codon in codon_list:
        current_codon = current_codon.upper()
        if current_codon in ['ATG', 'TGG', 'TAA', 'TAG', 'TGA']:
            recodnized_sequence += current_codon
        else:
            current_aa = codon_dict[current_codon]
            freq_dict = {}
            for codon in aa_to_codon[current_aa]:
                freq_dict[codon] = human_frequency_dict[codon.upper()]

            max_freq_codon = max(freq_dict, key=lambda x: freq_dict[x])
            #print(freq_dict, max_freq_codon)
            recodnized_sequence += max_freq_codon

    if not L % 3 == 0:
        p = L - (L // 3) * 3
        recodnized_sequence += dna_sequence[-p:]
    #print("Recodonized Sequence:", recodnized_sequence)
    return recodnized_sequence

def get_pam_positions(genomic_loci):
    genomic_loci = genomic_loci.upper()
    pam_positions = PAMposition(genomic_loci, "NGG")
    return pam_positions

def get_identity(seqA, seqB):
    identical = sum([1 for baseA, baseB in zip(seqA, seqB) if baseA == baseB])
    return (identical / len(seqA)) * 100

def get_sgrna_df(amino_acid_position, maximum_gRNA_distance, PAM,
                 genomic_loci, coding_sequence, reference_genome):
    # Run the cripr-tape specific function
    gRNA_df = specific_function(amino_acid_position, maximum_gRNA_distance, PAM,
                                   genomic_loci, coding_sequence, "", "", reference_genome)
    gRNA_df.to_csv("grna.csv", index=False)
    #print(gRNA_df)
    #gRNA_df = pd.read_csv(gRNA_filename, index_col=0, encoding='utf-8')
    # Drop rows with all NaN values
    #gRNA_df.dropna(how='all', axis=0, inplace=True, ignore_index=True)
    #gRNA_df = filter_grna(gRNA_df)  # filter gRNA_df
    gRNA_df = gRNA_df[gRNA_df['PAM'].str.contains("G", case=False)]
    gRNA_df['gRNA Sequence'] = gRNA_df['gRNA Sequence'].str.upper()
    return gRNA_df

def filter_grna(grna_df):
    # filter by off target
    #grna_df = grna_df[grna_df['Off Target Count'] == 0]
    # filter by GC%
    #print(grna_df)
    grna_df['G/C Content (%)'] = grna_df['G/C Content (%)'].astype(float)
    grna_df = grna_df[grna_df['G/C Content (%)'].between(40, 70)]
    # filter by tm
    grna_df['temp'] = grna_df['gRNA Sequence'].apply(lambda seq: MeltingTemp.Tm_NN(seq))
    grna_df = grna_df[grna_df['temp'].between(40, 70)]
    # add leading G
    #grna_df['gRNA Sequence'] = grna_df['gRNA Sequence'].apply(lambda seq: seq if seq.startswith("G") else "G"+seq)
    grna_df['gRNA Sequence'] = grna_df['gRNA Sequence'] + grna_df['PAM']
    return grna_df['gRNA Sequence'].to_list()

def get_core_mutational_template(genomic_loci, target_position):

    coding_seq = get_exons(genomic_loci)
    genome_seq = list(genomic_loci.upper())
    target_codon, indices = get_target_position(genomic_loci, get_exons(genomic_loci), target_position)
    target_codon = mutate_target_codon(target_codon)
    for i, index in enumerate(indices):
        genome_seq[index] = target_codon[i]
    genome_seq = "".join(genome_seq)
    pam_positions = PAMposition(genome_seq, "NGG")
    pam_dict = {pos : (grna, strand) for pos, grna, strand in zip(*pam_positions)}
    f0, fn = indices[0] - 40, indices[0] + 40
    pos_list = [pos for pos in pam_positions[0] if pos>f0 and pos<fn]

    posi = [pos for pos in pam_positions[0] if (pos > f0 and pos < indices[0]) and pam_dict[pos][1] == "reverse"]
    posf = [pos for pos in pam_positions[0] if pos > indices[0] and pos < fn]

    result = {}

    for start, end in product(posi, posf):
        grnas = (pam_dict[start][0], pam_dict[end][0])
        gcf, gcr = gc_fraction(grnas[0])*100, gc_fraction(grnas[1])*100
        #print(gcf, gcr)
        tempf, tempr = MeltingTemp.Tm_NN(grnas[0]), MeltingTemp.Tm_NN(grnas[1])

        if pam_dict[start][1] == "forward":
            e1 = start + 6
            s1 = e1 - 23
            grna1 = complement(genome_seq[s1:e1-3])
            fpam = mutate_pam(genomic_loci[e1-3 : e1])
        else:
            s1 = start - 5
            e1 = s1 + 23
            grna1 = reverse_complement(genome_seq[s1+3:e1])
            fpam = mutate_pam(genomic_loci[s1:s1+3])

        if pam_dict[end][1] == "forward":
            e2 = end + 6
            s2 = e2 - 23
            grna2 = complement(genome_seq[s2:e2-3])
            rpam = mutate_pam(genomic_loci[e2-3 : e2])
        else:
            s2 = end - 5
            e2 = s2 + 23
            grna2 = reverse_complement(genome_seq[s2+3:e2])
            rpam = mutate_pam(genomic_loci[s2 : s2+3])

        length = e2 - s1 - 6

        if length >= 50 and length <= 70:
            dna = f"{grna1}{genome_seq[e1:s2]}{grna2}"
            recodonized_dna = recodonize_sequence(dna)
            identity = get_identity(dna, recodonized_dna)
            result[identity] = (grnas, s1, e1, s2, e2)

    expected_grna, start, e1, s2, end = sorted(result.items(), key=lambda x: x[0])[0][1]
    #print(expected_grna)
    leftL = 50 - (indices[0] - start)
    rightL = 50 - (end - indices[0])

    leftS = genome_seq[start - leftL : start]
    leftS = recodonize_sequence(leftS)
    rightS = genome_seq[end : end + rightL]
    rightS = recodonize_sequence(rightS)

    donor_template = f"{fpam}{grna1}{leftS}{genome_seq[e1:s2]}{rightS}{grna2}{rpam}"

    homoL = start - leftL
    homoR = end + rightL
    up_homolog = genome_seq[homoL - 80 : homoL]
    #up_homolog = recodonize_sequence(up_homolog)
    down_homolog = genome_seq[homoR : homoR + 80]
    #down_homolog = recodonize_sequence(down_homolog)
    
    priL = homoL - 80
    priR = homoR + 80
    up_primer = genome_seq[priL - 20 : priL]
    #up_homolog = recodonize_sequence(up_homolog)
    down_primer = genome_seq[priR : priR + 20]
    #down_homolog = recodonize_sequence(down_homolog)

    core_mutational_template = f"{up_primer}{up_homolog}{donor_template}{down_homolog}{down_primer}"
    #print(core_mutational_template)
    return (core_mutational_template, donor_template, up_homolog, down_homolog,
            up_primer, down_primer, expected_grna[0], expected_grna[1], fpam, rpam)

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
    #coding_sequence = get_exons(genomic_loci)
    reference_genome = read_seq(args.reference_genome)

    core_mutational_template, donor_template, up_homolog, down_homolog, up_primer, down_primer, fgrna, rgrna, fpam, rpam = get_core_mutational_template(genomic_loci, args.cut_site)

    with open(args.outfile, "w") as out_obj:
        out_obj.write(f'Donor Template: {donor_template}\tLength: {len(donor_template)}\n')
        out_obj.write(f'Upstream primer: {up_primer}\tLength: {len(up_primer)}\n')
        out_obj.write(f'Downstream primer: {down_primer}\tLength: {len(down_primer)}\n')
        out_obj.write(f'Upstream homologous: {up_homolog}\tLength: {len(up_homolog)}\n')
        out_obj.write(f'Downstream homologous: {down_homolog}\tLength: {len(down_homolog)}\n')
        out_obj.write(f'Upstream gRNA: {fgrna}\tLength: {len(fgrna)}\n')
        out_obj.write(f'Downstream gRNA: {rgrna}\tLength: {len(rgrna)}\n')
        out_obj.write(f'Upstream PAM: {fpam}\tLength: {len(fpam)}\n')
        out_obj.write(f'Downstream PAM: {rpam}\tLength: {len(rpam)}\n')
        out_obj.write(f"Core Mutational Template: {core_mutational_template}\tLength: {len(core_mutational_template)}\n\n")

if __name__ == '__main__':
    # execute the whole script
    main()
    print("\nThe Script has Executed Successfully!\n")
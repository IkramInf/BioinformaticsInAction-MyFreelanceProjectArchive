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

codon_dict = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'B', 'AGT': 'B', 'AGA': 'D', 'AGG': 'D',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'O', 'CGC': 'O', 'CGG': 'O', 'CGT': 'O',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'J', 'TCC': 'J', 'TCG': 'J', 'TCT': 'J',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TGC': 'C', 'TGT': 'C',
        'TGG': 'W', 'TAA': 'ter', 'TAG': 'ter', 'TGA': 'ter'}
#print(codon_dict)
codon_table = CodonTable.CodonTable()
for codon, aa in codon_dict.items():
    if aa == 'ter':
        codon_table.forward_table[codon] = '*'
    else:
        codon_table.forward_table[codon] = aa

aa_to_codon = {
        'I':['ATA', 'ATT', 'ATC'], # s to l
        'M': 'ATG',
        'T':['ACT', 'ACC', 'ACA', 'ACG'], # s to l
        'N':['AAC', 'AAT'], 
        'K':['AAA', 'AAG'],
        'B':['AGC', 'AGT'],
        'D':['AGA', 'AGG'],
        'L':['CTT', 'CTG', 'CTC', 'CTA', 'TTA', 'TTG'], #s to l
        'P':['CCA', 'CCC', 'CCT', 'CCG'], #s to l
        'H':['CAC', 'CAT'],
        'Q':['CAA', 'CAG'],
        'O':['CGT', 'CGG', 'CGA', 'CGC'],
        'V':['GTA', 'GTT', 'GTG', 'GTC'],
        'A':['GCT', 'GCA', 'GCC', 'GCG'], #s to L
        'D':['GAC', 'GAT'],
        'E':['GAA', 'GAG'],
        'G':['GGG', 'GGT', 'GGA', 'GGC'],
        'J':['TCA', 'TCC','TCG', 'TCT'],
        'F':['TTC', 'TTT'],
        'Y':['TAC', 'TAT'],
        'C':['TGC', 'TGT'],
        'W': "TGG",
        'ter' : ['TAA', 'TAG', 'TGA']}

human_frequency_dict = {
        'ATA': 0.16, 'ATC': 0.48, 'ATT': 0.36, 
        'ATG': 1.00,
        'ACA': 0.28, 'ACC': 0.36, 'ACG': 0.12, 'ACT': 0.24,
        'AAC': 0.54, 'AAT': 0.46, 
        'AAA': 0.42, 'AAG': 0.58,
        'AGC': 0.24, 'AGT': 0.15,
        'AGA': 0.20, 'AGG': 0.20,
        'CTA': 0.07, 'CTC': 0.20, 'CTG': 0.41, 'CTT': 0.13,
        'CCA': 0.27, 'CCC': 0.33, 'CCG': 0.11, 'CCT': 0.28,
        'CAC': 0.59, 'CAT': 0.41, 
        'CAA': 0.25, 'CAG': 0.75,
        'CGA': 0.11, 'CGC': 0.19, 'CGG': 0.21, 'CGT': 0.08,
        'GTA': 0.11, 'GTC': 0.24, 'GTG': 0.47, 'GTT': 0.18,
        'GCA': 0.23, 'GCC': 0.40, 'GCG': 0.11, 'GCT': 0.26,
        'GAC': 0.54, 'GAT': 0.46,
        'GAA': 0.42, 'GAG': 0.58,
        'GGA': 0.25, 'GGC': 0.34, 'GGG': 0.25, 'GGT': 0.16,
        'TCA': 0.15, 'TCC': 0.22, 'TCG': 0.06, 'TCT': 0.18,
        'TTC': 0.55, 'TTT': 0.45, 
        'TTA': 0.07, 'TTG': 0.13,
        'TAC': 0.57, 'TAT': 0.43, 
        'TGC': 0.55, 'TGT': 0.45,
        'TGG': 1.00}

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
    # if pam in coding sequence
    synonymous_codons = []
    if pam.isupper():
        amino_acid = codon_dict[pam]
        synonymous_codons = aa_to_codon[amino_acid]
    second_base = random.choice(['C', 'T'])
    mutated_pam = f"{pam[0]}{second_base}{pam[2]}"
    if not mutated_pam in synonymous_codons:
        return mutated_pam.upper()
    else:
        bases = [base for base in "ACT" if base != mutated_pam[2]]
        third_base = random.choice(bases)
        return f"{mutated_pam[:2]}{third_base}".upper()

def get_target_position(genomic_loci, coding_sequence, aa_position):
    protein_dict = get_codon_index(genomic_loci, coding_sequence)
    target_dict = protein_dict[aa_position - 1]
    target_codon = f"{target_dict['base_1'][1]}{target_dict['base_2'][1]}{target_dict['base_3'][1]}"
    target_pos = (target_dict['base_1'][0], target_dict['base_2'][0], target_dict['base_3'][0])
    return (target_codon, target_pos)

def mutate_target_codon(target_codon):
    amino_acid = translate(target_codon, table=codon_table)

    abundant_codons = {}
    for codon in aa_to_codon[amino_acid]:
        abundant_codons[codon] = human_frequency_dict[codon]

    most_abundant_codon = max(abundant_codons, key=lambda x: abundant_codons[x])
    return most_abundant_codon

def recodonize_sequence(dna_seq):
    """
    Recodonise dna sequence to most abundant codons but not changing protein sequence
    """
    exon = get_exons(dna_seq)
    #changing codons to aa
    aa_sequence = ""
    exon_dict = {}
    for i in range(0, len(exon)-2, 3):
        codon = exon[i : i+3]
        aa = codon_dict[codon]
        aa_sequence += aa
        exon_dict[i // 3] = codon

    recodnized_sequence = ''
    for n, aa in enumerate(aa_sequence):
        current_codon = exon_dict[n]

        if aa == 'M' or aa == "W":
            recodnized_sequence += current_codon
        else:
            aa_dict = {codon : human_frequency_dict[codon] for codon in aa_to_codon[aa]}
            codon_freq_list = sorted(aa_dict.items(), key=lambda x: x[1], reverse=True)
            for i, codon_freq in enumerate(codon_freq_list):
                if (codon_freq[0] == current_codon) and (i == 0):
                    recodnized_sequence += codon_freq_list[1][0]  # lower down a rank
                elif codon_freq[0] == current_codon:
                    recodnized_sequence += codon_freq_list[i-1][0]  # upper a rank
    # return the recodonized sequence
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
    grna_df = grna_df[grna_df['Off Target Count'] == 0]
    # filter by GC%
    grna_df['G/C Content (%)'] = grna_df['G/C Content (%)'].astype(float)
    grna_df = grna_df[grna_df['G/C Content (%)'].between(50, 75)]
    # filter by tm
    grna_df['temp'] = grna_df['gRNA Sequence'].apply(lambda seq: MeltingTemp.Tm_NN(seq))
    grna_df = grna_df[grna_df['temp'].between(40, 70)]
    # add leading G
    #grna_df['gRNA Sequence'] = grna_df['gRNA Sequence'].apply(lambda seq: seq if seq.startswith("G") else "G"+seq)
    grna_df['gRNA Sequence'] = grna_df['gRNA Sequence'] + grna_df['PAM']
    return grna_df

def get_core_mutational_template(genomic_loci, grnas, target_position):

    coding_seq = get_exons(genomic_loci)
    genome_seq = genomic_loci.upper()
    #genome_seq = list(genomic_loci.upper())

    target_codon, indices = get_target_position(genomic_loci, get_exons(genomic_loci), target_position)
    target_codon = mutate_target_codon(target_codon)
    #for i, index in enumerate(indices):
    #    genome_seq[index] = target_codon[i]
    #genome_seq = "".join(genome_seq)
    pam_positions = PAMposition(genomic_loci, "NGG")

    pam_df = pd.DataFrame(zip(*pam_positions), columns=['pos', 'grna', 'strand'])
    pam_df['grna_seq'] = pam_df['grna'].apply(lambda seq: seq[:-3])
    pam_df = pam_df[pam_df['grna_seq'].isin(grnas)]
    if pam_df.empty:
        print("The Entered gRNA is not Found. Exitting Code!!!")
        sys.exit()
    pam_df.drop("grna_seq", axis=1, inplace=True)
    #print(pam_df)
    pam_dict = {pos : (grna, strand) for pos, grna, strand in pam_df.values}

    #f0, fn = indices[0] - 70, indices[0] + 70
    #pos_list = [pos for pos in pam_df['pos'] if pos>f0 and pos<fn]
    #posi = [pos for pos in pam_df['pos'] if pos > f0 and pos < indices[0]]
    #posf = [pos for pos in pam_df['pos'] if pos > indices[0] and pos < fn]
    #result = {}
    #for start, end in product(posi, posf):
    #grnas = (pam_dict[start][0], pam_dict[end][0])
    #gcf, gcr = gc_fraction(grnas[0])*100, gc_fraction(grnas[1])*100
    #print(gcf, gcr)
    #tempf, tempr = MeltingTemp.Tm_NN(grnas[0]), MeltingTemp.Tm_NN(grnas[1])

    for pos, grna, strand in pam_df.values:
        if grna.startswith(grnas[0]):
            start = pos
        elif grna.startswith(grnas[1]):
            end = pos
    try:
        if pam_dict[start][1] == "forward":
            #e1 = start + 7
            #s1 = e1 - 23
            s1 = genome_seq.index(pam_dict[start][0])
            e1 = s1 + 23
            grna1 = genome_seq[s1:e1-3]
            fpam = mutate_pam(genomic_loci[e1-3 : e1])
        else:
            #s1 = start - 5
            #e1 = s1 + 23
            s1 = genome_seq.index(reverse_complement(pam_dict[start][0]))
            e1 = s1 + 23
            grna1 = reverse_complement(genome_seq[s1+3:e1])
            fpam = mutate_pam(reverse_complement(genomic_loci[s1:s1+3]))

        if pam_dict[end][1] == "forward":
            #e2 = end + 7
            #s2 = e2 - 23
            s2 = genome_seq.index(pam_dict[end][0])
            e2 = s2 + 23
            grna2 = genome_seq[s2:e2-3]
            rpam = mutate_pam(genomic_loci[e2-3 : e2])
        else:
            #s2 = end - 5
            #e2 = s2 + 23
            s2 = genome_seq.index(reverse_complement(pam_dict[end][0]))
            e2 = s2 + 23
            grna2 = reverse_complement(genome_seq[s2+3:e2])
            rpam = mutate_pam(reverse_complement(genomic_loci[s2 : s2+3]))
    except ValueError:
        print("The Entered gRNA is not Found. Exitting Code!!!")
        sys.exit()

    #length = e2 - s1 - 6
    #if length >= 50 and length <= 70:
    #dna = f"{grna1}{genome_seq[e1:s2]}{grna2}"
    #recodonized_dna = recodonize_sequence(dna)
    #identity = get_identity(dna, recodonized_dna)
    #result[identity] = (grnas, s1, e1, s2, e2)
    #expected_grna, start, e1, s2, end = sorted(result.items(), key=lambda x: x[0])[0][1]
    #print(expected_grna)

    # Get Donor Template
    if (indices[0] > s1) and (indices[0] < e2):
        #leftL = 50 - (indices[0] - s1)
        #rightL = 50 - (e2 - indices[0])
        print("\nWarning! There may be an issue with the selection of the guide RNA (gRNA). As, excision sites directed by the two gRNAs doesn't contain target region.")

    d0 = 100 - (e2 - s1)
    leftL = d0 // 2
    rightL = d0 - leftL
    #rightL = leftL

    leftS = genome_seq[s1 - leftL : s1]
    #leftS = recodonize_sequence(leftS)
    rightS = genome_seq[e2 : e2 + rightL]
    #rightS = recodonize_sequence(rightS)
    donor_template = f"{grna1}{genome_seq[s1+23 : e2-23]}{grna2}"
    nd = len(donor_template)
    if nd % 3 != 0:
        r = nd % 3
        donor_template = f"{recodonize_sequence(donor_template[:-r])}{donor_template[-r:]}"
    else:
        donor_template = recodonize_sequence(donor_template)

    fstart, fend = start - (leftL + 80 + 20), end + rightL + 80 + 20
    # Get Upstream Homolog and Primer
    homoL = start - leftL
    if (homoL - 80) < 0:
        if homoL < 0:
            partial_seq = ""
            LEN = 80
        else:
            partial_seq = genome_seq[0:homoL]
            LEN = 80 - len(partial_seq)
        up_homolog = f"{genome_seq[fend : fend + LEN]}{partial_seq}{leftS}"
        fend = fend + LEN
    else:
        up_homolog = f"{genome_seq[homoL - 80 : homoL]}{leftS}"

    priL = homoL - 80
    if (priL - 20) < 0:
        if priL < 0:
            partial_seq = ""
            LEN = 20
        else:
            partial_seq = genome_seq[0:priL]
            LEN = 20 - len(partial_seq)
        up_primer = f"{genome_seq[fend : fend + LEN]}{partial_seq}"
    else:
        up_primer = genome_seq[priL - 20 : priL]

    N = len(genome_seq)
    # Get Downstream Homolog and Primer
    homoR = end + rightL
    if homoR + 80 > N:
        if homoR > N:
            partial_seq = ""
            LEN = 80
        else:
            partial_seq = genome_seq[homoR:]
            LEN = len(partial_seq)
        down_homolog = f"{rightS}{partial_seq}{genome_seq[fstart - LEN : fstart]}"
        fstart = fstart - LEN
    else:
        down_homolog = f"{rightS}{genome_seq[homoR : homoR + 80]}"

    priR = homoR + 80
    if priR + 20 > N:
        if priR > N:
            partial_seq = ""
            LEN = 20
        else:
            partial_seq = genome_seq[priR:]
            LEN = len(partial_seq)
        down_primer = f"{partial_seq}{genome_seq[fstart - LEN : fstart]}"
    else:
        down_primer = genome_seq[priR : priR + 20]

    core_mutational_template = f"{up_primer}{up_homolog}{fpam}{donor_template}{rpam}{down_homolog}{down_primer}"
    #print(core_mutational_template)
    return (core_mutational_template, donor_template, up_homolog, down_homolog,
            up_primer, down_primer, grnas[0], grnas[1], fpam, rpam)

def get_args():
    parser = argparse.ArgumentParser(description="Generate repair donor.py")
    parser.add_argument("-f", "--foreign_DNA", help="foreign DNA", required=True)
    #parser.add_argument("-r", "--reference_genome", help="reference genome", required=True)
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
    #reference_genome = read_seq(args.reference_genome)
    #amino_acid_position, maximum_gRNA_distance, PAM = args.cut_site, 1000, "NGG"
    #gRNA_df = get_sgrna_df(amino_acid_position, maximum_gRNA_distance, PAM,
    #                               genomic_loci, coding_sequence, reference_genome)
    # filter grna df
    #gRNA_df = filter_grna(gRNA_df)
    #gRNA1 = "GCAAAGAGCAGGACGGGTAC".upper()
    #gRNA2 = "GGTAAAACTGGGTGTCTTTG".upper()
    
    gRNA1 = input("Enter Upstream gRNA Sequence: ").upper()
    if len(gRNA1) == 23:
        gRNA1 = gRNA1[:-3]

    gRNA2 = input("Enter Downstream gRNA Sequence: ").upper()
    if len(gRNA2) == 23:
        gRNA2 = gRNA2[:-3]

    core_mutational_template, donor_template, up_homolog, down_homolog, up_primer, down_primer, fgrna, rgrna, fpam, rpam = get_core_mutational_template(genomic_loci, [gRNA1, gRNA2], args.cut_site)

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
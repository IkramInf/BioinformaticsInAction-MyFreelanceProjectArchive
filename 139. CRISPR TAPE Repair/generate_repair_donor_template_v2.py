import os
import re
import sys
import random
import argparse
import pandas as pd
from Bio.Data import CodonTable
from Bio.SeqUtils import MeltingTemp
from Bio.Seq import complement, reverse_complement, translate
from CRISPR_TAPE.shared_functions import get_codon_index, PAMposition
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

def mutate_pam(pam, gene_sequence, pam_positions):
    """
    Perform mutation in pam sequence
    """
    pam_dict = {seq[-3:] : (pos, strand) for pos, seq, strand in zip(*pam_positions)}
    position, strand = pam_dict[pam]
    if strand == "forward":
        pam = gene_sequence[position - 1 + 5 : position + 2 + 5]
    else:
        pam = gene_sequence[position - 1 - 4: position + 1 - 3]

    # True: coding sequence
    if pam.isupper():
        amino_acid = codon_dict[pam]
        synonymous_codons = aa_to_codon[amino_acid]
        synonymous_codons = [codon for codon in synonymous_codons if (not codon.endswith("AG")) and (codon != pam)]
        if synonymous_codons:
            mutated_pam = random.choice(synonymous_codons)
        else:
            mutated_pam = ""
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
            current_aa = codon_dict[current_codon.upper()]
            freq_dict = {}
            for codon in aa_to_codon[current_aa]:
                freq_dict[codon] = human_frequency_dict[codon.upper()]

            max_freq_codon = max(freq_dict, key=lambda x: freq_dict[x])
            #print(freq_dict, max_freq_codon)
            recodnized_sequence += max_freq_codon

    #print("Recodonized Sequence:", recodnized_sequence)
    return recodnized_sequence

def get_pam_positions(genomic_loci):
    genomic_loci = genomic_loci.upper()
    pam_positions = PAMposition(genomic_loci, "NGG")
    return pam_positions

def get_grna_cut_sites(genomic_loci, pam_positions, grna_guides):
    #print(grna_guides)
    grna_dict = {}
    for position, grna, strand in zip(*pam_positions):
        #print(grna)
        if grna in grna_guides:
            if strand == "forward":
                start, end = position - 21 + 5, position + 2 + 5
                cut_grna = genomic_loci[start : end]
                #print(cut_grna == grna)
                grna_dict.setdefault(strand, []).append((grna, start, end))
            else:
                start, end = position - 1 - 4, position + 21 - 3
                cut_grna = reverse_complement(genomic_loci[start : end])
                #print(cut_grna == grna)
                grna_dict.setdefault(strand, []).append((grna, start, end))

    return grna_dict

def get_excised_sequences(genomic_loci, grna_cut_sites):
    excised_sequences = []

    for fgrna in grna_cut_sites['forward']:
        fgrna, fstart, fend = fgrna
        f_grna, f_pam = fgrna[:-3], fgrna[-3:]
        for rgrna in grna_cut_sites['reverse']:
            rgrna, rstart, rend = rgrna
            r_grna, r_pam = rgrna[:-3], rgrna[-3:]
            #print(fgrna, fstart, fend)
            if fend < rstart:
                L = rstart - fend
                if (L <= 54) and (L >= 30):
                    start, end = fend, rstart
                    seq = genomic_loci[start : end]
                    excised_sequences.append((f_pam, f_grna, start, end, r_grna, r_pam, seq, "forward"))
            elif rend < fstart:
                L = fstart - rend
                if (L <= 54) and (L >= 30):
                    start, end = rend, fstart
                    seq = genomic_loci[start : end]
                    excised_sequences.append((f_pam, f_grna, start, end, r_grna, r_pam, seq, "reverse"))

    return excised_sequences

def get_donor_template(excised_sequence, genomic_loci, target_pos, mutated_target_codon, pam_positions):
    """
    Get the donor template for 300 bp sequence
    """
    fpam, fgrna, cutsite1, cutsite2, rgrna, rpam, internal_seq, strand = excised_sequence
    fgrna = reverse_complement(fgrna)
    fgrna = fgrna if fgrna.startswith("G") else "G" + fgrna[1:]
    rgrna = rgrna if rgrna.startswith("G") else "G" + rgrna[1:]
    fpam = reverse_complement(mutate_pam(fpam, genomic_loci, pam_positions))
    rpam = mutate_pam(rpam, genomic_loci, pam_positions)
    index = genomic_loci.index(internal_seq)
    internal_seq = recodonize_sequence(internal_seq)
    donor_template = f"{fpam}{fgrna}{internal_seq}{rgrna}{rpam}"
    
    target_base1 = target_pos[0] - index + 23
    if target_base1 < -23:
        return (None, fpam, rpam)
    target_base2 = target_pos[1] - index + 23
    target_base3 = target_pos[2] - index + 23
    
    positions = [target_base1, target_base2, target_base3]
    for base, position in zip(mutated_target_codon, positions):
        donor_template = donor_template[:position] + base + donor_template[position+1:]

    return (donor_template, fgrna, rgrna, fpam, rpam)

def get_homology_regions(excised_sequence, genomic_loci, reference_genome, homology_arm_length=80):
    
    try:
        genomic_loci_index = reference_genome.index(genomic_loci)
    except ValueError:
        genomic_loci_index = 0

    fpam, fgrna, cutsite_1, cutsite_2, rgrna, rpam, seq, strand = excised_sequence
    
    if strand == "forward":
        up_end = genomic_loci_index + cutsite_1
        up_start = up_end - homology_arm_length
        if up_start < 0:
            up_start = 0
        up_homolog = reference_genome[up_start : up_end]
        
        down_start = genomic_loci_index + cutsite_2
        down_end = down_start + homology_arm_length
        down_homolog = reference_genome[down_start : down_end]
    else:
        up_start = genomic_loci_index + cutsite_1
        up_end = up_start + homology_arm_length
        up_homolog = reference_genome[up_start : up_end]
        
        down_end = genomic_loci_index + cutsite_2
        down_start = down_end - homology_arm_length
        if down_start < 0:
            down_start = 0
        down_homolog = reference_genome[down_start : down_end]

    return (up_homolog, down_homolog)

def get_primer(reference_genome, genomic_loci, target_pos, aa_position=10, primer_length=20):

    try:
        genomic_loci_index = reference_genome.index(genomic_loci)
    except ValueError:
        genomic_loci_index = 0

    up_end = genomic_loci_index + target_pos[0]
    up_start = up_end - primer_length
    if up_start < 0:
        up_start = 0
    up_primer = reference_genome[up_start : up_end]

    down_start = genomic_loci_index + target_pos[2]
    down_end = down_start + primer_length
    down_primer = reference_genome[down_start : down_end]

    return (up_primer, down_primer)

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

def get_sgrna_df(amino_acid_position, maximum_gRNA_distance, PAM,
                 genomic_loci, coding_sequence, reference_genome):
    # Run the cripr-tape specific function
    gRNA_df = specific_function(amino_acid_position, maximum_gRNA_distance, PAM,
                                   genomic_loci, coding_sequence, "", "", reference_genome)
    #print(gRNA_df)
    #gRNA_df = pd.read_csv(gRNA_filename, index_col=0, encoding='utf-8')
    # Drop rows with all NaN values
    #gRNA_df.dropna(how='all', axis=0, inplace=True, ignore_index=True)
    #gRNA_df = filter_grna(gRNA_df)  # filter gRNA_df
    gRNA_df = gRNA_df[gRNA_df['PAM'].str.contains("G", case=False)]
    return gRNA_df

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
    #print(df)
    #check_grna(genomic_loci, coding_sequence, df, args.cut_site)
    #grna_list_pairs = get_sgrna_pairs(df)
    # filter grna df
    grna_guides = filter_grna(df)
    #print(grna_guides)
    pam_positions = get_pam_positions(genomic_loci)
    grna_cut_sites = get_grna_cut_sites(genomic_loci, pam_positions, grna_guides)
    #print(grna_cut_sites)
    if grna_cut_sites:
        excised_sequences = get_excised_sequences(genomic_loci, grna_cut_sites)
        #excised_sequences = get_excised_sequences(genomic_loci, df)
        target_codon, target_pos = get_target_position(genomic_loci, coding_sequence, args.cut_site)
        mutated_target_codon = mutate_target_codon(target_codon)
        up_primer, down_primer = get_primer(reference_genome, genomic_loci, target_pos)
        template_dict = {}  # store core mutational template info
        for excised_sequence in excised_sequences:
            donor_template, fgrna, rgrna, fpam, rpam = get_donor_template(excised_sequence, genomic_loci, target_pos, mutated_target_codon, pam_positions)
            if donor_template:
                up_homolog, down_homolog = get_homology_regions(excised_sequence, genomic_loci, reference_genome)
                core_mutational_template = f"{up_primer}{up_homolog}{donor_template}{down_homolog}{up_primer}"
                parts = [donor_template, up_primer, down_primer, up_homolog, down_homolog, fgrna, rgrna, fpam, rpam]
                template_dict[core_mutational_template] = parts

    core_mutational_template = max(template_dict, key=lambda k: len(k))
    #print(core_mutational_template, len(core_mutational_template))
    with open(args.outfile, "w") as out_obj:
        # output
        donor_template, up_primer, down_primer, up_homolog, down_homolog, fgrna, rgrna, fpam, rpam = template_dict[core_mutational_template]
        out_obj.write(f'Donor Template: {donor_template}\tLength: {len(donor_template)}\n')
        out_obj.write(f'Upstream primer: {up_primer}\tLength: {len(up_primer)}\n')
        out_obj.write(f'Downstream primer: {down_primer}\tLength: {len(down_primer)}\n')
        out_obj.write(f'Upstream homologous: {up_homolog}\tLength: {len(up_homolog)}\n')
        out_obj.write(f'Downstream homologous: {down_homolog}\tLength: {len(down_homolog)}\n')
        out_obj.write(f'Upstream gRNA: {fgrna}\tLength: {len(fgrna)}\n')
        out_obj.write(f'Downstream gRNA: {rgrna}\tLength: {len(rgrna)}\n')
        out_obj.write(f'Upstream PAM: {fpam}\tLength: {len(fpam)}\n')
        out_obj.write(f'Downstream PAM: {rpam}\tLength: {len(rpam)}\n')
        out_obj.write(f"Core Mutational Template: {core_mutational_template}\tLength: {len(core_mutational_template)}\n")

if __name__ == '__main__':
    # execute the whole script
    main()
    print("\nThe Script has Executed Successfully!\n")
translate_seq : Translate a dna or rna sequence into protein in six frames

## About
translate seq module contains several helper function with translate function to print similar output as EMBOSS Sixpack on screen.

## Functions explained
    [1] file_handler : Handle File opening errors 
    [2] extract_codon_table : Extract codon table from text file
    [3] read_fasta : Read fasta file
    [4] complement : Complement a dna sequence
    [5] triplets : Make list of codons from a dna sequence
    [6] grouping : Make list of sixty characters substring from a dna sequence
    [7] check_length : check the length of each protein and return all of them as a list
    [8] count_orf : Count orf in a protein sequence
    [9] translate : Translate dna sequence into protein in six frames
    [10] format_output : Format proteins and sequence to print on screen

## Usage
from translate_seq import read_fasta, translate
sequence = read_fasta("text_file_in_fasta_format")
translate(sequence) # print output on screen

please keep a codon_table.txt file in your directory.
Format of codon_table.txt : codon three_letter_amino_acid
e.g., aaa LYS
    
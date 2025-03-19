#!/usr/bin/env python
"""
Code takes a dna sequence as input, then:
 > transcribe dna into rna
 > find first open reading frame in the rna sequence
 > translate the orf into protein
 > finally print the result to screen
"""

# import required libraries
import re
import argparse
from Bio.Seq import Seq
from Bio import SeqIO

def get_args():
    """Return parsed command-line arguments."""

    parser = argparse.ArgumentParser(
        description="TODO: say what the script does.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # get the FASTA file of sequences
    parser.add_argument('filename',  # variable to access this data later: args.filename
                        metavar='FASTA', # shorthand to represent the input value
                        help='Provide name and path to FASTA file to process.', # message to the user, it goes into the help menu
                        type=str)
    parser.add_argument('-p', '--pattern',  # access with args.pattern
                        help='Provide a regex pattern for filtering FASTA entries',
                        default='^\d{1}\D*$')  # default works for Drosophila chromosomes

    return(parser.parse_args())


def find_first_orf(rna):
    """Return first open-reading frame of RNA sequence as a Bio.Seq object.

    Must start with AUG
    Must end with UAA, UAG, or UGA
    Must have even multiple of 3 RNA bases between
    """
    try:
        # regex to find the ORF
        orf = re.search("AUG(?:.{3})+?(?:UAA|UAG|UGA)", str(rna)).group(0)
    except AttributeError:  # if no match found, orf should be empty
        orf = ""
    return Seq(orf)


def translate_first_orf(dna):
    """
    Parameters:
            dna : a dna sequence as Bio.Seq object
                
    Returns:
            translate_first_orf transcribe dna into rna, find first orf with help of find_first_orf()
            and finally return translated orf
    """
    
    # transcribe dna into rna
    rna = dna.replace("T", "U")

    # find out first open reading frame using find_first_orf
    first_orf = find_first_orf(rna)

    # if orf is not empty translate it into protein
    if first_orf:
        translated_orf = first_orf.translate()
    else:
        translated_orf = ""

    # return the translated orf
    return translated_orf


if __name__ == "__main__":
    
    # get command-line arguments
    args = get_args()
    
    # get the fasta file records into a dictionary as {id : seq}
    # use the filename provided as command line argument
    sequences = {record.id : record.seq for record in SeqIO.parse(args.filename, "fastq")}
    
    # iterate over each sequence and translate the first orf into protein
    # then print the result to screen as like (Id:    translated_orf)
    for Id, sequence in sequences.items():
        translated_orf = translate_first_orf(sequence)
        if translated_orf:
            print(f"{Id}:\t{translated_orf}")

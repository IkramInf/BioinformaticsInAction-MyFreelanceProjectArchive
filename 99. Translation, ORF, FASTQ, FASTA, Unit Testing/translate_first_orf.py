#!/usr/bin/env python
"""
translate_first_orf.py : translate first orf of dna into protein
    get_args()
        get command line arguments (filename and pattern)
        
    find_first_orf()
        find out first orf from rna sequence through regex pattern
    
    translate_first_orf():
        transcribe dna into rna then call find_first_orf() to find out first orf and finally translate orf into protein
"""

import re  # handle regex
import argparse  # command line parsing module
from Bio import Seq, SeqIO  # to work with sequence

def get_args():
    """Return parsed command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Get filename and regex pattern to translate dna into protein",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # get the FASTA file of sequences
    parser.add_argument('filename',  # variable to access this data later: args.filename
                        metavar='FASTA', # shorthand to represent the input value
                        help='Provide name and path to FASTA file to process.', # message to the user, it goes into the help menu
                        type=str)
    parser.add_argument('-p', '--pattern',  # access with args.pattern
                        help='Provide a regex pattern for filtering FASTA entries',
                        default='aug(?:.{3})+?(?:uaa|uag|uga)')  # default works for Drosophila chromosomes

    return(parser.parse_args())


def find_first_orf(RNA):
    """Return first open-reading frame of RNA sequence as a Bio.Seq object.

    Must start with AUG
    Must end with UAA, UAG, or UGA
    Must have even multiple of 3 RNA bases between
    """
    try:
        # regex to find the ORF
        orf = re.search("aug(?:.{3})+?(?:uaa|uag|uga)", str(RNA), re.S|re.I).group(0)
    except AttributeError:  # if no match found, orf should be empty
        orf = ""
    return(Seq.Seq(orf))


def translate_first_orf(DNA):
    """
    Parameters:
            dna : a dna sequence as Bio.Seq object
                
    Returns:
            translate_first_orf transcribe dna into rna, find first orf with help of find_first_orf()
            and finally return translated orf
    """
    
    # transcribe dna into rna with transcribe()
    # find out first orf
    # and translate orf into protein and return results
    RNA = DNA.transcribe()
    ORF = find_first_orf(RNA)
    PROTEIN = ORF.translate()
    return PROTEIN


if __name__ == "__main__":
    # call get_args to access cmd arguments
    args = get_args()
    # args.filename to access cmd filename
    for record in SeqIO.parse(args.filename, "fastq"):
        PROTEIN = translate_first_orf(record.seq)
        # if the FASTA record's ID matches the regex pattern, print output
        if PROTEIN:
            print("{}:\t{}".format(record.id, PROTEIN))
            
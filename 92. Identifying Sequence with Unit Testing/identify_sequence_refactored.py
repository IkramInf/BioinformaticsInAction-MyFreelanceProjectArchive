#!/usr/bin/env python
"""identify_sequence_refactored.py"""

"""The code is identifing the given sequence using argparse"""

import argparse


def get_args():
    """Return parsed command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Find out sequence type from list of files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # get a list of text files to process
    parser.add_argument('file_list',  # variable to access this data later: args.file_list
                        metavar='FILE', # shorthand to represent the input value
                        help='Provide file name to process. For multiple files, separate their names with spaces.', # message to the user, it goes into the help menu
                        type=str,
                        nargs="+" # will combine multiple textfile inputs into a list 
                        )

    return parser.parse_args()


def identify_sequence(sequence):
    """TODO: Say what the function does"""
    dna = {'A', 'C', 'G', 'T'}
    rna = {'A', 'C', 'G', 'U'}
    protein = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

    # check if the sequence is dna, rna or protein
    # and store nucleic acid or amino acid in the variable sequence_type
    sequence = sequence.upper()
    if set(sequence).issubset(dna):
        sequence_type = "nucleic acid"
    elif set(sequence).issubset(rna):
        sequence_type = "rna"
    elif set(sequence).issubset(protein):
        sequence_type = "amino acid"
    else:
        sequence_type = "invalid sequence"

    return sequence_type


if __name__ == "__main__":
    # get all the command line arguments
    args = get_args()
    
    # loop through args.file_list to:
    #    1) Open each file
    #    2) Identify the sequence within the file (hint: call your indentify_sequence() function)
    #    3) Print the filename, and its identity to the Terminal (don't print the sequence itself)
    for filename in args.file_list:
        with open(filename, "r") as f:
            sequence_input = f.read().replace('\n','')
            sequence_value = identify_sequence(sequence_input)
            print(f"filename: {filename}\nIdentity: {sequence_value}")

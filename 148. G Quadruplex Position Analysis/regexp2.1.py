#!/usr/bin/env python
# coding: utf-8
import os  # handling os module
import re  # handling regex
import sys  # handling cmd arguments
import argparse  # handling cmd arguments

# Define a function to get the reverse complement of a DNA sequence
def reverse_complement(genome):
    genome = genome.upper()
    base_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_complemented_seq = ''.join([base_dict[base] for base in genome[::-1] if base in base_dict])
    return reversed_complemented_seq

# Define a function to calculate the GC percentage of a DNA sequence
def get_GC(genome):
    g_counts = genome.count("G")
    c_counts = genome.count("C")
    gc_percentage = ((g_counts + c_counts) / len(genome)) * 100
    return round(gc_percentage, 2)

def get_g_quadruplex_locations(sequence, strand = "+"):
    """
    Find out all G enriched positions in a genome
    """
    locations = []
    #Gpattern = r'G{2,}'  # look for 2 or more G's

    # Gtract (G{2,}) and refer back with group number-1
    Gpattern = r'(G{2,})\w*?(?:\1)\w*?(?:\1)\w*?(?:\1)'
    Gtracts = re.finditer(Gpattern, sequence, re.I|re.S)
    # iterate over all matches
    for match in Gtracts:
        start, end = match.span()
        #print(f'G tract found at positions: {start+1}, {end}, {strand}')
        locations.append(f'{start+1}, {end}, {strand}')
    # return locations as 1-based indexing
    return locations

def read_genome(genome):
    """
    Read genome fasta file
    """
    with open(genome, "r") as ifile:
        header = ifile.readline()
        genome = ifile.read().replace("\n", "")
    return genome

def get_args():
    # Create ArgumentParser object
    parser = argparse.ArgumentParser(description="Process genome data.")

    # Add command-line arguments
    parser.add_argument("genome", help="Path to the genome file")
    parser.add_argument("output", help="Path to the output file")

    # Parse the command-line arguments
    args = parser.parse_args()
    return args

# Main block of code to be executed when the script is run
if __name__ == "__main__":
    
    args = get_args()
    # Check if the specified genome file exists
    if not os.path.exists(args.genome):
        print(f"Error: File '{args.genome}' not found.")
        sys.exit()

    genome = read_genome(args.genome)
 
    print(f"GC% for the genome: {get_GC(genome)}%")
    print(genome)
    print("====================================================")
    reversed_genome = reverse_complement(genome)
    print(reversed_genome)
    
    # write output to a file
    with open(args.output, "w") as ofile:
        ofile.write("G tract found at positions:\n")
        plus_strand_locations = get_g_quadruplex_locations(sequence=genome, strand="+")
        plus_strand_locations = "\n".join(plus_strand_locations)
        minus_strand_locations = get_g_quadruplex_locations(sequence=reversed_genome, strand="-")
        minus_strand_locations = "\n".join(minus_strand_locations)
        # write locations into file
        ofile.write(plus_strand_locations)
        ofile.write(minus_strand_locations)

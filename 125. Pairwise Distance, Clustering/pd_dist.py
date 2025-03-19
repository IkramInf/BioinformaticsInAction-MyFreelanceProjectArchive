#!/usr/bin/env python
# coding: utf-8

import os
import re
import shutil
import argparse
import subprocess
from Bio import Phylo
import matplotlib.pyplot as plt

def calculate_pairwise_distances(filename, pthreads, slidingwindow):
    """
    Calculate pairwise distances for DNA sequences using FastTree and RAxML.

    This function takes a filename containing DNA sequences in FASTA format and calculates
    pairwise distances using FastTree to generate a newick tree file and RAxML to perform 
    phylogenetic analysis. The newick tree file is saved in the "output/" directory alongside
    other intermediate files generated during the process. The output filename is derived from
    the input filename.

    Parameters:
        filename (str): The path to the input file containing DNA sequences in FASTA format.
        pthreads (str): The number of threads to use for RAxML.
        slidingwindow (str): The sliding window size for RAxML.

    Returns:
        str: The path to the generated newick tree file.
    """
    out_dir = os.path.join(os.getcwd(), "output/")
    # Check if the output directory exists
    if os.path.exists(out_dir):
        # If it exists, delete the directory and its contents
        shutil.rmtree(out_dir)
    # Create the output directory
    os.makedirs(out_dir)
    outfile = os.path.splitext(os.path.basename(filename))[0]
    newick_path = os.path.join(out_dir, outfile + ".newick")

    # Define the FastTree command to generate a newick tree file
    fasttree_command = ['FastTree', '-nt', '-gtr', filename, '>', newick_path]
    # Execute the command
    p1 = subprocess.run(" ".join(fasttree_command), shell=True)

    # Define the RAxML command as a list of strings
    raxml_command = ["raxmlHPC", "-m", "GTRGAMMA", "-fx", "-s", filename, "-n", outfile,"-p", "1234",
                     "-t", newick_path, "-T", pthreads, "-w", out_dir, "-W", slidingwindow]
    # Execute the command and set the working directory to out_dir
    p2 = subprocess.run(raxml_command)
    return newick_path

def assign_subgroup(names):
    """
    Assign subgroups based on the sequence names.

    Parameters:
        names (list): List of sequence names.

    Returns:
        dict: A dictionary with subgroups as keys and lists of sequence names as values.
    """
    subgroups = {}
    for name in names:
        # Use regular expression to extract the subgroup identifier
        match = re.search(r'_([A-Za-z]+)_', name)
        if match:
            subgroup = match.group(1)
            matched_subgroup = None
            for key in subgroups:
                if subgroup[:-3] in key:
                    matched_subgroup = key
                    break
            if matched_subgroup:
                subgroups[matched_subgroup].append(name)
            else:
                subgroups.setdefault(subgroup, []).append(name)

    return subgroups

if __name__ == "__main__":
    # Add argument parser
    parser = argparse.ArgumentParser(description="Calculate pairwise distances, generate a tree, and assign subgroups.")
    parser.add_argument('-i', '--filename', required=True, help="Input FASTA file")
    parser.add_argument('-T', '--pthreads', default='2', help="Number of threads for RAxML (default: 2)")
    parser.add_argument('-W', '--slidingwindow', default='100', help="Sliding window size for RAxML (default: 100)")
    args = parser.parse_args()
    
    # Calculate pairwise distances
    newick_filepath = calculate_pairwise_distances(args.filename, args.pthreads, args.slidingwindow)
    
    # Create a matplotlib figure and draw the tree on it
    fig, ax = plt.subplots(figsize=(20, 20))  # You can adjust the figsize as needed
    tree = Phylo.read(newick_filepath, format="newick")
    Phylo.draw(tree, do_show=False, axes=ax)

    # Save the figure as an image
    output_file = os.path.splitext(os.path.basename(args.filename))[0] + "_tree.png"
    output_file = os.path.join("output", output_file)
    plt.savefig(output_file, format="png")
    plt.close()  # Close the figure to release resources
    print(f"Tree saved as an image: {output_file}")

    # Assign subgroups based on sequence names
    names = [clade.name for clade in tree.find_clades() if clade.name]
    subgroups = assign_subgroup(names)

    # Print the subgroups and their sequences
    print("\nPrinting Subgroups:")
    print("-"*20, "\n")
    for subgroup, names in subgroups.items():
        names = ",".join(names)
        print(f"{subgroup}: {names}\n")

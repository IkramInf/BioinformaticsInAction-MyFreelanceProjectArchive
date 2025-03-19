"""
viralGenomeFinder.py perform:
    [i] retrieve genbank record using accession id from NCBI
    [ii] convert genbank record to annotation file
    [iii] create and plot codon frequency for cds sequences
"""

import os
import sys
import argparse
import requests
from io import StringIO
from Bio import SeqIO
import matplotlib.pyplot as plt

def retrieve_genomic_sequence(accession_id):
    """
    Retrieve genbank file from NCBI
    """
    email = "your_email@example.com"  # Set your email for NCBI API
    # prepare url to search for accession id
    baseURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
    term = f"db=nuccore&id={accession_id}&strand=1,2&rettype=gb&retmode=text&email={email}"
    url = baseURL + term
    # Get genbank file as response object with requests.get method
    response = requests.get(url)
    # Parse the GenBank record
    record = SeqIO.read(StringIO(response.text), "genbank")

    # Save genomic sequence as a FASTA file
    with open(f"{accession_id}.fa", "w") as fasta_file:
        #write accession number and sequence
        fasta_file.write(f">{accession_id}\n{str(record.seq)}\n")

    # return the parsed genbank record
    return record

def create_gene_annotation_file(record):
    """
    Create annotation file from genbank record
    """
    # open file to write
    with open(f"{record.id}.gtf", "w") as ofile:
        # dictionary to get +/- strand from 1/-1
        get_strand = {1 : "+", -1 : "-"}
        cds_sequences = []  # store all cds sequences
        # iterate over feature table
        for feature in record.features:
            # check for gene feature
            if feature.type == 'gene':
                # find out all columns to write into gtf file
                ID = record.id
                feature_type = feature.type
                start, end = 1, 0
                # iterate over joined location
                for location in feature.location.parts:
                    start += location.start
                    end += location.end
                score = "."
                strand = get_strand[feature.strand]
                frame = "."
                attributes = f'gene_id"{feature.qualifiers["db_xref"][0]}";locus_tag"{feature.qualifiers["locus_tag"][0]}"'
                # write each column to file
                ofile.write(f"{ID}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}\n")

            # check for cds feature
            elif feature.type == 'CDS':
                # Extract the CDS sequence from the record and append to cds_sequences
                cds_sequences.append(feature.location.extract(record.seq))
                
                # find out all columns to write into gtf file
                ID = record.id
                feature_type = feature.type
                score = "."
                frame = (len(feature) + 1) % 3
                attributes = f'gene_id"{feature.qualifiers["db_xref"][0]}";locus_tag"{feature.qualifiers["locus_tag"][0]}"'
                # iterate over joined location
                for location in feature.location.parts:
                    start = location.start + 1
                    end = location.end
                    strand = get_strand[location.strand]
                    # write each column to file
                    ofile.write(f"{ID}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}\n")
                    
    # return extracted cds sequences
    return cds_sequences

def calculate_codon_frequency(cds_sequences, accession_id):
    """
    Calculate codon frequency
    """
    # Dictionary to store codon frequencies
    codon_frequency = {}
    # iterate over all CDS sequences
    for cds_sequence in cds_sequences:
        #iterate over each codon of cds sequence
        for i in range(0, len(cds_sequence) - 2, 3):
            # make codons from sequence
            codon = str(cds_sequence[i:i+3])
            # increase count 1 for codon
            codon_frequency[codon] = codon_frequency.get(codon, 0) + 1
            
    # calculate sum of counts
    total = sum(codon_frequency.values())
    # convert count value to frequency percentage
    codon_frequency = {codon: (count / total) * 100 for codon, count in codon_frequency.items()}

    #open txt file to write Codon table
    with open(f"{accession_id}.txt", "w") as codon_file:
        # write header in file
        codon_file.write("Codon\tFrequency(%)\n")
        for codon, frequency in codon_frequency.items():
            # write codon and frequency
            codon_file.write(f"{codon}\t{frequency:.4f}\n")

    # return the codon frequency dict
    return codon_frequency

def calculate_gc_percentage(sequence):
    #count 'G' and count 'C' then divide it by length of sequence and multiply by 100
    gc_content = (sequence.count("G") + sequence.count("C")) / len(sequence) * 100
    return gc_content

def plot_codon_frequency(codon_frequencies, accession_id):
    """
    Plot codon frequency as bar plot
    """
    # Extract codons and frequencies
    codons = list(codon_frequencies.keys())
    frequencies = list(codon_frequencies.values())

    # Create a bar plot
    plt.figure(figsize=(12, 8))
    plt.bar(codons, frequencies, color='skyblue', edgecolor='grey', alpha=0.7)  # bar plot
    plt.xlabel("Codon")  # x-label
    plt.ylabel("Frequency (%)")  # y-label
    plt.title(f"Codon Frequency Table for {accession_id}")  # figure title
    plt.xticks(rotation=90, fontsize=12)  # xticks fontsize
    plt.yticks(fontsize=12)  # yticks fontsize
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    # Adjust the figure layout
    plt.tight_layout()
    
    # Save the plot to a PDF file
    plt.savefig(f"{accession_id}.pdf", format="pdf", bbox_inches="tight")

    # Display the plot
    plt.show()

def main():
    # Create an argparse parser
    parser = argparse.ArgumentParser()
    # Add a positional argument for the GenBank file
    parser.add_argument("accession_id", help="NCBI accession id")
    # Parse the command-line arguments
    args = parser.parse_args()

    # retrieve genbank record
    record = retrieve_genomic_sequence(args.accession_id)
    # create annotation file and return cds sequences
    cds_sequences = create_gene_annotation_file(record)
    # Additional feature of codon frequency
    codon_frequencies = calculate_codon_frequency(cds_sequences, args.accession_id)
    # plot codon frequency
    plot_codon_frequency(codon_frequencies, args.accession_id)

if __name__ == "__main__":
    # executing the script
    main()

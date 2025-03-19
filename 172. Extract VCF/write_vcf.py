#!/usr/bin/env python
# coding: utf-8

# Required Libraries
import re  # handling regular expressions

# Read the list of genes from the file and store them in a set for fast lookup
with open("genes.txt", "r") as gene_reader:
    genes = set(gene_reader.read().splitlines())

# Open the output VCF file for writing
vcf_writer = open("SNPeff.selected.vcf", "w")

# Regular expression pattern to match gene identifiers
pattern = r'gene-[A-Za-z0-9_]+'

# Read the input VCF file
with open("SNPeff.all121.final.nonSyn.vcf", "r") as vcf_reader:
    for line in vcf_reader:
        if line.startswith("#"):
            # Write header lines directly to the output VCF file
            vcf_writer.write(line)
        else:
            # Extract the INFO field (8th column) from the VCF line
            info = line.split("\t")[7]
            # Search for the gene pattern in the INFO field
            match = re.search(pattern, info, re.I | re.S)
            if match:
                gene_id = match.group()

                # Write the line to the output VCF file if the gene is in the list
                if gene_id in genes:
                    vcf_writer.write(line)

# Close the output VCF file
vcf_writer.close()


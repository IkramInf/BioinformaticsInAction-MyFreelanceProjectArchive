#!/usr/bin/env python
# coding: utf-8

"""
The script read two files (Amigo_Mm_RNAproc.tsv & diffexp.csv) and create a output file
with three tab separated columns (gene name, gene descriptions, p-value)

make_testfiles: generate two test files for developing the code
read_go_annotations: read the go annotations file
read_diff_expressions: read the differential expressions file
write_to_file: write output into file

Last Date of Modification: 29 March, 2023
"""

# Function to make test files
def make_testfiles():
    """
    Parameters:
        No parameters
    Output:
        GO annotations test file
        Differential expressions test file
    """
    # make test file for go annotations
    with open("GOtestfile_studentid.tsv", "w") as ofile, open("Amigo_Mm_RNAproc.tsv", "r") as ifile:
        # take every tenth line to add into test file
        to_write = "".join(ifile.readlines()[0::10])
        ofile.write(to_write)
        
    # make test file for differential expressions
    with open("diffexptestfile_studentid.csv", "w") as ofile, open("diffexp.csv", "r") as ifile:
        # take every tenth line to add into test file
        to_write = "".join(ifile.readlines()[0::10])
        ofile.write(to_write)
        

# Function to read the go annotations file
def read_go_annotations(filename):
    """
    Parameters:
        filename: Input file location for go annotations
    Output:
        go_annot: A dictionary containing {gene_name: [(goterm, gene_desc),...], ...}
    """
    # GOterms in column five, gene names in column three and gene descriptions in column ten
    go_annot = {}
    with open(filename, "r") as f:
        # iterate over all lines
        for line in f.readlines():
            # split the line into columns (items)
            line = line.strip().split("\t")
            # python index start from 0
            go_annot.setdefault(line[2], []).append((line[4], line[9]))
            
    return go_annot
    

# Function to read the differential expressions file
def read_diff_expressions(filename): 
    """
    Parameters:
        filename: Input file location for differential expressions
    Output:
        diff_exp: A dictionary containing {gene_name: p_value, ...}
    """
    # gene names in column 1 and statistical significance reported as p-values in column 5
    diff_exp = {}
    with open(filename, "r") as f:
        # iterate over all lines
        for line in f.readlines():
            # split the line into columns (items)
            line = line.strip().split(",")
            # python index start from 0
            # str.strip('"') removes quotation mark from beginning and end
            diff_exp[line[0].strip('"')] = float(line[4])
    return diff_exp

# Function to write output into file
def write_to_file(data):
        """
        Parameters
            data: The data to write in output file
        Output
            Output file
        """
        with open("GOdiff_studentid.out", "w") as f:
            # iterate over all data
            for gene_name, desc_ps in data.items():
                for pair in desc_ps:
                    f.write(f"{gene_name}\t{pair[0]}\t{pair[1]}\n")


# calling the function make_testfiles
make_testfiles()
# calling the function read_go_annotations
go_annot = read_go_annotations("Amigo_Mm_RNAproc.tsv")
# calling the function read_diff_expressions
diff_exp = read_diff_expressions("diffexp.csv")

# create a dictionary for all genes that came up in the experiment as
# differentially expressed using a p-value cut-off of p<0.05
# and are associated with the GO-category GO:0003723
all_genes = {}
for gene_name, terms in go_annot.items():
    if gene_name in diff_exp:
        p_value = diff_exp[gene_name]
        for go_desc_pair in terms:
            go_category = go_desc_pair[0]
            # checking for conditions
            if p_value<0.5 and go_category=="GO:0003723":
                all_genes.setdefault(gene_name, []).append((go_desc_pair[1], p_value))

# sort same genes according to lowest p-value
for name in all_genes:
    all_genes[name] = sorted(all_genes[name], key=lambda x: x[1])
# calling the function write_to_file
write_to_file(all_genes)


"""
gffparser.py find list of genes within given QTL region coordinates.
It accepts five commandline options:
    “-i” => path to the GFF file
    “-c” => Chromosome name
    “-s” => Start coordinate of the region
    “-e” => End coordinate of the region
    “-o” => Where the output should be saved.
"""

import re  # handle regex
import argparse  # parse commandline arguments

def get_args():
    """
    Parse command-line arguments for the GFF parser.
    """
    description = 'Parse a GFF file to find list of genes within QTL region.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--input', required=True, help='Path to the GFF file')
    parser.add_argument('-c', '--chromosome', required=True, help='Chromosome name')
    parser.add_argument('-s', '--start', type=int, required=True, help='QTL start coordinate')
    parser.add_argument('-e', '--end', type=int, required=True, help='QTL end coordinate')
    parser.add_argument('-o', '--output', help='Output file path')
    args = parser.parse_args()
    return args

def get_genes_in_qtl_region(gff_file, chromosome, start, end):
    """
    Extract gene names within a specified genomic region from a GFF file.

    Args:
        gff_file (str): path to the GFF file
        chromosome (str): name of the chromosome to search
        start (int): QTL start coordinate
        end (int): QTL end coordinate

    Returns:
        list: list of gene names
    """
    genes_in_qtl = []  # store found genes

    # read gff file line by line
    with open(gff_file, "r", encoding="utf-8") as ifile:
        for annotation in ifile:
            if not annotation.startswith('#'):
                # separate tab delimited fields
                fields = annotation.strip().split('\t')
                # check for desired chromosome and feature
                if (fields[0] == chromosome) and (fields[2] == 'gene'):
                    gene_start, gene_end = int(fields[3]), int(fields[4])
                    # extract gene name
                    match = re.search(r"Name=([^;]+)", fields[-1], re.I|re.S)
                    if match:
                        gene_name = match.group(1)
                    # check coordinates to find out genes in qtl region
                    if (gene_start >= start) and (gene_end <= end):
                        # append gene names
                        genes_in_qtl.append(gene_name)

    # return the found genes
    return genes_in_qtl

def main():
    """
    Main function of the code.
    It parses command-line arguments, processes the GFF file,
    and outputs the list of genes within the QTL region.

    Usage:
        python gffparser.py -i TAIR10_GFF3_genes.gff -c Chr1 -s 1 -e 10000
    """
    # get the parsed cmd arguments
    args = get_args()
    # get genes in desired qtl region
    genes_in_qtl = get_genes_in_qtl_region(args.input, args.chromosome, args.start, args.end)
    genes_in_qtl = "\n".join(genes_in_qtl)

    # write output to file
    if args.output:
        with open(args.output, 'w', encoding="utf-8") as ofile:
            ofile.write(genes_in_qtl)
    else:
        # print on terminal if output path not given
        print(genes_in_qtl)

if __name__ == "__main__":
    # execute the script
    main()

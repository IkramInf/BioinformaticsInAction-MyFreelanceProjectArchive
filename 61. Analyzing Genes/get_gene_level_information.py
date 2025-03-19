import os
import sys
import re
from argparse import ArgumentParser
from assignment5 import config, io_utils

def update_host_name(host_name):
    """
    This function take the host name and checks for the available conversions from common to scientific names and return the scientific name.
    
    Parameters:
        host_name : host name
    Returns:
        scientific name of provided host name
    """
    host_keywords = config.get_keywords_for_hosts()
    if host_name in map(str.capitalize, host_keywords.keys()):
        return host_keywords[host_name]
    else:
        _print_directories_for_hosts()
    
    
def  _print_directories_for_hosts():
    """
    If the user asks for a directory that does not exist, the function prints out the host directories which do exist.
    
    Parameters:
        No parameters is needed
    Returns:
        None
    """
    print('\nEither the Host Name you are searching for is not in the database\n\nor If you are trying to use the scientific name please put the name in double quotes:\n\n"Scientific name"\n\nHere is a (non-case sensitive) list of available Hosts by scientific name\n\n')
    host_keywords = config.get_keywords_for_hosts()
    for i, v in enumerate(set(host_keywords.values())):
        print(f"{i+1}. {v}\n")
    print("\nHere is a (non-case sensitive) list of available Hosts by common name\n\n")
    for i, k in enumerate(sorted(host_keywords.keys(), key=lambda x:x.lower())):
        print(f"{i+1}. {k}\n")
    sys.exit(1)
    
def get_data_for_gene_file(file):
    """
    This function opens the file for the host and gene, extracts the list of tissues in which this gene is expressed and returns a sorted list of the tissues.
    
    Parameters:
        file : gene file name
    Returns:
        Sorted list of extracted tissues
    """
    ugene = io_utils.get_filehandle(file, "r")
    match = re.search("(?<=EXPRESS)(.*?)(?=\n)", ugene.read(), re.M | re.S)
    ugene.close()
    if match:
        tissue_strig = match.group(1)
        return sorted(list(map(str.strip, tissue_strig.strip().split("|"))))

def print_host_to_gene_name_output(host_name, gene_name, tissue_strig):
    """
    This function should print the tissue expression data for the gene
    
    Parameters:
        host_name : given host name
        gene_name : given gene name
        tissue_strig : List of tissues returned from get_data_for_gene_file
    Returns:
        None
    """
    print(f"In {host_name}, There are {len(tissue_strig)} tissues that {gene_name} is expressed in:\n\n")
    for i, tissue in enumerate(tissue_strig):
        print(f"{i+1}. {tissue.capitalize()}\n")
        
def main():
    # add command line options
    parser = ArgumentParser(description="Give the Host and Gene name")
    parser.add_argument('--host', help='Name of Host', default="Human")
    parser.add_argument('--gene', help='Name of Gene', default="TGM1")
    args = parser.parse_args()
    host = args.host.strip().replace("_", " ").capitalize()
    gene = args.gene.strip().upper()
    host = update_host_name(host)
    temp_host = host.replace("_", " ")
    
    file = os.path.join(config.get_directory_for_unigene(), host, gene + "." + config.get_extension_for_unigene())
    
    # check for the existence of file
    if io_utils.is_gene_file_valid(file):
        # using f-strings
        print(f"\nFound Gene {gene} for {temp_host}")
    else:
        print("Not found")
        print(f"Gene {gene} does not exist for {temp_host}. exiting now...", file=sys.stderr)
        sys.exit(1)
        
    tissue_strig = get_data_for_gene_file(file)
    print_host_to_gene_name_output(temp_host, args.gene, tissue_strig)

if __name__ == "__main__":
    main()
    
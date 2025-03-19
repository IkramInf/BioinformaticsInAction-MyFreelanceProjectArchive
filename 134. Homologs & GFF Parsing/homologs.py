"""
homologs.py finds out reciprocal best hit betwen two specieses
It accepts two comd arguments:
    --species1 : path to species1 blast xml output
    --species2 : path to species2 blast xml output
"""

import argparse  # handle cmd arguments
from Bio.Blast import NCBIXML  # parse blast xml output

def get_args():
    """
    Parse command-line arguments for the GFF parser.
    """
    parser = argparse.ArgumentParser(description="Identify reciprocal best hits")
    parser.add_argument("--species1", required=True, help="Blast XML output file for species 1")
    parser.add_argument("--species2", required=True, help="Blast XML output file for species 2")

    args = parser.parse_args()
    return args

def get_query_top_hit_pairs(xml_filepath):
    """
    Parses a BLAST XML file and extracts query-top hit pairs.

    Args:
        xml_filepath (str): path to the BLAST XML file

    Returns:
        dict: extracted query:top_hit pairs
    """
    # Parse the BLAST XML file
    result_handle = open(xml_filepath)
    blast_records = NCBIXML.parse(result_handle)

    query_tophit_pairs = {}
    # Iterate through the BLAST records and extract information
    for blast_record in blast_records:
        query = blast_record.query
        #query_gi = query.split("|")[1]
        if blast_record.alignments:
            top_hit = blast_record.alignments[0].title
            #top_hit_gi = top_hit.split("|")[1]
            query_tophit_pairs[query] = top_hit

    # Close the result handle
    result_handle.close()
    return query_tophit_pairs

def find_reciprocal_best_hits(species1_pairs, species2_pairs):
    """
    Finds reciprocal best hits

    Args:
        species1_pairs (dict): query-top hit pairs from species 1
        species2_pairs (dict): query-top hit pairs from species 2

    Returns:
        dict: reciprocal best hit pairs
    """
    reciprocal_pairs = {}  # store reciprocal pairs
    for query, tophit in species1_pairs.items():
        if tophit in species2_pairs:
            # check if both gene produces each other as top hit
            if query == species2_pairs[tophit]:
                # store into reciprocal_pairs
                reciprocal_pairs[query] = tophit
    return reciprocal_pairs

def write_to_file(reciprocal_pairs):
    """
    Writes reciprocal best hit pairs to a tab-separated text file
    """
    # open homologs.txt in write mode
    with open("homologs.txt", "w", encoding="utf-8")as ifile:
        for query, hit in reciprocal_pairs.items():
            # write each reciprocal pair into file
            ifile.write(f"{query}\t{hit}\n")

def main():
    """
    Main function to identify reciprocal best hits
    between two species and write them to a file
    
    Usage:
    python homologs.py --species1 drosophila_celegans.xml --species2 celegans_drosophila.xml
    """
    # get cmd arguments
    args = get_args()
    # get top hit for each query in species 1
    species1_pairs = get_query_top_hit_pairs(args.species1)
    # get top hit for each query in species 2
    species2_pairs = get_query_top_hit_pairs(args.species2)
    # find out reciprocal best blast pairs
    reciprocal_pairs = find_reciprocal_best_hits(species1_pairs, species2_pairs)
    # write results to a tab-separated file
    write_to_file(reciprocal_pairs)

if __name__ == "__main__":
    # execute the script
    main()
    
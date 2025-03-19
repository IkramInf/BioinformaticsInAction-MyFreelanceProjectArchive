#import required libraries
import requests
import argparse
from itertools import combinations

# adding command line arguments
parser = argparse.ArgumentParser(description='Extracting Phenotypes from BioLink API for Disease and Calculate Pairwise Jaccard Similarity for diseases...')
parser.add_argument('-i', '--pheno', required=True, help='An Input filename of disease.')
parser.add_argument('-o', '--output', default="similarity.tsv", help='An Output filename.')
args = parser.parse_args()

# read diseases into a list
with open(args.pheno, "r") as f:
    diseases = [disease.strip() for disease in f.readlines()]

# function to calculate jaccard similarity score
def jaccard_similarity(list1, list2):
    s1, s2 = set(list1), set(list2)
    return len(s1.intersection(s2)) / len(s1.union(s2))

# look into BioLink API and extract phenotypes
base_url = "https://api.monarchinitiative.org/api/bioentity"
phenotypes = {}
for disease in diseases:
    extension = f"/disease/{disease}/phenotypes"
    url = base_url + extension
    #print(url)
    response = requests.get(url).json()
    hp_nums = [block['object']['id'] for block in response['associations']]
    phenotypes[disease] = hp_nums

# write results into a tsv file
with open(args.output, "w") as f:
    f.write("disease_1\tdisease_2\tjaccard\n")
    for pair in combinations(phenotypes.keys(), 2):
        score = jaccard_similarity(phenotypes[pair[0]], phenotypes[pair[1]])
        f.write(f"{pair[0]}\t{pair[1]}\t{score}\n")

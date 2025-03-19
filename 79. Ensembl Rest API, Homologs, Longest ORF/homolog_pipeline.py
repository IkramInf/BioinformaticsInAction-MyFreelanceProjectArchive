# import required libraries
import sys
import re
import requests
from Bio.Seq import translate

# search ensembl id for gene 'MC1R'
gene = "MC1R"
server = "https://mygene.info/v3/query/"
endpoint = "?q={gene}&species=human&ensemblonly=1".format(gene=gene)
r = requests.get(server + endpoint).json()
entrez = r['hits'][0]['_id']

server = "http://mygene.info/v3/gene/"
endpoint = "{}?fields=ensembl".format(entrez)
r = requests.get(server + endpoint).json()
ensembl = r['ensembl']['gene']

# search dna sequence for the ensembl id
server = "https://rest.ensembl.org/"
endpoint = "sequence/id/{}?".format(ensembl)
r = requests.get(server+endpoint, headers={ "Content-Type" : "text/x-fasta"})

# write the dna sequence into a fasta file
texts = r.text.split("\n")
header = texts[0].strip()
seq = "".join(texts[1:]).replace("\n", "")

with open("mc1r.txt", "w") as f:
    f.write(f"{header}\n{seq}\n")


# find longest orf
orf = {}
N = len(seq)
for i in range(0, N, 3):
    if (i+3) <= N:
        if seq[i:i+3] == "ATG":
            for j in range(i, N, 3):
                if (j+3) <= N:
                    codon = seq[j:j+3]
                    if (codon == "TAA") or (codon == "TAG") or (codon == "TGA"):
                        orf.update({i : seq[i:j+3]})
                else:
                    orf.update({i : seq[i:j+3]})

longest_orf = sorted(orf.values(), key=len, reverse=True)[0]

# translate dna sequence into protein sequence
prot = translate(longest_orf)
with open("mc1r.txt", "a") as f:
    f.write(f"{header.split()[0]} protein sequence\n{prot}\n")


# find homolog species for the gene
server = "https://www.ncbi.nlm.nih.gov/homologene/"
endpoint = "?term={}".format(gene)
r = requests.get(server+endpoint)
texts = r.text

# write unique species into a text file
specieses = list(set(re.findall("<i>([A-Z]{1}\..+)</i>", texts, re.M)))

with open("mc1r_homology_list.txt", "w") as f:
    for species in specieses:
        f.write(f"{species}\n")




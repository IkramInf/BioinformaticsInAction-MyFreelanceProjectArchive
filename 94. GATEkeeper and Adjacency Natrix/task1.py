#!/usr/bin/env python
# coding: utf-8

# required python libraries
import os
import sys
import shutil
import subprocess
from Bio import SeqIO
import pandas as pd

# required input variables
filename = "Selected_Unique_COVID19_Genomes_Asia.fasta"
url = "https://github.com/chenyongrowan/GATEkeeper.git"
N = 100

# read all sequences into a dictionary as {id : sequence}
sequences = {record.id : str(record.seq) for record in SeqIO.parse(filename, "fasta")}


# if GATEkeeper isn't installed, install it
if not os.path.exists("GATEkeeper/"):
    #shutil.rmtree("GATEkeeper")
    p1 = subprocess.run(f"git clone {url} && cd GATEkeeper && make", shell=True)


# write 100 genomes into 100 separate fasta files in Output directory including NC_045512 genome
nc = [Id for Id in sequences.keys() if Id.startswith("NC_045512")]

if not os.path.exists("Output/"):
    os.mkdir("Output")

genomes = []
nc_id, nc_seq = nc[0], sequences[nc[0]]
with open(f"Output/{nc_id}.fasta", "w") as f:
    f.write(f">{nc_id}\n{nc_seq}")
genomes.append(f"Output/{nc_id}.fasta")

count = 1
for Id, sequence in sequences.items():
    if count == N:
        break
    if Id.startswith("NC_045512"):
        continue

    with open(f"Output/{Id}.fasta", "w") as f:
        f.write(f">{Id}\n{sequence}")
    genomes.append(f"Output/{Id}.fasta")
    count += 1


# calculate mutations for each pair and create adjacency matrix
adjacent_arr = [[0 for c in range(N)] for r in range(N)]

for r, genome1 in enumerate(genomes):
    for c, genome2 in enumerate(genomes):
        print(f"===== Running for ({genome1}, {genome2}) pair =====")
        if (adjacent_arr[r][c] == 0) and (genome1 != genome2):
            
            p2 = subprocess.run(f"GATEkeeper/bin/GATEkeeper -r {genome1} -q {genome2} -o Output/output", shell=True)
            #print(p2)
            if p2.returncode == 0:
                mutation_num, insert, subs = 0, 0, 0
                with open("Output/output.vcf", "r") as ifile:
                    mutation_num, insert, delete = 0, 0, 0
                    for line in ifile.readlines():
                        if not line.strip().startswith("#") and line:
                            mutation_num += 1

                            if "INSERT" in line.strip().split("\t")[-1].upper():
                                insert += 1

                            if "DELETE" in line.strip().split("\t")[-1].upper():
                                delete += 1

                subs = mutation_num - (insert + delete)
                #print(f"\nGATEkeeper identifies {subs} SNVs, {insert} insertions, and {delete} deletions [Output/output.vcf].\n")

                adjacent_arr[r][c] = mutation_num
                adjacent_arr[c][r] = mutation_num
                
            else:
                sys.exit(f"Error occurred for pair ({genome1}, {genome2})")
                
        else:
            print("Already Done or Both are Same Files!!!\n")
            
shutil.rmtree("Output")


# write adjacency matrix in a csv file
genomes = [name.split("/")[-1].split(".fa")[0] for name in genomes]
df = pd.DataFrame(adjacent_arr, columns=genomes, index = genomes)
df.to_csv("Adjacent Matrix.csv", index=True)

print("The program is successfully executed...")

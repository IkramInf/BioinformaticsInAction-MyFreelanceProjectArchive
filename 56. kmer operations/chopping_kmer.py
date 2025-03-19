import os
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", required=True, help="Name of reference sequence file")
parser.add_argument("-d", "--db", required=True, help="Name of database sequence file")
parser.add_argument("-k", "--length", type=int, default=15, help="Length of kmer")
args = parser.parse_args()
k = args.length

ref = str(SeqIO.read(args.ref, "fasta").seq)
kmers = [ref[i:i+k] for i in range(len(ref)-k+1)]
filename = os.path.splitext(args.ref)[0] + ".csv"

with open(filename, "w") as f:
    f.write("Index, kmer\n")
    for i, kmer in enumerate(kmers):    
        f.write(f"{i+1},{kmer}\n")

database = [str(record.seq) for record in SeqIO.parse(args.db, "fasta")]

counts = {}
for kmer in kmers:
    for seq in database:
        if kmer in seq:
            counts[kmer] = counts.get(kmer, 0) + 1

with open("counts_kmer.csv", "w") as w:
    w.write("kmer, counts\n")
    for k, v in counts.items():
        w.write(f"{k}, {v}\n")
        
print(f"Successfully extracted kmers. Check the following folder to see output: {os.getcwd()}")        
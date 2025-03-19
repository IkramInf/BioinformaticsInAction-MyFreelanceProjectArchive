# import required modules
from Bio import SeqIO
import os
import sys

alignment_sequence, threshold = sys.argv[1:]

try:
    # read all alignment sequences
    if os.path.isfile(alignment_sequence):
        records = [str(record.seq) for record in SeqIO.parse(alignment_sequence, "fasta")]
    else:
        print("The file is no longer exist! Enter a valid filename")
        sys.exit()
    
except Exception as e:
    sys.exit(e)
    
try:
    threshold = float(threshold)
except ValueError:
    sys.exit("Please enter a valid numeric value.")

# to store count of each position
a, c, g, t = [], [], [], []

# iterate over each position of all sequences
for position in zip(*records):
    a.append(position.count("A"))
    c.append(position.count("C"))
    g.append(position.count("G"))
    t.append(position.count("T"))

# find out rate of base for each position
# create a dictionary as like, position : {'A':rateOfA, 'C':rateOfC, 'G':rateOfG, 'T':rateOfT}
base_rate = {index : {"A":val[0]/sum(val), "C":val[1]/sum(val), "G":val[2]/sum(val), "T":val[3]/sum(val)}  for index, val in enumerate(zip(a, c, g, t))}


# iterate over base_rate dictionary
for index in base_rate:
    # sort the inner dictionary as descending order
    bases = sorted(base_rate[index].items(), key=lambda x: x[1], reverse = True)
    
    # take maximum base with rate
    max_base = bases[0]
    
    # take minimum base with rate
    min_base = bases[1]
    
    # to check if minimum base rate value is greater than threshold value
    if min_base[1] > threshold:
        print(f"Polymorphism at {index}, {max_base[0]} rate is : {max_base[1]}, {min_base[0]} rate is : {min_base[1]}")

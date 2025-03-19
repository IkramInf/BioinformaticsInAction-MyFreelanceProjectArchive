#!/usr/bin/env python
# coding: utf-8
# import required libraries
from math import log
from collections import defaultdict
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

## How to install dssp
# [1] install libcifpp firstly
#git clone https://github.com/PDB-REDO/libcifpp.git && cd libcifpp && mkdir build && cd build && cmake .. && cmake --build . --config Release && ctest -C Release && sudo cmake --install .
# [2] Then, install dssp
#git clone https://github.com/PDB-REDO/dssp.git && cd dssp && mkdir build && cd build && cmake .. && cmake --build . --config Release && ctest -C Release && sudo cmake --install .

# read dssp file and save informations in variables
with open("1shr.dssp", "r") as f:
    lines = f.readlines()
    index = [i for i, line in enumerate(lines) if line.strip().startswith("#")][0]
    # aa_types contains distribution frequencies of amino acid types
    # ss_types contains distribution frequencies of secondary acid types
    # probs contains conditional probability of P(aa | ss)
    aa_types, ss_types, probs = {}, {}, defaultdict(dict) 
    for line in lines[index+1:]:
        temp = line.strip().split()
        if temp[3].isalpha():
            aa_types[temp[3]] = aa_types.get(temp[3], 0) + 1
        if temp[4].isalpha():
            ss_types[temp[4]] = ss_types.get(temp[4], 0) + 1
            
        if temp[3].isalpha() and temp[4].isalpha():
            try:
                probs[temp[4]][temp[3]] += 1
            except KeyError:
                probs[temp[4]][temp[3]] = 0
                probs[temp[4]][temp[3]] += 1


# create pie chart of distribution plot for amino acid types
data, labels = [], []
for k, v in aa_types.items():
    data.append(v)
    labels.append(k)
#create pie chart
plt.figure(figsize=(12, 8))
plt.title("Distribution plot for amino acid types")
plt.pie(data, labels=labels, colors=sns.color_palette('pastel'), autopct='%.0f%%', textprops={'color':'k', 'fontsize': 12})
plt.show()

# create pie chart of distribution plot for secondary structure types
data, labels = [], []
for k, v in ss_types.items():
    data.append(v)
    labels.append(k)
#create pie chart
plt.figure(figsize=(12, 8))
plt.title("Distribution plot for secondary structure types")
plt.pie(data, labels=labels, colors=sns.color_palette('pastel'), autopct='%.0f%%', textprops={'color':'k', 'fontsize': 12})
plt.show()

# calculate propensity log ratio
# log( P(aa | ss) / P(aa) )
propensities = defaultdict(dict)
for k1, v1 in probs.items():
    for k2, v2 in v1.items():
        propensities[k1][k2] = abs(log(probs[k1][k2] / aa_types[k2]))

# visualize propensity value as bar plot
for k, v in propensities.items():
    y = list(v.keys())
    width = list(v.values())
    # creating the bar plot
    plt.title(f"Bar plot for {k}")
    plt.barh(y, width, height=0.8, color ='maroon')
    plt.show()

print("Program is Executed Successfully!")
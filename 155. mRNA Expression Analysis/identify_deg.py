# Required libraries
import os  # interacting with os
import sys  # system-related functionalities
import numpy as np  # arithmetic operations
import pandas as pd  # handling dataframe

# Read the mRNA file
mRNA_df = pd.read_csv("data_mrna_seq_v2_rsem.txt", sep="\t")

# Handling Sample IDs with more than just the “-01” sample
id_dict = {}  # store patiendID as key and sampleID as value
for sampleID in mRNA_df.columns:
    if sampleID.endswith("01"):
        patientID = "-".join(sampleID.split("-")[:-1])
        if f"{patientID}-02" in mRNA_df.columns:
            columns = [f"{patientID}-01", f"{patientID}-02"]
            mRNA_df[f"{patientID}-avg"] = mRNA_df[columns].mean(axis=1)
            mRNA_df.drop(columns, axis=1, inplace=True)
            id_dict[patientID] = sampleID
        else:
            id_dict[patientID] = sampleID

# Ask user to select p1 and p2 dataset
p1_dataset = input("Which p1 dataset would you like to select? ")
p1_dataset = f"{p1_dataset}.txt"
p2_dataset = input("Which p2 dataset would you like to select? ")
p2_dataset = f"{p2_dataset}.txt"

# Check if p1 and p2 dataset path exist
if not os.path.exists(p1_dataset):
    print(f"The specified path '{p1_dataset}' does not exist. Exiting...")
    sys.exit()
if not os.path.exists(p2_dataset):
    print(f"The specified path '{p2_dataset}' does not exist. Exiting...")
    sys.exit()

# Read the p1 and p2 dataset
p1_df = pd.read_csv(p1_dataset, skiprows=1)
p2_df = pd.read_csv(p2_dataset, skiprows=1)

## p1math_df ##
# Store mean and sd into a dict for each gene
# then convert the dict into dataframe
p1math_dict = {'sampleID':[], 'mean':[], 'sd':[]}
for Id in p1_df["Patient Identifier"].values:
    try:
        sampleID = id_dict[Id]
        expression_value = [np.log2(x) for x in mRNA_df[sampleID].values if x > 0]
        p1math_dict['sampleID'].append(sampleID)
        p1math_dict['mean'].append(np.mean(expression_value))
        p1math_dict['sd'].append(np.std(expression_value))
    except KeyError:
        continue

# convert dict into dataframe
p1math_df = pd.DataFrame(p1math_dict)
p1math_df.set_index('sampleID', inplace=True)
print("\np1math_df:")
print(p1math_df)

## p2math_df ##
# Store mean and sd into a dict for each gene
# then convert the dict into dataframe
p2math_dict = {'sampleID':[], 'mean':[], 'sd':[]}
for Id in p2_df["Patient Identifier"].values:
    try:
        sampleID = id_dict[Id]
        expression_value = [np.log2(x) for x in mRNA_df[sampleID].values if x > 0]
        p2math_dict['sampleID'].append(sampleID)
        p2math_dict['mean'].append(np.mean(expression_value))
        p2math_dict['sd'].append(np.std(expression_value))
    except KeyError:
        continue
        
# convert dict into dataframe
p2math_df = pd.DataFrame(p2math_dict)
p2math_df.set_index('sampleID', inplace=True)
print("\np2math_df:")
print(p2math_df)

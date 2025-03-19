#!/usr/bin/env python
# coding: utf-8

# importing required libraries and modules
import csv  # read/write csv file
import pandas as pd  # pip install pandas

# sra_table dictionary will store desired data
sra_table = {"SRA_Accession":[], "GEO_Accession":[], "Age":[], "Cell_type":[]}
# get SraRunTable.txt following the below url
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835
with open("SraRunTable.txt", "r") as f:
    # read the first line (headers)
    headers = f.readline().strip().split(",")
    # read the data with csv.reader(filename)
    for row in csv.reader(f):
        # add SRA, GEO accession, age and cell type into dict
        sra_table["SRA_Accession"].append(row[0])
        sra_table["GEO_Accession"].append(row[17])
        sra_table["Age"].append(row[1])
        sra_table["Cell_type"].append(row[9])
        
# make the dataframe from the dict
df = pd.DataFrame.from_dict(sra_table)
# convert dataframe into csv file
df.to_csv("output.csv", index=False)

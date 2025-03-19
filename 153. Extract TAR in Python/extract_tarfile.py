#!/usr/bin/env python
# coding: utf-8

# required libraries
import os
import wget
import tarfile
import pandas as pd

# Download data from cBioPortal using wget module
url = "https://cbioportal-datahub.s3.amazonaws.com/gbm_tcga.tar.gz"
filename = os.path.basename(url)
if not os.path.exists(filename):
    wget.download(url, out=filename)
else:
    print(f"{filename} already exists. Skipping download...")

# Extract the specific files from the tar.gz
specific_files = ["data_mrna_seq_v2_rsem.txt", "data_clinical_patient.txt"]
extracted_filepaths = []  # store path of specific file

# open tar.gz file as read mode
with tarfile.open(filename, "r:gz") as reader_obj:
    # check in specific files to extract them only
    for file in specific_files:
        for member in reader_obj.getmembers():
            if member.name.endswith(file):
                # extract the file
                reader_obj.extract(member)
                # add the file path into extracted_filepaths
                extracted_filepaths.append(member.name)


# Read the extracted files into pandas dataframes
mRNA_df = pd.read_csv(extracted_filepaths[0], sep="\t")
print("Dataframe for data_mrna_seq_v2_rsem.txt:")
print(mRNA_df.head())
clinical_df = pd.read_csv(extracted_filepaths[1], sep="\t", skiprows=range(1, 5), header=0)
print("\nDataframe for data_clinical_patient.txt:")
print(clinical_df.head())

# extract specific columns from data_mrna_seq_v2_rsem.txt
specific_columns = ['#Other Patient ID', 'Patient Identifier', 'Overall Survival Status',
                    'Overall Survival (Months)', 'Disease Free Status', 'Disease Free (Months)']
clinical_df = clinical_df[specific_columns]
print("\nDataframe for data_clinical_patient.txt with specific columns:")
print(clinical_df.head())
# save the new data as csv file
clinical_df.to_csv("new_data_clinical_patient.txt", index=False)

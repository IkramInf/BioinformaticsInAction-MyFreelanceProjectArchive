#!/usr/bin/env python
# coding: utf-8

# import required libraries
import re
import pandas as pd

# File location of List1 and List2
List1_filename = "SSCDR donor HLA type List 1.xlsx"
List2_filename = "Hypothetical donors-1-1-1.xlsx"

def get_allel(allel):
    """
    Extracts the HLA allele portion from a given allele string.

    Parameters:
    allel (str): The input allele string.

    Returns:
    str: The extracted HLA allele portion.

    """
    # Extract the first two components of the allele string
    allel = ":".join(allel.split(":")[:2])
    # Remove any non-alphanumeric characters from the allele string
    allel = re.sub("[a-zA-Z\s]", "", allel)
    return allel


def ishomo(alleles):
    """
    Counts the number of homozygous alleles in a given set of alleles.

    Parameters:
    alleles (dict): A dictionary containing the alleles for HLA-A1, HLA-A2, HLA-B1, HLA-B2, DRB1-1, and DRB1-2.

    Returns:
    int: The count of homozygous alleles in the given set.
    """
    count = 0
    if alleles['HLA-A1'] == alleles['HLA-A2']:
        count += 1
    if alleles['HLA-B1'] == alleles['HLA-B2']:
        count += 1
    if alleles['DRB1-1'] == alleles['DRB1-2']:
        count += 1
    return count

# read and filter List1 donor data
List1_df = pd.read_excel(List1_filename)
List1_df = List1_df.astype(str)
for column in List1_df.columns:
    List1_df[column] = List1_df[column].apply(get_allel)
List1_df['homozygous'] = List1_df.apply(lambda row: ishomo(row), axis=1)

# read List2 hypothetical donor data
List2_df = pd.read_excel(List2_filename)
List2_df['Haplotype'] = List2_df.apply(lambda row: list(map(str.strip, row)), axis=1)

# count haplotype matching and number of homozygous donors
count_lists = []  # saves like: [haplotype, number of matching donors, number of homozygous donors]
list2_columns = ['HLA-A1', 'HLA-A2', 'HLA-B1', 'HLA-B2', 'DRB1-1', 'DRB1-2']
for haplotype in List2_df['Haplotype']:
    dfs = []
    for column in list2_columns:
        dfs.append(List1_df[column].isin(haplotype))
    dfs = pd.concat(dfs, axis=1)
    count_hap = dfs.any(axis=1).sum()
    dfs = pd.concat([dfs.any(axis=1), List1_df['homozygous']], axis=1)
    dfs.columns = ['isdonor', 'homozygous']
    count_homo = dfs[(dfs['isdonor']==True) & (dfs['homozygous']>1)].shape[0]
    count_lists.append(["~".join(haplotype), count_hap, count_homo])

# make dataframe with the count_lists
columns = ['Haplotype', 'Number of Matching Donors', 'Number of Homozygous Donors']
df = pd.DataFrame(count_lists, columns=columns)

#sort values by the column 'Number of Matching Donors'
df = df.sort_values('Number of Matching Donors', ascending=False)

# find out number of samples
number_of_samples = df['Number of Matching Donors'].sum()
#number_of_samples = df.shape[0]

# calculate percentage matching donors
df['Percentage of Matching Donors'] = df['Number of Matching Donors'].apply(lambda n: round((n/number_of_samples)*100, 4))
# calculate cumulative sum
df['Cumulative Percentage'] = df['Percentage of Matching Donors'].cumsum()
df['Cumulative Percentage'] = df['Cumulative Percentage'].round(4)

# rearrange columns
homo_df = df['Number of Homozygous Donors']
df = pd.concat([df.drop('Number of Homozygous Donors', axis=1), homo_df], axis=1)

# save dataframe into a excel file
df.to_excel("output.xlsx", index=False)

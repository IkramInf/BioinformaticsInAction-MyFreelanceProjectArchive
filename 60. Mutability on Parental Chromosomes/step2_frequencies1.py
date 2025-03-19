
import os
import sys
import io
import gzip
import glob
import time
import pickle
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import LabelEncoder


def vcf_parser(filename):
    extension = os.path.splitext(filename)[1]
    if extension == '.gz':
        with gzip.open(filename, "rt") as f:
            lines = [l for l in f if not l.startswith('##')]
        cols = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        return pd.read_csv(io.StringIO(''.join(lines)), sep='\t', dtype={'#CHROM':float, 'POS':int}, usecols=lambda x: x not in cols).rename(columns={'#CHROM':'Chromosome', 'POS':'SNP', 'REF':'ReferenceVariant', 'ALT':'AlternativeVariant'})
    else:
        with open(filename, "r") as f:
            lines = [l for l in f if not l.startswith('##')]
        cols = ['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        return pd.read_csv(io.StringIO(''.join(lines)), sep='\t', dtype={'#CHROM':float, 'POS':int}, usecols=lambda x: x not in cols).rename(columns={'#CHROM':'Chromosome', 'POS':'SNP', 'REF':'ReferenceVariant', 'ALT':'AlternativeVariant'})

#VARIABLES USED THAT YOU CAN CHANGE/MODIFY

windows_range = 500000

#CODE

# In this part, we get all the names of the files we need to work. This part is meant to be modified according to your directories path

sample_location = "/TEST/paca/*/*/*/*/*/*/*/*/*.vcf.gz"
file_names = glob.glob(os.getcwd() + sample_location, recursive=True)

samples_files_names = []
for i, filename in enumerate(file_names):
    ch = filename.lower().split("_output_")[1].split(".txt")[0]
    code = filename.lower().split("/")[-1].split("_output_")[0]
    ID = code.split("-")[0]
    flag = ch.split("chr")[1]
    if flag.isdigit() and (flag != str(23)) and (ID == "3c86ba21"):
        samples_files_names.append(filename)

print("Total number of files in directory are : ", len(samples_files_names))
print(samples_files_names)

file_with_frequencies = "results_of_the_frequencies.txt"
table = pd.read_csv(file_with_frequencies, delimiter="\t")
table = table.dropna()
table['Chromosome'] = table['Chromosome'].astype(float)
print("All the files will be compared with the file : ", str(file_with_frequencies))
print(table.shape)
print(table.dtypes)
print(table.head())

if not os.path.exists("step2_images"):
    os.mkdir("step2_images")

# In this part, we take each chromosome file of the patient, and we generate the graphs with the list of all the SNPs we are evaluating

for l, file in enumerate(samples_files_names):
    print(f"===== {file} =====")
    df = vcf_parser(file)
    print(df.dtypes)
    df = df.merge(table, how='inner', on=['Chromosome', 'SNP', 'ReferenceVariant', 'AlternativeVariant'])
    print(df.shape)
    print(df.head())
    
    percentage_EUR = 0
    percentage_EAS = 0
    percentage_AFR = 0
    percentage_SAS = 0
    percentage_AMR = 0

    values_EUR = []
    values_EAS = []
    values_AFR = []
    values_SAS = []
    values_AMR = []
    
    positions = []
    counter_of_snps = df.shape[0]

    position_min = 0
    position_max = 0
    counter = 0

    for i in range(counter_of_snps):
        df1 = df.iloc[i]
        chromosome_sample, pos, genotype = df1['Chromosome'], df1['SNP'], df1['SAMP001']

        if position_min == 0:
            position_min = df1['SNP']
            position_max = df1['SNP'] + windows_range       

        if genotype == '0|1' or genotype == '1|0' or genotype == '1|1':
            counter += 1
            if genotype == '0|1' or genotype == '1|0':
                percentage_EUR += df1['FreqEURalt']
                percentage_EAS += df1['FreqEASalt']
                percentage_AFR += df1['FreqAFRalt']
                percentage_SAS += df1['FreqSASalt']
                percentage_AMR += df1['FreqAMRalt']

            elif genotype == '1|1':
                percentage_EUR += df1['FreqEURalt']
                percentage_EAS += df1['FreqEASalt']
                percentage_AFR += df1['FreqAFRalt']
                percentage_SAS += df1['FreqSASalt']
                percentage_AMR += df1['FreqAMRalt']

        if pos > position_max:
            values_EUR.append(percentage_EUR / counter)
            values_EAS.append(percentage_EAS / counter)
            values_AFR.append(percentage_AFR / counter)
            values_SAS.append(percentage_SAS / counter)
            values_AMR.append(percentage_AMR / counter)
            positions.append(str(position_min) + '-' + str(position_max))
            
                    
    #Graph for EUR
    plt.figure(figsize=(18,12)) 
    graph1_title = 'CHROMOSOME: ' + str(chromosome_sample) + '---> % EUR'
    graph1_x_tags = positions
    graph1_x_values = values_EUR
    plt.title(graph1_title)
    plt.xlabel('Position')
    plt.ylabel('SNP frequency')
    plt.bar(graph1_x_tags, graph1_x_values)
    plt.savefig(f'step2_images/graph_for_eur_{l+1}.png', bbox_inches='tight')
    plt.show()

    #Graph for EAS
    plt.figure(figsize=(18,12)) 
    graph1_title = 'CHROMOSOME: ' + str(chromosome_sample) + '---> % EAS'
    graph1_x_tags = positions
    graph1_x_values = values_EAS
    plt.title(graph1_title)
    plt.xlabel('Position')
    plt.ylabel('SNP frequency')
    plt.bar(graph1_x_tags, graph1_x_values)
    plt.savefig(f'step2_images/graph_for_eas_{l+1}.png', bbox_inches='tight')
    plt.show()

    #Graph for AFR
    plt.figure(figsize=(18,12)) 
    graph1_title = 'CHROMOSOME: ' + str(chromosome_sample) + '---> % AFR'
    graph1_x_tags = positions
    graph1_x_values = values_AFR
    plt.title(graph1_title)
    plt.xlabel('Position')
    plt.ylabel('SNP frequency')
    plt.bar(graph1_x_tags, graph1_x_values)
    plt.savefig(f'step2_images/graph_for_afr_{l+1}.png', bbox_inches='tight')
    plt.show()

    #Graph for SAS
    plt.figure(figsize=(18,12)) 
    graph1_title = 'CHROMOSOME: ' + str(chromosome_sample) + '---> % SAS'
    graph1_x_tags = positions
    graph1_x_values = values_SAS
    plt.title(graph1_title)
    plt.xlabel('Position')
    plt.ylabel('SNP frequency')
    plt.bar(graph1_x_tags, graph1_x_values)
    plt.savefig(f'step2_images/graph_for_sas_{l+1}.png', bbox_inches='tight')
    plt.show()

    #Graph for AMR
    plt.figure(figsize=(18,12)) 
    graph1_title = 'CHROMOSOME: ' + str(chromosome_sample) + '---> % AMR'
    graph1_x_tags = positions
    graph1_x_values = values_AMR
    plt.title(graph1_title)
    plt.xlabel('Position')
    plt.ylabel('SNP frequency')
    plt.bar(graph1_x_tags, graph1_x_values)
    plt.savefig(f'step2_images/graph_for_amr_{l+1}.png', bbox_inches='tight')
    plt.show()

    #Graph for all ancestors

    plt.figure(figsize=(18,12)) 
    graph1_title = 'CHROMOSOME: ' + str(chromosome_sample) + '---> % EUR/EAS/AFR/SAS/AMR'
    graph1_x_tags = positions
    eur = np.array(values_EUR)
    eas = np.array(values_EAS)
    afr = np.array(values_AFR)
    sas = np.array(values_SAS)
    amr = np.array(values_AMR)

    b_afr = eur + eas
    b_sas = eur + eas + afr
    b_amr = eur + eas + afr + sas

    plt.title(graph1_title)
    plt.xlabel('Position')
    plt.ylabel('SNP frequency')
    plt.bar(graph1_x_tags, eur, 0.4, label = "EUR")
    plt.bar(graph1_x_tags, eas, 0.4, bottom=eur, label = "EUR")
    plt.bar(graph1_x_tags, afr, 0.4, bottom=b_afr, label = "EUR")
    plt.bar(graph1_x_tags, sas, 0.4, bottom=b_sas, label = "EUR")
    plt.bar(graph1_x_tags, amr, 0.4, bottom=b_amr, label = "EUR")
    plt.savefig(f'step2_images/graph_for_all_ancestors_{l+1}.png', bbox_inches='tight')
    plt.show()


print("Take a look at the graphs now!!!")
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
from sklearn.decomposition import PCA

def vcf_parser(filename):
    extension = os.path.splitext(filename)[1]
    if extension == '.gz':
        with gzip.open(filename, "rt") as f:
            lines = [l for l in f if not l.startswith('##')]
        cols = ['ID','QUAL','FILTER','INFO','FORMAT']
        return pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int}, sep='\t', usecols=lambda x: x not in cols).rename(columns={'#CHROM': 'CHROM'})
    else:    
        with open(filename, "r") as f:
            lines = [l for l in f if not l.startswith('##')]
        cols = ['ID','QUAL','FILTER','INFO','FORMAT']
        return pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int}, sep='\t', usecols=lambda x: x not in cols).rename(columns={'#CHROM': 'CHROM'})
    

#CODE

#ANALYZE ALL FILES

if not os.path.exists("coincidences"):
    os.mkdir("coincidences")

patient = "3c86ba21"

chromosome_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

location_indels = '/mnt/iribhm/people/sagnetti/SNV/indel/3c86ba21-7b11-4ec7-9d20-a2325197c676.consensus.20161006.somatic.indel.vcf.gz'
location_snvs = '/mnt/iribhm/people/sagnetti/SNV/snv_mnv/3c86ba21-7b11-4ec7-9d20-a2325197c676.consensus.20160830.somatic.snv_mnv.vcf.gz'

file_names = glob.glob(os.getcwd()+"/TEST/paca/*/*/*/*/*/*/*/*/*.vcf.gz", recursive=True)
#Part to locate samples' data
location_samples = []
for i, filename in enumerate(file_names):
    ch = filename.lower().split("_output_")[1].split(".txt")[0]
    code = filename.lower().split("/")[-1].split("_output_")[0]
    ID = code.split("-")[0]
    flag = ch.split("chr")[1]
    if flag.isdigit() and flag != str(23):
        if ID == patient:
            location_samples.append(filename)


#Part to check 1vs1 each position of each chromosome

###### Comparision with indels #####
counter = 0
coincidences_indels = {}
#coincidences_txt = root_dir + "/coincidences/coincidences.txt"
loc_indels = vcf_parser(location_indels)
#print(loc_indels.head())

for i, file in enumerate(location_samples):
    if i == 1:
        break
        
    loc_sample = vcf_parser(file)
    #print(loc_sample.head())
    print(loc_sample.shape)
    print(loc_indels.shape)
    cols = list(loc_sample.columns)
    df = loc_sample.merge(loc_indels, how='inner', on=['CHROM', 'POS'], suffixes=('_1', '_2'))
    print(df.shape)
    print(df.head())
    chromosome = list(set(df['CHROM']))
    print(chromosome)
    pos_samp = pd.DataFrame(df[['POS', cols[-1]]].astype(str).drop_duplicates(keep='first').values, columns=['POS', 'SAMP'])
    print(pos_samp.shape)
    #pos_samp.columns = ['POS', 'SAMP']
    counter += len(pos_samp)
    coincidences_indels[chromosome[0]] = pos_samp


with open(f"coincidences/coincidences_indels_chr{chromosome[0]}.txt", 'w') as f:
    for k, v in coincidences_indels.items():
        v1, v2 = " ".join(v['POS']), " ".join(v['SAMP'])
        f.write(f"{k} : {v1}\n{v2}\n\n")
        
        
for k, v in coincidences_indels.items():
    val = list(v['SAMP'])
    ancestor1, ancestor2 = 0, 0
    for c in val:
        c1, c2 = c.split("|")
        if int(c1) == 1:
            ancestor1 += 1
        elif int(c2) == 1:
            ancestor2 += 1
            
    print(f"#Indels : Percentage of Ancestor 1 in Chr{k} is: {(ancestor1/(ancestor1 + ancestor2))*100}")
    print(f"#Indels : Percentage of Ancestor 2 in Chr{k} is: {(ancestor2/(ancestor1 + ancestor2))*100}")


###### Comparision with snvs #####
counter = 0
coincidences_snvs = {}    
loc_snvs = vcf_parser(location_snvs)
#print(loc_snvs.head())

for i, file in enumerate(location_samples):
    if i == 1:
        break
        
    loc_sample = vcf_parser(file)
    print(loc_sample.shape)
    print(loc_snvs.shape)
    cols = list(loc_sample.columns)
    df = loc_sample.merge(loc_snvs, how='inner', on=['CHROM', 'POS'], suffixes=('_1', '_2'))
    print(df.shape)
    print(df.head())
    chromosome = list(set(df['CHROM']))
    print(chromosome)
    pos_samp = pd.DataFrame(df[['POS', cols[-1]]].astype(str).drop_duplicates(keep='first').values, columns=['POS', 'SAMP'])
    print(pos_samp.shape)
    #pos_samp.columns = ['POS', 'SAMP']
    counter += len(pos_samp)
    coincidences_snvs[chromosome[0]] = pos_samp


with open(f"coincidences/coincidences_snvs_chr{chromosome[0]}.txt", 'w') as f:
    for k, v in coincidences_snvs.items():
        v1, v2 = " ".join(v['POS']), " ".join(v['SAMP'])
        f.write(f"{k} : {v1}\n{v2}\n\n")
        
        
for k, v in coincidences_snvs.items():
    val = list(v['SAMP'])
    ancestor3, ancestor4 = 0, 0
    for c in val:
        c1, c2 = c.split("|")
        if int(c1) == 1:
            ancestor3 += 1
        elif int(c2) == 1:
            ancestor4 += 1
            
    print(f"#Snvs : Percentage of Ancestor 1 in Chr{k} is: {(ancestor3/(ancestor3+ancestor4))*100}")
    print(f"#Snvs : Percentage of Ancestor 2 in Chr{k} is: {(ancestor4/(ancestor3+ancestor4))*100}")
    

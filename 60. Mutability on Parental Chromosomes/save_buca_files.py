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
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import plotly
import plotly.express as px
import plotly.graph_objects as go
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



print("===== Extracting filenames from BUCA directory =====")
file_names = glob.glob(os.getcwd()+"/*/*/*/*/*/*/*/*/*.vcf.gz", recursive=True)
print("Total number of files in BUCA directory are : ", len(file_names))


if not os.path.exists("buca"):
    os.mkdir("buca")
        
fnames = {}
for i, filename in enumerate(file_names):
    ch = filename.lower().split("_output_")[1].split(".txt")[0]
    code = filename.lower().split("/")[-1].split("_output_")[0]
    flag = ch.split("chr")[1]
    if flag.isdigit() and flag != str(23):
        #fname = "buca/" + ch + "_" + str(i) + ".csv"
        fname = "buca/" + code + "_" + ch + ".csv"
        df = vcf_parser(filename)
        df.to_csv(fname, index=False)
        fnames.setdefault(code, []).append(fname)        
        
    
with open('buca/fnames.pkl', 'wb') as dbfile:
    pickle.dump(fnames, dbfile, pickle.HIGHEST_PROTOCOL)

with open('buca/fnames.pkl', 'rb') as f:
    print(pickle.load(f))    
    
import os
import sys
import io
import re
import gzip
import glob
import time
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


# read the vcf file
df1 = pd.read_csv("chr22.csv")
print(df1.shape)
print(df1.head())

def vcf_parser(filename):
    extension = os.path.splitext(filename)[1]
    if extension == '.gz':
        with gzip.open(filename, "rt") as f:
            lines = [l for l in f if not l.startswith('##')]
        cols = ['ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
        return pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int}, sep='\t', usecols=lambda x: x not in cols).rename(columns={'#CHROM': 'CHROM'})
    else:    
        with open(filename, "r") as f:
            lines = [l for l in f if not l.startswith('##')]
        cols = ['ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
        return pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int}, sep='\t', usecols=lambda x: x not in cols).rename(columns={'#CHROM': 'CHROM'})
    
def numeric_form(df):
    # replace '0|0', '1|0' etc value in dataframe
    GT_conversion = {'./.':None, '0|0':0, '0|1':1, '1|0':1, '1|1':2}
    df = df.replace(GT_conversion)
    df = df.replace(to_replace=r'\d\|\d', value=3, regex=True)
    return df

# read igsr samples and make a list of ethnics corresponding of each sample ID of vcf file
igsr = pd.read_csv("/mnt/iribhm/people/nifernan/1000_genome/igsr_samples.tsv", delimiter="\t")
igsr_dict = {name:code for name, code in zip(igsr['Sample name'], igsr['Superpopulation code'])}

print("===== Extracting filenames from BUCA directory =====")
filenames = glob.glob(os.getcwd()+"/*/*/*/*/*/*/*/*/*.vcf.gz", recursive=True)
print("Total number of files in BUCA directory are : ", len(filenames))

# default value works fine for PCA
pca = PCA()
plt.figure(figsize=(18,12), dpi=300)
temp = pd.DataFrame()
for i, filename in enumerate(filenames):
    if ("chr22" in filename.lower()) or ("chr21" in filename.lower()):
        if i > 150:
            break
        df2 = vcf_parser(filename)
        common = set(df1['POS']).intersection(df2['POS'])
        df1 = df1.loc[df1['POS'].isin(common)]        
        df2 = df2.loc[df2['POS'].isin(common)]
        #df3 = df1[['CHROM', 'POS']]
        df2 = df2.set_index('POS')
        df2 = df2.reindex(index=df1['POS'])
        df2 = df2.reset_index()
        df2 = df2.dropna(how='any', axis=0)
        df1 = numeric_form(df1)        
        print(pd.concat([temp, df1], axis=0).drop_duplicates(keep=False))
        temp = df1
        print(df1.shape, df2.shape)
        cols = list(df1.columns)[2:-1]
        cvec = [igsr_dict[Id] for Id in cols]
        df2 = numeric_form(df2)
        X1 = df1.iloc[:, 2:-1].T
        X2 = df2.iloc[:, 2:].T
        samples = list(df2.columns)[-1]
        #X2.append(list(df2.iloc[:, 2:].T.values))
        # pca for X
        exon_pca = pd.DataFrame(pca.fit_transform(X1))
        # transform the dataframe
        proj = pd.DataFrame(pca.transform(X2))
        
        if "chr22" in filename.lower():
            inp = "22"
        else:
            inp = "21"
      
        ax = sns.scatterplot(x=0, y=1, data=proj, marker='X', s=300)
        plt.setp(ax.lines, zorder=100)
        plt.setp(ax.collections, zorder=100, label=samples+inp)
        G = sns.scatterplot(x=0, y=1, data=exon_pca, hue=cvec, marker='o', linewidth=0.5, s=150, alpha=0.8, ax=ax)
        #ax.legend(markerscale=1)
        G.legend(markerscale=1.5)       
        
#X2 = pd.DataFrame(np.squeeze(np.array(X2)))
#print(len(X1), len(X2))
#print(X2.shape)
#print(X1)
#print(X2)

#print("===== Performing PCA =====")
# color for projection
#cvec2 = ['PCA_samples']*len(proj)

if not os.path.exists("images"):
    os.mkdir("images")

handles, labels = G.get_legend_handles_labels()
unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
G.legend(*zip(*unique))
plt.savefig('images/foo.png', bbox_inches='tight')
plt.show()


print("Check the images folder now!!!")
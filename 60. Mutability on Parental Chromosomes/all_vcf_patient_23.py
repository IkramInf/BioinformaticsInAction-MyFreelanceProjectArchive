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

def numeric_form(df):
    # replace '0|0', '1|0' etc value in dataframe
    GT_conversion = {'./.':None, '0|0':0, '0|1':1, '1|0':1, '1|1':2}
    df = df.replace(GT_conversion)
    df = df.replace(to_replace=r'\d\|\d', value=None, regex=True)
    return df

def legend_without_duplicate_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique))
    
file_names = glob.glob(os.getcwd()+"/*/*/*/*/*/*/*/*/*.vcf.gz", recursive=True)
print("Total number of files in BUCA directory are : ", len(file_names))

fnames = {}
for i, filename in enumerate(file_names):
    ch = filename.lower().split("_output_")[1].split(".txt")[0]
    code = filename.lower().split("/")[-1].split("_output_")[0]
    flag = ch.split("chr")[1]
    if flag.isdigit() and flag != str(23):
        #fname = "buca/" + ch + "_" + str(i) + ".csv"
        #fname = "buca/" + ch + "_" + code + "_" + str(i) + ".csv"
        #df = vcf_parser(filename)
        #df.to_csv(fname, index=False)
        fnames.setdefault(code, []).append(filename)
        
        
iterable_files = glob.glob(os.getcwd()+"/vcf_csvs/*")
all_vcf = []

for i, iterable in enumerate(iterable_files):     
    #if i == 8:
    #    break
    # read the vcf file
    print(f"=== Extracting for {iterable} ===")
    ext = os.path.splitext(iterable)[1]
    if ext == ".csv":
        df1 = pd.read_csv(iterable)
    else:
        df1 = pd.read_csv(iterable, delimiter="\t")
        
    all_vcf.append(df1)
    
df = pd.concat(all_vcf, axis=0)
#print(df.shape)
#df = df.drop_duplicates(subset=['CHROM', 'POS'])
print(df.shape)
print(df.head())

# read igsr samples and make a list of ethnics corresponding of each sample ID of vcf file
igsr = pd.read_csv("/mnt/iribhm/people/nifernan/1000_genome/igsr_samples.tsv", delimiter="\t")
igsr_dict = {name:code for name, code in zip(igsr['Sample name'], igsr['Superpopulation code'])}

#with open('buca/fnames.pkl', 'rb') as f:
#    file_names = pickle.load(f)

plt.figure(figsize=(18,12), dpi=300)    
cmap = {'EUR': 'red', 'EAS': 'green', 'AMR': 'blue', 'SAS': 'orange', 'AFR': 'cyan'}  

with open("coordinates.txt", "w") as f:
    f.write("Co-ordinate : Patient ID\n")
    for i, pair in enumerate(fnames.items()):
        ID = pair[0].split("-")[0]
        #if i == 5:
        #    break

        one_buca = []    
        for j, filename in enumerate(pair[1]):
            #if j == 10:
            #    break
            df1 = vcf_parser(filename)
            one_buca.append(df1)

        dfp = pd.concat(one_buca, axis=0)    
        #print(dfp.shape)
        #dfp = dfp.drop_duplicates(subset=['CHROM', 'POS'])
        print(dfp.shape)
        print(dfp.head())

        df = numeric_form(df)
        dfp = numeric_form(dfp)
        df = df.dropna(how='any', axis=0) 
        dfp = dfp.dropna(how='any', axis=0) 
        common = set(df['POS']).intersection(dfp['POS'])
        print(len(common))
        df = df.loc[df['POS'].isin(common)]       
        print(df.shape)
        dfp = dfp.loc[dfp['POS'].isin(common)]   
        print(dfp.shape)

        df = df.set_index('POS')
        dfp = dfp.set_index('POS')
        df = df[~df.index.duplicated(keep='first')]
        dfp = dfp[~dfp.index.duplicated(keep='first')]
        print(df.shape)
        print(dfp.shape)
        df = df.reset_index()

        #dfp = dfp.set_index('POS')
        dfp = dfp.reindex(index=df['POS'])
        dfp = dfp.reset_index()
        dfp = dfp.dropna(how='any', axis=0) 

        if "REF" in df.columns:       
            X1 = df.iloc[:, 4:-1].T.values
            cols = list(df.columns)[4:-1]
        else:
            X1 = df.iloc[:, 2:-1].T.values
            cols = list(df.columns)[2:-1]

        X2 = dfp.iloc[:, 4:].T.values

        cvec = [igsr_dict[Id] for Id in cols]
        #colors = [cmap[e] for e in cvec]
        le = LabelEncoder()
        cn = le.fit_transform(cvec)
        #colors = ListedColormap(colors)

        print("===== Performing PCA =====")
        # default value works fine for PCA
        pca = PCA()
        # pca for X
        exon_pca = pd.DataFrame(pca.fit_transform(X1))
        # transform the dataframe
        proj = pd.DataFrame(pca.transform(X2))

        # color for projection
        #cvec2 = ['PCA_samples']*len(proj)
        #print(proj.iloc[0,0], proj.iloc[0,1])
        x, y = proj.iloc[0,2], proj.iloc[0,3]
        ax = plt.scatter(x=x, y=y, c=[0], cmap=ListedColormap(['k']), marker='X', s=200, edgecolor='k', label="PCA samples", zorder=2)
        plt.text(x, y, ID, fontsize=10)
        G = plt.scatter(x=exon_pca.iloc[:,2], y=exon_pca.iloc[:,3], c=cn, marker='o', linewidth=0.5, s=100, alpha=0.8, zorder=1)
        f.write(f"({x}, {y}) : {ID}\n")
        #plt.legend(handles=G.legend_elements()[0], labels=cvec)
        #plt.legend()
        #ax.legend(markerscale=1)
        #G.legend(markerscale=1.5)

if not os.path.exists("images"):
    os.mkdir("images")

#plt.legend(markerscale=1)
#legend_without_duplicate_labels(G)
handle1, label1 = ax.legend_elements()    
handle, label = G.legend_elements()
label = list(le.inverse_transform([int(lb.split("mathdefault{")[1].split("}")[0]) for lb in label]))
plt.legend(handles=handle+handle1, labels=label+['PCA samples'])
plt.savefig('images/chr_proj_buca_23.png', bbox_inches='tight')
plt.show()


print("Check the images folder now!!!")
    


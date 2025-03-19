import os
import sys
import requests
import io
import gzip
import glob
import time
import asyncio
import nest_asyncio
nest_asyncio.apply()
import argparse
import numpy as np
import pandas as pd
import plotly
import plotly.express as px
import plotly.graph_objects as go
from sklearn.decomposition import PCA
from pyensembl import EnsemblRelease
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_extraction.text import TfidfVectorizer

start = time.time()

# add command line options
parser = argparse.ArgumentParser(description="Provide input vcf file name.")
parser.add_argument('-i', '--input', help='Input vcf filename', default="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
args = parser.parse_args()


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

# release 75 uses human reference genome GRCh38
ensembl = EnsemblRelease(75)
ensembl.download()
ensembl.index()

chr22 = pd.read_csv("/mnt/iribhm/people/nifernan/biasSomaticMutationsParents/list_SNPs_pcawg/chr22.tsv", delimiter="\t")
chr22_pos = np.array(chr22['POS'])

def get_exon(contig, pos):
    # get exon id
    exon_id = ensembl.exon_ids_at_locus(contig=contig, position=pos)
    if exon_id:
        return exon_id[0]
    else:
        return None

def numeric_form(df):
    # replace '0|0', '1|0' etc value in dataframe
    GT_conversion = {'./.':None, '0|0':0, '0|1':1, '1|0':1, '1|1':2}
    df = df.replace(GT_conversion)
    df = df.replace(to_replace=r'\d\|\d', value=3, regex=True)
    return df

async def do_exon_pca(filename):    
    print("===== Extracting all files =====")
    df = vcf_parser(filename)
    print("===== Total Rows and columns of the file : ", df.shape)
    allele = pd.read_csv("/mnt/iribhm/people/nifernan/biasSomaticMutationsParents/hg19_alleles/chr22.1kg.phase3.v5a_GRCh37nounref_allele_index.txt", delimiter="\t")
    
    print("===== Keeping same positions as allele and chr22 files from the vcf file =====")
    df = df[df['POS'].apply(lambda pos: pos in np.array(allele['position']))]
    df = df[df['POS'].apply(lambda pos: pos in np.array(chr22['POS']))]
    print("===== Total Rows and columns of the file after reducing : ", df.shape)
    L1 = df.shape[0]
    print("===== get exons for each contig and position =====")
    df['exons'] = df.apply(lambda x: get_exon(x.CHROM, x.POS), axis=1)
    # drop rows with NaN values
    df = df.dropna(how = 'any', axis = 0)
    L2 = df.shape[0]
    print(f"Percentage of exons are {(L2/L1)*100}%")
    df = numeric_form(df)
    df = df.iloc[:, 2:].T
    df.columns = df.iloc[-1]
    df = df.drop(df.index[-1])
    df['ID'] = df.index
    df.index = range(df.shape[0])
    igsr = pd.read_csv("/mnt/iribhm/people/nifernan/1000_genome/igsr_samples.tsv", delimiter="\t")
    igsr_dict = {name:code for name, code in zip(igsr['Sample name'], igsr['Superpopulation code'])}
    cvec = [igsr_dict[Id] for Id in df['ID']]
    label_encoder = LabelEncoder()
    df['ID'] = label_encoder.fit_transform(df['ID'])
    # scale data in [0, 1]
    scaler = MinMaxScaler()
    df = scaler.fit_transform(df)
    return (df, cvec)


async def do_projection():
    print("===== Extracting filenames from BUCA directory =====")
    filenames = glob.glob(os.getcwd()+"/*/*/*/*/*/*/*/*/*.vcf.gz", recursive=True)
    print("Total number of files in BUCA directory are : ", len(filenames))
    infos = {}
    for filename in filenames:
        if "chr22" in filename.lower():
            df1 = vcf_parser(filename)
            df1 = df1[df1['POS'].apply(lambda pos: pos in chr22_pos)]
            df1 = df1.drop(['CHROM', 'POS'], axis=1)
            col = df1.columns[0]
            infos.setdefault(col, []).extend(df1[col])
    df = pd.DataFrame(infos)
    print(df.head())
    df = numeric_form(df)
    #samples = df.columns
    df = df.T
    df['samples'] = df.index
    df.index = range(df.shape[0])
    color = df['samples']
    df = df.drop(['samples'], axis=1)
    # scale data in [0, 1]
    scaler = MinMaxScaler()
    df = scaler.fit_transform(df)
    return (df, color)

async def async_main():
    res = await asyncio.gather(do_exon_pca(args.input), do_projection())
    return res

def main():
    res = asyncio.run(async_main())
    return res
    
df_cvec1, df_cvec2 = main()

df1, cvec1 = df_cvec1
df2, cvec2 = df_cvec2

print("===== Performing PCA =====")
# default value works fine for PCA
pca = PCA()
exon_pca = pca.fit_transform(df1)
proj = pca.transform(np.array(df2).reshape(1, -1))

print("===== Plotting =====")
fig1 = px.scatter(exon_pca, x=0, y=1, color=cvec1)
fig2 = px.scatter(proj, color=cvec2, symbol_sequence='x')
fig2.update_traces(marker={'color':'black', 'size': 8})

fig3 = go.Figure(data=fig1.data + fig2.data)
fig3.show()

if not os.path.exists("images"):
    os.mkdir("images")

# html file
plotly.offline.plot(fig3, filename='images/exon_pca_projection.html')    
# png file
fig3.write_image("images/exon_pca_projection.png")
# pdf file
fig3.write_image("images/exon_pca_projection.pdf")

elapsed = (time.time() - start)
print("Total execution time : ", time.strftime("%Hh%Mm%Ss", time.gmtime(elapsed)))

import os
import sys
import requests
import io
import gzip
import argparse
import numpy as np
import pandas as pd
import plotly
import plotly.express as px
from sklearn.decomposition import PCA
from pyensembl import EnsemblRelease
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_extraction.text import TfidfVectorizer


def vcf_parser(filename):
    with open(filename, "r") as f:
        lines = [l for l in f if not l.startswith('##')]
    cols = ['ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    return pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int}, sep='\t', usecols=lambda x: x not in cols).rename(columns={'#CHROM': 'CHROM'})

# release 75 uses human reference genome GRCh38
ensembl = EnsemblRelease(75)
ensembl.download()
ensembl.index()

def get_exon(contig, pos):
    # get exon id
    exon_id = ensembl.exon_ids_at_locus(contig=contig, position=pos)
    if exon_id:
        return exon_id[0]
    else:
        return None

if __name__ == "__main__":
    # add command line options
    parser = argparse.ArgumentParser(description="Provide input vcf file name.")
    parser.add_argument('-i', '--input', help='Input vcf filename', default="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
    args = parser.parse_args()
 
    print("===== Extracting all files =====")
    df = vcf_parser(args.input)
    
    print("===== Total Rows and columns of the file : ", df.shape)
    allele = pd.read_csv("/mnt/iribhm/people/nifernan/biasSomaticMutationsParents/hg19_alleles/chr22.1kg.phase3.v5a_GRCh37nounref_allele_index.txt", delimiter="\t")
    chr22 = pd.read_csv("/mnt/iribhm/people/nifernan/biasSomaticMutationsParents/list_SNPs_pcawg/chr22.tsv", delimiter="\t")
    
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
    
    # replace '0|0', '1|0' etc value in dataframe
    GT_conversion = {'./.':None, '0|0':0, '0|1':1, '1|0':1, '1|1':2}
    df = df.replace(GT_conversion)
    df = df.replace(to_replace=r'\d\|\d', value=3, regex=True)
    
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

    print("===== Performing PCA =====")
    # scale data in [0, 1]
    scaler = MinMaxScaler()
    df = scaler.fit_transform(df)

    # default value works fine for PCA
    pca = PCA()
    exon_pca = pca.fit_transform(df)
    #print(pca.explained_variance_ratio_)

    print("===== Plotting =====")
    fig = px.scatter(exon_pca, x=0, y=1, color=cvec)
    fig.show()
    
    if not os.path.exists("images"):
        os.mkdir("images")
    
    # html file
    plotly.offline.plot(fig, filename='images/exon_pca.html')    
    # png file
    fig.write_image("images/vcf_pca.png")
    # pdf file
    fig.write_image("images/vcf_pca.pdf")

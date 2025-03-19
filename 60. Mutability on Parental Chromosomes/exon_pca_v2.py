import os
import sys
import requests
import gzip
import argparse
import pandas as pd
from pyensembl import EnsemblRelease
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA
import plotly.express as px
from sklearn.feature_extraction.text import TfidfVectorizer


def vcf_parser(filename):
    extension = os.path.splitext(filename)[1]
    if extension == '.gz':
        ifile = gzip.open(filename , "rt")
    else:
        ifile = open(filename, "r")
    with open("vcf.tsv", "w") as ofile:
        for line in ifile:
            if not line.strip().startswith("##"):
                ofile.write("\t".join(line.strip().split("\t")))
                ofile.write("\n")
    ifile.close()    
    

def get_exon(contig, pos):
    # release 77 uses human reference genome GRCh38
    ensembl = EnsemblRelease(75)
    gene_name = ensembl.gene_names_at_locus(contig=contig, position=pos)
    #print(gene_name)
    if gene_name:
        # get all exons associated gene_id[0]
        exon_id  = ensembl.exon_ids_of_gene_name(gene_name[0])
        if exon_id:
            try:
                server = "http://rest.ensembl.org"
                query = "/sequence/id/" + exon_id[0] + "?content-type=text/plain"
                r = requests.get(server+query, headers={ "Content-Type" : "text/plain"})
                if r.ok:
                    return r.text
                else:
                    return None
            except (ValueError, Exception):
                return None
    return None

if __name__ == "__main__":
    # add command line options
    parser = argparse.ArgumentParser(description="Provide input vcf file name.")
    parser.add_argument('-i', '--input', help='Input vcf filename', default="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
    args = parser.parse_args()
 
    print("=====Extracting vcf file into tsv file=====")
    vcf_parser(args.input)

    print("======Loading tsv file and extracting exons=====")
    test = pd.read_csv("vcf.tsv", delimiter="\t").rename(columns={'#CHROM':'CHROM'})

    # replace '0|0', '1|0' etc value in dataframe
    GT_conversion = {'./.':None, '0|0':0, '0|1':1, '1|0':1, '1|1':2}
    test = test.replace(GT_conversion)
    test = test.replace(to_replace=r'\d\|\d', value=3, regex=True)

    # get exons for each contig and position
    test['exons'] = test.apply(lambda x: get_exon(x.CHROM, int(x.POS)), axis=1)

    # unique exons in dataframe
    unique = set(test['exons'])
    print("Number of unique exons are ", len(unique))
    #print(unique)

    # drop rows with NaN values
    test = test.dropna(how = 'any', axis = 0)

    # Represent exons into numeric form
    vectorizer = TfidfVectorizer()
    E = vectorizer.fit_transform(test['exons'])
    E = pd.DataFrame(E.toarray())

    # scale data in [0, 1]
    scaler = MinMaxScaler()
    X = scaler.fit_transform(E)

    # default value works fine for PCA
    pca = PCA()
    exon_pca = pca.fit_transform(X)
    print(pca.explained_variance_ratio_)

    # plotting
    fig = px.scatter(exon_pca, x=0, y=1)
    fig.show()
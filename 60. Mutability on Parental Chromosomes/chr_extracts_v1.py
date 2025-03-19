import os
import sys
import requests
import io
import re
import gzip
import glob
import argparse
import numpy as np
import pandas as pd
from pyensembl import EnsemblRelease
from multiprocessing import Pool

# release 75 uses human reference genome GRCh38
ensembl = EnsemblRelease(75)
ensembl.download()
ensembl.index()
allele = pd.read_csv("/mnt/iribhm/people/nifernan/biasSomaticMutationsParents/hg19_alleles/chr22.1kg.phase3.v5a_GRCh37nounref_allele_index.txt", delimiter="\t")
chr22 = pd.read_csv("/mnt/iribhm/people/nifernan/biasSomaticMutationsParents/list_SNPs_pcawg/chr22.tsv", delimiter="\t")

def vcf_parser(filename):
    extension = os.path.splitext(filename)[1]
    if extension == '.gz':
        with gzip.open(filename, "rt") as f:
            lines = [l for l in f if not l.startswith('##')]
        cols = ['ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
        df = pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int}, sep='\t', usecols=lambda x: x not in cols).rename(columns={'#CHROM': 'CHROM'})
        
        df = df[df['POS'].apply(lambda pos: pos in np.array(allele['position']))]
        df = df[df['POS'].apply(lambda pos: pos in np.array(chr22['POS']))]
        df['exons'] = df.apply(lambda x: get_exon(str(x.CHROM), int(x.POS)), axis=1)
        # drop rows with NaN values
        df = df.dropna(how = 'any', axis = 0)
        m = re.search("ALL\.(.+?)\.", filename, re.I)
        inp = m.group(1)
        df.to_csv(inp+".csv", index=False)
    else:    
        with open(filename, "r") as f:
            lines = [l for l in f if not l.startswith('##')]
        cols = ['ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
        df = pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int}, sep='\t', usecols=lambda x: x not in cols).rename(columns={'#CHROM': 'CHROM'})
        df = df[df['POS'].apply(lambda pos: pos in np.array(allele['position']))]
        df = df[df['POS'].apply(lambda pos: pos in np.array(chr22['POS']))]
        df['exons'] = df.apply(lambda x: get_exon(str(x.CHROM), int(x.POS)), axis=1)
        # drop rows with NaN values
        df = df.dropna(how = 'any', axis = 0)
        m = re.search("ALL\.(.+?)\.", filename, re.I)
        inp = m.group(1)
        df.to_csv(inp+".csv", index=False)

def get_exon(contig, pos):
    # get exon id
    exon_id = ensembl.exon_ids_at_locus(contig=contig, position=pos)
    if exon_id:
        return 1
    else:
        return None
    
if __name__ == "__main__":
    chr_filenames = glob.glob("/mnt/iribhm/people/nifernan/1000_genome/*.vcf", recursive=True)
    chr_n = glob.glob("/mnt/iribhm/people/nifernan/biasSomaticMutationsParents/list_SNPs_pcawg/*.tsv", recursive=True)
    #filenames = list(*zip(chr_filenames, chr_n))
    print(chr_filenames, chr_n)
    fn = ["/mnt/iribhm/people/nifernan/1000_genome/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf"]              
    
    pool = Pool()
    pool.map_async(vcf_parser, fn)
    pool.close()
    pool.join()

    
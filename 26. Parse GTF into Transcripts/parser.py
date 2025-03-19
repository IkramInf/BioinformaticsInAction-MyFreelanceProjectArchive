#!/usr/bin/env python
# coding: utf-8

# to install gtfparse : pip install gtfparse

from Bio.SeqIO import parse
from gtfparse import read_gtf
from collections import defaultdict

# read the fasta file to a dictionary as 'geneId : sequence' pair
f1 = {record.id : str(record.seq) for record in parse("mouse_w200_SuperTrans.fasta", "fasta")}
# read the gtf file
table = read_gtf("mouse_w200.SuperTrans.gtf")

# only keep 'gene_id', 'transcript_id', 'start', 'end' columns in table
table = table[['gene_id', 'transcript_id', 'start', 'end']]
# take exons from sequence by start and end index
table['exons'] = [f1[g][s:e] for g, s, e in zip(table.gene_id, table.start, table.end)]
# drop start and end columns
table = table.drop(['start', 'end'], axis=1)

# group the table by the column 'gene_id'
df1 = table.groupby('gene_id')
# save the groups into a variable
group1 = df1.groups

# create a dictionary to store list value
genes = defaultdict(list)

# save 'gene_id : [transcript_id]' into genes dictionary
for Id in table.gene_id:
    genes[Id] = list(set([table.at[i, 'transcript_id'] for i in group1[Id]]))

# group the table by the column 'transcript_id'
df2 = table.groupby('transcript_id')
# save the groups into a variable
group2 = df2.groups

# create a dictionary to store transcripts
transcripts = defaultdict(str)
# save 'transcript_id : transcripts' into transcripts dictionary
for Id in table.transcript_id:
    seq = "".join([table.at[i, 'exons'] for i in group2[Id]])
    transcripts[Id] = seq

# write to a fasta file
with open("output.fasta", "w") as f:
    for g, t in genes.items():
        f.write(f"Gene ID : {g}\n")
        for k in t:
            f.write(f"Transcript_{k} : {transcripts[k]}\n")
        f.write('-' * 80 + "\n")


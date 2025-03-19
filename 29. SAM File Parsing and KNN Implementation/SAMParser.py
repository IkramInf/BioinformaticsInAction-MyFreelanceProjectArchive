# import libraries
import pysam
from collections import defaultdict

# read bam file wth pysam
samfile = pysam.AlignmentFile("sample_sort.bam", "rb")

# read gene annotations from gtf file and store choromosome, start, end and gene name into annotations list
annotations = []
with open("Arabidopsis.gtf", "r") as f:
    lines = f.readlines()
    for line in lines:
        line = line.split("\t")
        if line[2] == "exon":
            annotations.append((line[0],line[3],line[4],line[-1].split(' gene_name "')[1].replace('";\n',"")))


# count the number of reads with pysam            
gene_counts = defaultdict(int)
for annot in annotations:
    if annot[0] not in ['ChrC', 'ChrM']:
        count = samfile.count(annot[0], int(annot[1]), int(annot[2]))
        gene_counts[annot[3]] += count
    else:
        continue


# sorted the dictionary in descending order
gene_counts = dict(sorted(gene_counts.items(), key=lambda x: x[1], reverse=True))

# show output as matrix
gene_counts = [[k, v] for k, v in gene_counts.items()]
print(gene_counts)        
#import required libraries
# to install pysam : pip install pysam
import pysam

# read vcf file
vcf = pysam.VariantFile("Final.recode.outlier.vcf")
# read fasta file
genome = pysam.FastaFile("GCF_000214255.1_Bter_1.0_genomic.fna")
# define by how many bases the variant should be flanked
flank = 500

# store all snps into a list
snps = []

# iterate over each variant
for record in vcf:
    # python index starts from 0 but vcf position starts from 1...
    # so, substract 1 from position to use as python index
    pos = record.pos - 1
    chrom, ref, alt = record.chrom, record.ref, record.alts[0]
    # length of ref
    L = len(ref)
    if pos < flank:
        flank = pos
        seq = genome.fetch(reference=chrom, start=pos-flank, end=pos+L+flank)
        seq_format = f"{seq[:flank]}[{ref}/{alt}]{seq[flank+L:]}"
    else:
        seq = genome.fetch(reference=chrom, start=pos-flank, end=pos+L+flank)
        seq_format = f"{seq[:flank]}[{ref}/{alt}]{seq[flank+L:]}"
    snps.append(seq_format)


# write results into a fasta file
with open("output.fasta", "w") as f:
    for i, seq in enumerate(snps):
        f.write(f">snp{i+1}\n{seq}\n")

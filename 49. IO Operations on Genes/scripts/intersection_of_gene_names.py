from assignment4.io_utils import get_filehandle
from argparse import ArgumentParser

# add command line options
parser = ArgumentParser(description="Provide two gene list (ignore header line), find intersection")
parser.add_argument('-i1', '--infile1', help='Gene list 1 to open', required=True)
parser.add_argument('-i2', '--infile2', help='Gene list 2 to open', required=True)
args = parser.parse_args()

chr21 = get_filehandle(args.infile1, "r")
hugo = get_filehandle(args.infile2, "r")
iw = get_filehandle("OUTPUT/intersection_output.txt", "w")

chr_genes = set([line.strip().split("\t")[0] for line in chr21.readlines()[1:]])
hugo_genes = set([line.strip().split("\t")[0] for line in hugo.readlines()[1:]])
common = chr_genes.intersection(hugo_genes)
iw.write("\n".join(sorted(common)))

print(f"Number of unique gene names in {args.infile1}: {len(chr_genes)}")
print(f"Number of unique gene names in {args.infile2}: {len(hugo_genes)}")
print("Number of common gene symbols found: ", len(common))
print("Output stored in OUTPUT/intersection_output.txt")
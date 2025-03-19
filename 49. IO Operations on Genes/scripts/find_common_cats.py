from assignment4.io_utils import get_filehandle
from argparse import ArgumentParser

# add command line options
parser = ArgumentParser(description="Combine on gene name and count the category occurrence")
parser.add_argument('-i1', '--infile1', help='Path to the gene description file to open', required=True)
parser.add_argument('-i2', '--infile2', help='Path to the gene category to open', required=True)
args = parser.parse_args()

chr21 = get_filehandle(args.infile1, "r")
desp = get_filehandle(args.infile2, "r")
cw = get_filehandle("OUTPUT/categories.txt", "w")

categories = {}
for line in chr21.readlines()[1:]:
    elements = line.strip().split("\t")
    categories.setdefault(elements[-1], []).append(elements[0])
    
category_desp = {}
for line in desp.readlines():
    data = line.strip().split("\t")
    category_desp.update({data[0] : data[1]})
    
cw.write("Category\tOccurrence\tDescription\n")
categories = dict(sorted(categories.items(), key=lambda x:x[0].lower()))
for k, v in categories.items():
    if k in category_desp.keys():
        cw.write(f"{k}\t{len(v)}\t{category_desp[k]}\n")
        
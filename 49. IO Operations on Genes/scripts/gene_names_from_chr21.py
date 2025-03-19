from assignment4.io_utils import get_filehandle
from argparse import ArgumentParser

# add command line options
parser = ArgumentParser(description="Open chr21_genes.txt, and ask user for a gene name")
parser.add_argument('-i', '--infile', help='Path to file to open', required=True)
args = parser.parse_args()

chr21 = get_filehandle(args.infile, "r")
chr21_genes = {}
for line in chr21.readlines()[1:]:
    elements = line.strip().split("\t")
    chr21_genes.update({elements[0].lower() : "".join(elements[1:])})
    
while True:
    name = input("Enter gene name of interest. Type quit to exit: ")
    
    if name.lower() in ["quit", "exit"]:
        print("Thanks for querying the data.")
        break
        
    if name.lower() in chr21_genes.keys():
        print(f"{name} found! Here is the description:\n{chr21_genes[name.lower()]}\n")
    else:
        print("Not a valid gene name.")
        
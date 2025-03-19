import sys
# pip install bcbio-gff
from BCBio import GFF

# separate command line options i.e. -i, -c, -s, -e and -o
options = [opt for opt in sys.argv[1:] if opt.startswith("-")]
# separate command line arguments
args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]

# create a dictionary of commandline options and commandline arguments
arguments = dict(zip(options, args))

# list all chromosomes
with open("TAIR10_GFF3_genes.gff", "r") as handle:
    chromosomes = [record.id for record in GFF.parse(handle)]

# create a list to store all gene names    
genes = []

# to check -c option is given or not in command line
if '-c' in options:
    chromosomes = arguments['-c'].split()
    
limit_info = {"gff_id" : chromosomes}

# read input file as read mode
with open("TAIR10_GFF3_genes.gff", "r") as handle:
    # parse gff file
    for record in GFF.parse(handle, limit_info=limit_info):
        for feature in record.features:
            # only take the gene feature
            if feature.type == "gene":
                # to check -s and -e option are given or not in command line
                if '-s' in options and '-e' in options:
                    if feature.location.start > int(arguments['-s']) and feature.location.end < int(arguments['-e']):
                        genes.append(feature.id)
                else:
                    # adding gene names to genes list
                    genes.append(feature.id)
                        
# to check -o option is given or not in command line                        
if '-o' in options:
    # open a file as write mode as given name in command line
    with open(arguments['-o'], "w") as f:
        for name in genes:
            # write to output file
            f.write(f"{name}\n")
else:
    # if -o is not given, directly print gene names
    print("\n".join(genes))
    #print(len(genes))
       
# read file and save gene and ib numbers into a dictionary
genes = {}
with open("blat_output.psl", "r") as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip().split("\t")
        gene = line[13].split('|')[-1].split('-')[0]
        genes.setdefault(gene, []).append(line[9])       

# write the output into a tsv file
with open("score_genes.tsv", "w") as writer:
    for gene, ib in genes.items():
        writer.write(f"{gene}\t{len(ib)}\n")


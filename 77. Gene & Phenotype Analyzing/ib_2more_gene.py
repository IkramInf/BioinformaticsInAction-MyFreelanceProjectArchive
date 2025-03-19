# extract iBs and genes from Off-Target-Gene-Liste.tsv
ibs = {}
with open("Off-Target-Gene-Liste.tsv", "r") as f:
    for line in f.readlines():
        data = line.strip().split("\t")
        if len(data) > 2:
            ibs.setdefault(data[0], data[1:])
            
# write output into a tsv file
with open("ib_2more_genes.tsv", "w") as f:
    f.write("iBs (with more than 1 targeted gene)\tGenes\n")
    for ib, gene in ibs.items():
        f.write(f"{ib}\t{gene}\n")



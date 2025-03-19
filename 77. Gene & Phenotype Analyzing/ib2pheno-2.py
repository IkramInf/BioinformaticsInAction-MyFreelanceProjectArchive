# import required libraries
import ast
import requests

# extract iBs and genes from Off-Target-Gene-Liste.tsv
ibs = {}
with open("Off-Target-Gene-Liste.tsv", "r") as f:
    for line in f.readlines():
        data = line.strip().split("\t")
        if len(data) == 2:
            ibs.setdefault(data[0], data[1])
            
# extract phenotype over internet
base_url = "http://ibb-test.vm19002.virt.gwdg.de/ibb/api/lethality/v1/lethality/"
phenotypes = {ib : requests.get(base_url+ib).text for ib in ibs.keys()}

# write output into a tsv file
with open("ib_phenotype.tsv", "w") as f:
    f.write("iBs (with only 1 targeted gene)\tGenes\tPhenotype\n")
    for ib, gene in ibs.items():
        pheno = ast.literal_eval(phenotypes[ib])
        flag = []
        for v in pheno.values():
            if not v == -1:
                flag.append(True)
            else:
                flag.append(False)
        if any(flag):
            f.write(f"{ib}\t{gene}\t{pheno}\n")


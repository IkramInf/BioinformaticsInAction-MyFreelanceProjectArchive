# read query text file into a dictionary as chrom:[pos] pair
with open("queryRDA.txt", "r") as f:
    queries = {}
    for line in f.readlines()[1:]:
        values = line.strip().split("\t")
        if len(values) == 2:
            chrom, pos = values
            queries.setdefault(chrom, []).append(pos)

# read the vcf file and write the subset file        
with open("Final.recode.outlier01.vcf", "r") as ifile, open("Final.recode.outlier01_subset.vcf", "w") as ofile:
    # read all lines as a list
    vcf = ifile.readlines()
    # find index of #CHROM line
    index = [i for i, line in enumerate(vcf) if line.strip().startswith("#CHROM")][0]
    # separate all header informations
    headers = "".join(vcf[:index+1])
    ofile.write(headers)
    # check if the chrom and pos of vcf file present in query text file
    for line in vcf[index+1:]:
        infos = line.strip().split("\t")
        chrom, pos = infos[0].strip(), infos[1].strip()
        if chrom in queries.keys():
            if pos in queries[chrom]:
                ofile.write(line)


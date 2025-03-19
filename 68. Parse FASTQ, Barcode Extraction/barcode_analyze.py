# import required libraries
import os

# extract barcodes from harrington_clinical_data.txt
barcodes = {}
with open("harrington_clinical_data.txt", "r") as f:
    for line in f.readlines()[1:]:
        data = line.strip().split("\t")
        barcodes.setdefault(data[-1], []).append(data[0])

def fd_index(qual):
    i = len(qual)
    if "DD" in qual:
        i1 = qual.index("DD")
        if i1 < i:
            i = i1
    if "FF" in qual:
        i1 = qual.index("FF")
        if i1 < i:
            i = i1
    if "DF" in qual:
        i1 = qual.index("DF")
        if i1 < i:
            i = i1
    if "FD" in qual:
        i1 = qual.index("FD")
        if i1 < i:
            i = i1

    if i == len(qual):
        return 0
    else:
        return i
        
# parsing the fastq file        
with open("hawkins_pooled_sequences.fastq", "r") as fq:
    lines = fq.readlines()
    index = [d for d in range(len(lines)) if lines[d].strip().startswith("@")]
    # extract title and sequence
    seqs = {}
    for i in range(len(index)):
        if index[i] != index[-1]:
            info = lines[index[i]:index[i+1]]
            seq = info[1].strip()
            qual = info[-1].strip()
            ind = fd_index(qual)
            if ind == 0:
                seqs.setdefault(seq, "".join(info).strip())
            else:
                info = info[0] + seq[0:ind] + "\n" + info[2] + qual[0:ind]
                seqs.setdefault(seq, info.strip())
        else:
            info = lines[index[i]:]
            seq = info[1].strip()
            qual = info[-1].strip()
            ind = fd_index(qual)
            print("index : ", ind)
            if ind == 0:
                seqs.setdefault(seq, "".join(info).strip())
            else:
                info = info[0] + seq[0:ind] + "\n" + info[2] + qual[0:ind]
                print(seq)
                print(seq[0:ind+1])
                print(info)
                seqs.setdefault(seq, info.strip())
            
# matching between genetic barcodes with fastq sequence             
results = {}
for seq in seqs.keys():
    for barcode in barcodes.keys():
        if barcode in seq:
            results.setdefault(barcode, []).append(seq)

# create the output directory            
if not os.path.exists("fastqs"):
    os.mkdir("fastqs")

# write the results    
for barcode, seq in results.items():
    dir_name = "fastqs/"+barcode+"/"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    with open(dir_name+barcodes[barcode][0]+".fastq", "w") as f:
        for s in seq:
            f.write(f"{seqs[s]}\n")    
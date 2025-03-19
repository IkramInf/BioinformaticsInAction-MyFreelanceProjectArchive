# import required libraries
import os
import re
import glob
import shutil
import pysam
import subprocess

##### Demultiplexing pooled fastq data #####
# import required libraries

# extract barcodes from harrington_clinical_data.txt
barcodes = {}
with open("harrington_clinical_data.txt", "r") as f:
    for line in f.readlines()[1:]:
        data = line.strip().split("\t")
        barcodes.setdefault(data[-1], []).append(data[0])

def trim_barcode(barcode, x, z):
    m = re.search(barcode, x, re.M|re.S)
    s, e = m.start(), m.end()
    sqs = x[0:s] + x[e:]
    quals = z[0:s] + z[e:]
    return (sqs, quals)
                
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
        
# parsing the fastq file and trimming each sequence       
with open("hawkins_pooled_sequences.fastq", "r") as fq:
    lines = fq.readlines()
    index = [d for d in range(len(lines)) if lines[d].strip().startswith("@")]
    # extract title and sequence
    seqs = {}
    for i in range(len(index)):
        if index[i] != index[-1]:
            w, x, y, z = lines[index[i]:index[i+1]]
            x, z = x.strip(), z.strip()
            for barcode in barcodes.keys():
                if x.startswith(barcode):
                    xx, zz = trim_barcode(barcode, x, z)
                    ind = fd_index(zz)
                    if ind == 0:
                        wxyz = "".join(w + xx + "\n" + y + zz)
                        seqs.setdefault(barcode, []).append(wxyz)
                    else:
                        wxyz = "".join(w + xx[0:ind] + "\n" + y + zz[0:ind])
                        seqs.setdefault(barcode, []).append(wxyz)
        else:
            w, x, y, z = lines[index[i]:]
            x, z = x.strip(), z.strip()
            for barcode in barcodes.keys():
                if barcode in x:
                    xx, zz = trim_barcode(barcode, x, z)
                    ind = fd_index(zz)
                    if ind == 0:
                        wxyz = "".join(w + xx + "\n" + y + zz)
                        seqs.setdefault(barcode, []).append(wxyz)
                    else:
                        wxyz = "".join(w + xx[0:ind] + "\n" + y + zz[0:ind])
                        seqs.setdefault(barcode, []).append(wxyz)

# create the output directory            
if not os.path.exists("fastqs"):
    os.mkdir("fastqs")

# write the results    
for barcode, seq in seqs.items():
    dir_name = "fastqs/"+barcode+"/"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    with open(dir_name+barcodes[barcode][0]+"_trimmed.fastq", "w") as f:
        for data in seq:
            f.write(f"{data}\n")    
            
##### Perform alignment on each FASTQ to the reference sequence #####
##### Convert sam files into bam files #####
##### Discover variants #####

with open("dgorgon_reference.fasta", "r") as f:
    header = f.readline()[1:]
    reference = f.read().replace("\n", "")

indexed_ref = subprocess.run(['bwa', 'index', 'dgorgon_reference.fasta'], stdout=subprocess.PIPE)

# extract colors from harrington_clinical_data.txt
colors = {}
with open("harrington_clinical_data.txt", "r") as f:
    for line in f.readlines()[1:]:
        data = line.strip().split("\t")
        colors.setdefault(data[0], data[1])            

def pileup(sorted_bam_file):
    #test file, replaced with the sorted.bam you are using. Make sure it is indexed! (Use samtools index yourbam.sorted.bam)
    samfile = pysam.AlignmentFile(sorted_bam_file, "rb")
    name = sorted_bam_file.split("/")[-1].split(".sorted")[0]
    #use a dictionary to count up the bases at each position
    ntdict = {}
    counts = {}
    #Since our reference only has a single sequence, we're going to pile up ALL of the reads. Usually you would do it in a specific region (such as chromosome 1, position 1023 to 1050 for example)
    for pileupcolumn in samfile.pileup():       
        #print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        N = pileupcolumn.n
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                #print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))                
                pos = pileupread.query_position
                base = pileupread.alignment.query_sequence[pos]
                # hold all of the base read counts per nucletoide per position into a dictionary
                if len(reference) > pos:
                    ref = reference[pos]
                    if base != ref:
                        ntdict.setdefault(pos, []).append((base, ref))
                        counts[name] = counts.get(name, 0) + 1
                        
    return ntdict, counts, N
    samfile.close()
    
#nt = pileup("bams/Shelby.sorted.bam")

if not os.path.exists("temp"):
    os.mkdir("temp")
dir_names = glob.glob(os.getcwd()+"/fastqs/*", recursive=True)
if not os.path.exists("bams"):
    os.mkdir("bams")

pileup_results = {}    
for dir_name in dir_names:
        
    sequences = glob.glob(dir_name+"/*.fastq")   
    #print(sequences)
    name = sequences[0].split("/")[-1].split("_trimmed")[0]
    sorted_bam_file = "bams/" + name + ".sorted.bam"
    
    p1 = subprocess.run(f"bwa mem dgorgon_reference.fasta {sequences[0]} > temp/{name}.sam", stdout=subprocess.PIPE, text=True, shell=True)
    p2 = subprocess.run(f"samtools view -bS temp/{name}.sam > temp/{name}.bam", stdout=subprocess.PIPE, text=True, shell=True)
    p3 = subprocess.run(f"samtools sort -m 100M -o {sorted_bam_file} temp/{name}.bam", stdout=subprocess.PIPE, text=True, shell=True)
    p4 = subprocess.run(f"samtools index {sorted_bam_file}", stdout=subprocess.PIPE, text=True, shell=True)
    
    pileup_results[name] = pileup(sorted_bam_file)
    
    
with open("report.txt", "w") as f:    
    for k, v in pileup_results.items():
        for k1, v1 in v[0].items():
            b = max(v1, key=v1.count)
            f.write(f"The {colors[k]} mold was caused by a mutation in position {k1}. The wildtype base was {b[1]} and the mutation was {b[0]}.\n")
    f.write("\n")
    for k, v in pileup_results.items():
        for k1, v1 in v[0].items():
            f.write(f"Sample {k} had a {colors[k]} mold, {v[1][k]} reads, and had {int((v[1][k]/v[2])*100)}% of the reads at position {k1} had the mutation {v1[0][1]}.\n")    
            
# remove temp directory
shutil.rmtree("temp")
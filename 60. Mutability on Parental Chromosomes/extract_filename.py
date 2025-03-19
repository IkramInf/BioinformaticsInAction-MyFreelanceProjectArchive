import os
import sys
import glob
import pickle


print("===== Extracting filenames =====")
file_names = glob.glob(os.getcwd()+"/*/*/*/*/*/*/*/*/*.vcf.gz", recursive=True)
print("Total number of files in the directory are : ", len(file_names))

fnames = {}
for i, filename in enumerate(file_names):
    ch = filename.lower().split("_output_")[1].split(".txt")[0]
    code = filename.lower().split("/")[-1].split("_output_")[0]
    flag = ch.split("chr")[1]
    if flag.isdigit():
        #fname = "buca/" + ch + "_" + str(i) + ".csv"
        #fname = "buca/" + ch + "_" + code + "_" + str(i) + ".csv"
        #df = vcf_parser(filename)
        #df.to_csv(fname, index=False)
        fnames.setdefault(code, []).append(filename)
    
    
with open('names/fnames_buca.pkl', 'wb') as dbfile:
    pickle.dump(fnames, dbfile, pickle.HIGHEST_PROTOCOL)

with open('names/fnames_buca.pkl', 'rb') as f:
    files = pickle.load(f)
    print([len(v) for k, v in files])
    
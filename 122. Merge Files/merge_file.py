#!/usr/bin/env python
# coding: utf-8

# required libraries
import os
import re
import glob
import gzip
import subprocess
from collections import defaultdict

# extract .fastq.gz filenames
folder_location = input("Enter the input folder name: ")
if folder_location.endswith("/"):
    folder_location = folder_location[:-1]
filenames = glob.glob(f"{folder_location}/**/*.fastq.gz", recursive=True)
#print(filenames)


# keep all file as [individual][R1/2][filenames]
file_dict = defaultdict(dict)
s_dict = defaultdict(dict)
for filename in filenames:
    name = filename.split("/")[-1]
    individual = name.split("_S")[0]
    try:
        S = int(re.findall("_S(\d+)_", name)[0])
        R = re.findall("_(R\d+)_", name)[0]
    except:
        print(f"FileNameError: Invalid filename of {filename}.")
        continue
    if R not in file_dict[individual]:
        file_dict[individual][R] = [filename]
        s_dict[individual][R] = S
    else:
        file_dict[individual][R].append(filename)
        s_dict[individual][R] += S

# merge files
print("Wait! until finishing merging files. You will have a Output folder in your current path.")
if not os.path.exists("Outputs/"):
    os.mkdir("Outputs/")

for individual, infos in file_dict.items():
    for R, fnames in infos.items():
        outfile = "Outputs/" + re.sub("_S\d+_", f"_S{s_dict[individual][R]}_", fnames[0]).split("/")[-1]
        fnames = " ".join(fnames)
        p = subprocess.run(f"zcat {fnames} | pigz -p {os.cpu_count()} > {outfile}", shell=True)

print("Successfully merged all files.")
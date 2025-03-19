#!/usr/bin/env python
# coding: utf-8

# read HWE.txt file in read mode
with open("HWE.txt", "r") as f:
    # add chr and pos to results for entry whose obs het > 60
    results = []
    # iterate over each line except first header line
    for line in f.readlines()[1:]:
        # extract chr, pos and obs entry
        Chr, pos, obs = line.strip().split("\t")[:3]
        # extract obs het from obs
        obs_het = int(obs.split("/")[1])
        
        # check if obs het > 60
        if obs_het > 60:
            # write it into results list
            results.append(f"{Chr}\t{pos}")

# create a file output.txt and write results into it
with open("output.txt", "w") as f:
    # write header CHR and POS
    print("CHR\tPOS", file=f)
    # write each entry of the results list
    for entry in results:
        print(entry, file=f)

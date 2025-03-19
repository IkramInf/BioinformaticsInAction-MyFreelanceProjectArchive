#!/bin/bash

# Specify Input Filepath
vcf_path=${1:-"trios.vcf"}

# Print headers to stdout
grep "^#" "$vcf_path"

awk -F'\t' '
!/^#/ {
    if ($7 != "PASS") next
    
    split($9, format, ":")
    for(i=1; i<=length(format); i++) {
        if(format[i] == "GT") gt_idx = i
        if(format[i] == "GQ") gq_idx = i
        if(format[i] == "DP") dp_idx = i
    }
    
    split($10, child, ":")
    split($11, father, ":")
    split($12, mother, ":")
    
    if(child[gt_idx] ~ /0\/0|0\|0/) next
    if(father[gt_idx] !~ /0\/0|0\|0/) next
    if(mother[gt_idx] !~ /0\/0|0\|0/) next
    if(child[gq_idx] < 20 || father[gq_idx] < 20 || mother[gq_idx] < 20) next
    if(child[dp_idx] < 10 || father[dp_idx] < 10 || mother[dp_idx] < 10) next
    
    # if de novo mutations found, print the entire line to stdout
    print $0
}' "$vcf_path"


#!/usr/bin/env bash

# Specify Input Filepath
vcf_path=${1:-"trios.vcf"}
echo "Input VCF path: $vcf_path"

awk -F'\t' '
!/^#/ {
	# Initialize variables
    ref_len = length($4)
    max_alt_len = 0
    is_snp = 1
    chrom = $1
    
    # Process ALT field
    split($5, alts, ",")
    for (i in alts) {
        alt_len = length(alts[i])
        if (alt_len > max_alt_len) max_alt_len = alt_len
        if (alt_len != 1) is_snp = 0
    }
    
    # Assign variants into SNP, INS, DEL and Other
    if (ref_len == 1 && is_snp) {
        counts[chrom,"SNP"]++
    } else if (ref_len > max_alt_len) {
        counts[chrom,"Deletion"]++
    } else if (ref_len < max_alt_len) {
        counts[chrom,"Insertion"]++
    } else {
        counts[chrom,"Other"]++
    }
    
    # Track unique chromosome once
    if (!(chrom in seen)) {
        chroms[chrom] = chrom
        seen[chrom] = 1
    }
}

END {
    # Print results for each chromosome
    for (chrom in chroms) {
        printf "Chromosome %s:\n", chrom
        printf "  Deletion: %d\n", (counts[chrom,"Deletion"] ? counts[chrom,"Deletion"] : 0)
        printf "  Insertion: %d\n", (counts[chrom,"Insertion"] ? counts[chrom,"Insertion"] : 0)
        printf "  SNP: %d\n", (counts[chrom,"SNP"] ? counts[chrom,"SNP"] : 0)
        printf "  Other: %d\n", (counts[chrom,"Other"] ? counts[chrom,"Other"] : 0)
    }
}' "$vcf_path"


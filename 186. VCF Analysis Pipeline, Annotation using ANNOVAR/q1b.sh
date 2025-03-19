#!/usr/bin/env bash

# Specify Input Filepath
vcf_path=${1:-"trios.vcf"}
echo "Input VCF path: $vcf_path"

# Summarize INFO fields
awk -F'\t' '
BEGIN {
    # Initialize variant counters
    total_variants = 0
    
    # Specify key metrics to track
    split("AC,AF,AN,DP,MQ,QD,FS,VQSLOD", key_metrics, ",")
    for (i in key_metrics) {
        metric = key_metrics[i]
        min[metric] = 999999
        max[metric] = -999999
        sum[metric] = 0
        count[metric] = 0
    }
}

# Processing INFO field
!/^#/ && $8 {
    total_variants++
    
    # Split INFO field
    n = split($8, info_fields, ";")
    
    # Process each field
    for (i = 1; i <= n; i++) {
        # Split each one into key-value pairs
        split(info_fields[i], kv, "=")
        key = kv[1]
        value = kv[2]
        
        # Store count values into key metrics
        for (m in key_metrics) {
            metric = key_metrics[m]
            if (key == metric && value != "") {
                count[metric]++
                sum[metric] += value + 0  # Convert to number
                
                # Update min and max
                if (value + 0 < min[metric]) min[metric] = value + 0
                if (value + 0 > max[metric]) max[metric] = value + 0
            }
        }
        
        # Track unique combinations
        info_combinations[info_fields[i]]++
    }
}

END {
    print "Total variants processed:", total_variants

    for (m in key_metrics) {
        metric = key_metrics[m]
        if (count[metric] > 0) {
			printf "%s:\n", metric
			printf "  Count: %d\n", count[metric]
			printf "  Range: %.2f to %.2f\n", min[metric], max[metric]
			printf "  Average: %.3f\n", sum[metric] / count[metric]
		}
    }
    
    print "\nQuality Metrics in INFO Field Used:"
    if (count["QD"] > 0) print "- Quality by Depth (QD)"
    if (count["FS"] > 0) print "- Fisher Strand (FS)"
    if (count["MQ"] > 0) print "- RMS Mapping Quality (MQ)"
    if (count["VQSLOD"] > 0) print "- VQSR LOD Score (VQSLOD)"
    
    print "\nDepth Statistics:"
    if (count["DP"] > 0) {
        printf "Average depth: %.1f\n", sum["DP"]/count["DP"]
        printf "Depth range: %d to %d\n", min["DP"], max["DP"]
    }
    
    print "\nAllele Statistics:"
    if (count["AF"] > 0) {
        printf "Average allele frequency: %.3f\n", sum["AF"]/count["AF"]
        printf "AF range: %.2f to %.2f\n", min["AF"], max["AF"]
    }
}' "$vcf_path"


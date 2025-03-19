# Trimming Reads
filtlong --assembly rsv.fa --min_length 1000 --min_mean_q 15 --target_bases 500000000 --keep_percent 50 --mean_q_weight 10 --trim --split 1000 barcode08.filtered.fastq.gz | gzip > barcode08.filtered.final.fastq.gz
# Index Reference genome 
samtools faidx rsv.fa
# Mapping and sorting 
minimap2 -ax map-ont rsv.fa barcode08.filtered.fastq.gz | samtools sort -o barcode08.sorted.bam
samtools index barcode08.sorted.bam
# Variant calling and extracting consensus seq 
medaka_consensus -i barcode08.filtered.fastq.gz -d rsv.fa -o medaka_output/


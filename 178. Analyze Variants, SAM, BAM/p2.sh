#!/bin/bash

# Handling cmd arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <reference_genome> <fastq1> <fastq2> <threads> <vcf_filepath>"
    exit 1
fi

# Assign cmd arguments into variables
REFERENCE=$1
FASTQ1=$2
FASTQ2=$3
THREADS=$4
VCF_PATH=$5

# Output files
REFERENCE_BASE=$(basename "$REFERENCE" .fa)
RAW_VCF="file.vcf"
SAM_FILE="aligned.sam"
UNSORTED_BAM="aligned.bam"
SORTED_BAM="aligned.sorted.bam"

# A. Index the reference genome with bwa
echo "A. Indexing the reference genome with bwa..."
bwa index -p "$REFERENCE_BASE" "$REFERENCE"

# B & C. Align paired-end reads to the reference index with bwa mem and sort the output
echo -e "\nB & C. Aligning reads with BWA and sorting output..."
bwa mem -t "$THREADS" "$REFERENCE_BASE" "$FASTQ1" "$FASTQ2" > "$SAM_FILE"
samtools view -@ "$THREADS" -Sb "$SAM_FILE" > "$UNSORTED_BAM"
samtools sort -@ "$THREADS" "$UNSORTED_BAM" -o "$SORTED_BAM"
samtools index "$SORTED_BAM"
#rm "$SAM_FILE" "$UNSORTED_BAM"

# D. Call variants and generate raw VCF file
echo -e "\nD. Calling variants with bcftools..."
bcftools mpileup -Ou -f "$REFERENCE" "$SORTED_BAM" | bcftools call -mv -Ov -o "$RAW_VCF"

# E. Final vcf file by filtering low-confidence variants
echo -e "\nE. Filtering low-confidence variants..."
bcftools filter -i 'QUAL>20 && DP>20' "$RAW_VCF" -Ov -o "$VCF_PATH"

# Count the variant calls
CALLS=$(grep -v "^#" "$VCF_PATH" | wc -l)
echo -e "Total number of high-confidence variant calls: $CALLS"


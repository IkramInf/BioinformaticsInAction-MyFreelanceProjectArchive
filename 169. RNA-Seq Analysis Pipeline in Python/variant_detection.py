## Required Libraries
import os
import re
import sys
import glob
import argparse
import subprocess

### Extract command line arguments ###
parser = argparse.ArgumentParser(description='RNA-SEQ Analysis.')
parser.add_argument('-m', '--mode', type=str, choices=['SE', 'PE'], required=True, help='Mode of sequencing')
parser.add_argument('-r1', '--r1_fastq', type=str, required=True, help='Path to the R1 FASTQ file')
parser.add_argument('-r2', '--r2_fastq', type=str, required=False, help='Path to the R2 FASTQ file')
parser.add_argument('-r', '--reference_genome', type=str, required=True, help='Path to the reference genome')
parser.add_argument('-g', '--gtf', type=str, required=True, help='Path to the GTF file')
args = parser.parse_args()

# Create Output Directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# Extract cmd args
mode = args.mode
genome_path = args.reference_genome
gtf_path = args.gtf
numthreads = "4"

# Define genome index file
index_prefix = f"{output_dir}/genome_index"
## Build HISAT2 index ##
# Check if the index already exists
index_files = [f"{index_prefix}.{i+1}.ht2" for i in range(8)]
if all(os.path.exists(file) for file in index_files):
    print("HISAT2 index already exists. Skipping index building.")
else:
    subprocess.run(["hisat2-build", "-p", numthreads, genome_path, index_prefix])

# Define output SAM and BAM file
output_sam = f"{output_dir}/aligned_reads.sam"
output_bam = f"{output_dir}/aligned_reads.bam"
sorted_bam = f"{output_dir}/aligned_reads_sorted.bam"

### ANALYSIS ###
if mode == "PE":
    r1_fastq_path = args.r1_fastq
    r2_fastq_path = args.r2_fastq
    TRUSEQ = "TruSeq3-PE-2.fa"

    ## Quality Control
    subprocess.run(["fastqc", r1_fastq_path, r2_fastq_path, "-o", output_dir])
    
    ## Trimming
    fname1 = os.path.basename(r1_fastq_path).split("R1")[0]
    fname2 = os.path.basename(r2_fastq_path).split("R2")[0]
    
    r1_trim_path = f"{output_dir}/{fname1}trimmed_R1.fastq"
    r1_untrim_path = f"{output_dir}/{fname1}untrimmed_R1.fastq"
    r2_trim_path = f"{output_dir}/{fname2}trimmed_R2.fastq"
    r2_untrim_path = f"{output_dir}/{fname2}untrimmed_R2.fastq"
    
    subprocess.run(["TrimmomaticPE", "-threads", numthreads, "-phred33",
                    r1_fastq_path, r2_fastq_path,
                    r1_trim_path, r1_untrim_path, r2_trim_path, r2_untrim_path,
                    f"ILLUMINACLIP:{TRUSEQ}:2:30:10", "LEADING:1",
                    "TRAILING:1","SLIDINGWINDOW:4:15", "MINLEN:15"])
    
    ## Mapping
    # Align reads to the reference genome with HISAT2
    subprocess.run(["hisat2", "-p", numthreads, "-x", index_prefix,
                    "-1", r1_trim_path, "-2", r2_trim_path,
                    "-S", output_sam])
    
else:
    fastq_path = args.r1_fastq
    TRUSEQ = "TruSeq3-SE.fa"

    ## Quality Control
    subprocess.run(["fastqc", fastq_path, "-o", output_dir])

    ## Trimming
    fname = os.path.basename(fastq_path).split(".fastq")[0]
    trim_path = f"{output_dir}/{fname}_trimmed.fastq"

    subprocess.run(["TrimmomaticSE", "-threads", numthreads, "-phred33",
                    fastq_path, trim_path,
                    f"ILLUMINACLIP:{TRUSEQ}:2:30:10", "LEADING:1",
                    "TRAILING:1","SLIDINGWINDOW:4:15", "MINLEN:15"])
    
    ## Quality control for trimmed fastq
    subprocess.run(["fastqc", trim_path, "-o", output_dir])
    
    ## Mapping
    # Align reads to the reference genome with HISAT2
    subprocess.run(["hisat2", "-p", numthreads, "-x", index_prefix,
                    "-U", trim_path, "-S", output_sam])


# Convert SAM to BAM
subprocess.run(["samtools", "view", "-bS", output_sam, "-o", output_bam])
# Sort BAM file
subprocess.run(["samtools", "sort", output_bam, "-o", sorted_bam])
# Index sorted BAM file
subprocess.run(["samtools", "index", sorted_bam])

## Count output filepath
htcount_file1 = f'{output_dir}/htcount.txt'
htcount_file2 = f'{output_dir}/htcount2.txt'

# Define the commands as lists
htseq_count_cmd = ['htseq-count', '-m', 'union', '-f', 'bam', '--additional-attr=transcript_id',
                   '-s', 'yes', sorted_bam, gtf_path]
# Redirect the output to a file
with open(htcount_file1, 'w') as htcount_output:
    subprocess.run(htseq_count_cmd, stdout=htcount_output)

# Define the sed command as a list
sed_cmd = ['sed', '/^__/d']
# Run the sed command to process the htcount.txt
with open(htcount_file2, 'w') as htcount2_output:
    with open(htcount_file1, 'r') as htcount_input:
        subprocess.run(sed_cmd, stdin=htcount_input, stdout=htcount2_output)

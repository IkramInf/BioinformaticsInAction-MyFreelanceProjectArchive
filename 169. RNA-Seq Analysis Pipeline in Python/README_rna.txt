# RNA-SEQ Analysis Pipeline

This Python script performs RNA-SEQ analysis including quality control, trimming, alignment, and feature counting.

## Required Libraries

Make sure the following libraries are installed:
    - os
    - re
    - sys
    - glob
    - argparse
    - subprocess

You can install any missing libraries using `pip`.

## Usage

### Command Line Arguments
The script requires the following command line arguments:

    - -m or --mode: Mode of sequencing (`SE` for Single-End or `PE` for Paired-End). (Required)
    - -r1 or --r1_fastq: Path to the R1 FASTQ file. (Required)
    - -r2 or --r2_fastq: Path to the R2 FASTQ file (only required for Paired-End mode).
    - -r or --reference_genome: Path to the reference genome. (Required)
    - -g or --gtf: Path to the GTF file. (Required)

### Running the Script
    1. Ensure all required files and libraries are available.
    2. Open a terminal and navigate to the directory containing the script.
    3. Run the script with the necessary arguments.

#### Example Commands
For Single-End mode:
python3 script.py -m SE -r1 path/to/r1.fastq -r path/to/reference_genome.fa -g path/to/annotations.gtf

For Paired-End mode:
python3 script.py -m PE -r1 path/to/r1.fastq -r2 path/to/r2.fastq -r path/to/reference_genome.fa -g path/to/annotations.gtf

i.e., python3 variant_detection.py -m SE -r1 SRR29234182.fastq.gz -r hg38.fa -g hg38.ncbiRefSeq.gtf.gz

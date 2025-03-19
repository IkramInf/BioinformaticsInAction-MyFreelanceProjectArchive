#!/usr/bin/env python3
"""interleaved.py writes fastq files into fasta file.
    get_args():
        args.r1 : path for R1 fastq file
        args.r2 : path for R2 fastq file
        args.output : output fasta filename
        
    pathLogFile():
        return log file path
        
    interleave():
        interleave two fastq files into a list
    
    logInterleave():
        create log of all process running
"""
import os  # to create and remove directory
import argparse  # for command line argument parsing
from datetime import datetime  # for getting current timestamp
from Bio import SeqIO  # for reading/writing FASTQ/A files

def get_args():
    """Return parsed command-line arguments."""

    parser = argparse.ArgumentParser(
             description="Interleave mate-pair FASTQ sequences into a single FASTA file.",
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Get the first mate FASTQ file name (or path)
    parser.add_argument('-l', '--left',  # variable to access this data later: args.r1
                        metavar='FASTQ', # shorthand to represent the input value
                        help='Path and filename for the left or R1 FASTQ mate-pair file.',
                        type=str,
                        required=True)

    # Get the second mate FASTQ file name
    parser.add_argument('-r', '--right',  # variable to access this data later: args.r2
                        metavar='FASTQ', # shorthand to represent the input value
                        help='Path and filename for the left or R2 FASTQ mate-pair file.',
                        type=str,
                        required=True)

    # Get output FASTA file name
    parser.add_argument('-o', '--output',  # variable to access this data later: args.output
                        metavar='FASTA', # shorthand to represent the input value
                        help='Provide the path for the output FASTA file.',
                        type=str,
                        required=True)

    # extra arguments to help us format our log file output
    parser.add_argument('--logFolder',  # variable to access this data later: args.logFolder
                        help='Provide the folder for log files.', # message to the user, it goes into the help menu
                        type=str,
                        default="../results/logs/")
    parser.add_argument('--logBase',  # variable to access this data later: args.logBase
                        help='Provide the base for the log file name',
                        type=str,
                        default=parser.prog)  # get the name of the script

    return(parser.parse_args())


def pathLogFile(logFolder, logBase):
    """Return a log file path and name using the current time and script name."""
    timestamp = datetime.now().strftime("%Y-%m-%d-%H%M")  # get current time in YYYY-MM-DD-HHMM
    return(f"{logFolder}{timestamp}_{logBase}.log")


def interleave(mate1, mate2):
    """Return list of interleaved SeqRecords.

    Assumes mate1 and mate2 inputs are SeqIO.parse iterator objects.
    """
    interleaved = [(x, y) for x, y in zip(mate1, mate2)]
    return list(sum(interleaved, ()))


def logInterleave(args):
    """Create log of Interleave progress."""
    logFile = pathLogFile(args.logFolder, args.logBase)

    with open(logFile, 'w') as log:
        log.write(f"Running interleaved.py on {datetime.now()}\n")

        log.write("\n**** Summary of arguments ****")
        # log the two mate files and the output file
        log.write("\nLeft FASTQ Filename: {}".format(args.left))
        log.write("\nRight FASTQ Filename: {}".format(args.right))
        log.write("\nOutput FASTA Filename: {}".format(args.output))
        log.write("\n\n")

        #  1. Get the FASTQ sequences with SeqIO.parse
        #  2. Get the interleaved list of SeqRecord objects
        #  3. Write the interleaved list of SeqRecord objects to our FASTA file with SeqIO.write
        
        left, right = [], []
        for record in SeqIO.parse(args.left, "fastq"):
            left.append(record)
        for record in SeqIO.parse(args.right, "fastq"):
            right.append(record)
            
        interleaved = interleave(left, right)
        total_records = SeqIO.write(sequences=interleaved, handle=args.output, format="fasta")
        log.write("Total records found in {} is {}.\n".format(args.output, total_records))
        log.write(f"\nScript has finished at {datetime.now()}")


if __name__ == "__main__":
    logInterleave(get_args())  # pass arguments directly into the primary function

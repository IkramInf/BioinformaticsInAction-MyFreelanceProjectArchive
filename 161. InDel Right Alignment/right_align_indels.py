# Required Libraries
import pysam
import argparse
from Bio import SeqIO

def expand_cigar(reference, query, cigartuples):
    """
    Expands the CIGAR string representation of sequence alignment into the actual alignment strings.

    Parameters:
        reference (str): The reference sequence.
        query (str): The query sequence.
        cigartuples (list): List of tuples representing the CIGAR string.
                            Each tuple contains two elements: operation and length.
                            Operation can be one of the following:
                            - 0: Match or Mismatch
                            - 1: Insertion
                            - 2: Deletion
                            - 3: Skipped region from reference
                            - 7: Match
                            - 8: Mismatch
    Returns:
        tuple: A tuple containing two strings: (new_ref, new_query).
               new_ref: The aligned reference sequence after expansion.
               new_query: The aligned query sequence after expansion.
    """
    new_ref = ""
    new_query = ""
    ref_pos = 0
    query_pos = 0

    for operation, length in cigartuples:
        if operation == 0:  # Match or Mismatch
            ref_pos += length
            query_pos += length
            new_ref += reference[ref_pos - length : ref_pos]
            new_query += query[query_pos - length : query_pos]

        elif operation == 1:  # Insertion
            query_pos += length
            new_ref += "-" * length
            new_query += query[query_pos - length : query_pos]
            
        elif operation == 2:  # Deletion
            ref_pos += length
            new_ref += reference[ref_pos - length : ref_pos]
            new_query += "-" * length
            
        elif operation == 3:  # Skipped region from reference
            ref_pos += length

        elif operation == 7:  # Match
            ref_pos += length
            query_pos += length
            new_ref += reference[ref_pos - length : ref_pos]
            new_query += query[query_pos - length : query_pos]
            
        elif operation == 8:  # Mismatch
            ref_pos += length
            query_pos += length
            new_ref += reference[ref_pos - length : ref_pos]
            new_query += query[query_pos - length : query_pos]
            
        else:
            pass

    return (new_ref, new_query)

def do_right_alignment_of_indels(reference, query, cigartuples):
    """
    Performs right alignment of indels in a sequence alignment.

    Parameters:
        reference (str): The reference sequence.
        query (str): The query sequence.
        cigartuples (list): List of tuples representing the CIGAR string.

    Returns:
        tuple: A tuple containing two strings: (reference, query) after right alignment of indels.
    """
    ref_pos, query_pos = 0, 0

    for operation, length in cigartuples:
        if operation in [0, 7, 8]:
            ref_pos += length
            query_pos += length

        elif operation == 1:  # Insertion
            pos1, pos2 = ref_pos, query_pos
            ref_pos += length
            query_pos += length
            while True:
                ref1 = reference[pos1 : pos1 + length]
                query1 = query[pos2 : pos2 + length]
                ref2 = reference[pos1 + length : pos1 + 2*length]
                query2 = query[pos2 + length : pos2 + 2*length]
                if ref2 == query1:
                    reference = f"{reference[:pos1]}{ref2}{ref1}{reference[pos1 + 2*length:]}"
                    query = f"{query[:pos2]}{query2}{query1}{query[pos2 + 2*length:]}"
                    pos1 += length
                    pos2 += length
                else:
                    break

        elif operation == 2:  # Deletion
            pos1, pos2 = ref_pos, query_pos
            ref_pos += length
            query_pos += length
            while True:
                ref1 = reference[pos1 : pos1 + length]
                query1 = query[pos2 : pos2 + length]
                ref2 = reference[pos1 + length : pos1 + 2*length]
                query2 = query[pos2 + length : pos2 + 2*length]
                if ref1 == query2:
                    reference = f"{reference[:pos1]}{ref2}{ref1}{reference[pos1 + 2*length:]}"
                    query = f"{query[:pos2]}{query2}{query1}{query[pos2 + 2*length:]}"
                    pos1 += length
                    pos2 += length
                else:
                    break
        else:
            pass

    return (reference, query)

def generate_cigar(seq1, seq2):
    """
    Generates the CIGAR string representation of sequence alignment between two sequences.

    Parameters:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        tuple: A tuple containing two elements:
               - cigarstring (str): The CIGAR string representing the alignment.
               - cigartuples (list): List of tuples representing the CIGAR string.
    """
    cigarstring = ""
    cigartuples = []
    current_operation = None
    current_count = 0

    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            operation = 7
        elif seq1[i] == "-":
            operation = 1
        elif seq2[i] == "-":
            operation = 2
        else:
            operation = 8  # Mismatch (or S in some tools)

        if operation == current_operation:
            current_count += 1
        else:
            if current_operation:
                cigartuples.append((current_operation, current_count))
            current_operation = operation
            current_count = 1

    # Add the last operation and count
    if current_operation:
        cigartuples.append((current_operation, current_count))
    
    for operation, length in cigartuples:
        if operation == 1:
            cigarstring += f"{length}I"
        elif operation == 2:
            cigarstring += f"{length}D"
        elif operation == 7:
            cigarstring += f"{length}="
        else:
            cigarstring += f"{length}X"

    return (cigarstring, cigartuples)

if __name__ == "__main__":
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description="Perform RIGHT alignment of indels for BAM file")
    parser.add_argument("-i", "--input_bam", help="Path to aligned BAM file")
    parser.add_argument("-r", "--ref_fasta", help="Path to reference FASTA file")
    parser.add_argument("-o", "--output_bam", help="Path to new BAM file")
    args = parser.parse_args()

    #ref_fasta = pysam.FastaFile(args.ref_fasta)
    ref_fasta = {record.id : str(record.seq) for record in SeqIO.parse(args.ref_fasta, "fasta")}
    #print(ref_fasta)
    input_bam = pysam.AlignmentFile(args.input_bam, "rb")
    # open a output bam file to write
    output_bam = pysam.AlignmentFile(args.output_bam, "wb", header=input_bam.header)

    # Iterate through reads
    for read in input_bam:
        #reference = read.get_reference_sequence()
        #print(read.reference_name)
        try:
            reference = ref_fasta[read.reference_name]
            query = read.query_sequence
            cigartuples = read.cigartuples
            new_ref, new_cigar = expand_cigar(reference, query, cigartuples)
            right_aligned_ref, right_aligned_query = do_right_alignment_of_indels(new_ref, new_cigar, cigartuples)
            new_cigarstring, new_cigartuples = generate_cigar(right_aligned_ref, right_aligned_query)
            read.cigarstring = new_cigarstring
            read.cigartuples = new_cigartuples
            output_bam.write(read)
        except KeyError:
            output_bam.write(read)
            continue

    input_bam.close()
    output_bam.close()

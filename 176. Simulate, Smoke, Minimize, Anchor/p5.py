import re
import sys

def print_alignment(seq1, seq2, cigar):
    """
    Aligns two sequences based on a CIGAR string and returns the formatted alignment.

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        cigar (str): CIGAR string describing the alignment (e.g., "5M1I3M1D2M").

    Returns:
        tuple: A tuple containing three strings:
            - alignment1: Aligned version of the first sequence.
            - symbols: Symbols indicating matches ('|'), mismatches ('*'), and gaps.
            - alignment2: Aligned version of the second sequence.
    """
    
    alignment1 = ""  # String to store aligned first sequence
    symbols = ""     # String to store matching symbols (| for match, * for mismatch/gap)
    alignment2 = ""  # String to store aligned second sequence
    n1 = 0  # Pointer for seq1
    n2 = 0  # Pointer for seq2
    
    # Parse the CIGAR string into tuples of (length, operation)
    cigar_tuples = re.findall(r'(\d+)([MX=DI])', cigar)

    # Iterate through each (length, operation) pair in the CIGAR string
    for length, operator in cigar_tuples:
        length = int(length)  # Convert the length from string to integer

        # Match operation (M): Check for matches/mismatches and append to alignments
        if operator == "M":
            for i in range(1, length+1):
                symbol = "|" if seq1[n1] == seq2[n2] else "*"
                alignment1 += seq1[n1]
                symbols += symbol
                alignment2 += seq2[n2]
                n1 += 1
                n2 += 1
        
        # Exact match operation (=): Perfect match, use '|' for all positions
        elif operator == "=":
            alignment1 += seq1[n1 : n1 + length]
            symbols += "|" * length
            alignment2 += seq2[n2 : n2 + length]
            n1 += length
            n2 += length
        
        # Mismatch operation (X): Mark mismatches with '*'
        elif operator == "X":
            alignment1 += seq1[n1 : n1 + length]
            symbols += "*" * length
            alignment2 += seq2[n2 : n2 + length]
            n1 += length
            n2 += length

        # Deletion operation (D): Sequence 1 has a gap ('-'), seq2 has characters
        elif operator == "D":
            alignment1 += "-" * length
            symbols += "*" * length  # Mark gaps with '*'
            alignment2 += seq2[n2 : n2 + length]
            n2 += length

        # Insertion operation (I): Sequence 2 has a gap ('-'), seq1 has characters
        elif operator == "I":
            alignment1 += seq1[n1 : n1 + length]
            symbols += "*" * length  # Mark gaps with '*'
            alignment2 += "-" * length
            n1 += length
        
        # Invalid operator: Print error message and exit
        else:
            print("Exiting!!! Invalid Operator.")
            sys.exit()

        # Add space between each block in the alignment for readability
        alignment1 += " "
        symbols += " "
        alignment2 += " "

    return (alignment1, symbols, alignment2)

if __name__ == "__main__":
    # Ensure correct usage of command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python p5.py <sequence 1> <sequence 2> <cigar>")
        sys.exit()

    # Read command-line arguments
    seq1 = sys.argv[1]  # First sequence
    seq2 = sys.argv[2]  # Second sequence
    cigar = sys.argv[3]  # CIGAR string
    
    # Perform the alignment and capture the result
    alignment1, symbols, alignment2 = print_alignment(seq1, seq2, cigar)
    
    # Print the formatted alignment
    print(f"{alignment1}\n{symbols}\n{alignment2}")

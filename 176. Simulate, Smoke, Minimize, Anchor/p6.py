import sys
import re

# Function to get the reverse complement of a k-mer
def reverse_complement(seq):
    """
    Returns the reverse complement of the given nucleotide sequence.

    Args:
        seq (str): A nucleotide sequence (e.g., "ATCG").

    Returns:
        str: The reverse complement of the input sequence (e.g., "CGAT" for "ATCG").
    """
    # Translate nucleotides and reverse the sequence
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

# Function to find the k-minimizer, canonical k-minimizer, and ORCOM k-minimizer
def minimize(read, k):
    """
    Finds the regular k-minimizer, canonical k-minimizer, and ORCOM k-minimizer from a nucleotide read.

    Args:
        read (str): A nucleotide sequence read (e.g., from a FASTQ file).
        k (int): Length of the k-mer.

    Returns:
        tuple: A tuple containing:
            - k_minimizer (str): Lexicographically smallest k-mer from the read.
            - canonical_k_minimizer (str): Lexicographically smallest k-mer or its reverse complement.
            - orcom_k_minimizer (str): ORCOM k-minimizer (smallest k-mer/reverse complement with no 3+ identical nucleotides in a row).
    """
    # Find all k-mers in the read
    kmers = [read[i:i + k] for i in range(len(read) - k + 1)]

    # Regular k-minimizer (lexicographically smallest k-mer)
    k_minimizer = min(kmers)

    # Canonical k-minimizer (smallest k-mer or its reverse complement)
    canonical_kmers = [min(kmer, reverse_complement(kmer)) for kmer in kmers]
    canonical_k_minimizer = min(canonical_kmers)

    # ORCOM k-minimizer (smallest k-mer or reverse complement with no more than two identical nucleotides in a row)
    pattern = r"(A{3,}|C{3,}|G{3,}|T{3,})"  # Regular expression to find k-mers with 3 or more identical bases in a row
    orcom_kmers = [min(kmer, reverse_complement(kmer)) for kmer in kmers if not re.search(pattern, kmer)]
    
    # If no ORCOM minimizers found (all k-mers have repetitive nucleotides), return k-dashes
    if orcom_kmers:
        orcom_k_minimizer = min(orcom_kmers)
    else:
        orcom_k_minimizer = '-' * k  # No valid ORCOM k-mer found

    return k_minimizer, canonical_k_minimizer, orcom_k_minimizer


if __name__ == "__main__":
    """
    This script reads a FASTQ file, extracts the nucleotide sequence from each read,
    and computes the k-minimizer, canonical k-minimizer, and ORCOM k-minimizer.
    
    Usage:
        python script.py <fastq_file> <k>
    
    Arguments:
        fastq_file (str): Path to the FASTQ file.
        k (int): Length of the k-mers to find minimizers.
    """
    # Get command-line arguments: FASTQ file and k-mer length
    fastq_file = sys.argv[1]  # Path to FASTQ file
    k = int(sys.argv[2])      # Length of k-mers

    # Process the FASTQ file
    with open(fastq_file, 'r') as fastq_reader:
        while True:
            # Read 4 lines for each record in the FASTQ file
            read_id = fastq_reader.readline().strip()  # Line 1: Sequence identifier
            if not read_id.startswith("@"):
                break  # End of file or invalid FASTQ format

            sequence = fastq_reader.readline().strip()  # Line 2: Nucleotide sequence
            separator = fastq_reader.readline().strip()  # Line 3: Separator line (often '+')
            quality = fastq_reader.readline().strip()  # Line 4: Quality scores in ASCII format

            # Strip '@' from the read name to get the unique identifier
            read_id = read_id.split()[0][1:]

            # Find the minimizers for the current sequence
            k_minimizer, canonical_k_minimizer, orcom_k_minimizer = minimize(sequence, k)

            # Print the read ID along with its k-minimizers
            print(f"{read_id}\t{k_minimizer}\t{canonical_k_minimizer}\t{orcom_k_minimizer}")

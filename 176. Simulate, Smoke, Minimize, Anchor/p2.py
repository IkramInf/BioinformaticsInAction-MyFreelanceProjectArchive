import sys
from Bio import SeqIO

def simulate(reference_genome, sequence_name, l, t, d, output_file="sim.sam"):
    """
    Simulates paired-end reads from a reference genome and outputs them in SAM format.

    Args:
        reference_genome (str): Path to the reference genome in FASTA format.
        sequence_name (str): The name of the chromosome or sequence to simulate reads from.
        l (int): Length of each read.
        t (int): Template length (distance between forward and reverse reads).
        d (int): Depth of coverage (number of reads per base).
        output_file (str): The SAM output file to write. Defaults to 'sim.sam'.
    """
    
    # Read the reference genome and store the sequences in a dictionary
    sequences = {record.id : str(record.seq) for record in SeqIO.parse(reference_genome, "fasta")}
    sequence = sequences[sequence_name]  # Extract the specific sequence for the chromosome

    # Calculate the number of fragments to generate based on the coverage, template, and sequence length
    # num_fragments = int(d * len(sequence) / (2 * l))  # Old formula
    num_fragments = int(((len(sequence) - t + 1) * d) / t)  # New formula considering fragment overlap

    # Open the output file for writing
    with open(output_file, "w") as f:
        # Write the SAM header
        f.write("@HD\tVN:1.6\tSO:coordinate\n")  # HD: Header with version and sorting info
        for chrom, seq in sequences.items():
            f.write(f"@SQ\tSN:{chrom}\tLN:{len(seq)}\n")  # SQ: Sequence dictionary entry
        
        # Define FLAG values for SAM (representing read properties)
        PAIRED = 0x1  # Paired-end read
        PROPER_PAIR = 0x2  # Properly paired
        READ_UNMAPPED = 0x4  # Read is unmapped
        MATE_UNMAPPED = 0x8  # Mate is unmapped
        READ_REVERSE_STRAND = 0x10  # Read is on reverse strand
        MATE_REVERSE_STRAND = 0x20  # Mate is on reverse strand
        FIRST_IN_PAIR = 0x40  # First read in pair
        SECOND_IN_PAIR = 0x80  # Second read in pair

        # Forward read FLAG: paired, properly paired, first in pair, mate on reverse strand
        fFLAG = PAIRED | PROPER_PAIR | FIRST_IN_PAIR | MATE_REVERSE_STRAND
        # Reverse read FLAG: paired, properly paired, second in pair, read on reverse strand
        rFLAG = PAIRED | PROPER_PAIR | SECOND_IN_PAIR | READ_REVERSE_STRAND

        # Generate paired-end reads
        for i in range(num_fragments):
            # Calculate the starting position of the fragment in the sequence
            fragment_pos = (i * t) // d

            # Read identifier (QNAME) for the read pair
            QNAME = f"read_{sequence_name}_{fragment_pos + 1}"

            # Extract the template sequence of length t from the reference sequence
            template = sequence[fragment_pos : fragment_pos + t]
            
            # If the template is shorter than the expected length (near the end), skip the fragment
            if len(template) < t:
                continue  # Avoid incomplete reads
            
            # Generate forward and reverse reads
            fSEQ = template[:l]  # First 'l' bases for forward read
            rSEQ = template[-l:]  # Last 'l' bases for reverse read

            # Calculate positions for the forward and reverse reads on the reference
            fPOS = fragment_pos + 1  # Position of forward read (1-based)
            rPOS = fragment_pos + t - l + 1  # Position of reverse read

            # SAM fields
            RNAME = sequence_name  # Reference sequence name (chromosome)
            MAPQ = 0  # Mapping quality (0 as no actual mapping is done here)
            CIGAR = f"{l}M"  # CIGAR string (l matches)
            RNEXT = "="  # Same reference for the next read
            fPNEXT = rPOS  # Position of the mate (reverse read) for forward read
            rPNEXT = fPOS  # Position of the mate (forward read) for reverse read
            TLEN = t  # Template length (distance between forward and reverse reads)
            QUAL = "*"  # Placeholder for quality scores

            # Write forward read to SAM file
            f.write(f"{QNAME}\t{fFLAG}\t{RNAME}\t{fPOS}\t{MAPQ}\t{CIGAR}\t{RNEXT}\t{fPNEXT}\t{TLEN}\t{fSEQ}\t{QUAL}\n")

            # Write reverse read to SAM file
            f.write(f"{QNAME}\t{rFLAG}\t{RNAME}\t{rPOS}\t{MAPQ}\t{CIGAR}\t{RNEXT}\t{rPNEXT}\t{TLEN}\t{rSEQ}\t{QUAL}\n")

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 6:
        print("Usage: python p2.py <reference genome> <chromosome name> <read length> <template length> <coverage>")
        sys.exit()

    # Parse command-line arguments
    reference_genome = sys.argv[1]  # Path to reference genome file, e.g., "sacCer3.fa"
    sequence_name = sys.argv[2]  # Chromosome name, e.g., "chrI"
    l = int(sys.argv[3])  # Read length, e.g., 100
    t = int(sys.argv[4])  # Template length, e.g., 500
    d = int(sys.argv[5])  # Depth of coverage, e.g., 20

    # Call the simulate function to generate paired-end reads and write them to a SAM file
    simulate(reference_genome, sequence_name, l, t, d, "sim.sam")

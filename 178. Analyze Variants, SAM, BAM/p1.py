import sys
import pysam
from Bio import SeqIO
from collections import defaultdict

# Handling command line arguments
if len(sys.argv) != 3:
    print("Usage: python3 p1.py <reference_filepath> <sam_filepath>")
    sys.exit(1)

# Assign cmd arguments into variables
reference_filepath = sys.argv[1]
sam_filepath = sys.argv[2]

# Read Reference genome into a dict
ref_genome = {record.id: str(record.seq).upper() for record in SeqIO.parse(reference_filepath, "fasta")}
# Read Sam File using pysam
samfile = pysam.AlignmentFile(sam_filepath, "r")
# Create a dict to store position wise alleles
alleles = defaultdict(lambda: defaultdict(int))

# Iterate over each read in the file
for read in samfile.fetch():
    # Get reference details
    ref_chr = read.reference_name
    ref_pos = read.reference_start
    # Get read details
    read_seq = read.query_sequence
    read_pos = 0

    # Iterate over cigar string
    for op, length in read.cigartuples:
        if op in [0, 7, 8]:  # [M, =, X]
            ref_seq = ref_genome[f"chr{ref_chr}"][ref_pos : ref_pos + length]
            query_seq = read_seq[read_pos : read_pos + length]
            # Count the observed alleles
            for ref, query, pos in zip(ref_seq, query_seq, range(ref_pos, ref_pos + length)):
                alleles[pos][query] += 1
                alleles[pos]['ref'] = (f"chr{ref_chr}", pos, ref_genome[f"chr{ref_chr}"][pos])

            # Update positions after matching sequence
            read_pos += length
            ref_pos += length
        elif op in [1, 4]:  # [I, S]
            read_pos += length
        elif op in [2, 3]:  # [D, N]
            ref_pos += length
        else:
            continue

# Close the samfile object
samfile.close()

# Analyze variants
for pos, counts in alleles.items():
    if 'ref' not in counts:
        continue
    
    # Retrieve from counts['ref']
    chrom, location, ref_base = counts['ref']
    # Total coverage excluding 'ref' tuple
    total_coverage = sum(count for base, count in counts.items() if isinstance(count, int))
    
    if total_coverage == 0:
        continue

    # Calculate frequencies for each allele
    frequencies = {base: count/total_coverage for base, count in counts.items() if isinstance(count, int)}
    # Sort alleles by frequency
    sorted_alleles = sorted(frequencies.items(), key=lambda x: x[1], reverse=True)

    # Check for non-reference homozygous variant (>80% coverage)
    if sorted_alleles[0][1] >= 0.8:
        frequent_allele = sorted_alleles[0][0]
        if frequent_allele != ref_base:
            print(f"{chrom} {location + 1} HOM {frequent_allele}/{frequent_allele}")

    # Check for heterozygous variant (>40% coverage of both alleles)
    elif sorted_alleles[0][1] >= 0.4 and sorted_alleles[1][1] >= 0.4:
        allele1, allele2 = sorted([sorted_alleles[0][0], sorted_alleles[1][0]])
        print(f"{chrom} {location + 1} HET {allele1}/{allele2}")


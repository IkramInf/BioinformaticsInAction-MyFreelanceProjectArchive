
def valid_sequence(sequence):
    """Verifies that a nucleotide sequence contains only DNA letters.
    
    Valid characters include A, C, G, T and their lower-case variants.
    
    Args:
        sequence (str): Nucleotide sequence.
            
    Returns:
        bool: `True` if `sequence` contains only valid characters,
            `False` otherwise.
    """
    valid_letters = ['A', 'C', 'G', 'T']
    
    for i in sequence.upper():
        if not i in valid_letters:
            return False
            
    return True


def find_longest_orf(sequence):
    """Reports longest open reading frame in a DNA sequence.

    Returns the start position (0-based) and length (in nucleotides) of
    the longest open reading frame (ORF) in the input sequence. Each
    ORF is delimited by a start (ATG) and a stop (TAA, TAG, TGA)
    codon. If there are multiple equivalent solutions, the one with the
    smallest start position/index is reported.

    Args:
        sequence (str): DNA sequence in lower and/or uppercase.
    
    Returns:
        tuple[int, int]: Start position/index and length of the longest
            open reading frame. If not ORF was found, return -1 for
            each.

    Raises:
        ValueError: `sequence` is not a valid DNA sequence.
    """
    
    stop_codons = ['TAA', 'TAG', 'TGA']
    orf = {}
    if valid_sequence(sequence):
        seq = sequence.upper()
        for i in range(len(seq)):
            if seq[i:i+3] == "ATG":
                for j in range(i, len(seq)):
                    if seq[j:j+3] in stop_codons:
                        orf.update({i : j-i+3})
        
    else:
        raise ValueError("The sequence is not a valid DNA sequence")
    
    if orf:
        start, length = max(orf.items(), key=lambda x:x[1])
        return [start, length]
    else:
        return [-1, -1]


def translate_dna(sequence, start, length):
    """Translates DNA to the corresponding peptide sequence.

    Translates the DNA sequence starting at a specified position 
    (0-based) and extending a specified number of nucleotides 
    (length) into the corresponding peptide sequence, given the 
    universal genetic code.

    The genetic code is implemented as a dictionary 
    {“codon” : “AA”} (AA = amino acid in 1-letter code) and 
    includes entries for the stop codons, which are translated 
    to “*”, a non-amino acid character. 

    Args:
        sequence (str): DNA sequence in lower and/or uppercase.
        start (int): Sequence position/index from which to start
            nucleotide to amino acid translation.
        length (int): Number of nucleotides to translate.

    Returns:
        str: Peptide sequence corresponding to the indicated DNA
            sequence region. If the DNA sequence is an empty string,
            return an empty string.

        Raises:
            ValueError: `sequence` is not a valid DNA sequence.
            ValueError: `start` is out of bounds.
            ValueError: `length` is not a multiple of 3.
            ValueError: The sum of `start` and `length` - 1 goes
                beyond the end of the sequence.
    """
    
    CodonTable = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    
    if valid_sequence(sequence):
        seq = sequence[start:start+length]
        protein = "".join([CodonTable[seq[i:i + 3]] for i in range(0, len(seq), 3)])
        return protein
        
    else:
        raise ValueError("The sequence is not a valid DNA sequence")


def predict_orf_seq(sequence):
    """Predicts peptide sequence for longest ORF in DNA sequence.

    Finds the longest ORF within a DNA sequence and, if available,
    translates it to the peptide sequence (1-letter code) it encodes.
    If multiple ORFs of the same length are found, only the first
    (smaller start position) is reported.

    Args:
        sequence (str): DNA sequence in lower and/or uppercase.

    Returns:
        tuple[int, int, str]: Tuple of start position, length and
            peptide sequence. If no ORF was found, start position
            and length return values are -1. If `sequence` is an
            empty string, return an empty string for the peptide
            sequence.
    """
    
    if sequence:
        start, length = find_longest_orf(sequence)

        if start != -1 and length != -1:
            protein = translate_dna(sequence, start, length)
            return [start, length, protein]
  
    return ""


def main(fasta):
    """Predicts peptide sequences for open reading frames in input file.
    
    Given a FASTA file of DNA sequences, for each sequence, the longest
    open reading frame (ORF) is identified and the corresponding
    peptide sequence predicted.
    
    Results are printed to the screen in the following comma-separated
    format:
    
    sequence_name,index_start_orf,length_orf,peptide_sequence
    
    where
    
    sequence_name: is taken from the identifier line of each record in
        the input FASTA file
    index_start_orf: defines the 0-based index where the longest ORF
        was found
    length_orf: defines the length of the longest ORF in nucleotides
    peptide_sequence: is the peptide sequence (1-letter code) predicted
        to be translated from that ORF according to the universal
        genetic code

    An example output record might look like this:
    
    ABC1,0,18,MSEYQP*
    
    Args:
        fasta(str): Path to a FASTA file of DNA sequences.
    """
    
    with open(fasta, "r") as f:
        lines = f.readlines()
        index = [d for d in range(len(lines)) if lines[d].startswith(">")]
        titles = [t.replace("\n", "") for t in lines if t.startswith(">")]
        sequences = []
        for i in range(len(index)):
            if index[i] == index[-1]:
                sequences.append("".join(lines[index[i]+1:]).replace("\n", ""))
            else:
                sequences.append("".join(lines[index[i]+1:index[i+1]]).replace("\n", ""))
            
    for title, seq in zip(titles, sequences):
        result = predict_orf_seq(seq)
        if result:
            start, length, protein = result
            print(f"{title.split()[0][1:]}, {start}, {length}, {protein}")


# calling the main function to test on fasta file
main("seqs.fa")


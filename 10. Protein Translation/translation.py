# Your code here
def translate(sequence):
    # dictionary of codon : amino acid pair
    table = {
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
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }
    
    # create a dictionary to store proteins
    proteins = {}
    
    # for iterating over forward sequence
    for frame in range(3):
        length = 3 * ((len(sequence)-frame) // 3) # Multiple of three
        seq = sequence[frame:frame+length]
        proteins["f"+str(frame+1)] = "".join([table[seq[i:i + 3]] for i in range(0, len(seq), 3)])
    
    # reverse the sequence
    rev_seq = sequence[::-1]
    
    # for iterating over reverse sequence
    for frame in range(3):
        length = 3 * ((len(rev_seq)-frame) // 3) # Multiple of three
        rev_seq_framed = rev_seq[frame:frame+length]
        proteins["r"+str(frame+1)] = "".join([table[rev_seq_framed[j:j + 3]] \
                                              for j in range(0, len(rev_seq_framed), 3)])
        
    return proteins


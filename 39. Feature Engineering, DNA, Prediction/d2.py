import csv
import d1


def frequencyVector(seq, seq_type):
    
    dna_letters = {'A','C','G','T'}
    rna_letters = {'A','C','G','U'}
    protein_letters = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}

    # initialize dictionary with 0
    dna = {b:0 for b in dna_letters}
    rna = {b:0 for b in rna_letters}
    protein = {aa:0 for aa in protein_letters}

    if seq_type == 1:
        for base in seq:
            dna[base] += 1
        return list(dict(sorted(dna.items())).values())
    elif seq_type == 2:
        for base in seq:
            rna[base] += 1
        return list(dict(sorted(rna.items())).values())
    elif seq_type == 3:
        for base in seq:
            protein[base] += 1
        return list(dict(sorted(protein.items())).values())
    else:
        print("Invalid Sequence Type!!!")


def compositionVector(seq, seq_type):
    
    dna_letters = {'A','C','G','T'}
    rna_letters = {'A','C','G','U'}
    protein_letters = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}

    # initialize dictionary with 0
    dna = {b:0 for b in dna_letters}
    rna = {b:0 for b in rna_letters}
    protein = {aa:0 for aa in protein_letters}

    if seq_type == 1:
        for base in seq:
            dna[base] += 1
        comp_vec = list(dict(sorted(dna.items())).values())
        comp_vec = [round(v/len(seq), 3) for v in comp_vec]
        return comp_vec
    elif seq_type == 2:
        for base in seq:
            rna[base] += 1
        comp_vec = list(dict(sorted(rna.items())).values())
        comp_vec = [round(v/len(seq), 3) for v in comp_vec]
        return comp_vec
    elif seq_type == 3:
        for base in seq:
            protein[base] += 1
        comp_vec = list(dict(sorted(protein.items())).values())
        comp_vec = [round(v/len(seq), 3) for v in comp_vec]
        return comp_vec
    else:
        print("Invalid Sequence Type!!!")


def absolutePositioningVector(seq, seq_type):
    
    dna_letters = {'A','C','G','T'}
    rna_letters = {'A','C','G','U'}
    protein_letters = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}

    # initialize dictionary with 0
    dna = {b:0 for b in dna_letters}
    rna = {b:0 for b in rna_letters}
    protein = {aa:0 for aa in protein_letters}

    if seq_type == 1:
        for i, base in enumerate(seq):
            dna[base] += i+1
        return list(dict(sorted(dna.items())).values())
    elif seq_type == 2:
        for i, base in enumerate(seq):
            rna[base] += i+1
        return list(dict(sorted(rna.items())).values())
    elif seq_type == 3:
        for i, base in enumerate(seq):
            protein[base] += i+1
        return list(dict(sorted(protein.items())).values())
    else:
        print("Invalid Sequence Type!!!")


def relativePositioningVector(seq, seq_type):
    
    dna_letters = {'A','C','G','T'}
    rna_letters = {'A','C','G','U'}
    protein_letters = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}

    # initialize dictionary with 0
    dna = {i : {j:0 for j in dna_letters} for i in dna_letters}
    rna = {i : {j:0 for j in rna_letters} for i in rna_letters}
    protein = {i : {j:0 for j in protein_letters} for i in protein_letters}
    
    rev_vec = []
    
    if seq_type == 1:
        for i, occur in enumerate(dna_letters):
            if occur in seq:
                ind = seq.index(occur)
                for j, base in enumerate(dna_letters):
                    if base in seq and i!=j:
                        indices = [k-ind for k, x in enumerate(seq) if x == base]
                        dna[occur][base] =  sum(indices)
        
        dna = dict(sorted(dna.items(), key=lambda x:x[0].lower()))
        for value in dna.values():
            rev_vec.extend(dict(sorted(value.items(), key=lambda x:x[0].lower())).values())
            
    elif seq_type == 2:
        for i, occur in enumerate(rna_letters):
            if occur in seq:
                ind = seq.index(occur)
                for j, base in enumerate(rna_letters):
                    if base in seq and i!=j:
                        indices = [k-ind for k, x in enumerate(seq) if x == base]
                        rna[occur][base] =  sum(indices)
        
        rna = dict(sorted(rna.items(), key=lambda x:x[0].lower()))
        for value in rna.values():
            rev_vec.extend(dict(sorted(value.items(), key=lambda x:x[0].lower())).values())

    elif seq_type == 3:
        for i, occur in enumerate(protein_letters):
            if occur in seq:
                ind = seq.index(occur)
                for j, base in enumerate(protein_letters):
                    if base in seq and i!=j:
                        indices = [k-ind for k, x in enumerate(seq) if x == base]
                        protein[occur][base] =  sum(indices)
        
        protein = dict(sorted(protein.items(), key=lambda x:x[0].lower()))
        for value in protein.values():
            rev_vec.extend(dict(sorted(value.items(), key=lambda x:x[0].lower())).values())
            
    else:
        print("Invalid Sequence Type!!!")
        
    return rev_vec


def calculate(filename, seq_type):
    
    # read data with parser of project 1
    data = d1.parser(filename)
    
    outputfilename = "".join(filename.split(".")[:-1]).split("/")[-1] + ".csv"
    print(outputfilename)

    with open(outputfilename, "w") as csvfile:
        writer = csv.writer(csvfile)

        for k, seq in data.items():
            row = frequencyVector(seq, seq_type) + compositionVector(seq, seq_type) + \
                absolutePositioningVector(seq, seq_type) + relativePositioningVector(seq, seq_type)
            writer.writerow(row)
            








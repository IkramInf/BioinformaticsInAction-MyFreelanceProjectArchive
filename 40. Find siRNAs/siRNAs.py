import re
from Bio import SeqIO


ib = {record.id : str(record.seq) for record in SeqIO.parse("ib.fasta", "fasta")}


def isAsymmetric(seq):
    l = len(seq)
    if l % 2 == 0:
        if seq[0:int(l/2)] == seq[int(l/2):] or seq == seq[::-1]:
            return False
    return True


def checkGC(seq):
    seq = seq.upper()
    gc = float(seq.count('G')+seq.count('C'))/len(seq)*100
    if gc >= 20 and gc <= 50:
        return True
    else:
        return False


def notConsecutive(seq):
    seq = seq.upper()
    matches = re.findall("A{4,}|C{4,}|G{4,}|T{4,}", seq)
    if matches:
        return False
    else:
        return True


def findsiRNAs(sequences, length=21):
    #siRNAs = {}
    with open("updated1.fasta", "w") as f:
        for Id, seq in sequences.items():
            for i in range(len(seq)-length+1):
                kmer = seq[i:i+length]
                if isAsymmetric(kmer) and checkGC(kmer) and notConsecutive(kmer):
                    #siRNAs[str(Id)+'.'+str(i)] = kmer
                    f.write(f">{str(Id)+'.'+str(i)}\n{kmer}\n")
                
    #return siRNAs


# calling the function                            
findsiRNAs(ib)


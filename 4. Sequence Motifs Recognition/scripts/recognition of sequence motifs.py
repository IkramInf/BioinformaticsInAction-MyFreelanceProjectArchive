from Bio import Seq, SeqIO
from Bio import motifs
from Bio.SeqRecord import SeqRecord

record = SeqIO.read('Human_Chr19_all.gb', "gb") # Read chromosome 19
record_rev_seq = record.seq.reverse_complement() # reverse complement of the seq

# to save the positions of exon end and exon start
exon_start = []
exon_end = []

for feature in record.features: #for testing only the first features 
    if feature.type =='mRNA':
        for part in feature.location.parts:
            #the mRNA consists of several exons whose coordinates in location. parts are stored
            if part.strand == 1: # this is only the forward beach!
                exon_start.append(part.start.position)
                exon_end.append(part.end.position)
                
                if part.start.position > part.end.position:
                    raise Exception('Incorrect coordinates in {}'.format(feature)) #this should never happen


# to save the positions of exon end and exon start
exon_start_rev = []
exon_end_rev = []

for feature in record.features: #for testing only the first features 
    if feature.type =='mRNA':
        for part in feature.location.parts:
            #the mRNA consists of several exons whose coordinates in location. parts are stored
            if part.strand == -1: # this is only the forward beach!
                exon_start_rev.append(part.start.position)
                exon_end_rev.append(part.end.position)
                
                if part.start.position > part.end.position:
                    raise Exception('Incorrect coordinates in {}'.format(feature)) #this should never happen

# remove the first position from the exon_start and exon_start_rev list
exon_start = exon_start[1:]
exon_start_rev = exon_start_rev[1:]

# to store donor sites and acceptor sites
donorsites = []
acceptorsites = []

# to ensure that the sequences obtained are as expected: GT--Intron- AG
# with GT at the beginning of an intron (donor site) and AG at the end of an intron (acceptor site)
for ee, es in zip(exon_end, exon_start):
    #remember, [i:j] indicates, start point (i) is inclusive and end point (j) is exclusive in python indexing
    # ee was exclusive and es was inclusive while adding to exon_end and exon_start lists
    if record.seq[ee:ee+2] == "GT" and record.seq[es-2:es] == "AG":
        # Last 7 nucleotides of the exon and the first 7 nucleotides of the intron
        donorsites.append(record.seq[ee-7 : ee+7])
        
        # Last 7 nucleotides of the intron and the first 7 nucleotides of the exon
        acceptorsites.append(record.seq[es-7 : es+7])

# to store donor sites and acceptor sites
donorsites_rev = []
acceptorsites_rev = []

# to ensure that the sequences obtained are as expected: GT--Intron- AG
# with GT at the beginning of an intron (donor site) and AG at the end of an intron (acceptor site)
for eer, esr in zip(exon_end_rev, exon_start_rev):
    #remember, [i:j] indicates, start point (i) is inclusive and end point (j) is exclusive in python indexing
    # ee was exclusive and es was inclusive while adding to exon_end and exon_start lists
    if record_rev_seq[eer:eer+2] == "GT" and record_rev_seq[esr-2:esr] == "AG":
        # Last 7 nucleotides of the exon and the first 7 nucleotides of the intron
        donorsites_rev.append(record_rev_seq[eer-7 : eer+7])
        
        # Last 7 nucleotides of the intron and the first 7 nucleotides of the exon
        acceptorsites_rev.append(record_rev_seq[esr-7 : esr+7])

# accumulating donor sites of forward strand and reverse strand and removing duplicates
donorsites = list(set(donorsites + donorsites_rev))

# accumulating acceptor sites of forward strand and reverse strand and removing duplicates
acceptorsites = list(set(acceptorsites + acceptorsites_rev))

# donorsites: List of sequences of the motif
donor = motifs.create(donorsites) #create motif
print(donor.consensus) #consensus sequence
print(donor.pwm) #relative frequencies
print(len(donor.instances)) #number of sequences for motif creation

#donor sites write to mfasta file
from Bio.SeqRecord import SeqRecord
records=[SeqRecord(e, id='donor', description='') for e in donor.instances]
SeqIO.write(records, 'donor.mfasta', 'fasta')

# acceptorsites: List of sequences of the motif
acceptor = motifs.create(acceptorsites) #create motif
print(acceptor.consensus) #consensus sequence
print(acceptor.pwm) #relative frequencies
print(len(acceptor.instances)) #number of sequences for motif creation

#acceptor sites write to mfasta file
records=[SeqRecord(e, id='acceptor', description='') for e in acceptor.instances]
SeqIO.write(records, 'acceptor.mfasta', 'fasta')

# creating weblogo of donorsites and acceptorsites
donor.weblogo('donorsites.png', format='PNG', size='large', logo_title= 'Donor Site (Exon/Intron)')
acceptor.weblogo('acceptorsites.png', format='PNG', size='large', logo_title= 'Acceptor Site (Intron/Exon)')


#---------------------------------------------------------------------------------
#BIOPYTHON EXERCISE - Biol568 BMI students
#---------------------------------------------------------------------------------
#
# You will need to install Biopython and dependencies, or use a computer
#   where they are installed.  https://biopython.org
#
# To write the functions, you can learn how to use Biopython using the
#   Biopython Tutorial http://biopython.org/DIST/docs/tutorial/Tutorial.html
#   I recommend trying out Biopython code on the command line first.

# To use Biopython, you will need to import the Biopython libraries needed.
# Follow the Tutorial. You will certainly need to read about how Biopython
# stores sequence data.
#
#For example:
#
#from Bio import SeqIO

#
#PART 1 - You will need a Fasta file. One with MULTIPLE SEQUENCE RECORDS.
#   You can use the one in the tutorial (ls_orchid.fasta). 
#   See section 2.4.1  Simple FASTA parsing example

#FUNCTION NAME: ParseFastaFile
#PARAMETERS: 1 (File Name -a string of the name of a fasta sequence file)
#PURPOSE: The function should:
#           (1) Use Biopython to parse a fasta file with multiple sequences.
#           (2) Make a dictionary: The keys will be the gi number 
#               and the values the sequences.
#
# The gi number is the first number of the title.
# Example title sequence:
# >gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA
#
# For this entry you will want '2765657' as the key.
# This is the dictionary entry {'2765657':'CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAGAATATATGATCGAGTG...'}
#
#RETURN VALUES: The dictionary with ALL the entries.

from Bio import SeqIO
#== FUNCTION 1 ==
def ParseFastaFile(filename):
    parsed_fasta = {r.id.split("|")[1] : str(r.seq) for r in SeqIO.parse("ls_orchid.fasta", "fasta")}
    return parsed_fasta

# Example function call
#print(ParseFastaFile("ls_orchid.fasta"))

#---------------------------------------------------------------------------------
#PART 2 - You will need a Genbank file. One with MULTIPLE SEQUENCE RECORDS.
#   You can use the one in the tutorial (ls_orchid.gbk). 
#   See section 2.4.2  Simple GenBank parsing example
#   and section 4.2.3  SeqRecord objects from GenBank files


#FUNCTION NAME: GetGenbankInfo
#PARAMETERS: 1 (File Name -a string of the name of a GenBank sequence file)
#PURPOSE: The function should:
#           (1) Use Biopython to parse the GenBank file with multiple sequences.
#           (2) Make a dictionary -
#                   (i)  The key will be the record.id
#                   (ii) The value will be a LIST with two items:
#                       1. the description
#                       2. the sequence

# For instance, from this record (only top part shown):
"""
LOCUS       NC_005816               9609 bp    DNA     circular BCT 21-JUL-2008
DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete
            sequence.
ACCESSION   NC_005816
VERSION     NC_005816.1  GI:45478711
PROJECT     GenomeProject:10638
"""
# The dictionary entry should be:
# {'NC_005816.1': [''Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence.', ''TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG']
# 
#
#RETURN VALUES: The dictionary with ALL the entries.
#

from Bio import SeqIO
#== FUNCTION 2 ==
def ParseGenbankFile(filename):
    parsed_gb = {}
    for r in SeqIO.parse(filename, "genbank"):
        parsed_gb.setdefault(r.id, []).append(r.description)
        parsed_gb.setdefault(r.id, []).append(str(r.seq))
    
    return parsed_gb

# Example function call
#print(ParseGenbankFile("ls_orchid.gbk"))

#---------------------------------------------------------------------------------
#PART 3 -   See sectons 3.8 Transcription AND 3.9  Translation
#   Note: Assume the sequence is a CODING STRAND (don't reverse complement).

#FUNCTION NAME: TranscribeTranslate
#PARAMETERS: 1 (A list of sequence nmers)
#PURPOSE: The function should:
#           (1) Iterate through the LIST of DNA sequence nmers. You make them up.
#           (2) Using Biopython, it should transcribe and then translate each sequence.
#           (3) These amino acid sequences should be stored in a list.
  
#RETURN VALUES: A LIST of all the amino acid sequences (strings)

from Bio.Seq import Seq
#== FUNCTION 3 ==
def TranscribeTranslate(seq_lists):
    proteins = []
    for seq in seq_lists:
        seq = Seq(seq)
        rna = seq.transcribe()
        protein = rna.translate()
        proteins.append(str(protein))
    return proteins

# Example function call
#print(TranscribeTranslate(['CGTAACAAGGTT', 'TCCGTAGGTGAACCT', 'TGCGGAAGG', 'ATCATTGATGAGACCGTG']))
      
#---------------------------------------------------------------------------------
#PART 4 -  See section 9.6  EFetch: Downloading full records from Entrez
#       This function will use a pubmed ID number to fetch a file from the internet
#       and save it to a file on your computer, under the pubmed id name.

#FUNCTION NAME: FetchGenbankFile
#PARAMETERS: 1 (A string - Pubmed ID number)
#PURPOSE: The function should:
#           (1) Get the genbank file as text from the nucleotide database.
#           (2) Write it to a file. (See below)
#RETURN VALUES: 1 (True)

from Bio import Entrez
#== FUNCTION 4 ==
def FetchGenbankFile(pmid):
    Entrez.email = "example@example.org"
    handle = Entrez.efetch(db="nucleotide", id=pmid, rettype="gb", retmode="text")
    with open(pmid+".gbk.txt", "w") as f:
        f.write(handle.read())
        
    return 1

# Example function call
#print(FetchGenbankFile('186972394'))
"""
If you were to call the FetchGenbankFile with the following id:

x=FetchGenbankFile('186972394')

You should get a file named 186972394.gbk.txt in your computer with the Genbank file below.

From the tutorial:

>>> from Bio import Entrez
>>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
>>> handle = Entrez.efetch(db="nucleotide", id="186972394", rettype="gb", retmode="text")
>>> print handle.read()
LOCUS       EU490707                1302 bp    DNA     linear   PLN 05-MAY-2008
DEFINITION  Selenipedium aequinoctiale maturase K (matK) gene, partial cds;
            chloroplast.
ACCESSION   EU490707
VERSION     EU490707.1  GI:186972394
KEYWORDS    .
SOURCE      chloroplast Selenipedium aequinoctiale
  ORGANISM  Selenipedium aequinoctiale
            Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
            Spermatophyta; Magnoliophyta; Liliopsida; Asparagales; Orchidaceae;
            Cypripedioideae; Selenipedium.
REFERENCE   1  (bases 1 to 1302)
  AUTHORS   Neubig,K.M., Whitten,W.M., Carlsward,B.S., Blanco,M.A.,
            Endara,C.L., Williams,N.H. and Moore,M.J.
  TITLE     Phylogenetic utility of ycf1 in orchids
  JOURNAL   Unpublished
REFERENCE   2  (bases 1 to 1302)
  AUTHORS   Neubig,K.M., Whitten,W.M., Carlsward,B.S., Blanco,M.A.,
            Endara,C.L., Williams,N.H. and Moore,M.J.
  TITLE     Direct Submission
  JOURNAL   Submitted (14-FEB-2008) Department of Botany, University of
            Florida, 220 Bartram Hall, Gainesville, FL 32611-8526, USA
FEATURES             Location/Qualifiers
     source          1..1302
                     /organism="Selenipedium aequinoctiale"
                     /organelle="plastid:chloroplast"
                     /mol_type="genomic DNA"
                     /specimen_voucher="FLAS:Blanco 2475"
                     /db_xref="taxon:256374"
     gene            <1..>1302
                     /gene="matK"
     CDS             <1..>1302
                     /gene="matK"
                     /codon_start=1
                     /transl_table=11
                     /product="maturase K"
                     /protein_id="ACC99456.1"
                     /db_xref="GI:186972395"
                     /translation="IFYEPVEIFGYDNKSSLVLVKRLITRMYQQNFLISSVNDSNQKG
                     FWGHKHFFSSHFSSQMVSEGFGVILEIPFSSQLVSSLEEKKIPKYQNLRSIHSIFPFL
                     EDKFLHLNYVSDLLIPHPIHLEILVQILQCRIKDVPSLHLLRLLFHEYHNLNSLITSK
                     KFIYAFSKRKKRFLWLLYNSYVYECEYLFQFLRKQSSYLRSTSSGVFLERTHLYVKIE
                     HLLVVCCNSFQRILCFLKDPFMHYVRYQGKAILASKGTLILMKKWKFHLVNFWQSYFH
                     FWSQPYRIHIKQLSNYSFSFLGYFSSVLENHLVVRNQMLENSFIINLLTKKFDTIAPV
                     ISLIGSLSKAQFCTVLGHPISKPIWTDFSDSDILDRFCRICRNLCRYHSGSSKKQVLY
                     RIKYILRLSCARTLARKHKSTVRTFMRRLGSGLLEEFFMEEE"
ORIGIN      
        1 attttttacg aacctgtgga aatttttggt tatgacaata aatctagttt agtacttgtg
       61 aaacgtttaa ttactcgaat gtatcaacag aattttttga tttcttcggt taatgattct
      121 aaccaaaaag gattttgggg gcacaagcat tttttttctt ctcatttttc ttctcaaatg
      181 gtatcagaag gttttggagt cattctggaa attccattct cgtcgcaatt agtatcttct
      241 cttgaagaaa aaaaaatacc aaaatatcag aatttacgat ctattcattc aatatttccc
      301 tttttagaag acaaattttt acatttgaat tatgtgtcag atctactaat accccatccc
      361 atccatctgg aaatcttggt tcaaatcctt caatgccgga tcaaggatgt tccttctttg
      421 catttattgc gattgctttt ccacgaatat cataatttga atagtctcat tacttcaaag
      481 aaattcattt acgccttttc aaaaagaaag aaaagattcc tttggttact atataattct
      541 tatgtatatg aatgcgaata tctattccag tttcttcgta aacagtcttc ttatttacga
      601 tcaacatctt ctggagtctt tcttgagcga acacatttat atgtaaaaat agaacatctt
      661 ctagtagtgt gttgtaattc ttttcagagg atcctatgct ttctcaagga tcctttcatg
      721 cattatgttc gatatcaagg aaaagcaatt ctggcttcaa agggaactct tattctgatg
      781 aagaaatgga aatttcatct tgtgaatttt tggcaatctt attttcactt ttggtctcaa
      841 ccgtatagga ttcatataaa gcaattatcc aactattcct tctcttttct ggggtatttt
      901 tcaagtgtac tagaaaatca tttggtagta agaaatcaaa tgctagagaa ttcatttata
      961 ataaatcttc tgactaagaa attcgatacc atagccccag ttatttctct tattggatca
     1021 ttgtcgaaag ctcaattttg tactgtattg ggtcatccta ttagtaaacc gatctggacc
     1081 gatttctcgg attctgatat tcttgatcga ttttgccgga tatgtagaaa tctttgtcgt
     1141 tatcacagcg gatcctcaaa aaaacaggtt ttgtatcgta taaaatatat acttcgactt
     1201 tcgtgtgcta gaactttggc acggaaacat aaaagtacag tacgcacttt tatgcgaaga
     1261 ttaggttcgg gattattaga agaattcttt atggaagaag aa
//"""

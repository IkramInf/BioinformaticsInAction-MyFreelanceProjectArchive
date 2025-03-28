# import required libraries
import sys
import argparse

def get_filehandle(filename, mode):
    try:
        f = open(filename, mode)
        return f
    except (OSError, ValueError, Exception) as error:
        print(error)
        
def _verify_lists(list_headers, list_seqs):
    if len(list_headers) != len(list_seqs):
        sys.exit("The size of the sequence list and the size of the header list is not same!")
    else:
        return True
    
def get_fasta_lists(fh_in):
    lines = fh_in.readlines()
    index = [d for d in range(len(lines)) if lines[d].startswith(">")]
    headers = [t.strip()[1:] for t in lines if t.startswith(">")]
    sequences = []
    for i in range(len(index)):
        if index[i] == index[-1]:
            sequences.append("".join(lines[index[i]+1:]).replace("\n", ""))
        else:
            sequences.append("".join(lines[index[i]+1:index[i+1]]).replace("\n", ""))
            
    if not sequences[-1]:
        sequences.remove(sequences[-1])
        
    flag = _verify_lists(headers, sequences)
    
    if flag:
        return headers, sequences

if __name__ == "__main__":
    # add command line options
    parser = argparse.ArgumentParser(description="Provide a FASTA file to perform splitting on sequence and secondary structure")
    parser.add_argument('-i', '--infile', help='Path to file to open', required=True)
    args = parser.parse_args()
    
    # get file handle of sequence file and split headers and sequences
    fh_in = get_filehandle(args.infile, 'r')
    list_headers, list_seqs = get_fasta_lists(fh_in)
    # get a file handle for write output
    fh_out1 = get_filehandle("pdb_protein.fasta", "w")
    fh_out2 = get_filehandle("pdb_ss.fasta", "w")
    
    for header, seq in zip(list_headers, list_seqs):
        if "sequence" in header:
            fh_out1.write(f">{header}\n{seq}\n")
        else:
            fh_out2.write(f">{header}\n{seq}\n")
    
    # close all the file
    fh_in.close()
    fh_out1.close()
    fh_out2.close()
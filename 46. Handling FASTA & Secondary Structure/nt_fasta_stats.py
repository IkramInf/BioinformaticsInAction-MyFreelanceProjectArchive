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
    
def _get_ncbi_accession(header_string):
    return header_string.strip().split()[0]    
    
def _get_num_nucleotides(base, seq):
    if base in ['A', 'C', 'G', 'T', 'N']:
        return seq.count(base)
    else:
        sys.exit("Did not code this condition")
        
def _calculate_gc(seq):
    return float(seq.count('G')+seq.count('C'))/len(seq)*100

def output_seq_statistics(list_headers, list_seqs, fh_out):
    fh_out.write(f"Number\tAccession\tA's\tG's\tC's\tT's\tN's\tLength\tGC%\n")
    
    for i, seq in enumerate(zip(list_headers, list_seqs)):
        a_nt = _get_num_nucleotides('A', seq[1])
        c_nt = _get_num_nucleotides('C', seq[1])
        g_nt = _get_num_nucleotides('G', seq[1])
        t_nt = _get_num_nucleotides('T', seq[1])
        n_nt = _get_num_nucleotides('N', seq[1])
        accession_string = _get_ncbi_accession(seq[0])
        gc = round(_calculate_gc(seq[1]), 1)
        fh_out.write(f"{i+1}\t{accession_string}\t{a_nt}\t{g_nt}\t{c_nt}\t{t_nt}\t{n_nt}\t{len(seq[1])}\t{gc}\n")
        

if __name__ == "__main__":
    # add command line options
    parser = argparse.ArgumentParser(description="Provide a FASTA file to generate nucleotide statistics")
    parser.add_argument('-i', '--infile', help='Path to file to open', required=True)
    parser.add_argument('-o', '--outfile', help='Path to file to write', required=True)
    args = parser.parse_args()
    # get file handle of sequence file and split headers and sequences
    fh_in = get_filehandle(args.infile, 'r')
    list_headers, list_seqs = get_fasta_lists(fh_in)
    
    # get a file handle for write output
    fh_out = get_filehandle(args.outfile, "w")
    output_seq_statistics(list_headers, list_seqs, fh_out)
    
    # close all the file
    fh_in.close()
    fh_out.close()
    
import argparse

def read_seqs(filename):
    """
    Read fasta format sequence after '## FASTA' in gff file
    
    Parameters:
        filename: GFF filename
    Returns:
        A dictionary containing fasta file as {header:sequence} pair
    """
    with open(filename, "r") as f:
        lines = f.readlines()
        index = [d for d in range(len(lines)) if lines[d].startswith(">")]
        titles = [t.strip()[1:] for t in lines if t.startswith(">")]
        sequences = []
        for i in range(len(index)):
            if index[i] == index[-1]:
                sequences.append("".join(lines[index[i]+1:]).replace("\n", ""))
            else:
                sequences.append("".join(lines[index[i]+1:index[i+1]]).replace("\n", ""))

    return {title:seq for title, seq in zip(titles, sequences)}

def extract_info(filename):
    """
    Extract the annotations from gff file
    
    Parameters:
        filename: GFF filename
    Returns:
        A dictionary containing the desired annotations
    """
    infos = {}
    with open(filename, "r") as f:
        lines = f.readlines()
        onehash = [i for i, l in enumerate(lines) if l.startswith("#") and not l.startswith("##")][-1]
        threehash = [i for i, l in enumerate(lines) if l.startswith("###")][-1]
        lines = lines[onehash+1:threehash]

        for line in lines:
            l = line.split("\t")
            attrs = [l[2], l[0], l[3], l[4], l[6]]
            for x in l[-1].split(";"):
                d = x.split("=")
                infos.setdefault(d[1].strip(), []).append([d[0].strip()]+attrs)
    return infos

def format_print(fastas, item):
    """
    Format output to print on screen
    
    Parameters:
        fastas: a dictionary of header:sequence pairs
        item : a list --> [seqid, type, start, end, strand] for corresponding attribute value
    Returns:
        None
    """
    seqid, s, e, strand = item[2:]
    s, e = int(s), int(e)
    seq = fastas[seqid]
    if strand == "+":
        print(f">{args.type}:{args.attribute}:{args.value}\n{seq[s:e]}")
    elif strand == "-":
        seq = seq[::-1]
        print(f">{args.type}:{args.attribute}:{args.value}\n{seq[s:e]}")
    else:
        print(f">{args.type}:{args.attribute}:{args.value}\n{seq}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Provide a gff file to parse it")
    parser.add_argument('-i', '--source_gff', required=True, help="Enter GFF Filename...")
    parser.add_argument('-t','--type', required=True, help="Enter Type of Sequence...")
    parser.add_argument('-a','--attribute', required=True, help="Enter Attribute Name...")
    parser.add_argument('-v','--value', required=True, help="Enter Attribute Value...")
    args = parser.parse_args()
    
    fastas = read_seqs(args.source_gff)
    infos = extract_info(args.source_gff)

    items = []
    if args.value in infos.keys():
        for val in infos[args.value]:
            if (args.attribute in val) and (args.type in val):
                items.append(val)                                            

        if len(items) == 1:
            format_print(fastas, items[0])

        elif len(items) == 2:
            print("Warning: More than one match found!!!")
            for item in items:
                format_print(fastas, item)
        else:
            print("Warning: No match found!!!")
    else:
        print("Warning: No match found!!!")
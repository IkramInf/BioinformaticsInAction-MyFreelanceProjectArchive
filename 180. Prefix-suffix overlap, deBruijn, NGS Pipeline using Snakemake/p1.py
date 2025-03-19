# Import Libraries
import sys
from Bio import SeqIO
import graphviz

def find_overlap(read1, read2, min_overlap=10):
    """
    Compute prefix-suffix overlap between two reads
    """
    max_overlap = 0
    N = min(len(read1), len(read2))
    for index in range(min_overlap, N + 1):
        if read1[-index:] == read2[:index]:
            max_overlap = index
    return max_overlap

def overlap_graph(reads, min_overlap=10):
    """
    Build overlap graph
    """
    graph = {}
    for i, read1 in enumerate(reads):
        for j, read2 in enumerate(reads):
            if (i != j):
                overlap_length = find_overlap(read1, read2, min_overlap)
                if overlap_length >= min_overlap:
                    graph[(read1, read2)] = overlap_length
 
    return graph

def plot_graph(overlap_graph):
    """
    Plot the graph using Graphviz
    """
    dot = graphviz.Digraph(format="png")
    for (read1, read2), weight in overlap_graph.items():
        dot.edge(read1, read2, label=str(weight))
    dot.render("overlap", format="png", cleanup=True)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 p1.py <fastq_filename>")
        sys.exit(1)

    # Extract input filename
    input_filename = sys.argv[1]
    
    reads = [str(record.seq) for record in SeqIO.parse(input_filename, "fastq")]
    overlap_graph = overlap_graph(reads)
    plot_graph(overlap_graph)

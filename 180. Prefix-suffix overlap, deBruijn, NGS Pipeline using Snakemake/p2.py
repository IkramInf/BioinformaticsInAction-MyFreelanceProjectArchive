# Import Libraries
import sys
from Bio import SeqIO
import graphviz

def debruijn_graph(sequence, k=10):
    """
    Construct De Bruijn graph from a fasta sequence
    """
    nodes = {}
    edges = []

    for i in range(len(sequence) - (k - 1)):
        # Generate k-mer and (k-1)-mers
        kmer = sequence[i:i+k]
        prefix = kmer[:-1]
        suffix = kmer[1:]

        # Add edge and node
        edges.append((prefix, suffix))
        nodes.setdefault(prefix, []).append(suffix)

    return nodes, edges

def plot_graph(edges):
    """
    Plot the De Bruijn graph using Graphviz
    """
    dot = graphviz.Digraph(format="png")
    for source, target in edges:
        dot.edge(source, target)
    dot.render("debruijn", format="png", cleanup=True)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 p2.py <fasta_filename>")
        sys.exit(1)

    # Extract input filename
    input_filename = sys.argv[1]
    
    # Read the FASTA file
    sequence = str(SeqIO.read(input_filename, "fasta").seq)
    # Build the De Bruijn graph
    k = 15
    nodes, edges = debruijn_graph(sequence, k)
    # Plot the De Bruijn graph
    plot_graph(edges)

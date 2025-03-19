"""
Implementation of the graph-theoretic Shortest Superstring problem
"""
class ShortestSuperstring:
    def __init__(self, words):
        """
        initialize class with list of word
        """
        # run discard function to discard proper substring of any other fragment
        self.words = self.discard(words)
        
    def discard(self, words):
        """
        discard proper substring of any other fragment
        """
        new = []
        for i, w1 in enumerate(words):
            # if word1 isn't substring of all word2, add word1 into new list
            if all([w1 not in w2 for j, w2 in enumerate(words) if i != j]):
                new.append(w1)
        return new

    def w(self, u, v):
        """
        implement w(uv) = max{|t| : t suffix of u, t prefix of v}
        to find out maximum overlap length
        """
        for l in range(min(len(u), len(v)), 0, -1):
            # if suffix of u matches with prefix of v
            if u.endswith(v[:l]):
                return l
        else:
            return 0

    def overlapGraphs(self):
        """
        iterate over each read assumming them as node
        {node: {edge1: w(node, edge1), ...}, ...}
        """
        graph = {}
        for w1 in (self.words):
            for w2 in (self.words):
                if w1 != w2:
                    weight = self.w(w1, w2)
                    if weight != 0:
                        # iterate over each pair of words where, word1 != word2
                        # calculate their overlap length, if length != 0, add the pair into dict
                        graph.setdefault(w1, {})[w2] = weight
        return graph

    def hamiltonianPath(self):
        """
        Find out Hamiltonian path: a path that visits all vertices exactly once
        """
        N = len(self.words) - 1
        graph = self.overlapGraphs()
        # check each node and associated edges for hamiltonian path
        # return a sorted path according to their weight value
        ham_graph = {k : v for k, v in graph.items() if len(v.keys()) == N}
        if ham_graph:
            ham_path = sorted(ham_graph.items(), key=lambda x: sum(x[1].values()), reverse=True)[0]
            ham, edges = ham_path[0], ham_path[1]
            edges = dict(sorted(edges.items(), key=lambda x: x[1], reverse=True))
            return [ham] + list(edges.keys())
        else:
            path = sorted(graph.items(), key=lambda x: len(x[1]) and sum(x[1].values()), reverse=True)[0]
            node, edges = path[0], path[1]
            edges = dict(sorted(edges.items(), key=lambda x: x[1], reverse=True))
            return [node] + list(edges.keys())

    def reassembledSuperstring(self):
        """
        merge each read in hamiltonian path to construct shortest superstring
        """
        ham_path = self.hamiltonianPath()
        superstring, p0 = [], ""
        for p1 in ham_path:
            # merge each read
            superstring.append(p1[self.w(p0, p1):])
            p0 = p1
        return ''.join(superstring)

if __name__ == "__main__":
    # read first file and find out shortest superstring
    with open("fragex1.txt", "r") as f:
        words1 = f.read().splitlines()
    ss = ShortestSuperstring(words1)
    superstring = ss.reassembledSuperstring()
    print(superstring)
    
    # read first file and find out shortest superstring
    with open("fragtext162.txt", "r") as f:
        words2 = f.read().splitlines()
    ss = ShortestSuperstring(words2)
    superstring = ss.reassembledSuperstring()
    print(superstring)
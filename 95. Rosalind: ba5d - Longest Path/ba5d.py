# import networkx module to create and investigate graph
import networkx as nx

class DirectedAcyclicGraph():
    """
    DirectedAcyclicGraph class find longest path between source and sink and calculate the path distance
    __init__():
        Parameters:
            None
        Return:
            None
            Create directed graph using networkx module and save weight value in a dictionary during initialization
    add_edge():
        Parameters:
            edge: edge of network
            weight: weight value corresponding the edge
        Return:
            None
    add_weight():
        Parameters:
            edge: edge of network
            weight: weight value corresponding the edge
        Return:
            None
    calculate_distance():
        Parameters:
            path: path to calculate distance
        Return:
            distance: calculated distance of given path 
    add_edge():
        Parameters:
            source: start node
            sink: end node
        Return:
            longest_path: longest path among all possible path
            distance: distance of longest path
        
    """
    def __init__(self):
        self.graph = nx.DiGraph()  # directed graph using networkx
        self.weight = {}  # dictionary to save weight as {edge : weight}
        
    def add_edge(self, edge, weight):
        # add edge to self.graph graph
        self.graph.add_edge(edge[0], edge[1], weight=weight)
        
    def add_weight(self, edge, weight):
        # add weight to self.weight dictionary
        self.weight[edge] = weight
        
    def calculate_distance(self, path):
        # iterate over each edge over the path and sum their weight value
        distance = sum([self.weight[(path[i], path[i+1])] for i in range(len(path)-1)])
        return distance

    def calculate_longest_path(self, source, sink):
        # find out all possible paths 
        all_possible_paths = nx.all_simple_paths(self.graph, source=source, target=sink)

        # dictionary to save path and their corresponding distance value
        distance_dicts = {}
        
        for path in all_possible_paths:
            distance = self.calculate_distance(path)
            path = list(map(str, path))
            path = "->".join(path)
            distance_dicts[path] = distance
            
        # find out longest path and distance from the distance_dicts dictionary
        longest_path, distance = sorted(distance_dicts.items(), key=lambda x:x[1], reverse=True)[0]
        # return the results
        return longest_path, distance

if __name__ == '__main__':
    
    # create object for the DirectedAcyclicGraph class
    dag = DirectedAcyclicGraph()
    
    # read the input file
    with open("rosalind_ba5d.txt", "r") as f:
        # read first line as source
        source = int(f.readline().strip())
        # read second line as sink
        sink = int(f.readline().strip())
        # read edges and weights from third line
        for val in f.readlines():
            edge = tuple(map(int, val.split(":")[0].split("->")))
            weight = int(val.split(":")[-1])
            
            # add edge and weight
            dag.add_edge(edge, weight)
            dag.add_weight(edge, weight)

    # print the result to screen
    longest_path, distance = dag.calculate_longest_path(source, sink)
    print(f"{distance}\n{longest_path}")
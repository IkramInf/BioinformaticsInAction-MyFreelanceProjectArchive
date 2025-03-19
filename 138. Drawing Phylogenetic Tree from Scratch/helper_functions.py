"""Helper functions for HW3"""
import numpy as np
from copy import deepcopy
from matplotlib.axes import Axes
from typing import Optional
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import networkx as nx
from networkx.readwrite import json_graph
from networkx.drawing.nx_agraph import graphviz_layout

class Node:
    def __init__(
        self,
        name: str,
        left: "Node",
        left_distance: float,
        right: "Node",
        right_distance: float,
        confidence: float = None,
    ):
        """A node in a binary tree produced by neighbor joining algorithm.

        Parameters
        ----------
        name: str
            Name of the node.
        left: Node
            Left child.
        left_distance: float
            The distance to the left child.
        right: Node
            Right child.
        right_distance: float
            The distance to the right child.
        confidence: float
            The confidence level of the split determined by the bootstrap method.
            Only used if you implement Bonus Problem 1.

        Notes
        -----
        The current public API needs to remain as it is, i.e., don't change the
        names of the properties in the template, as the tests expect this kind
        of structure. However, feel free to add any methods/properties/attributes
        that you might need in your tree construction.

        """
        self.name = name
        self.left = left
        self.left_distance = left_distance
        self.right = right
        self.right_distance = right_distance
        self.confidence = confidence
        self.children = []

def build_nj_matrix(D, totalDistanceD):
    #totalDistanceD = D.sum(axis=0)
    D_new = np.zeros(D.shape)
    n = D.shape[0]
    for i in range(n):
        for j in range(i+1):
            if i == j:
                D_new[i, j] = np.inf
            else:
                D_new[i, j] = (n - 2) * D[i, j] - totalDistanceD[i] - totalDistanceD[j]
                D_new[j, i] = D_new[i, j]
    return D_new

def update_matrix(dm, i, j):
    # calculate distance from m to other leaves
    m_to_other_leaves = 0.5 * (dm[i, :] + dm[j, :] - dm[i, j])
    # add distances at first row
    dm = np.vstack([m_to_other_leaves, dm])
    # add distances at first column
    dm = np.insert(dm, 0, np.append(0, m_to_other_leaves), axis=1)
    # delete i and j column
    dm = np.delete(dm, [i+1, j+1], axis=0)
    dm = np.delete(dm, [i+1, j+1], axis=1)
    # return updated matrix
    return dm

def neighbor_joining(distances: np.ndarray, labels: list) -> Node:
    """The Neighbor-Joining algorithm.

    For the same results as in the later test dendrograms;
    add new nodes to the end of the list/matrix and
    in case of ties, use np.argmin to choose the joining pair.

    Parameters
    ----------
    distances: np.ndarray
        A 2d square, symmetric distance matrix containing distances between
        data points. The diagonal entries should always be zero; d(x, x) = 0.
    labels: list
        A list of labels corresponding to entries in the distances matrix.
        Use them to set names of nodes.

    Returns
    -------
    Node
        A root node of the neighbor joining tree.

    """
    n = distances.shape[0]
    totalDistanceD = distances.sum(axis=0)

    if n == 2:
        left, right = labels[0], labels[1]
        if not isinstance(left, Node):
            left = Node(left, None, 0, None, 0)
        if not isinstance(right, Node):
            right = Node(right, None, 0, None, 0)
        clustered_label = f"({left.name},{right.name})"
        new_node = Node(clustered_label, left, distances[0, 1] / 2, right, distances[0, 1] / 2)
        new_node.children.extend([left, right])
        return new_node

    D_prime = build_nj_matrix(distances, totalDistanceD)
    i, j = np.unravel_index(np.argmin(D_prime), D_prime.shape)

    deltaD = (totalDistanceD[i] - totalDistanceD[j]) / (n - 2)
    limbLengthi = 0.5 * (distances[i, j] + deltaD)
    limbLengthj = 0.5 * (distances[i, j] - deltaD)
    distances = update_matrix(distances, i, j)

    left, right = labels[i], labels[j]
    if not isinstance(left, Node):
        left = Node(left, None, 0, None, 0)
    if not isinstance(right, Node):
        right = Node(right, None, 0, None, 0)
    clustered_label = f"({left.name},{right.name})"
    new_node = Node(clustered_label, left, limbLengthi, right, limbLengthj)
    new_node.children.extend([left, right])

    labels.insert(0, new_node)
    labels = [label for _, label in enumerate(labels) if _ not in [i + 1, j + 1]]

    return neighbor_joining(distances, labels)

def create_tree_edges(node, position, edges):
    if node.left is not None:
        edges.append((node.name, node.left.name, node.left_distance))
        create_tree_edges(node.left, pos, edges)

    if node.right is not None:
        edges.append((node.name, node.right.name, node.right_distance))
        create_tree_edges(node.right, position, edges)

def plot_nj_tree(tree: Node, ax: Axes = None, fig_name: str = None,
                 prog: str = "twopi", groups: dict = None, **kwargs) -> None:
    """A function for plotting neighbor joining phylogeny dendrogram.

    Parameters
    ----------
    tree: Node
        The root of the phylogenetic tree produced by `neighbor_joining(...)`.
    ax: Axes
        A matplotlib Axes object which should be used for plotting.
    kwargs
        Feel free to replace/use these with any additional arguments you need.
        But make sure your function can work without them, for testing purposes.

    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>>
    >>> tree = neighbor_joining(distances)
    >>> fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    >>> plot_nj_tree(tree=tree, ax=ax)
    >>> fig.savefig("example.png")

    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
        
    G = nx.DiGraph()
    position = {}
    edges = []

    def traverse_tree(node, x, y, distance):
        position[node.name] = (x, y)
        if node.left is not None:
            edges.append((node.name, node.left.name, node.left_distance))
            x_left = x - distance / 2
            traverse_tree(node.left, x_left, y - 1, distance / 2)

        if node.right is not None:
            edges.append((node.name, node.right.name, node.right_distance))
            x_right = x + distance / 2
            traverse_tree(node.right, x_right, y - 1, distance / 2)

    traverse_tree(tree, 0, 0, 100)
    G.add_weighted_edges_from(edges)
    
    position = graphviz_layout(G, prog=prog)

    leaf_nodes = [node for node in G.nodes if G.out_degree(node) == 0]
    leaf_pos = {node : position[node] for node in leaf_nodes}
    
    if groups:
        colors = cm.tab10.colors
        group_colors = {group : colors[i] for i, group in enumerate(groups)}

        node_colors = []
        for node in G.nodes:
            if node in groups['Alpha']:
                node_colors.append(group_colors['Alpha'])
            elif node in groups['Beta']:
                node_colors.append(group_colors['Beta'])
            elif node in groups['Gamma']:
                node_colors.append(group_colors['Gamma'])
            elif node in groups['Delta']:
                node_colors.append(group_colors['Delta'])
            elif node in groups['Unclassified']:
                node_colors.append(group_colors['Unclassified'])
            else:
                node_colors.append("skyblue")

        # Draw edges
        nx.draw_networkx_edges(G, position, edgelist=G.edges, edge_color='grey', arrows=False, width=2.0, alpha=0.7)

        # Draw nodes with colors
        nx.draw_networkx_nodes(G, position, nodelist=G.nodes, node_size=300,
                               node_shape='o', node_color=node_colors, alpha=0.8, ax=ax)

        # Draw labels for leaf nodes
        for leaf, pos in leaf_pos.items():
            leaf_name = getattr(leaf, 'name', str(leaf))
            ax.text(pos[0], pos[1], leaf_name, color='black', ha='center', va='center', fontsize=10, rotation=90)

        patch_list = []
        for group, color in group_colors.items():
            patch_list.append(mpatches.Patch(color=color, label=group))
        ax.legend(handles=patch_list, loc= 'upper left')
        
    else:
        nx.draw(G, position, with_labels=False, arrows=False,
                node_size=300, node_shape='o', node_color='skyblue',
                width=2.0, edge_color='lightgrey')

        for leaf, pos in leaf_pos.items():
            leaf_name = getattr(leaf, 'name', str(leaf))
            ax.text(pos[0], pos[1], leaf_name, color='black', ha='center', va='center', fontsize=10, rotation=90)
        
    plt.axis('off')
    plt.tight_layout()
    if fig_name:
        plt.savefig(fig_name, dpi=300, format='png')
    plt.show()
    
    return ax

def _find_a_parent_to_node(tree: Node, node: Node) -> tuple:
    """Utility function for reroot_tree"""
    stack = [tree]

    while len(stack) > 0:

        current_node = stack.pop()
        if node.name == current_node.left.name:
            return current_node, "left"
        elif node.name == current_node.right.name:
            return current_node, "right"

        stack += [
            n for n in [current_node.left, current_node.right] if n.left is not None
        ]

    return None


def _remove_child_from_parent(parent_node: Node, child_location: str) -> None:
    """Utility function for reroot_tree"""
    setattr(parent_node, child_location, None)
    setattr(parent_node, f"{child_location}_distance", 0.0)


def reroot_tree(original_tree: Node, outgroup_node: Node) -> Node:
    """A function to create a new root and invert a tree accordingly.

    This function reroots tree with nodes in original format. If you
    added any other relational parameters to your nodes, these parameters
    will not be inverted! You can modify this implementation or create
    additional functions to fix them.

    Parameters
    ----------
    original_tree: Node
        A root node of the original tree.
    outgroup_node: Node
        A Node to set as an outgroup (already included in a tree).
        Find it by it's name and then use it as parameter.

    Returns
    -------
    Node
        Inverted tree with a new root node.
    """
    tree = deepcopy(original_tree)

    parent, child_loc = _find_a_parent_to_node(tree, outgroup_node)
    distance = getattr(parent, f"{child_loc}_distance")
    _remove_child_from_parent(parent, child_loc)

    new_root = Node("new_root", parent, distance / 2, outgroup_node, distance / 2)
    child = parent

    while tree != child:
        parent, child_loc = _find_a_parent_to_node(tree, child)

        distance = getattr(parent, f"{child_loc}_distance")
        _remove_child_from_parent(parent, child_loc)

        empty_side = "left" if child.left is None else "right"
        setattr(child, f"{empty_side}_distance", distance)
        setattr(child, empty_side, parent)

        if tree.name == parent.name:
            break
        child = parent

    other_child_loc = "right" if child_loc == "left" else "left"
    other_child_distance = getattr(parent, f"{other_child_loc}_distance")

    setattr(child, f"{empty_side}_distance", other_child_distance + distance)
    setattr(child, empty_side, getattr(parent, other_child_loc))

    return new_root

def count_leaves(tree: Node) -> int:
    """
    Count the number of leaves in the subtree rooted at the given node.
    """
    if not tree.left and not tree.right:
        return 1
    left_leaves = count_leaves(tree.left) if tree.left else 0
    right_leaves = count_leaves(tree.right) if tree.right else 0
    return left_leaves + right_leaves

def sort_children(tree: Node) -> None:
    """Sort the children of a tree by their corresponding number of leaves.

    The tree can be changed inplace.

    Paramteres
    ----------
    tree: Node
        The root node of the tree.

    """
    tree.children.sort(key=lambda x: count_leaves(x), reverse=True)
    for child in tree.children:
        sort_children(child)

def plot_nj_tree_radial(tree: Node, ax: Axes = None, fig_name: str = None, 
                        prog: str = "twopi", groups: dict = None, **kwargs) -> None:
    """A function for plotting neighbor joining phylogeny dendrogram
    with a radial layout.

    Parameters
    ----------
    tree: Node
        The root of the phylogenetic tree produced by `neighbor_joining(...)`.
    ax: Axes
        A matplotlib Axes object which should be used for plotting.
    kwargs
        Feel free to replace/use these with any additional arguments you need.

    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>>
    >>> tree = neighbor_joining(distances)
    >>> fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    >>> plot_nj_tree_radial(tree=tree, ax=ax)
    >>> fig.savefig("example_radial.png")

    """
    plot_nj_tree(tree, ax, fig_name, prog, groups)

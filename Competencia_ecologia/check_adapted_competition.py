import sys
import pandas as pd
import networkx as nx

# read network file from command line argument:
network_file = sys.argv[1]

# read adjacency matrix
mat = pd.read_csv(network_file, index_col=0)

# build directed graph from adjacency matrix:
G = nx.from_pandas_adjacency(mat, create_using=nx.DiGraph())

# convert to undirected (collapse A to B and B to A into one edge):
G_undir = G.to_undirected()

# compute complement:
comp = nx.complement(G_undir)

# print edge list without header and index, space-separated, for bash consumption:
pair = pd.DataFrame(comp.edges(), columns=["node_1", "node_2"])
# print only if non-empty (prints nothing when empty so bash -z check works):
if not pair.empty:
    print(pair.to_string(index=False, header=False))
    
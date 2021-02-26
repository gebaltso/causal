
from collections import defaultdict

import itertools
import time
import sys
import copy
import csv
from itertools import combinations
from loadGraph import loadGraph
from Louvain import Louvain
from FindEndogenousViaRandomWalks import endogenous
from rankMetricM import rankMetricM
from rankMetricRC import rankMetricRC
from ResAndDis import resAndDis
import igraph as ig
from initPartition import initPartition
import numpy as np
import networkx as nx
import netwulf as nw


edgeFile = "datasets/coraEdgelist.csv"
#edgeFile = "datasets/email.csv"
    

# Load the graph
G = loadGraph(edgeFile)

# Number of nodes of the initial G 
nodesLen = G.number_of_nodes()
#nodesLen = g.vcount()
#for edges: G.ecount()

g = ig.Graph.TupleList(G.edges(), directed=True)

# Find the initial partitioning from the given comm file
InitialPartition, lisInitialPartition = initPartition(g)

partitionI = InitialPartition



# Create a network
#G = nx.random_partition_graph([10, 10, 10], .25, .01)

# Change 'block' node attribute to 'group'
#for k, v in G.nodes(data=True):
#    v['group'] = v['block']; del v['block']

# Or detect communities and encode them in 'group' attribute
# import community
# bb = community.best_partition(G)
nx.set_node_attributes(G,partitionI, 'group')

# Set node 'size' attributes
for n, data in G.nodes(data=True):
    data['size'] = np.random.random()

# Set link 'weight' attributes
#for n1, n2, data in G.edges(data=True):
#    data['weight'] = np.random.random()

nw.visualize(G)












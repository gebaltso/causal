#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 17:32:21 2020

@author: georgia
"""

import community as community_louvain
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from collections import defaultdict
import copy
import sys
import random
from loadGraph import loadGraph
from Louvain import Louvain
from FindEndogenousViaRandomWalks import endogenous
from metric import metric



# Load the graph
G = loadGraph()

# Compute the best partition
partition, Q = Louvain(G)
#print("Modularity of Louvain partioning is:", Q)

# communities_dict contains the communities and nodes of each community
communities = set(partition.values())
communities_dict = {c: [k for k, v in partition.items() if v == c] for c in communities}


# Compute the endogenous set (via random walks). Give as parameters the graph, 
# the number of random walks, the length of the random walks and the initial node
# returns E which is the list of all the random walks and endoEdges which is the set of edges in these walks
E, endoEdges = endogenous(G, 10, 3, 9)

print(E)
print(endoEdges)


# Keep the nodes of the above endoEdges set in endoNodes list
flat_E = [item for sublist in E for item in sublist]
endoNodes = list((set(flat_E)))


# Compute the metric M for the above results and node as 1rst parameter
# the 1rst parameter should be changed to every(?) node of endoEdges
M = metric(17, partition, communities_dict, G)


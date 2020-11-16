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
import itertools
from loadGraph import loadGraph
from Louvain import Louvain
from FindEndogenousViaRandomWalks import endogenous
from rankMetricM import rankMetricM
from rankMetricRC import rankMetricRC



# Load the graph
G = loadGraph()

# Compute the best partition
partition, Q = Louvain(G)
#partition = {0: 0, 1: 0, 2: 0, 3: 0, 4: 1, 5: 1, 6: 1, 7: 0, 8: 2, 9: 2, 10: 1, 11: 0, 12: 0, 13: 0, 14: 2, 15: 2, 16: 1, 17: 0, 18: 2, 19: 0, 20: 2, 21: 0, 22: 2, 23: 3, 24: 3, 25: 3, 26: 2, 27: 3, 28: 3, 29: 2, 30: 2, 31: 3, 32: 2, 33: 2}
print("Modularity of Louvain partioning is:", Q)
#print(partition)

# communities_dict contains the communities and nodes of each community
communities = set(partition.values())
communities_dict = {c: [k for k, v in partition.items() if v == c] for c in communities}

# The initial query node
initnode = 8 

# The number of items to examine as causal (selected from endogenous)
nOfCausal = 3 

# select how to compute the endogenous set
method = 'incident' 
rankMetric = 'rc'

##########################################################################################

if method == 'exo':
    # the user gives the endogenous set exogenously
    endoEdges = {(4,5),(5,6),(3,6)}
    endoNodes = set(itertools.chain.from_iterable(list(endoEdges)))
        
elif method == 'random':
    # Compute the endogenous set via random walks. Give as parameters the graph, 
    # the number of random walks, the length of the random walks and the initial node
    # returns E which is the list of all the random walks and endoEdges which is the set of edges in these walks
    nrw = 10 #number of random walks
    l = 3 #length of random walks
    E, endoEdges = endogenous(G, nrw, l, initnode)
        
    # Keep the nodes of the above endoEdges set in endoNodes list
    flat_E = [item for sublist in E for item in sublist]
    endoNodes = list((set(flat_E)))
    
elif method == 'incident':
    endoEdges = list(G.edges(initnode))
    endoNodes = set(itertools.chain.from_iterable(list(endoEdges)))
    endoNodes.remove(initnode) # remove the initial node as G.edges include it

elif method == 'incomm':
    # Find the community C that node n belongs
    C = partition.get(initnode)    
    # Find all the nodes of community C (n belongs to C)
    communityTmp =  [v for k, v in communities_dict.items() if k == C]    
    # Create flat list to make one list from a list with lists    
    community = [item for sublist in communityTmp for item in sublist]
        
    endoEdges = []
    for i in community:
        for j in community:
            if G.has_edge(i, j):
                if (j,i) not in endoEdges:
                    endoEdges.append((i,j))
    
    endoNodes = set(itertools.chain.from_iterable(list(endoEdges)))
    endoNodes.remove(initnode) # remove the initial node as G.edges include it

###########################################################################################

# Dictionary to keep nodes of endoNodes and their corresponding M values
# in order to select which edges to examine as causal
possibleCausal = defaultdict(dict)
finalCausalNodes = []

# The edges that will be examined as causal
CausalEdges = []

# Compute the metric M for the endoNodes found before
if rankMetric == 'emb':
    for i in endoNodes:
        possibleCausal[i] = rankMetricM(i, partition, communities_dict, G)
    sortedPossibleCausal = {k: v for k, v in sorted(possibleCausal.items(), key=lambda item: item[1], reverse=True)}
    finalCausalNodes = list(sortedPossibleCausal)[0:nOfCausal]
    for e in list(endoEdges):
        if (e[0] in finalCausalNodes) or (e[1] in finalCausalNodes):
            CausalEdges.append(e)

    
elif rankMetric == 'rc':
    for i in endoNodes:
        possibleCausal[i] = rankMetricRC(i, partition, communities_dict, G)
    sortedPossibleCausal = {k: v for k, v in sorted(possibleCausal.items(), key=lambda item: item[1], reverse=True)}
#    print(sortedPossibleCausal)
    finalCausalNodes = list(sortedPossibleCausal)[0:nOfCausal]
    for e in list(endoEdges):
        if (e[0] in finalCausalNodes) or (e[1] in finalCausalNodes):
            CausalEdges.append(e)

###########################################################################################
 
print("The edges that will be examined as causal are: ", CausalEdges)
                           
# Remove from G the edges proposed as causal and re-run community detection algorithm
G.remove_edge(CausalEdges[1][0], CausalEdges[1][1])
G.remove_edge(CausalEdges[2][0], CausalEdges[2][1])
partition2, Q2 = Louvain(G)  

print("Modularity of Louvain partioning is:", Q2)
#if partition2[initnode] != partition[initnode]:
#    print("The query node changed community and gone to: ", partition2[initnode])
    
# Find all the nodes that have changed communities
nodesChanged = []
for i in partition:
    if partition2[i] != partition[i]:
        nodesChanged.append(i)
            
print("The nodes that due to causals have changed communities: ", nodesChanged)
    
# Compute 
    
    
    
    
    

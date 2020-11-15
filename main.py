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
#print("Modularity of Louvain partioning is:", Q)

# communities_dict contains the communities and nodes of each community
communities = set(partition.values())
communities_dict = {c: [k for k, v in partition.items() if v == c] for c in communities}

initnode = 9 #the initial query node

# select how to compute the endogenous set
method = 'incomm' 
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
    
#    print("Random walks of length ", l, ": ", E)
#    print("Edges induced by the previous random walks: ", endoEdges)
    
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

if rankMetric == 'emb':
    # Compute the metric M for the above results and node as 1rst parameter
    # the 1rst parameter should be changed to every(?) node of endoEdges
    M = rankMetricM(17, partition, communities_dict, G)
    
elif rankMetric == 'rc':
    M = rankMetricRC(33, partition, communities_dict, G)
    print(M)


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 16:32:37 2020

@author: georgiabaltsou
"""


import networkx as nx
from collections import defaultdict
import operator


def rankMetricM(n, partition, communities_dict, G):
    
    # Find the community C that node n belongs
    C = partition.get(n)
    
    # Find all the nodes of community C (n belongs to C)
    communityTmp =  [v for k, v in communities_dict.items() if k == C]    
    # Create flat list to make one list from a list with lists    
    community = [item for sublist in communityTmp for item in sublist]
    
    # Find the degree of n inside C
#    degInCofn = 0
#    for i in community:
#        if G.has_edge(n, i):
#            degInCofn +=1
    
    # Find the degree of each node inside its community
    degInC = defaultdict(dict)
    for i in community:
        degInC[i] = 0
        
    for i in degInC:
        for j in community:
            if G.has_edge(i, j):
                degInC[i] +=1
            
    # Keep the max value of the above dict
    maxDegInC = max(degInC.items(), key=operator.itemgetter(1))[1]
    # If I wanted the index aka the node with the maximum in degree in community
    # nodeWithMaxDegInC = max(degInC.items(), key=operator.itemgetter(1))[0]

    # Compute embededdness as the ratio between the degree of node n inside community C
    # and the total degree of node n which is deg
    deg = len(list(G.neighbors(n)))
    emb = (degInC[n]/deg)
    
    # Compute metric M
    M = emb*degInC[n]/maxDegInC
    
    return M
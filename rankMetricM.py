#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 16:32:37 2020

@author: georgiabaltsou
"""


from collections import defaultdict
import operator


def rankMetricM(n, partition, memb, G):
    
    # Find the community C that node n belongs
    C = memb[n]
    # Find all the nodes of community C (n belongs to C)
    community = [i for i, e in enumerate(memb) if e == C]
    
#    print("comm=", community, n)

    
    # Find the degree of each node inside its community
    degInC = defaultdict(dict)
    for i in community:
        degInC[i] = 0
#    print(degInC)
        
    for i in degInC:
        for j in community:
            if i in G.neighbors(j):
                degInC[i] +=1
#    print(degInC)
        
    # Keep the max value of the above dict
    maxDegInC = max(degInC.items(), key=operator.itemgetter(1))[1]
    # If I wanted the index aka the node with the maximum in degree in community
    # nodeWithMaxDegInC = max(degInC.items(), key=operator.itemgetter(1))[0]

    # Compute embededdness as the ratio between the degree of node n inside community C
    # and the total degree of node n which is deg
#    deg = len(list(G.neighbors(n)))
    deg = G.degree(n) #HERE IT TAKES INTO ACCOUNT SELF EDGES AS DOUBLE COUNT
    emb = (degInC[n]/deg)
    
    # Compute metric M
    M = emb*degInC[n]/maxDegInC
    
    return M
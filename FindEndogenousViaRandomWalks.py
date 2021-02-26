#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 12:20:25 2020

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
import igraph as ig

def messagePropagation(n, N, k, G):
    
    
    # initialize the probability dictionary
    Pr = defaultdict(dict)
    # initialize the walk of each iteration
    walk = []
    
    addingEdge = []
    
    
    while N < k:
    
        walk.append(n)
        
        # s is the source node of the edge
        s = n
         
        # find the adjacent nodes to node n 
        neiOfn = G.neighbors(n)
        incomNei = G.neighbors(n, mode="in")
        outcomNei = G.neighbors(n, mode="out")

#        neiOfn = G.neighbors(n, mode="out") # for incoming neis or out for outgoing neis
    
        if neiOfn == []: # If a node does not have any successors stop the procedure
            break
    
#        neiOfn = list(G.neighbors(n))
#        deg = len(neiOfn)
#        deg = G.outdegree(n) 
        deg = G.degree(n) #deg = G.indegree(n) for in degree of n, outdegree for out degree of n. 


#        print("n=", n, "nei=", neiOfn, "deg=", G.degree[n], "adj=", list(G[n]))
        
        #calculate the possibility of visiting a neighbour node
        Pr = defaultdict(dict)
        for i in neiOfn:
            Pr[i] = 1/deg
            
        
        # find the next node in random way considering the probabilities of edges    
        n = (np.random.choice(list(Pr.keys()), p=list(Pr.values())))
    
        # t is the target node of the edge
        t = n

        N += 1
        
        # addingEdge list has the endogenous edges
        if N != k:
            if t in incomNei:
                addingEdge.append((t,s))
            else:
                addingEdge.append((s,t))
    
    return walk, addingEdge
    


def endogenous(G, p, k, n):
    
    allTheWalks = []  
    allEndogenous = []

    for i in range(1, p+1):
        N = 0 #counter to check the length of the random walk                  
        walk, addingEdge = messagePropagation(n, N, k, G)
        allTheWalks.append(walk)
        allEndogenous.append(addingEdge)
     
        
    # Create flat list to make one list from a list with lists    
    flat_list = [item for sublist in allEndogenous for item in sublist]

    
    # Return the walks and the set of endogenous edges without double items
    # If I want a list again should use list(set(flat_list))
    return allTheWalks, set(flat_list)


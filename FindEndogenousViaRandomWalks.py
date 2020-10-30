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
         
        # find the adjacent nodes to node un        
        neiOfn = list(G.neighbors(n))
        
        
        #calculate the possibility of visiting a neighbour node
        Pr = defaultdict(dict)
        for i in neiOfn:
            Pr[i] = 1/G.degree[n]
            
#        print(list(Pr.keys()))
                
        # find the next node in random way considering the probabilities of edges    
        n = (np.random.choice(list(Pr.keys()), p=list(Pr.values())))
        
        # t is the target node of the edge
        t = n

        N += 1
        
        # addingEdge list has the endogenous edges
        if N != k:
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


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
    
    
    while N < k:
    
        walk.append(n)
         
        # find the adjacent nodes to node un        
        neiOfn = list(G.neighbors(n))
        
        #calculate the possibility of visiting a neighbour node
        Pr = defaultdict(dict)
        for i in neiOfn:
            Pr[i] = 1/G.degree[n]
                
        # find the next node in random way considering the probabilities of edges    
        n = (np.random.choice(list(Pr.keys()), p=list(Pr.values())))
    
        N += 1
    
    return walk
    

G = nx.karate_club_graph()

#how many random walks will be held
p = 10 
#what length these random walks will have     
k = 3
#initial node
n = 9               
#the result of the random walks after p iterations
allTheWalks = []    

for i in range(1, p+1):
    N = 0 #counter to check the length of the random walk                  
    walk = messagePropagation(n, N, k, G)
    allTheWalks.append(walk)
    
print(allTheWalks)
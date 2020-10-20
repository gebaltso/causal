##!/usr/bin/env python3
## -*- coding: utf-8 -*-
#"""
#Created on Mon Oct 19 14:27:29 2020
#
#@author: georgia
#"""

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
    
    # initialize the T containing the nodes already visited
    T = []
    
#    # initialize the weight dict
#    weightDict = defaultdict(dict)
    
    # initialize the probability dictionary
    Pr = defaultdict(dict)
    
    while N < k and G.number_of_nodes() > len(T) :
         
        # find the adjacent nodes to node un        
        neiOfn = list(G.neighbors(n))
                
        # Ihat is the list of neighbour nodes that haven't been visited yet
        Ihat = neiOfn
        for i in T:
            Ihat.remove(i)
        
#        # weight = 1 for all edges in Ihat
#        for i in Ihat:
#            weightDict[i] = 1                             
        
  
        Pr = defaultdict(dict)
        for i in Ihat:
            Pr[i] = 1/G.degree[i]
            
        print(Pr)
        
        
            
        # find the index of the edge at random way considering the probabilities of edges    
        newnIndex = (np.random.choice(list(Pr.keys()), p=list(Pr.values())))
        
        print(newnIndex)
        
        # find the exact edge from the above index
        newn = list(Pr)[newnIndex]
        
        # find the node reached by the above edge
        newNode = newn[1]  
        
        
#        # increase by 1 the weight of the above edge        
#        Graph[un][newNode] += 1
#        Graph[newNode][un] += 1
#        weightDict[newUn] += 1
        
        # add the edge to the visited edges        
        T.add(newn)        
        
        # un is now the new edge
        n = newNode
        
        #increase N by 1
        N += 1
        
        
    return G
    

    
G = nx.karate_club_graph()

# depends on the preferred community size and depicts how many random walks will be held
p = 10 
      
k = 10
 
n = 9               
    
for i in range(1, p+1):
    N = 0 #counter to check the length of the k-path
                    
    Graph = messagePropagation(n, N, k, G)
                
            
for key in Graph:
    for i in Graph[key]:
        Graph[key][i] = Graph[key][i]/p
                        


###############################################################################################################
#import networkx as nx
#import random
#import matplotlib.pyplot as plt
#import operator
#import numpy as np
#
##select random graph using gnp_random_graph() function of networkx
#Graph = nx.karate_club_graph()
#nx.draw(Graph, with_labels=True, node_color='green') #draw the network graph 
#plt.figure(figsize=(15,10))
#plt.show() #to show the graph by plotting it
#
## random_node is the start node selected randomly
##random_node = random.choice([i for i in range(Graph.number_of_nodes())])
#random_node = 29
#dict_counter = {} #initialise the value for all nodes as 0
#for i in range(Graph.number_of_nodes()):
#    dict_counter[i] = 0
## increment by traversing through all neighbor nodes
#dict_counter[random_node] = dict_counter[random_node]+1
#
#
##Traversing through the neighbors of start node
#for i in range(10):
#    list_for_nodes = list(Graph.neighbors(random_node))
#    if len(list_for_nodes)==0:# if random_node having no outgoing edges
#        random_node = random.choice([i for i in range(Graph.number_of_nodes())])
#        dict_counter[random_node] = dict_counter[random_node]+1
#    else:
#        random_node = random.choice(list_for_nodes) #choose a node randomly from neighbors
#        dict_counter[random_node] = dict_counter[random_node]+1
# 
#
## using pagerank() method to provide ranks for the nodes        
#rank_node = nx.pagerank(Graph)
#
#
##sorting the values of rank and random walk of respective nodes
#sorted_rank = sorted(rank_node.items(), key=operator.itemgetter(1))
#sorted_random_walk = sorted(dict_counter.items(), key=operator.itemgetter(1))
#
#print(sorted_rank)
#print(sorted_random_walk)
#############################################################################################################

#import networkx as nx
#import numpy as np
#import random
#
#
#G = nx.karate_club_graph()
#
#pr = nx.pagerank_numpy(G, personalization={29: 1})
#
#next_node = np.random.choice(list(pr.keys()), p=list(pr.values())) 
#
#print(next_node)

#num_walks=10 #number of walks (iteration times)
#num_steps=2 #length of walks
#
#walks = list()
#i = 10 #the starting node
#for walk in range(num_walks):
#    curr_walk = [i]
#    curr = i
#    for step in range(num_steps):
#        neighbors = list(G.neighbors(curr))
#        if len(neighbors) > 0:
#            curr = random.choice(neighbors)
#        curr_walk.append(curr)
#    walks.append(curr_walk)
#print(walks)











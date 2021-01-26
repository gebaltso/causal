#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 17:32:21 2020

@author: georgia
"""

import numpy as np
from collections import defaultdict
from gephistreamer import graph
from gephistreamer import streamer
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
from igraph import *



# Load the graph
G, p = loadGraph()

# Number of nodes of the initial G 
nodesLen = G.number_of_nodes()

# If p = False I do not have the partition
if p == False:
    # Compute the best partition
    partition = Louvain(G)
#    partition, Q = Louvain(G)
#    print("Modularity of Louvain partioning is:", Q, "\n")
else:
    # I have already the initial partition
    partition = defaultdict(dict)
    
    with open("datasets/aminerEdgelistWithComms.csv", 'r') as f:
        csv_reader = csv.reader(f, delimiter=',')        
        for row in csv_reader:
            if int(row[0]) not in partition:
                partition[int(row[0])] = int(row[2])
            if int(row[1]) not in partition:
                partition[int(row[1])] = int(row[3])            
    partition = dict(partition)




##########################################################################################
## For connecting and visualizing with Gephi api
#stream = streamer.Streamer(streamer.GephiWS(hostname="localhost",port=8080,workspace="workspace5"))
#for source, target in G.edges():   
#    node_source = graph.Node(source, size=50, community=partition[source]) #for Aminer
#    node_target = graph.Node(target, size=50, community=partition[target]) #for Aminer
##    node_source = graph.Node(source, size=50, community=partition[source], label=G.nodes[source]['old_labels'])
##    node_target = graph.Node(target, size=50, community=partition[target], label=G.nodes[target]['old_labels'])
#    stream.add_node(node_source,node_target)
#    stream.add_edge(graph.Edge(node_source,node_target, directed=False))
#time.sleep(1) #It might be possible the script runs too fast and last action arn't sent properly
##########################################################################################

# communities_dict contains the communities and nodes of each community
communities = set(partition.values())
communities_dict = {c: [k for k, v in partition.items() if v == c] for c in communities}

# The initial query node
initnode = 78100
 
# select how to compute the endogenous set
method = 'incomm' # Select among 'exo' for exogenously given edges, 'random' for random walk, 
#'incident' for all edges incident to query node, 'incomm' for edges inside the same community as the query node
# random is the most time consuming

# The number of items to examine as causal (selected from endogenous)
nOfCausal = 6  # This is not used when I choose to check all the edges as causals
rankMetric = 'none' # Select among 'emb' for embeddedness, 'rc' for relative commitment, 'none' for considering all the endogenous as possible causal ones

# After the computation of ρ and γ keep the top k as causals to be removed
#k = 2 

##########################################################################################

if method == 'exo':
    # the user gives the endogenous set exogenously
    endoEdges = {(4,5),(5,6),(3,6)}
    endoNodes = set(itertools.chain.from_iterable(list(endoEdges)))
        
elif method == 'random':
    # Compute the endogenous set via random walks. Give as parameters the graph, 
    # the number of random walks, the length of the random walks and the initial node
    # returns E which is the list of all the random walks and endoEdges which is the set of edges in these walks
    nrw = 5 #number of random walks
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
    

    
    if initnode in endoNodes:
        endoNodes.remove(initnode) # remove the initial node as G.edges include it

        
print("Endogenous edges found.\n")
sys.exit()
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
            if (e[1],e[0]) not in CausalEdges: # this is to avoid adding both (a,b) & (b,a)
                CausalEdges.append(e)

    
elif rankMetric == 'rc':
    for i in endoNodes:
        possibleCausal[i] = rankMetricRC(i, partition, communities_dict, G)
    sortedPossibleCausal = {k: v for k, v in sorted(possibleCausal.items(), key=lambda item: item[1], reverse=True)}
#    print(sortedPossibleCausal)
    finalCausalNodes = list(sortedPossibleCausal)[0:nOfCausal]
    for e in list(endoEdges):
        if (e[0] in finalCausalNodes) or (e[1] in finalCausalNodes):
            if (e[1],e[0]) not in CausalEdges:
                CausalEdges.append(e)
                
                
# When I want to examine all the endogenous edges as causal ones
elif rankMetric == 'none':
    for e in list(endoEdges):
        if (e[1],e[0]) not in CausalEdges:
            CausalEdges.append(e)

print("Causal edges found.\n")        
###########################################################################################
            
# CausalEdges contains the edges of endogenous that has been selected to be examined
# as causal due to rankMetric (emb, rc or none)
#print("The edges that will be examined as causal are: ", CausalEdges, "\n")

# Dictionary that contains the causal edges ranked based on their ρ and γ values
rankOfCausals = defaultdict(dict)

for cardinality in range(1,4):
    contigencyLen = cardinality-1
    
    if cardinality == 1:        
        # Use a deep copy of the initial network to test the edge removals    
        testG = copy.deepcopy(G)          
        for e in CausalEdges:
            testG.remove_edge(e[0], e[1])
        
            partition2 = Louvain(testG)
#            partition2, Q2 = Louvain(testG)
            communities2 = set(partition.values())
            communities_dict2 = {c: [k for k, v in partition.items() if v == c] for c in communities}
#
#            print("1=",communities_dict[partition.get(initnode)])
#            print("2=",communities_dict2[partition2.get(initnode)])
        

            if partition2[initnode] != partition[initnode]:
                # Means that the query node changed community      
                rankOfCausals = resAndDis(partition, partition2, initnode, contigencyLen, nodesLen, rankOfCausals, e)
                testG.add_edge(e[0], e[1])
            else:
                testG.add_edge(e[0], e[1])
                continue                          
#        print("The edges found with cardinality = 1 are: ", rankOfCausals, "\n")  

        
    if cardinality == 2:
        testG = copy.deepcopy(G)          
#        rankOfCausals = defaultdict(dict)        
        # Keep in comb all the possible combinations of CausalEdges tuples with max length = cardinality
        comb = list(combinations(CausalEdges, cardinality))
        for e in comb:
            testG.remove_edge(e[0][0], e[0][1])
            testG.remove_edge(e[1][0], e[1][1])
            partition2 = Louvain(testG)            
#            partition2, Q2 = Louvain(testG)
            if partition2[initnode] != partition[initnode]:
                # Means that the query node changed community
                rankOfCausals = resAndDis(partition, partition2, initnode, contigencyLen, nodesLen, rankOfCausals, e)
                testG.add_edge(e[0][0], e[0][1])
                testG.add_edge(e[1][0], e[1][1])
            else:
                testG.add_edge(e[0][0], e[0][1])
                testG.add_edge(e[1][0], e[1][1])
                continue   
#        print("The edges found with cardinality = 2 are: ", rankOfCausals, "\n")

        
    if cardinality == 3:
            testG = copy.deepcopy(G)               
#            rankOfCausals = defaultdict(dict)  
            # Keep in comb all the possible combinations of CausalEdges tuples with max length = cardinality
            comb = list(combinations(CausalEdges, cardinality))
            for e in comb:
                testG.remove_edge(e[0][0], e[0][1])
                testG.remove_edge(e[1][0], e[1][1])
                testG.remove_edge(e[2][0], e[2][1]) 
                partition2 = Louvain(testG)
#                partition2, Q2 = Louvain(testG)
                if partition2[initnode] != partition[initnode]:
                    # Means that the query node changed community
                    rankOfCausals = resAndDis(partition, partition2, initnode, contigencyLen, nodesLen, rankOfCausals, e)
                    testG.add_edge(e[0][0], e[0][1])
                    testG.add_edge(e[1][0], e[1][1])
                    testG.add_edge(e[2][0], e[2][1])    
                else:
                    testG.add_edge(e[0][0], e[0][1])
                    testG.add_edge(e[1][0], e[1][1])
                    testG.add_edge(e[2][0], e[2][1])
                    continue                       
#            print("The edges found with cardinality = 3 are: ", rankOfCausals, "\n")

# Case when query node does not change community
if len(rankOfCausals) == 0:
    print("###################################################")
    print("No change with this query node. The process stops")
    print("###################################################\n")
    sys.exit()
# Sort rankOfCausals based on responsibility values
sortedRankOfCausals = sorted(rankOfCausals.items(), key=lambda x: (x[1], x[1][1]), reverse=True)
print("The causal edges ranked based on their ρ and γ values are: ", sortedRankOfCausals, "\n")

###########################################################################################

## Keep the top k ranked causalEdges in a list
#causalToBeOut= sortedRankOfCausals[:k] 
#print("The top", k , "causal edges based on their ρ and γ values to be removed: ", causalToBeOut, "\n")

###########################################################################################

# Keep in remList the first edge/s of the above and compute the modularity of the new partitioning. Only 
# the first item of the above (it may contain 1, 2 or 3 edges depending on the contigency)
remList = []
for e in (list(sortedRankOfCausals[0][0])):
    remList.append(e) 
#print("The edge/s to be removed is/are:", remList, "\n")

for e in remList:
    if type(e)== int or isinstance(e, (int, np.integer)):
        G.remove_edge(remList[0], remList[1])
        break
    else:
        G.remove_edge(e[0],e[1])

partitionF = Louvain(G)
#partitionF, QF = Louvain(G)
#print("Modularity of new Louvain partioning is:", QF, "\n")

###########################################################################################

# For connecting and visualizing with Gephi api
#stream = streamer.Streamer(streamer.GephiWS(hostname="localhost",port=8080,workspace="workspace13"))
#for source, target in G.edges():   
#    node_source = graph.Node(source, size=50, community=partitionF[source], label=G.nodes[source]['old_labels'])
#    node_target = graph.Node(target, size=50, community=partitionF[target], label=G.nodes[target]['old_labels'])
#    stream.add_node(node_source,node_target)
#    # time.sleep(0.5) # Make it slower
#    stream.add_edge(graph.Edge(node_source,node_target, directed=False))
#time.sleep(1) #It might be possible the script runs too fast and last action arn't sent properly

###########################################################################################









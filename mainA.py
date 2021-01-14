#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 16:24:26 2020

@author: georgiabaltsou
"""

import numpy as np
from collections import defaultdict
from gephistreamer import graph
from gephistreamer import streamer
import itertools
import operator
import time
import sys
import copy
import csv
from itertools import combinations
from collections import Counter
from loadGraph import loadGraph
from Louvain import Louvain
from FindEndogenousViaRandomWalks import endogenous
from rankMetricM import rankMetricM
from rankMetricRC import rankMetricRC
from ResAndDis import resAndDis

##########################################################################################
def inComm(node, comm, G):
    count = 0
    for i in comm:
        if G.has_edge(node, i):
            count += 1
    return count
##########################################################################################
    
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
# For connecting and visualizing with Gephi api
#stream = streamer.Streamer(streamer.GephiWS(hostname="localhost",port=8080,workspace="workspace1"))
#for source, target in G.edges():   
#    node_source = graph.Node(source, size=50, community=partition[source], label=G.nodes[source]['old_labels'])
#    node_target = graph.Node(target, size=50, community=partition[target], label=G.nodes[target]['old_labels'])
#    stream.add_node(node_source,node_target)
#    stream.add_edge(graph.Edge(node_source,node_target, directed=False))
#time.sleep(1) #It might be possible the script runs too fast and last action arn't sent properly
##########################################################################################

# communities_dict contains the communities and nodes of each community
communities = set(partition.values())
communities_dict = {c: [k for k, v in partition.items() if v == c] for c in communities}

# The initial query node
initnode = 78100

# The number of highest in-community degree values that I'll keep the nodes of to check for causals. E.g. if highNodes=3 means that I'll keep the 3 nodes with highest in-community degree to check for edge adding.
highNodes = 4 

# Find the community C that node n belongs
C = partition.get(initnode)

##########################################################################################

# Find the communities that initial node has heigbours in
neighboringComms = []
for n in G.neighbors(initnode):
    if partition.get(n) not in neighboringComms:
        neighboringComms.append(partition.get(n))
if C in neighboringComms: # remove the community of initial node
    neighboringComms.remove(C)
#print("neighboringComms=", neighboringComms)
lenComms = len(neighboringComms) # Find how many are the neighbouring communities

##########################################################################################
        
# Find the highest degree nodes of the neighbouring to initial node communities
highest_degree = defaultdict(dict)
allNodesComms = defaultdict(dict)       
for k in range(lenComms):
    temp = defaultdict(dict)
    maxVal = []
    for i in communities_dict[neighboringComms[k]]:  
        temp[i] = inComm(i, communities_dict[neighboringComms[k]], G)        
#    highest_degree[neighboringComms[k]] = max(temp.items(), key=operator.itemgetter(1))[0] # this is to return only one key whose value is max, not all with max values   
#        if temp[i] not in maxVal:
#            maxVal.append(temp[i])
#    maxVal = sorted(maxVal, reverse=True)
#    print(maxVal)
#    allNodesComms[neighboringComms[k]] = sorted(temp.items(), key=lambda x: (x[1]), reverse=True)
    
    # keep in mcf the highNodes top ranked nodes that have highest in-community degree
    c = Counter(temp)
    mc = c.most_common(highNodes)
    mcf = [x[0] for x in mc]
#    print("mc=", mcf)
    
    highest_degree[neighboringComms[k]] = mcf
#    highest_degree[neighboringComms[k]] = [keys for keys,values in temp.items() if values == max(temp.values())]    

print("The neighbouring communities with their corresponding highest in-community degree nodes are:", dict(highest_degree))


if len(dict(highest_degree))> 1:
    # User selects which community to check
    commToCheck = input("Choose from the communities above which to check for causalities: ") 
else:
    commToCheck = list(dict(highest_degree).keys())[0]

# Keep the nodes found before from the selected community
nodesMax = highest_degree[int(commToCheck)] # keep only the max degree in-community nodes

##########################################################################################
# Check if there are edges connecting initial node to nodesMax and add them in finalCausalNodes
finalCausalNodes = []
CausalEdges = []
for node in nodesMax:
    if G.has_edge(initnode, node):
        continue
    else:
        finalCausalNodes.append(node)
        CausalEdges.append((initnode, node))

##########################################################################################
           
# CausalEdges contains the edges of endogenous that has been selected to be examined
# as causal.
print("The edges that will be examined as causal are: ", CausalEdges, "\n")

# Find causals only as unique, douples or triples
cardinality = lenComms
if cardinality > 3:
    cardinality = 3

# Dictionary that contains the causal edges ranked based on their ρ and γ values
rankOfCausals = defaultdict(dict)

for cardinality in range(1,4):
    contigencyLen = cardinality-1
    
    if cardinality == 1:        
        # Use a deep copy of the initial network to test the edge removals    
        testG = copy.deepcopy(G)          
        for e in CausalEdges:
            testG.add_edge(e[0], e[1])
        
            partition2 = Louvain(testG)
#            partition2, Q2 = Louvain(testG)
            if partition2[initnode] != partition[initnode]:
                # Means that the query node changed community      
                rankOfCausals = resAndDis(partition, partition2, initnode, contigencyLen, nodesLen, rankOfCausals, e)
                testG.remove_edge(e[0], e[1])
            else:
                testG.remove_edge(e[0], e[1])
                continue                          
#        print("The edges found with cardinality = 1 are: ", rankOfCausals, "\n")  

        
    if cardinality == 2:
        testG = copy.deepcopy(G)          
#        rankOfCausals = defaultdict(dict)        
        # Keep in comb all the possible combinations of CausalEdges tuples with max length = cardinality
        comb = list(combinations(CausalEdges, cardinality))
        for e in comb:
            testG.add_edge(e[0][0], e[0][1])
            testG.add_edge(e[1][0], e[1][1])
            partition2 = Louvain(testG)            
#            partition2, Q2 = Louvain(testG)
            if partition2[initnode] != partition[initnode]:
                # Means that the query node changed community
                rankOfCausals = resAndDis(partition, partition2, initnode, contigencyLen, nodesLen, rankOfCausals, e)
                testG.remove_edge(e[0][0], e[0][1])
                testG.remove_edge(e[1][0], e[1][1])
            else:
                testG.remove_edge(e[0][0], e[0][1])
                testG.remove_edge(e[1][0], e[1][1])
                continue   
#        print("The edges found with cardinality = 2 are: ", rankOfCausals, "\n")

        
    if cardinality == 3:
            testG = copy.deepcopy(G)               
#            rankOfCausals = defaultdict(dict)  
            # Keep in comb all the possible combinations of CausalEdges tuples with max length = cardinality
            comb = list(combinations(CausalEdges, cardinality))
            for e in comb:
                testG.add_edge(e[0][0], e[0][1])
                testG.add_edge(e[1][0], e[1][1])
                testG.add_edge(e[2][0], e[2][1])
                partition2 = Louvain(testG)                
#                partition2, Q2 = Louvain(testG)
                if partition2[initnode] != partition[initnode]:
                    # Means that the query node changed community
                    rankOfCausals = resAndDis(partition, partition2, initnode, contigencyLen, nodesLen, rankOfCausals, e)
                    testG.remove_edge(e[0][0], e[0][1])
                    testG.remove_edge(e[1][0], e[1][1])
                    testG.remove_edge(e[2][0], e[2][1])    
                else:
                    testG.remove_edge(e[0][0], e[0][1])
                    testG.remove_edge(e[1][0], e[1][1])
                    testG.remove_edge(e[2][0], e[2][1])
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
addList = []
for e in (list(sortedRankOfCausals[0][0])):
    addList.append(e) 
#print("The edge/s to be addeed is/are:", addList, "\n")

for e in addList:
    if type(e)== int or isinstance(e, (int, np.integer)):
        G.add_edge(addList[0], addList[1])
        break
    else:
        G.add_edge(e[0],e[1])

partitionF = Louvain(G)
#partitionF, QF = Louvain(G)
#print("Modularity of new Louvain partioning is:", QF, "\n")

###########################################################################################

## For connecting and visualizing with Gephi api
#stream = streamer.Streamer(streamer.GephiWS(hostname="localhost",port=8080,workspace="workspace3"))
#for source, target in G.edges():   
#    node_source = graph.Node(source, size=50, community=partitionF[source], label=G.nodes[source]['old_labels'])
#    node_target = graph.Node(target, size=50, community=partitionF[target], label=G.nodes[target]['old_labels'])
#    stream.add_node(node_source,node_target)
#    # time.sleep(0.5) # Make it slower
#    stream.add_edge(graph.Edge(node_source,node_target, directed=False))
#time.sleep(1) #It might be possible the script runs too fast and last action arn't sent properly

###########################################################################################









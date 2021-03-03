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
import igraph as ig
from initPartition import initPartition
from jaccard import jaccard
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score

##########################################################################################
def inComm(node, comm, G):
    count = 0
    for i in comm:
        if i in G.neighbors(node):
            count += 1
    return count
##########################################################################################
# The edgelist file  
#edgeFile = "datasets/publEdgelist.csv"
edgeFile = "datasets/polbooksEdgelist.csv" #105 nodes 441 edges

# The number of highest in-community degree values that I'll keep the nodes of to check for causals. E.g. if highNodes=3 means that I'll keep the 3 nodes with highest in-community degree to check for edge adding.
highNodes = 4    

# Load the graph
G = loadGraph(edgeFile)
ig.summary(G)

# Number of nodes of the initial G 
nodesLen = G.vcount()
#print(nodesLen)


# Find the initial partitioning from the given comm file
partitionI, memb, changedPartI = initPartition(G)


Gx = G.to_networkx()
 
#print(type(partitionI[0]))
###########################################################################################
## For connecting and visualizing with Gephi api
#stream = streamer.Streamer(streamer.GephiWS(hostname="localhost",port=8080,workspace="workspace3"))
#for source, target in Gx.edges():   
##    node_source = graph.Node(source, size=50, community=partitionI[source]) #for Aminer
##    node_target = graph.Node(target, size=50, community=partitionI[target]) #for Aminer
#    node_source = graph.Node(source, size=20, community=changedPartI[source])
#    node_target = graph.Node(target, size=20, community=changedPartI[target]) 
##    node_source = graph.Node(source, size=50, community=partition[source], label=G.nodes[source]['old_labels'])
##    node_target = graph.Node(target, size=50, community=partition[target], label=G.nodes[target]['old_labels'])
#    stream.add_node(node_source,node_target)
#    stream.add_edge(graph.Edge(node_source,node_target, directed=True))
#time.sleep(1) #It might be possible the script runs too fast and last action arn't sent properly
##########################################################################################

# The initial query node
#150, 1030, 1312 works for cora
initnode = 5
 
# Find the community C that node n belongs
C = memb[initnode]
# Find all the nodes of community C (n belongs to C)
community = [i for i, e in enumerate(memb) if e == C]


# communities_dict contains the communities and nodes of each community
communities = set(changedPartI.values())
communities_dict = {c: [k for k, v in changedPartI.items() if v == c] for c in communities}


##########################################################################################

# Find the communities that initial node has neigbours in
neighboringComms = []
for n in G.neighbors(initnode):
    if changedPartI.get(n) not in neighboringComms:
        neighboringComms.append(changedPartI.get(n))
if C in neighboringComms: # remove the community of initial node
    neighboringComms.remove(C)
#print("neighboringComms=", neighboringComms)
lenComms = len(neighboringComms) # Find how many are the neighbouring communities

#sys.exit()
##########################################################################################
        
# Find the highest degree nodes of the neighbouring to initial node communities
highest_degree = defaultdict(dict)
allNodesComms = defaultdict(dict)       
for k in range(lenComms):
    temp = defaultdict(dict)
    maxVal = []
    for i in communities_dict[neighboringComms[k]]:  
        temp[i] = inComm(i, communities_dict[neighboringComms[k]], G)        
#    print(temp)
    # keep in mcf the highNodes top ranked nodes that have highest in-community degree
    c = Counter(temp)
    mc = c.most_common(highNodes)
    mcf = [x[0] for x in mc]
#    print("mc=", mcf)
    
    highest_degree[neighboringComms[k]] = mcf
#    highest_degree[neighboringComms[k]] = [keys for keys,values in temp.items() if values == max(temp.values())]    

print("The neighbouring communities with their corresponding highest in-community degree nodes are:", dict(highest_degree))


if len(dict(highest_degree)) == 0:
    print("###################################################")
    print("No neighoring communities with this query node. The process stops")
    print("###################################################\n")
    sys.exit()


if len(dict(highest_degree))> 1:
    # User selects which community to check
    commToCheck = input("Choose from the communities above which to check for causalities: ") 
else:
    commToCheck = list(dict(highest_degree).keys())[0]

# Keep the nodes found before from the selected community
nodesMax = highest_degree[int(commToCheck)] # keep only the max degree in-community nodes

#sys.exit()
##########################################################################################
# Check if there are edges connecting initial node to nodesMax and add them in finalCausalNodes
finalCausalNodes = []
CausalEdges = []
for node in nodesMax:
    if node in G.neighbors(initnode): # if it already has an edge I continue as I cannot add edge
        continue
    else: # if it doesn't have edge, I add the corresponding node to finalCausalNodes and edge in CausalEdges
        finalCausalNodes.append(node)
        CausalEdges.append((node, initnode)) # I add the edge that points the initial node
print("Causal edges found.\n")
#sys.exit()
##########################################################################################
           
# CausalEdges contains the edges of endogenous that has been selected to be examined
# as causal.
print("The edges that will be examined as causal are: ", CausalEdges, "\n")
#sys.exit()
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


#        testG.add_edges([(5, 8)])
#        testG.add_edges([(1, 8)])
#        testG.add_edges([(5, 12)])
#        testG.add_edges([(1, 12)])
#        testG.add_edges([(811, 4037)])
#        testG.add_edges([(811, 1222)])

#        testG.add_edges([(1323, 811)])
#        testG.add_edges([(4037, 811)])
#        testG.add_edges([(1222, 811)])

#        testG.delete_edges([(0, 1)])
#        testG.delete_edges([(1, 5)])


#        partition2, memb2, changedPart2 = Louvain(testG, memb)
#        # Find the community C that node n belongs
#        C2 = memb2[initnode]
##        C3 = memb2[3942]
#        # Find all the nodes of community C (n belongs to C)
#        community2 = [i for i, e in enumerate(memb2) if e == C2]
##        community3 = [i for i, e in enumerate(memb2) if e == C3]
#         
#        for k in range(0,165):
#            communityI = [i for i, e in enumerate(memb) if e == k]
#            if jaccard(communityI,community2) > 0.5:
#                print(k)
##        
#        sys.exit()






        
        for e in CausalEdges:
            testG.add_edges([(e[0], e[1])])            
            partition2, memb2, changedPart2 = Louvain(testG, memb)
            # Find the community C that node n belongs
            C2 = memb2[initnode]
            # Find all the nodes of community C (n belongs to C)
            community2 = [i for i, e in enumerate(memb2) if e == C2]
            
            if (jaccard(community,community2) < 0.5):
                rankOfCausals[e] = 1/(1+contigencyLen)
                testG.delete_edges([(e[0], e[1])])
            else:
                testG.delete_edges([(e[0], e[1])])
                continue 
            
                          
#        print("The edges found with cardinality = 1 are: ", rankOfCausals, "\n")  

        
    if cardinality == 2:
        testG = copy.deepcopy(G)          
#        rankOfCausals = defaultdict(dict)        
        # Keep in comb all the possible combinations of CausalEdges tuples with max length = cardinality
        comb = list(combinations(CausalEdges, cardinality))
        for e in comb:
            testG.add_edges([(e[0][0], e[0][1])])
            testG.add_edges([(e[1][0], e[1][1])])
            
            partition2, memb2, changedPart2 = Louvain(testG,memb) 
            # Find the community C that node n belongs
            C2 = memb2[initnode]
            # Find all the nodes of community C (n belongs to C)
            community2 = [i for i, e in enumerate(memb2) if e == C2]
#            partition2, Q2 = Louvain(testG)
#            if partition2[initnode] != partitionI[initnode]:
##            if memb[initnode] != memb2[initnode]:
            if (jaccard(community,community2) < 0.5):
                # Means that the query node changed community
#                rankOfCausals = resAndDis(memb, memb2, initnode, contigencyLen, nodesLen, rankOfCausals, e)
                rankOfCausals[e] = 1/(1+contigencyLen)
                testG.delete_edges([(e[0][0], e[0][1])])
                testG.delete_edges([(e[1][0], e[1][1])])
            else:
                testG.delete_edges([(e[0][0], e[0][1])])
                testG.delete_edges([(e[1][0], e[1][1])])
                continue   
#        print("The edges found with cardinality = 2 are: ", rankOfCausals, "\n")

        
    if cardinality == 3:
            testG = copy.deepcopy(G)               
#            rankOfCausals = defaultdict(dict)  
            # Keep in comb all the possible combinations of CausalEdges tuples with max length = cardinality
            comb = list(combinations(CausalEdges, cardinality))
            for e in comb:
                testG.add_edges([(e[0][0], e[0][1])])
                testG.add_edges([(e[1][0], e[1][1])])
                testG.add_edges([(e[2][0], e[2][1])]) 
                partition2, memb2, changedPart2 = Louvain(testG,memb)
                # Find the community C that node n belongs
                C2 = memb2[initnode]
                # Find all the nodes of community C (n belongs to C)
                community2 = [i for i, e in enumerate(memb2) if e == C2]
#                partition2, Q2 = Louvain(testG)
#                if partition2[initnode] != partitionI[initnode]:
##                if memb[initnode] != memb2[initnode]:
                if (jaccard(community,community2) < 0.5):
                    # Means that the query node changed community
#                    rankOfCausals = resAndDis(memb, memb2, initnode, contigencyLen, nodesLen, rankOfCausals, e)
                    rankOfCausals[e] = 1/(1+contigencyLen)
                    testG.delete_edges([(e[0][0], e[0][1])])
                    testG.delete_edges([(e[1][0], e[1][1])])
                    testG.delete_edges([(e[2][0], e[2][1])])    
                else:
                    testG.delete_edges([(e[0][0], e[0][1])])
                    testG.delete_edges([(e[1][0], e[1][1])])
                    testG.delete_edges([(e[2][0], e[2][1])])
                    continue                                          
#            print("The edges found with cardinality = 3 are: ", rankOfCausals, "\n")

# Case when query node does not change community
if len(rankOfCausals) == 0:
    print("###################################################")
    print("No change with this query node. The process stops")
    print("###################################################\n")
    sys.exit()
# Sort rankOfCausals based on responsibility values
sortedRankOfCausals = sorted(rankOfCausals.items(), key=lambda x: (x[1], x[1]), reverse=True)
print("The causal edges ranked based on their ρ values are: ", sortedRankOfCausals, "\n")
#sys.exit()

###########################################################################################

## Keep the top k ranked causalEdges in a list
#causalToBeOut= sortedRankOfCausals[:k] 
#print("The top", k , "causal edges based on their ρ and γ values to be removed: ", causalToBeOut, "\n")

###########################################################################################

# Keep in remList the first edge/s of the above and compute the new partitioning. Only 
# the first item of the above (it may contain 1, 2 or 3 edges depending on the contigency)
addList = []
for e in (list(sortedRankOfCausals[0][0])):
    addList.append(e) 
#print("The edge/s to be addeed is/are:", addList, "\n")

flag = False # flag used to understand if I have empty contigency or not
for e in addList:
    if type(e)== int or isinstance(e, (int, np.integer)):
        G.add_edges([(addList[0], addList[1])])
        ed = (addList[0], addList[1])
        cont = []
        print("The actual/counterfactual cause is:", ed, "with contigency set: ", cont)

        break
    else:
        flag = True
        G.add_edges([(e[0], e[1])])

if flag == True:
    cont = addList
    if len(cont) == 2:
        print("The actual/counterfactual cause is:",addList[0], "with contigency set: ", addList[1])
    else:
        print("The actual cause is:",addList[0], "with contigency set[:", addList[1],",", addList[2],"]")

partitionF, membF, changedPartF = Louvain(G, memb)
ig.summary(G)
GxF = G.to_networkx()

###########################################################################################

## For connecting and visualizing with Gephi api
#stream = streamer.Streamer(streamer.GephiWS(hostname="localhost",port=8080,workspace="workspace5"))
#for source, target in GxF.edges(): 
#    node_source = graph.Node(source, size=30, community=changedPartF[source]) 
#    node_target = graph.Node(target, size=30, community=changedPartF[target])
#    #node_source = graph.Node(source, size=50, community=partitionF[source], label=G.nodes[source]['old_labels'])
#    #node_target = graph.Node(target, size=50, community=partitionF[target], label=G.nodes[target]['old_labels'])
#    stream.add_node(node_source,node_target)
#    stream.add_edge(graph.Edge(node_source,node_target, directed=True))
#time.sleep(1) #It might be possible the script runs too fast and last action arn't sent properly

###########################################################################################

## In order to find out if the community that the initial node has moved to is the same as the one before it's movement.
CF = membF[initnode]
communityF = [i for i, e in enumerate(membF) if e == CF]

for k in range(0,4942):
    communityI = [i for i, e in enumerate(memb) if e == k]
    if jaccard(communityI,communityF) > 0.5:
        print(k)


print("NMI=", normalized_mutual_info_score(memb,membF))
print("ARI=", adjusted_rand_score(memb,membF))
print("RI=", rand_score(memb,membF))
print("AMI=", adjusted_mutual_info_score(memb,membF))
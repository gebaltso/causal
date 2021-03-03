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
import igraph as ig
from initPartition import initPartition
import cdlib
import random
from jaccard import jaccard
import networkx as nx
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score


edgeFile = "datasets/polbooksEdgelist.csv" #105 nodes 441 edges
#edgeFile = "datasets/aminer.csv" # 19720 nodes 26210 edges

    
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
#stream = streamer.Streamer(streamer.GephiWS(hostname="localhost",port=8080,workspace="workspace54"))
#for source, target in Gx.edges():   
##    node_source = graph.Node(source, size=50, community=partitionI[source]) #for Aminer
##    node_target = graph.Node(target, size=50, community=partitionI[target]) #for Aminer
#    node_source = graph.Node(source, size=20, community=changedPartI[source])
#    node_target = graph.Node(target, size=20, community=changedPartI[target]) 
##    node_source = graph.Node(source, size=50, community=partition[source], label=G.nodes[source]['old_labels'])
##    node_target = graph.Node(target, size=50, community=partition[target], label=G.nodes[target]['old_labels'])
#    stream.add_node(node_source,node_target)
#    stream.add_edge(graph.Edge(node_source,node_target, directed=False))
#time.sleep(1) #It might be possible the script runs too fast and last action arn't sent properly
##########################################################################################

# The initial query node
#150, 1030, 1312 works for cora
initnode = 68
 
# select how to compute the endogenous set
method = 'incident' # Select among 'exo' for exogenously given edges, 'random' for random walk, 
#'incident' for all edges incident to query node, 'incomm' for edges inside the same community as the query node
# random is the most time consuming

# The number of items to examine as causal (selected from endogenous)
nOfCausal = 6  # This is not used when I choose to check all the edges as causals
rankMetric = 'none' # Select among 'emb' for embeddedness, 'rc' for relative commitment, 'none' for considering all the endogenous as possible causal ones

# After the computation of ρ and γ keep the top k as causals to be removed
#k = 2 

#sys.exit()

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
    if initnode in endoNodes:
        endoNodes.remove(initnode) # remove the initial node as G.edges include it

    
elif method == 'incident':
#    endoEdges = list(G.in_edges(initnode))
#    endoEdgesInit = G.incident(initnode, mode ='out') # change mode to in for in coming edges
    endoEdgesInit = G.incident(initnode)
    endoEdges = [edge.tuple for edge in G.es[endoEdgesInit]]
    endoNodes = set(itertools.chain.from_iterable(list(endoEdges)))
    if initnode in endoNodes:
        endoNodes.remove(initnode) # remove the initial node as G.edges include it
    

elif method == 'incomm':
    # Find the community C that node n belongs
    C = memb[initnode]
    # Find all the nodes of community C (n belongs to C)
    community = [i for i, e in enumerate(memb) if e == C]
    
    endoNodes = []
    for i in community:
        if i in G.neighbors(initnode):
            endoNodes.append(i)
#    endoEdgesInit = G.incident(initnode, mode ='out') # change mode to in for in coming edges
    endoEdgesInit = G.incident(initnode)
    endoEdges = [edge.tuple for edge in G.es[endoEdgesInit]]

    if initnode in endoNodes:
        endoNodes.remove(initnode) # remove the initial node as G.edges include it

        
print("Endogenous edges found.\n")
print("Endogenous nodes:", endoNodes)
#sys.exit()
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
        possibleCausal[i] = rankMetricM(i, partitionI, memb, G)
    sortedPossibleCausal = {k: v for k, v in sorted(possibleCausal.items(), key=lambda item: item[1], reverse=True)}
    finalCausalNodes = list(sortedPossibleCausal)[0:nOfCausal]
    for e in list(endoEdges):
        if (e[0] in finalCausalNodes) or (e[1] in finalCausalNodes):
            if (e[1],e[0]) not in CausalEdges: # this is to avoid adding both (a,b) & (b,a)
                CausalEdges.append(e)

    
elif rankMetric == 'rc':
    for i in endoNodes:
        possibleCausal[i] = rankMetricRC(i, partitionI, memb, G)
#    print("possibleCausal", possibleCausal)
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

# CausalEdges contains the edges of endogenous that has been selected to be examined
# as causal due to rankMetric (emb, rc or none)
print("Causal edges found.\n")
print("The edges that will be examined as causal are: ", CausalEdges, "\n")
#sys.exit()        
###########################################################################################
 
# Find the community C that node n belongs
C = memb[initnode]
# Find all the nodes of community C (n belongs to C)
community = [i for i, e in enumerate(memb) if e == C]


# Dictionary that contains the causal edges ranked based on their ρ (and γ) values
rankOfCausals = defaultdict(dict)


for cardinality in range(1,4):
    contigencyLen = cardinality-1
    
    if cardinality == 1:        
        # Use a deep copy of the initial network to test the edge removals    
        testG = copy.deepcopy(G)
               
#        testG.delete_edges([(1, 0)])
#        testG.delete_edges([(1, 5)])
#        testG.delete_edges([(6, 7)])
#        testG.delete_edges([(5, 2)])
#        testG.delete_edges([(1, 6)])
#        testG.delete_edges([(145, 811)])
##        testG.delete_edges([(144, 811)])
##        testG.delete_edges([(603, 6543)])#            testG.delete_edges([(717,738)])
#
#        partition2, memb2, changedPart2 = Louvain(testG, memb)
#        # Find the community C that node n belongs
#        C2 = memb2[initnode]
##        C3 = memb2[3942]
#        # Find all the nodes of community C (n belongs to C)
#        community2 = [i for i, e in enumerate(memb2) if e == C2]
##        community3 = [i for i, e in enumerate(memb2) if e == C3]
#         
#        for k in range(0,4942):
#            communityI = [i for i, e in enumerate(memb) if e == k]
#            if jaccard(communityI,community2) > 0.5:
#                print(k)
##        
#        sys.exit()
#        
        
        
       
        
        
        for e in CausalEdges:

            testG.delete_edges([(e[0], e[1])])


            partition2, memb2, changedPart2 = Louvain(testG, memb)
            # Find the community C that node n belongs
            C2 = memb2[initnode]
            # Find all the nodes of community C (n belongs to C)
            community2 = [i for i, e in enumerate(memb2) if e == C2]

#            myList = []
#            for i in community:
#                if i not in community2:
#                    myList.append(i)
#            print(myList)
#            print("cont1", jaccard(community,community2))
            
            if (jaccard(community,community2) < 0.5):
                rankOfCausals[e] = 1/(1+contigencyLen)
                testG.add_edges([(e[0], e[1])])
            else:
                testG.add_edges([(e[0], e[1])])
                continue          
            
#            print("Sim=", jaccard(community,community2))
#            
#            sys.exit()

#            if partition2[initnode] != partitionI[initnode]:
##            if memb[initnode] != memb2[initnode]:
##                # Means that the query node changed community      
###                rankOfCausals = resAndDis(memb, memb2, initnode, contigencyLen, nodesLen, rankOfCausals, e)
# #               rankOfCausals[e] = 1/(1+contigencyLen)
##                testG.add_edges([(e[0], e[1])])
# #           else:
# #               testG.add_edges([(e[0], e[1])])
##                continue                          
#        print("The edges found with cardinality = 1 are: ", dict(rankOfCausals), "\n")  

#    sys.exit()
        
    if cardinality == 2:
        testG = copy.deepcopy(G) 
                 
#        rankOfCausals = defaultdict(dict)        
        # Keep in comb all the possible combinations of CausalEdges tuples with max length = cardinality
        comb = list(combinations(CausalEdges, cardinality))
        for e in comb:
            testG.delete_edges([(e[0][0], e[0][1])])
            testG.delete_edges([(e[1][0], e[1][1])])
#            testG.remove_edge(e[0][0], e[0][1])
#            testG.remove_edge(e[1][0], e[1][1])
#            testg = ig.Graph.from_networkx(testG)
#            testg = cdlib.utils.convert_graph_formats(testG, ig.Graph, directed=True)
            
#            testg.vs['id'] = list(range(nodesLen))
            
            
            partition2, memb2, changedPart2 = Louvain(testG,memb) 
            # Find the community C that node n belongs
            C2 = memb2[initnode]
            # Find all the nodes of community C (n belongs to C)
            community2 = [i for i, e in enumerate(memb2) if e == C2]
#            partition2, Q2 = Louvain(testG)
#            if partition2[initnode] != partitionI[initnode]:
##            if memb[initnode] != memb2[initnode]:
            
#            print("cont2", jaccard(community,community2))
            
            if (jaccard(community,community2) < 0.5):
                # Means that the query node changed community
#                rankOfCausals = resAndDis(memb, memb2, initnode, contigencyLen, nodesLen, rankOfCausals, e)
                rankOfCausals[e] = 1/(1+contigencyLen)
                testG.add_edges([(e[0][0], e[0][1])])
                testG.add_edges([(e[1][0], e[1][1])])
            else:
                testG.add_edges([(e[0][0], e[0][1])])
                testG.add_edges([(e[1][0], e[1][1])])
                continue   
#        print("The edges found with cardinality = 2 are: ", rankOfCausals, "\n")

        
    if cardinality == 3:
            testG = copy.deepcopy(G)
                           
#            rankOfCausals = defaultdict(dict)  
            # Keep in comb all the possible combinations of CausalEdges tuples with max length = cardinality
            comb = list(combinations(CausalEdges, cardinality))
            for e in comb:
                testG.delete_edges([(e[0][0], e[0][1])])
                testG.delete_edges([(e[1][0], e[1][1])])
                testG.delete_edges([(e[2][0], e[2][1])])
#                testg = ig.Graph.from_networkx(testG)
#                testg = cdlib.utils.convert_graph_formats(testG, ig.Graph, directed=True)
#                
#                testg.vs['id'] = list(range(nodesLen))
                
                
                partition2, memb2, changedPart2 = Louvain(testG,memb)
                # Find the community C that node n belongs
                C2 = memb2[initnode]
                # Find all the nodes of community C (n belongs to C)
                community2 = [i for i, e in enumerate(memb2) if e == C2]
#                partition2, Q2 = Louvain(testG)
#                if partition2[initnode] != partitionI[initnode]:
##                if memb[initnode] != memb2[initnode]:
                
#                print("cont3", jaccard(community,community2))
                
                if (jaccard(community,community2) < 0.5):
                    # Means that the query node changed community
#                    rankOfCausals = resAndDis(memb, memb2, initnode, contigencyLen, nodesLen, rankOfCausals, e)
                    rankOfCausals[e] = 1/(1+contigencyLen)
                    testG.add_edges([(e[0][0], e[0][1])])
                    testG.add_edges([(e[1][0], e[1][1])])
                    testG.add_edges([(e[2][0], e[2][1])])    
                else:
                    testG.add_edges([(e[0][0], e[0][1])])
                    testG.add_edges([(e[1][0], e[1][1])])
                    testG.add_edges([(e[2][0], e[2][1])])
                    continue                       
#            print("The edges found with cardinality = 3 are: ", rankOfCausals, "\n")

# Case when query node does not change community
if len(rankOfCausals) == 0:
    print("###################################################")
    print("No change with this query node. The process stops")
    print("###################################################\n")
    sys.exit()
    
# Sort rankOfCausals based on responsibility values
#sortedRankOfCausals = sorted(rankOfCausals.items(), key=lambda x: (x[1], x[1][1]), reverse=True)
sortedRankOfCausals = sorted(rankOfCausals.items(), key=lambda x: (x[1], x[1]), reverse=True)
print("The causal edges ranked based on their ρ value are: ", sortedRankOfCausals, "\n")

###########################################################################################

## Keep the top k ranked causalEdges in a list
#causalToBeOut= sortedRankOfCausals[:k] 
#print("The top", k , "causal edges based on their ρ and γ values to be removed: ", causalToBeOut, "\n")

###########################################################################################

#k = copy.deepcopy(G)
#testk = cdlib.utils.convert_graph_formats(k, ig.Graph, directed=True)
#partitionK = Louvain(testk,memb)


# Keep in remList the first edge/s of the above and compute the new partitioning. Only 
# the first item of the above (it may contain 1, 2 or 3 edges depending on the contigency)
remList = []
for e in (list(sortedRankOfCausals[0][0])):
    remList.append(e) 
#print("The edge/s to be removed is/are:", remList, "\n")

flag = False # flag used to understand if I have empty contigency or not
for e in remList:
    if type(e)== int or isinstance(e, (int, np.integer)):
#        G.remove_edge(remList[0], remList[1])
        G.delete_edges([(remList[0], remList[1])])
        ed = (remList[0], remList[1])
        cont = []
        print("The actual/counterfactual cause is:", ed, "with contigency set: ", cont)
        break
    else:
        flag = True
#        G.remove_edge(e[0],e[1])
        G.delete_edges([(e[0], e[1])])
 
if flag == True:
    cont = remList
    if len(cont) == 2:
        print("The actual/counterfactual cause is:",remList[0], "with contigency set: ", remList[1])
    else:
        print("The actual cause is:",remList[0], "with contigency set[:", remList[1],",", remList[2],"]")
        
#g = ig.Graph.from_networkx(G)
#finalg = cdlib.utils.convert_graph_formats(G, ig.Graph, directed=True)
#
#finalg.vs['id'] = list(range(nodesLen))



partitionF, membF, changedPartF = Louvain(G, memb)
ig.summary(G)
GxF = G.to_networkx()
#partitionF, QF = Louvain(G)
#print("Modularity of new Louvain partioning is:", QF, "\n")

###########################################################################################

## For connecting and visualizing with Gephi api
#stream = streamer.Streamer(streamer.GephiWS(hostname="localhost",port=8080,workspace="workspace52"))
#for source, target in GxF.edges(): 
#    node_source = graph.Node(source, size=30, community=changedPartF[source]) 
#    node_target = graph.Node(target, size=30, community=changedPartF[target])
#    #node_source = graph.Node(source, size=50, community=partitionF[source], label=G.nodes[source]['old_labels'])
#    #node_target = graph.Node(target, size=50, community=partitionF[target], label=G.nodes[target]['old_labels'])
#    stream.add_node(node_source,node_target)
#    stream.add_edge(graph.Edge(node_source,node_target, directed=False))
#time.sleep(1) #It might be possible the script runs too fast and last action arn't sent properly
###########################################################################################

## In order to find out if the community that the initial node has moved to is the same as the one before it's movement.
CF = membF[initnode]
communityF = [i for i, e in enumerate(membF) if e == CF]

for k in range(0,4942):
    communityI = [i for i, e in enumerate(memb) if e == k]
    if jaccard(communityI,communityF) > 0.5:
        print(k)


#for v in nx.optimize_graph_edit_distance(Gx,GxF):
#    minv = v
#
#print("Edit distance:", minv)

print("NMI=", normalized_mutual_info_score(memb,membF))
print("ARI=", adjusted_rand_score(memb,membF))
print("RI=", rand_score(memb,membF))
print("AMI=", adjusted_mutual_info_score(memb,membF))
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 11:06:41 2020

@author: georgia
"""

import community as community_louvain
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import sys
import csv
import igraph as ig
import louvain
from collections import defaultdict
import leidenalg



def Louvain(G, memb):
        
    
#    optimiser = leidenalg.Optimiser()

#    g0 = ig.Graph.GRG(100, 0.2)
#    g0.vs['id'] = list(range(100))
#    g1 = ig.Graph.GRG(100, 0.1)
#    g1.vs['id'] = list(range(100))
    
#    membership, improv = leidenalg.find_partition_temporal([G0, G1], leidenalg.CPMVertexPartition)    
    
    
#    for i in (membership[0]):
#        if membership[0][i] != membership[1][i]:
#            print("Nodes that changed comm:", i)
    
#    print(membership)
            
            
#    layers, interslice_layer, G_full = leidenalg.time_slices_to_layers([G0, G1])
#    partitions = [leidenalg.CPMVertexPartition(H, node_sizes='node_size')for H in layers]
#    
#    interslice_partition = leidenalg.CPMVertexPartition(interslice_layer, resolution_parameter=0, node_sizes='node_size', weights='weight')
#    diff = optimiser.optimise_partition_multiplex(partitions + [interslice_partition])
#    
# 
#    for i in range(len(partitions)):
#        if partitions[0][i] != partitions[1][i]:
#            print(i)
#    print((partitions[0][10]))
#    print((partitions[1][10]))
    
#    print(membership[0][1312])
#    print(membership[1][1312])



        
#    partition = louvain.find_partition(G0, louvain.ModularityVertexPartition, initial_membership=memb, seed=0)
    
    
    
#    partition = louvain.find_partition(G, louvain.CPMVertexPartition, initial_membership=lisInitialPartition, seed=0, resolution_parameter=0.00001)
 
    partition = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition,initial_membership=memb, seed=0)
#    partition = leidenalg.find_partition(G, leidenalg.CPMVertexPartition,initial_membership=lisInitialPartition, seed=0)

    memb2 = partition.membership 
    
#    membership, improvement = leidenalg.find_partition_temporal(
#                             [G0, G],leidenalg.CPMVertexPartition,interslice_weight=0.1)
    
    
#    layers, interslice_layer, G_full = leidenalg.time_slices_to_layers([G0, G])
#    partitions = [leidenalg.ModularityVertexPartition(H) for H in layers]
#    interslice_partition = leidenalg.ModularityVertexPartition(interslice_layer)
#    
#    
#    print(partitions)
    
#    print(len(partition))
#    print("********************************************")
    
#    optimiser = louvain.Optimiser()
#    diff = optimiser.optimise_partition(partition)
    
#    print(len(partition))    
#    sys.exit()
#    print(partition.quality()) #for modularity
    
    changedPart = defaultdict(dict)
    for i in range(0,len(partition)):
        for j in partition[i]:
            changedPart[j] = i    
    changedPartF = dict(sorted(iter(changedPart.items())))

    
    return  partition, memb2, changedPartF



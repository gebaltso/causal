#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 17:12:51 2021

@author: georgiabaltsou
"""

from collections import defaultdict
import csv
import community as community_louvain
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import sys
import igraph as ig
import louvain
from collections import OrderedDict
import numpy as np
import leidenalg
from cdlib import algorithms, readwrite, viz
import cdlib


def initPartition(G):

#    G = cdlib.utils.convert_graph_formats(g, nx.Graph, directed=True)
#    
#    coms = algorithms.leiden(G)
#    viz.plot_community_graph(G, coms)
    
#    readwrite.write_community_csv(coms, "communities.csv", ",")

 
#    InitialPartition = louvain.find_partition(G, louvain.ModularityVertexPartition, initial_membership=None, seed=0)
#    InitialPartition = louvain.find_partition(G, louvain.CPMVertexPartition, initial_membership=None, seed=0, resolution_parameter=0.001)

    InitialPartition = leidenalg.find_partition(G, leidenalg.ModularityVertexPartition,seed=0)
#    InitialPartition = leidenalg.find_partition(G, leidenalg.CPMVertexPartition,seed=0,resolution_parameter=1) #modularity with resolution parameter


#    InitialPartition = leidenalg.find_partition(G, leidenalg.RBERVertexPartition,seed=0,max_comm_size=290,resolution_parameter=0.002)

    print(InitialPartition.quality()) #for modularity
    
#    partition.set_membership(lisInitialPartition)
    
#    print(InitialPartition)
#    print(InitialPartition[78])
#    for i in range(len(InitialPartition)):
#      print(InitialPartition[i])

#    optimiser = leidenalg.Optimiser()
#    diff = optimiser.optimise_partition(InitialPartition)
#
#    print(diff)
    
    memb = InitialPartition.membership
    

    
    changedPart = defaultdict(dict)
    for i in range(len(InitialPartition)):
        for j in InitialPartition[i]:
            changedPart[j] = i             
    changedPartI = dict(sorted(iter(changedPart.items())))
#    changedPartI =OrderedDict(sorted(changedPart.items()))

#    print(InitialPartition)
    
    
    return InitialPartition, memb, changedPartI

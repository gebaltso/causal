#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 19:10:02 2021

@author: georgiabaltsou
"""

from cdlib import algorithms
import dynetx as dn
import networkx as nx
from collections import defaultdict

dg = dn.DynDiGraph(edge_removal=True)
#for x in range(4):
#    g = nx.erdos_renyi_graph(150, 0.05)
#    dg.add_interactions_from(list(g.edges()), t=x)


edgeFile = "datasets/coraEdgelist.csv"
G = nx.read_edgelist(edgeFile, create_using=nx.DiGraph(), delimiter=" ", encoding='utf-8-sig')
dg.add_interactions_from(list(G.edges()), t=1)
dg.add_interactions_from(list(G.edges()), t=2)
dg.add_interaction(u=1113831, v=35,t=1, e=2)
dg.add_interaction(u=1113831, v=248425,t=1, e=2)


    
coms = algorithms.tiles(dg)

#print(coms.get_observation_ids()) #the timestamps I have

c = coms.get_clustering_at(1)
comunity_ids_m = c.named_communities.values()
communities = defaultdict(dict)

# For one tid only
for i in c.named_communities:
    l = []
    for j in c.named_communities[i]:
        l.append(j)        
    communities[i] = l    
print("comms at 1:",dict(communities), "\n\n\n")

c = coms.get_clustering_at(2)
comunity_ids_m = c.named_communities.values()
communities = defaultdict(dict)

# For one tid only
for i in c.named_communities:
    l = []
    for j in c.named_communities[i]:
        l.append(j)        
    communities[i] = l    
print("comms at 1:",dict(communities), "\n\n\n")


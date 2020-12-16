#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 16:41:34 2020

@author: georgiabaltsou
"""

import community as community_louvain
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from collections import defaultdict
import csv
import sys
import os


def loadGraph():
    
    myFile = "netScEdgelist.csv"
    
#    G = nx.read_edgelist("hdnEdgelist.csv", create_using=nx.Graph(), delimiter=",", encoding='utf-8-sig') #human disease network with 1419 nodes and 2738 edges
#    G = nx.read_edgelist("hdnDiseaseEdgelist.csv", create_using=nx.Graph(), delimiter=",", encoding='utf-8-sig') #human disease network with 516 nodes and 1188 edges
#    G = nx.read_edgelist("yeastEdgelist.csv", create_using=nx.Graph(), delimiter=",", encoding='utf-8-sig') #protein protein interaction network with 2361 nodes and 7182 edges/6646 without self loops
    G = nx.read_edgelist(myFile, create_using=nx.Graph(), delimiter=",", encoding='utf-8-sig') #network scientists coauthorship network with 1461 nodes and 2742 edges    
#    G = nx.karate_club_graph() # 34 nodes 78 edges
#    G = nx.les_miserables_graph() # 77 nodes 254 edges
#    G = nx.read_adjlist('ca-GrQc.txt', delimiter='\t') # 5242 nodes 14496 edges
    
    # To remove self loops
#    G.remove_edges_from(nx.selfloop_edges(G))
    
#    print(G.number_of_nodes())
#    print(G.number_of_edges())
#    
#    sys.exit()
    

    # Check if node ids are integers. If not convert them to integers
    if all(isinstance(n, int) for n in list(G.nodes)):
        print(" ")      
    else:
        start = 0
        G = nx.convert_node_labels_to_integers(G,first_label=start,ordering='sorted', label_attribute='old_labels')

    # Write new and old labels in a file
    with open('Labels/labels_'+str(myFile) + '.csv', 'a') as out_file:
        writer = csv.writer(out_file, delimiter=';')            
        if os.stat('Labels/labels_'+str(myFile) + '.csv').st_size == 0:
            writer.writerow(["Old Label", "New Label"])                
            for n in list(G.nodes()):
                writer.writerow([G.nodes[n]['old_labels'], n]) 


    return G









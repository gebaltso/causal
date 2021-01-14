#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 16:41:34 2020

@author: georgiabaltsou
"""

import networkx as nx
import os
import csv
from collections import defaultdict


def loadGraph():
    
    
    myFile = "datasets/aminerEdgelist.csv"
#    myFile = "datasets/finalB.csv"
    
    if myFile == "datasets/aminerEdgelist.csv":
        p = True
    else:
        p = False
    
#    G = nx.read_edgelist("hdnEdgelist.csv", create_using=nx.Graph(), delimiter=",", encoding='utf-8-sig') #human disease network with 1419 nodes and 2738 edges
#    G = nx.read_edgelist("hdnDiseaseEdgelist.csv", create_using=nx.Graph(), delimiter=",", encoding='utf-8-sig') #human disease network with 516 nodes and 1188 edges
#    G = nx.read_edgelist("yeastEdgelist.csv", create_using=nx.Graph(), delimiter=",", encoding='utf-8-sig') #protein protein interaction network with 2361 nodes and 7182 edges/6646 without self loops
    G = nx.read_edgelist(myFile, create_using=nx.Graph(), delimiter=",", encoding='utf-8-sig') #network scientists coauthorship network with 1461 nodes and 2742 edges    / Largest component 379 nodes 914 edges
#    G = nx.karate_club_graph() # 34 nodes 78 edges
#    G = nx.les_miserables_graph() # 77 nodes 254 edges
#    G = nx.read_adjlist('datasets/ca-GrQc.txt', delimiter='\t') # 5242 nodes 14496 edges
    
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
        if p == True:
            mapping = defaultdict(dict)
            with open(myFile, 'r') as f:
                csv_reader = csv.reader(f, delimiter=',')        
                for row in csv_reader:
                    mapping[row[0]] = int(row[0])
                    mapping[row[1]] = int(row[1]) 
            G = nx.relabel_nodes(G, mapping, copy=False)
        else:
            start = 0
            G = nx.convert_node_labels_to_integers(G,first_label=start,ordering='sorted', label_attribute='old_labels')

#    # Write new and old labels in a file
#    with open('datasets/labels/labels'+ '.csv', 'a') as out_file:
#        writer = csv.writer(out_file, delimiter=';')            
#        if os.stat('datasets/labels/labels' + '.csv').st_size == 0:
#            writer.writerow(["Old Label", "New Label"])                
#            for n in list(G.nodes()):
#                writer.writerow([G.nodes[n]['old_labels'], n]) 
            
                     

    return G, p









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
import copy
import sys
import random


def loadGraph():
    
    G = nx.karate_club_graph() # 34 nodes 78 edges
#    G = nx.les_miserables_graph() # 77 nodes 254 edges
#    G = nx.read_adjlist('ca-GrQc.txt', delimiter='\t') # 5242 nodes 14496 edges
    
#    print(G.number_of_nodes())
#    print(G.number_of_edges())
#    
#    sys.exit()
    
    # Check if node ids are integers. If not convert them to integers
    if all(isinstance(n, int) for n in list(G.nodes)):
        print(" ")
    else:
        start = 0
        G = nx.convert_node_labels_to_integers(G,first_label=start,ordering='sorted')

    
    
    return G
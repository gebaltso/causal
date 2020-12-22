#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 13:23:52 2020

@author: georgiabaltsou
"""

import networkx as nx
from netwulf import visualize

#G = nx.barabasi_albert_graph(100,m=1)
myFile = "hdnEdgelist.csv"
G = nx.read_edgelist(myFile, create_using=nx.Graph(), delimiter=",", encoding='utf-8-sig')
visualize(G)
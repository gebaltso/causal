#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 16:41:34 2020

@author: georgiabaltsou
"""

import networkx as nx
import sys
import os
import csv
from collections import defaultdict
from igraph import *


def loadGraph(edgeFile):
    
    # igraph load edgelist from csv
    G = Graph.Read_Ncol(edgeFile, directed=False)
    
    G.vs["id"] = G.vs["name"]
    
#    print(G.vcount())
#    print(G.ecount())
    
#    nodes = G.vs.indices
#    print([G.vs[4707]['name']])
#    print([G.vs[1674]['name']])
#    print([G.vs[5107]['name']])
#    print([G.vs[4372]['name']])
#    print([G.vs[3870]['name']])
#    print([G.vs[537]['name']])
#    print([G.vs[536]['name']])
#    print([G.vs[533]['name']])
#    print([G.vs[4879]['name']])
#    print([G.vs[3686]['name']])
#    print([G.vs[603]['name']])
#    print([G.vs[11909]['name']])
    
#    print(G.vs.find(name="2801"))
#    print(G.vs.find(name="47957"))

#    sys.exit()
    
#    # Write new and old labels in a file
#    with open('datasets/labels/labels'+ '.csv', 'a') as out_file:
#        writer = csv.writer(out_file, delimiter=';')            
#        if os.stat('datasets/labels/labels' + '.csv').st_size == 0:
#            writer.writerow(["Old Label", "New Label"])                
#            for n in list(G.nodes()):
#                writer.writerow([G.nodes[n]['old_labels'], n]) 
    

    return G









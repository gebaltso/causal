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



def Louvain(G):
    
    partition = community_louvain.best_partition(G)
    modularity = community_louvain.modularity(partition, G)
    
    
    return partition, modularity

#    # draw the graph
#    pos = nx.spring_layout(G)
#    # color the nodes according to their partition
#    cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
#    nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=80,
#                           cmap=cmap, node_color=list(partition.values()))
#    nx.draw_networkx_edges(G, pos, alpha=0.5)
#    nx.draw_networkx_labels(G,pos,font_size=7)
#    #plt.savefig("plot.png", dpi=1000)
#    plt.show()
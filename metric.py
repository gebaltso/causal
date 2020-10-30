#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 16:32:37 2020

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


def metric(n, partition, communities_dict, G):
    
    # Find the community C than node n belongs
    C = partition.get(n)
    
    # Find all the nodes of community C (n belongs to C)
    community =  [v for k, v in communities_dict.items() if k == C]
    
    print(community)
    
    
    
    return 0
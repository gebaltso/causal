#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:29:31 2020

@author: georgia
"""

import networkx as nx
from random_walks import RandomWalk


karate_g = nx.karate_club_graph()

random_walk = RandomWalk(karate_g, walk_length=3, num_walks=10, p=1, q=1, workers=6)

walklist = random_walk.walks

for w in walklist:
    print(w)
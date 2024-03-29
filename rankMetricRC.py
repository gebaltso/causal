#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 17:53:19 2020

@author: georgiabaltsou
"""

import networkx as nx
from collections import defaultdict
import operator
import math


def FindComm(n, partition, memb, G):
    
    # Find the community C that node n belongs
    C = memb[n]
    # Find all the nodes of community C (n belongs to C)
    community = [i for i, e in enumerate(memb) if e == C]
    
    return community
    
def FinddegInAndMaxDeg(community, G):
    
    # Find the degree of each node inside its community
    degInC = defaultdict(dict)
    for i in community:
        degInC[i] = 0
        
    for i in degInC:
        for j in community:
            if i in G.neighbors(j):
                degInC[i] +=1
            
    # Keep the max value of the above dict
    maxDegInC = max(degInC.items(), key=operator.itemgetter(1))[1]
    # If I wanted the index aka the node with the maximum in degree in community
    # nodeWithMaxDegInC = max(degInC.items(), key=operator.itemgetter(1))[0]
    
    return degInC, maxDegInC


def Findkin(community, G):
    
    kin = defaultdict(dict)
    for i in community:
        kin[i] = []
        
    for i in kin:
        for j in community:
            if i in G.neighbors(j):
                kin[i].append(j) 
    
    return kin


def FindBasicForAll(n, partition, memb, G):
    
    # Find the community C that node n belongs
    C = memb[n]
    # Find all the nodes of community C (n belongs to C)
    community = [i for i, e in enumerate(memb) if e == C]
    
    # The degree into community of n and
    # the maximum degree inside this community
    degInC, maxDegInC = FinddegInAndMaxDeg(community, G)
#    print(degInC)
    
    # The dictionary kin contains the nodes of community of
    # n and their neighbours only inside community
    kin = Findkin(community, G)
    
    return community, degInC, maxDegInC, kin


def Rccomputation(n, degInC, maxDegInC):
    # Compute metric RC for node = n
    RCn = math.log(degInC[n]+1)/math.log(maxDegInC+1)
#    print("RCn=", RCn, "degInC=", degInC[n],"maxDegInC=", maxDegInC )
    return RCn



###############################################################################

def rankMetricRC(n, partition, memb, G):
    
    # Find all the basic values
    community, degInC, maxDegInC, kin = FindBasicForAll(n, partition, memb, G)
    
    # Compute RC of node n
    RCn = Rccomputation(n, degInC, maxDegInC)
    
    
    # Compute the nominator of RC
    nom = 0
    # for each neighbor of n compute the RC and sum all the RCs
    for i in kin[n]:
        nom = nom + Rccomputation(i, degInC, maxDegInC)
 
        
    # Find all the neighbours of node n
#    nei1 = [i for i in G[n]]
#    nei2 = [i for i in G.predecessors(n)]
#    nei = nei1 + nei2
    nei = G.neighbors(n)
#    print("Neigbhours= ", nei)

    
    # Find the external neighbours of node n (that are outside the community of n)
    extNei = []
    for i in nei:
        if i not in community:
            extNei.append(i) 
#    print("External neighbours= ", extNei)
    
    dem2 = 0
    for i in extNei:
        community, degInC, maxDegInC, kin = FindBasicForAll(i, partition, memb, G)
        for j in kin[i]:
            dem2 = dem2 + Rccomputation(j, degInC, maxDegInC)
    
    
    # Compute total metric RC
    RC = RCn * ((nom)/(nom + dem2))
    
    return RC
    
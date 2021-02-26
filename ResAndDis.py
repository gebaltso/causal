#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 12:06:28 2020

@author: georgiabaltsou
"""



def resAndDis(memb, memb2, initnode, contigencyLen, nodesLen, rankOfCausals, e):
                 
   res = 1/(1+contigencyLen)
   # Find all the nodes that have changed communities
   nodesChanged = []
   for i in range(len(memb)):
       if memb2[i] != memb[i]:
           nodesChanged.append(i) 
   if initnode in nodesChanged:
        nodesChanged.remove(initnode)
        dis = len(nodesChanged)/((nodesLen)-1)
   rankOfCausals[e] = (res, dis)
   
   return rankOfCausals
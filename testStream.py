#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 13:17:13 2020

@author: georgiabaltsou
"""

#from gephistreamer import graph
#from gephistreamer import streamer
#import random as rn
#
#stream = streamer.Streamer(streamer.GephiWS(hostname="localhost",port=8080,workspace="workspace1"))
#
#szfak = 100  # this scales up everything - somehow it is needed
#cdfak = 3000
#
#nodedict = {}
#def addfnode(fname):
#  # grab the node out of the dictionary if it is there, otherwise make a newone
#  if (fname in nodedict):
#    nnode = nodedict[fname]
#  else:
#    nnode = graph.Node(fname,size=szfak,x=cdfak*rn.random(),y=cdfak*rn.random(),color="#8080ff",type="f")
#    nodedict[fname] = nnode # new node into the dictionary
#  return nnode
#
#def addnodes(pname,fnodenamelist):
#  pnode = graph.Node(pname,size=szfak,x=cdfak*rn.random(),y=cdfak*rn.random(),color="#ff8080",type="p")
#  stream.add_node(pnode)
#  for fname in fnodenamelist:
#    print(pname+"-"+fname)
#    fnode = addfnode(fname)
#    stream.add_node(fnode)
#    pfedge = graph.Edge(pnode,fnode,weight=rn.random())
#    stream.add_edge(pfedge)
#
#person1friends = ['mike','alex','arker','locke','dave','david','ross','rachel','anna','ann','darl','carl','karle']
#person2friends = ['mika','adlex','parker','ocke','ave','david','rosse','rachel','anna','ann','darla','carla','karle']
#person3friends = ['mika','alex','parker','ocke','ave','david','rosse','ross','anna','ann','darla','carla','karle','sasha','daria']
#
#addnodes("p1",person1friends)
#addnodes("p2",person2friends)
#addnodes("p3",person3friends)

##########################################################################################################################



from gephistreamer import graph
from gephistreamer import streamer
import itertools
import random
import time

stream = streamer.Streamer(streamer.GephiWS(hostname="localhost",port=8080,workspace="workspace2"))
test =  [x for x in itertools.permutations('ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890', 2) ]
random.shuffle(test)
for source, target in test:
    
    node_source = graph.Node(source)
    node_target = graph.Node(target)

    stream.add_node(node_source,node_target)
    # time.sleep(0.5) # Make it slower :D
    stream.add_edge(graph.Edge(node_source,node_target))

time.sleep(1) #It might be possible the script runs too fast and last action anr't sent properly


#stream = streamer.Streamer(streamer.GephiWS(hostname="localhost",port=8080,workspace="workspace1"))
#
## Create a node with a custom_property
#node_a = graph.Node("A",custom_property=1)
#
## Create a node and then add the custom_property
#node_b = graph.Node("B")
#node_b.property['custom_property']=2
#
## Add the node to the stream
## you can also do it one by one or via a list
## l = [node_a,node_b]
## stream.add_node(*l)
#stream.add_node(node_a,node_b)
#
## Create edge 
## You can also use the id of the node :  graph.Edge("A","B",custom_property="hello")
#edge_ab = graph.Edge(node_a,node_b,custom_property="hello")
#stream.add_edge(edge_ab)























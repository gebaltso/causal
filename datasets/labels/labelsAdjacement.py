#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 17:45:00 2021

@author: georgiabaltsou
"""

import csv
import os
from collections import defaultdict


edgeFile = "publEdgelist.csv"
labelFile = "labels.csv"
outFile = "new.csv"


with open(edgeFile, 'r') as edges:
    readerEdges = csv.reader(edges, delimiter=',')
    
    for rowEdge in readerEdges: # 0 30
        line = []
        for node in rowEdge: # 0
            with open(labelFile, 'r') as csvInfile:
                readerNames = csv.reader(csvInfile, delimiter=',')                                
                for row in readerNames:                    
                    if node == row[0]:                         
                        line.append(row[1])
                        
        with open(outFile, 'a') as out_file:
            writer = csv.writer(out_file, delimiter=',')            
            writer.writerow([line])
                
with open(outFile, 'r') as infile, \
     open("test.txt", 'a') as outfile:
    data = infile.read()
    data = data.replace("'", "")
    data = data.replace("[", "")
    data = data.replace(", ", ",")
    data = data.replace("]", "")
    outfile.write(data)  

with open("test.txt", 'r') as in_file:
    stripped = (line.strip() for line in in_file)
    lines = (line.split(",") for line in stripped if line)  #change split to "\t" for tabbed txt file
    with open("test.csv", 'w') as out_file:
        writer = csv.writer(out_file, delimiter = ',')
        writer.writerows(lines)        
        
with open("test.csv", 'r') as infile, \
     open("final.csv", 'a') as outfile:
    data = infile.read()
    data = data.replace('"', "")
    outfile.write(data) 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 16:56:10 2020

@author: georgiabaltsou
"""

import csv


csv_Infile = 'netScLCEdgeList.csv'
csv_Outfile = 'numsToNames.csv'
csvNames = 'netScLCNodes.csv'
finalNames_txt = 'final.txt' # the edgeList file with edges as Newman, M Park L in txt format


with open(csv_Infile, 'r') as edges:
    readerEdges = csv.reader(edges, delimiter=',')
    
    for rowEdge in readerEdges: # 0 30
        line = []
        for node in rowEdge: # 0
            with open(csvNames, 'r') as csvInfile:
                readerNames = csv.reader(csvInfile, delimiter=',')                                
                for row in readerNames:                    
                    if node == row[0]:                         
                        line.append(row[1])
                        
        with open(csv_Outfile, 'a') as out_file:
            writer = csv.writer(out_file, delimiter=',')            
            writer.writerow([line])
            
with open(csv_Outfile, 'r') as infile, \
     open(finalNames_txt, 'a') as outfile:
    data = infile.read()
    data = data.replace("'", "")
    data = data.replace("[", "")
    data = data.replace(", ", ",")
    data = data.replace("]", "")
    outfile.write(data)
    

csv_file = 'finalA.csv' # the edgeList csv file with edges as "Newman M Park L"
csv_fileB = 'finalB.csv' # the final edgeList csv file with edges as Newman M Park L
                                
                
with open(finalNames_txt, 'r') as in_file:
    stripped = (line.strip() for line in in_file)
    lines = (line.split(",") for line in stripped if line)  #change split to "\t" for tabbed txt file
    with open(csv_file, 'w') as out_file:
        writer = csv.writer(out_file, delimiter = ',')
        writer.writerows(lines)        
        
with open(csv_file, 'r') as infile, \
     open(csv_fileB, 'a') as outfile:
    data = infile.read()
    data = data.replace('"', "")
    outfile.write(data)        
        
     
        
        
        
        
        
        
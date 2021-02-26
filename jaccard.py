#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:23:50 2021

@author: georgiabaltsou
"""


### Calculates Jaccard Similarity. If I want Jaccard distance this is 1-jaccard similarity.

def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union
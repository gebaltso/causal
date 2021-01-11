#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 18:37:24 2021

@author: georgiabaltsou
"""

import pandas as pd
import os

os.chdir('/Users/georgiabaltsou/Downloads/')

data = pd.io.stata.read_stata('dryad.dta')
data.to_csv('dryad.csv')
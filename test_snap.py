#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 15:24:52 2020

@author: georgia
"""

status = False
try:
    import snap
    version = snap.Version
    i = snap.TInt(5)
    if i == 5:
        status = True
except:
    pass

if status:
    print("SUCCESS, your version of Snap.py is %s" % (version))
else:
    print("*** ERROR, no working Snap.py was found on your computer")

print("####")
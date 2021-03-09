#!/usr/bin/env python3

"""
Greedy Homology Basis Genarator

1. calculate MST T according to edge length
2. (G\T)* dual graph construction
   - with edge weight | \sigma(e) |
3. for e in ((G\T)* - maximal spanning tree of (G\T)*):
   - calculate shortest loop with e
"""


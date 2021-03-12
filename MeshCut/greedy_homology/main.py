#!/usr/bin/env python3

"""
Greedy Homology Basis Genarator

1. calculate MST T according to edge length
2. (G\T)* dual graph construction
   - with edge weight | \sigma(e) |
3. for e in ((G\T)* - maximal spanning tree of (G\T)*):
   - calculate shortest loop with e
"""

from .graph import Graph
import openmesh as om
import sys
import logging

def main(options):
   logging.basicConfig(level=logging.INFO)

   if len(options) != 1:
      print(f"Options: obj_file")
      sys.exit(1)
   
   mesh = om.read_trimesh(options[0])
   graph = Graph.from_openmesh(mesh)

   graph.build_mst()
   graph.get_cycle_basis()
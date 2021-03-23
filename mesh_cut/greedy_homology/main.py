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
import numpy as np
import sys
import logging

import pyvista as pv

def om_to_vis_polydata(mesh: om.TriMesh):
   """ Plot (triangulated) OpenMesh via pyvista """
   points = mesh.points()
   # print(points)

   orig_indices = mesh.fv_indices()
   # print(orig_indices)

   trans_indices = np.zeros((orig_indices.shape[0], orig_indices.shape[1] + 1), dtype=np.int32)
   trans_indices[:,1:4] = orig_indices
   trans_indices[:,0] = 3

   trans_indices = np.hstack(trans_indices)
   # print(trans_indices)

   surf = pv.PolyData(points, trans_indices)
   return surf

def lines_to_vis_polydata(points: np.ndarray, edges: list):
   """edges: [(vs, vd), ...]"""
   poly = pv.PolyData()
   poly.points = np.copy(points)
   num_edges = len(edges)

   lines = np.full((num_edges, 3), 2, dtype=np.int_)
   idx = 0
   for v1, v2 in edges:
      lines[idx, 1] = v1
      lines[idx, 2] = v2
      idx += 1

   assert(num_edges == idx)
   poly.lines = lines

   return poly

def combine_plot(*args):
    assert(len(args) >= 1)
    p = pv.Plotter()
    for idx, arg in enumerate(args):
        p.add_mesh(arg[0], **arg[1])
    p.show()

def offscreen_combine_plot(filename, *args):
    assert(len(args) >= 1)
    p = pv.Plotter(off_screen=True, window_size=[1920, 1080])
    p.set_position(np.array([4.0, 6.4, 5.6]))
    for idx, arg in enumerate(args):
        p.add_mesh(arg[0], **arg[1])
    
    p.screenshot(filename)

def main(options):
   logging.basicConfig(level=logging.INFO)

   if len(options) != 1:
      print(f"Options: obj_file")
      sys.exit(1)
   
   mesh = om.read_trimesh(options[0])
   graphInst = Graph.from_openmesh(mesh)

   graphInst.build_mst()
   graphInst.get_cycle_basis()

   hbasis, hpivot, bc_tilde = graphInst.get_h1_basis()

   dim_h1 = hbasis.shape[1]
   if dim_h1 == 0:
      print(f"Expected nonzero H1 dimension, abort.")
      sys.exit(1)
   
   base_data = om_to_vis_polydata(mesh)
   for i in range(0, dim_h1):
      edge_list = graphInst.edge_list_from_vector(hbasis[:, i])
      edge_data = lines_to_vis_polydata(mesh.points(), edge_list)
      offscreen_combine_plot(f"genus1_{i}.png",
         (
            edge_data,
            {
               'color': 'red',
               'line_width': 3.0
            }
         ),
         (
            base_data,
            {
               'color': 'tan',
               'opacity': 0.5,
               'style': 'surface',
               'show_edges': True
            }
         )
      )

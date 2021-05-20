#!/usr/bin/env python3

"""
Greedy Homology Basis Genarator

1. calculate MST T according to edge length
2. (G\T)* dual graph construction
   - with edge weight | \sigma(e) |
3. for e in ((G\T)* - maximal spanning tree of (G\T)*):
   - calculate shortest loop with e
"""

from mesh_cut.handle_loop.homology_opt import HomologyBasisOptimizer
from mesh_cut.handle_loop.annotator import Annotator
from mesh_cut.handle_loop.graphbase import GraphBase
import openmesh as om
import numpy as np
import sys, os
import logging

import pyvista as pv

logger = logging.getLogger(__name__)

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
    #p.set_position(np.array([4.0, 6.4, 5.6]))
    for idx, arg in enumerate(args):
        p.add_mesh(arg[0], **arg[1])
    
    p.screenshot(filename)

def main(options):
   logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)40s - %(levelname)s - %(message)s')

   if len(options) != 1:
      print(f"Options: obj_file")
      sys.exit(1)

   logger.info(f"Reading {options[0]}")
   mesh = om.read_trimesh(options[0])

   logger.info("Constructing GraphBase..")
   graphBase = GraphBase.from_openmesh(mesh)

   logger.info("Constructing volumetric GraphBase...")
   volumetricGraphBase = GraphBase.volumetric_from_openmesh(mesh)
   annotator = Annotator(volumetricGraphBase)

   logger.info("Calculating annotation..")
   annotation, null_vector = annotator.compute_annotation()
   graphBase.set_annotation(annotation, null_vector)

   optim = HomologyBasisOptimizer(graphBase)

   logger.info("Computing optimal basis..")
   cycles = optim.compute_optimal_basis()
   logger.info("Optimal basis computation finished.")

   base_data = om_to_vis_polydata(mesh)
   for i in range(0, len(cycles)):
      edge_data = lines_to_vis_polydata(
            mesh.points(),
            graphBase.edge_list_from_vector(graphBase.get_path_vector(cycles[i][1]))
         )
      resname = os.path.split(options[0])[-1].split(".")[0]
      offscreen_combine_plot(f"{resname}_{i}_optim.png",
      #combine_plot(
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

#!/usr/bin/env python3

import openmesh as om
import numpy as np
from om_plot.plot import om_to_vis_polydata, combine_plot, offscreen_combine_plot
from cutmesh.mesh import CutMesh
import random

def unit_test():
    pass

cur_shot_id = 0

def visualize(cmesh, mesh, plot_func):
    new_plot = cmesh.to_vis_polydata()
    plot_func(
        (
            new_plot,
            {
                "color": "red",
                "line_width": 3.0
            }
        ),
        (
            om_to_vis_polydata(mesh),
            {
                "color": "tan",
                "opacity": 0.5,
                #"style": "wireframe"
                "style": "surface",
                "show_edges": True
            }
        )
    )

def get_plot_func():
    # global cur_shot_id
    # this_shot = cur_shot_id
    # def pfunc(*args):
    #     offscreen_combine_plot(f'scr_{this_shot}', *args)
    
    # cur_shot_id += 1
    # return pfunc
    return combine_plot

def make_initial_cut(mesh, seed=None):
    """
    Remove seed triangle.
    while there remains an edge e adjacent to only one triangle t
        Remove e and t.
    while there remains a vertex v adjacent to only one edge e
        Remove v and e.
    """
    # remove seed triangle
    if seed is None:
        seed = random.randint(0, mesh.n_faces() - 1)
    v1, v2, v3 = mesh.fv_indices()[seed].tolist()

    print(f"[make_initial_cut] Seed: {seed}, {v1} {v2} {v3}")
    cmesh = CutMesh(mesh)

    # for va, vb in cmesh.get_edge_iterator():
    #     print(f"{va} {vb}")

    #visualize(cmesh, mesh, get_plot_func())
    cmesh.remove_face(v1, v2, v3)
    #visualize(cmesh, mesh, get_plot_func())

    while True:
        removed = False
        for va, vb in cmesh.get_edge_iterator():
            vc_set = cmesh.edge_adj_faces(va, vb)
            if len(vc_set) == 1:
                # an edge adjacent to only one triangle
                cmesh.remove_edge(va, vb)
                #visualize(cmesh, mesh, get_plot_func())
                removed = True
                break

        if not removed:
            break

    while True:
        removed = False
        for v in range(0, cmesh.n_points):
            adj_set = cmesh.adj_verts(v)
            if len(adj_set) == 1:
                (adj_v, ) = adj_set
                cmesh.remove_edge(v, adj_v)
                #visualize(cmesh, mesh, get_plot_func())
                removed = True
                break
        
        if not removed:
            break
    
    # for va, vb in cmesh.get_edge_iterator():
    #     print(f"{va} {vb}")
    visualize(cmesh, mesh, get_plot_func())

if __name__ == '__main__':
    mesh = om.read_trimesh('./meshes/Genus2.obj')
    # print(mesh.points())

    make_initial_cut(mesh, 27)
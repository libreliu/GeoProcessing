import unittest
import openmesh as om
from mesh_cut.greedy_homology import graph
from mesh_cut.greedy_homology.linalg import *

import numpy as np
class GraphTest(unittest.TestCase):
    def setUp(self):
        MESH_BASEPATH = "./meshes"

        self.meshes = {
            'genus0': om.read_trimesh(f"{MESH_BASEPATH}/Genus0.obj"),
            'genus1': om.read_trimesh(f"{MESH_BASEPATH}/Genus1.obj")
        }
    
    def test_boundary_basis(self):
        graphInst = graph.Graph.from_openmesh(self.meshes['genus0'])
        bbasis, bpivot = graphInst.get_boundary_basis()
        self.assertEqual(bbasis.shape, (18, 12))

        print(bbasis)
        print(get_Bopt_column(bbasis))
    
    def test_cycle_basis(self):
        graphInst = graph.Graph.from_openmesh(self.meshes['genus0'])
        cbasis = graphInst.get_cycle_basis()
        
        print(cbasis)
        print(get_Bopt_column(cbasis))

    def test_h1_basis(self):
        graphInst = graph.Graph.from_openmesh(self.meshes['genus0'])
        hbasis, hpivot, bc_tilde = graphInst.get_h1_basis()
        
        assert(len(hpivot) == 0)


        graphInst = graph.Graph.from_openmesh(self.meshes['genus1'])
        hbasis, hpivot, bc_tilde = graphInst.get_h1_basis()
        
        print(hbasis)
        assert(len(hpivot) == 2)

    def test_cycle_annotations(self):
        graphInst = graph.Graph.from_openmesh(self.meshes['genus0'])
        graphInst.get_cycle_annotations()
        
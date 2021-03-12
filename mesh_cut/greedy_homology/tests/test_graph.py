import unittest
import openmesh as om
from mesh_cut.greedy_homology import graph
import numpy as np
class GraphTest(unittest.TestCase):
    def setUp(self):
        MESH_BASEPATH = "./meshes"

        self.meshes = {
            'genus0': om.read_trimesh(f"{MESH_BASEPATH}/Genus0.obj")
        }

    def test_Bopt(self):
        a = np.array([
            [1, 0, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 0, 1]
        ], dtype=np.int8)
        b = np.array([
            [1, 0, 1, 0],
            [0, 1, 0, 1],
            [0, 0, 1, 1]
        ], dtype=np.int8)
        c = np.array([
            [1, 0, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 1]
        ], dtype=np.int8)
        d = np.array([
            [1, 0, 1, 1],
            [0, 0, 0, 0],
            [0, 0, 0, 0]
        ], dtype=np.int8)
        e = np.array([
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]
        ], dtype=np.int8)
        f = np.array([
            [1, 0, 0, 1],
            [0, 0, 0, 1],
            [0, 1, 0, 0]
        ], dtype=np.int8)

        graphInst = graph.Graph.from_openmesh(self.meshes['genus0'])
        self.assertEqual(graphInst.get_Bopt_column(a), [0, 2, 3])
        self.assertEqual(graphInst.get_Bopt_column(b), [0, 1, 2])
        self.assertEqual(graphInst.get_Bopt_column(c), [0, 2])
        self.assertEqual(graphInst.get_Bopt_column(d), [0])
        self.assertEqual(graphInst.get_Bopt_column(e), [])
        self.assertEqual(graphInst.get_Bopt_column(f), [0, 1, 3])
    
    def test_boundary_basis(self):
        graphInst = graph.Graph.from_openmesh(self.meshes['genus0'])
        bbasis, bpivot = graphInst.get_boundary_basis()
        self.assertEqual(bbasis.shape, (18, 12))

        print(bbasis)
        print(graphInst.get_Bopt_column(bbasis))
    
    def test_cycle_basis(self):
        graphInst = graph.Graph.from_openmesh(self.meshes['genus0'])
        cbasis = graphInst.get_cycle_basis()
        
        print(cbasis)
        print(graphInst.get_Bopt_column(cbasis))

    def test_h1_basis(self):
        graphInst = graph.Graph.from_openmesh(self.meshes['genus0'])
        hbasis, hpivot = graphInst.get_h1_basis()
        
        print(hbasis)
        print(hpivot)
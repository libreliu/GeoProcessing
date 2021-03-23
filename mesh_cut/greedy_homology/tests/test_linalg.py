import unittest
from mesh_cut.greedy_homology.linalg import *

import numpy as np

class LinAlgTest(unittest.TestCase):
    def test_solve_lstsq(self):

        c_basis = np.array([
                [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1],
                [1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0],
                [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1],
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0],
                [0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1],
                [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1]
            ], dtype=np.int8)

        x = solve_z2(c_basis, [
                1,
                1,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                1,
                1,
                0,
                0,
                0,
                0,
            ])

        # two possible? [1,0,0,1,0,...] or [0, 1, 1, ... ,1, 0] ?
        # however [0, 1, 1, ... ,1, 0] is not in Ax = b, only AtAx = Atb

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

        self.assertEqual(get_Bopt_column(a), [0, 2, 3])
        self.assertEqual(get_Bopt_column(b), [0, 1, 2])
        self.assertEqual(get_Bopt_column(c), [0, 2])
        self.assertEqual(get_Bopt_column(d), [0])
        self.assertEqual(get_Bopt_column(e), [])
        self.assertEqual(get_Bopt_column(f), [0, 1, 3])
    
    def test_solve(self):
        a = np.array([
            [1, 1],
            [0, 0],
            [0, 1]
        ], dtype=np.int8)

        a_b = np.array([
            [1],
            [0],
            [0]
        ], dtype=np.int8)

        res = solve_z2(a, a_b)
        self.assertTrue(((a @ res - a_b) % 2 == 0).all())
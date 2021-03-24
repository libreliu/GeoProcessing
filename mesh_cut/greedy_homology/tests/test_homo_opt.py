from mesh_cut.greedy_homology.graphbase import *
from mesh_cut.greedy_homology.annotator import Annotator
from mesh_cut.greedy_homology.homology_opt import HomologyBasisOptimizer
import unittest
import openmesh as om

class HomologyOptimizerTest(unittest.TestCase):
    def setUp(self) -> None:
        MESH_BASEPATH = "./meshes"

        self.meshes = {
            'genus0': om.read_trimesh(f"{MESH_BASEPATH}/Genus0.obj"),
            'genus1': om.read_trimesh(f"{MESH_BASEPATH}/Genus1.obj")
        }

        logging.basicConfig(level=logging.DEBUG)

    def test_optim(self):
        graphBase = GraphBase.from_openmesh(self.meshes['genus1'])
        annotator = Annotator(graphBase)

        annotation, null_vector = annotator.compute_annotation()
        graphBase.set_annotation(annotation, null_vector)

        optim = HomologyBasisOptimizer(graphBase)

        optim.compute_optimal_basis()
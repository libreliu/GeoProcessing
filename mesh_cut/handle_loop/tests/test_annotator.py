from mesh_cut.handle_loop.sp_tree import SpanningTree
from mesh_cut.handle_loop.graphbase import *
from mesh_cut.handle_loop.annotator import Annotator
from mesh_cut.handle_loop.main import offscreen_combine_plot, combine_plot, om_to_vis_polydata, lines_to_vis_polydata
import unittest
import openmesh as om

class AnnotatorTest(unittest.TestCase):
    def setUp(self) -> None:
        MESH_BASEPATH = "./meshes"

        self.meshes = {
            'genus0': om.read_trimesh(f"{MESH_BASEPATH}/Genus0.obj"),
            'genus1': om.read_trimesh(f"{MESH_BASEPATH}/Genus1.obj")
        }

        logging.basicConfig(level=logging.DEBUG)

    def test_annotation(self):
        graphBase = GraphBase.from_openmesh(self.meshes['genus1'])
        annotator = Annotator(graphBase)
        annotation, null_vector = annotator.compute_annotation()
        
        graphBase = GraphBase.volumetric_from_openmesh(self.meshes['genus1'])
        annotator = Annotator(graphBase)

        annotation, null_vector = annotator.compute_annotation()

        print(annotation)
    
    def test_vis_homology_basis(self):
        voluGraphBase = GraphBase.volumetric_from_openmesh(self.meshes['genus1'])
        annotator = Annotator(voluGraphBase)

        sp_tree = SpanningTree(voluGraphBase)
        sp_tree.build_mst()
        
        bbasis, bpivot = annotator.get_boundary_basis()
        cbasis = annotator.get_cycle_basis(sp_tree)
        dim_boundary = len(bpivot)
        
        z_tilde, dim_h1 = annotator.get_h1_basis(bbasis, cbasis, bpivot)
        h1_vec = z_tilde[:, (z_tilde.shape[1] - dim_h1):].reshape(voluGraphBase.n_edges)

        print(h1_vec)

        base_data = om_to_vis_polydata(self.meshes['genus1'])
        edge_data = lines_to_vis_polydata(
            voluGraphBase._points,
            voluGraphBase.edge_list_from_vector(h1_vec)
        )
        combine_plot(
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


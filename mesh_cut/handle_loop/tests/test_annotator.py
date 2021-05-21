from mesh_cut.handle_loop.homology_opt import HomologyBasisOptimizer, OptimizedHomologyBasisOptimizer
from mesh_cut.handle_loop.sp_tree import SpanningTree
from mesh_cut.handle_loop.graphbase import *
from mesh_cut.handle_loop.annotator import Annotator
from mesh_cut.handle_loop.main import offscreen_combine_plot, combine_plot, om_to_vis_polydata, lines_to_vis_polydata
import unittest
import openmesh as om
import sys

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

    def test_annotation_for_known_h1(self):
        mesh = self.meshes['genus1']

        print("Constructing GraphBase..")
        graphBase = GraphBase.from_openmesh(mesh)

        print("Constructing volumetric GraphBase...")
        volumetricGraphBase = GraphBase.volumetric_from_openmesh(mesh)
        annotator = Annotator(volumetricGraphBase)

        print("Calculating annotation..")
        annotation, null_vector = annotator.compute_annotation()

        sp_tree = SpanningTree(volumetricGraphBase)
        sp_tree.build_mst()
        
        bbasis, bpivot = annotator.get_boundary_basis()
        cbasis = annotator.get_cycle_basis(sp_tree)
        dim_boundary = len(bpivot)
        
        z_tilde, dim_h1 = annotator.get_h1_basis(bbasis, cbasis, bpivot)
        h1_vec = z_tilde[:, (z_tilde.shape[1] - dim_h1):].reshape(volumetricGraphBase.n_edges)

        h1_edge_list = volumetricGraphBase.edge_list_from_vector(h1_vec)

        res = null_vector
        for vpair in h1_edge_list:
            try:
                res += annotation[vpair]
            except KeyError:
                pass
        
        print(res)

    def test_vis_annotated_edges(self):
        mesh = self.meshes['genus1']

        print("Constructing GraphBase..")
        graphBase = GraphBase.from_openmesh(mesh)

        print("Constructing volumetric GraphBase...")
        volumetricGraphBase = GraphBase.volumetric_from_openmesh(mesh)
        annotator = Annotator(volumetricGraphBase)

        print("Calculating annotation..")
        annotation, null_vector = annotator.compute_annotation()
        graphBase.set_annotation(annotation, null_vector)

        optim = HomologyBasisOptimizer(graphBase)

        print("Computing optimal basis..")
        cycles = optim.compute_optimal_basis()
        print("Optimal basis computation finished.")
        print(f"cycles: {len(cycles)}")

        base_data = om_to_vis_polydata(mesh)
        for i in range(0, len(cycles)):
            edge_data = lines_to_vis_polydata(
                    volumetricGraphBase._points,
                    graphBase.edge_list_from_vector(graphBase.get_path_vector(cycles[i][1]))
                )
            # offscreen_combine_plot(f"testnew_{i}.png",

            cycle_edge_list = graphBase.edge_list_from_vector(graphBase.get_path_vector(cycles[i][1]))
            res = null_vector
            for vpair in cycle_edge_list:
                try:
                    res += annotation[vpair]
                except KeyError:
                    pass
            print(res)
        

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


            
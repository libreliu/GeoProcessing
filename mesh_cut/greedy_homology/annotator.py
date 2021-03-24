from .linalg import check_z2, get_Bopt_column, solve_z2_sequential
from .sp_tree import SpanningTree
from .graphbase import GraphBase
import numpy as np
import logging

logger = logging.getLogger(__name__)

class Annotator:
    def __init__(self, graphBase: GraphBase):
        self.graphBase = graphBase

    # Boundary group mechanism, not cached
    def get_boundary_basis(self):
        """Get boundary group basis (n_edges x n_faces)"""
        boundary_basis = np.zeros(
                            (self.graphBase.n_edges, self.graphBase.n_faces),
                            dtype=np.int8
                        )

        for idx, face in enumerate(self.graphBase._fv_indices):
            e0_idx = self.graphBase.edge_lookup[tuple(sorted((face[0], face[1])))]
            e1_idx = self.graphBase.edge_lookup[tuple(sorted((face[1], face[2])))]
            e2_idx = self.graphBase.edge_lookup[tuple(sorted((face[0], face[2])))]

            # boundary(f_idx) w.r.t self.coeff_field
            boundary_basis[:, idx][e0_idx] = 1
            boundary_basis[:, idx][e1_idx] = 1
            boundary_basis[:, idx][e2_idx] = 1

        bpivot = get_Bopt_column(boundary_basis)
        return boundary_basis, bpivot

    # Cycle group mechanism, not cached
    def get_cycle_basis(self, sp_tree: SpanningTree):
        residual_edges = sp_tree.get_residual_edges()
        dim_cycle = len(residual_edges)

        cycle_basis = np.ndarray((self.graphBase.n_edges, dim_cycle), dtype=np.int8)
        for idx, (vs, vd) in enumerate(residual_edges):
            path = sp_tree.get_path(vs, vd)
            path_vector = self.graphBase.get_path_vector(path)
            residual_edge_vector = self.graphBase.get_path_vector([vs, vd])
            
            cycle_basis[..., idx] = path_vector + residual_edge_vector
            assert(check_z2(cycle_basis[..., idx]))

        return cycle_basis

    def get_h1_basis(self, bbasis: np.ndarray, cbasis: np.ndarray, bpivot: list):
        bcmatrix = np.ndarray(
                        (self.graphBase.n_edges, bbasis.shape[1] + cbasis.shape[1]),
                        dtype=np.int8
                    )
        bcmatrix[:, 0:bbasis.shape[1]] = bbasis
        bcmatrix[:, bbasis.shape[1]:] = cbasis

        bc_pivot = get_Bopt_column(bcmatrix)
        assert(bc_pivot[:len(bpivot)] == bpivot)
        dim_h1 = len(bc_pivot) - len(bpivot)
        z_tilde = bcmatrix[:, bc_pivot]

        return z_tilde, dim_h1

    def compute_annotation(self):
        """Returns (annotation, annotation_null_vector) tuple
        annotation[(vs, vd)] is not None if have non-null annotation
        (where vs < vd)
        """

        sp_tree = SpanningTree(self.graphBase)
        sp_tree.build_mst()

        residual_edges = sp_tree.get_residual_edges()
        dim_cycle = len(residual_edges)
        logger.info(f"Cycle basis dimension: {dim_cycle}")

        boundary_basis, bpivot = self.get_boundary_basis()
        dim_boundary = len(bpivot)
        logger.info(f"Boundary basis dimension: {dim_boundary}")

        cycle_basis = self.get_cycle_basis(sp_tree)
        z_tilde, dim_h1 = self.get_h1_basis(boundary_basis, cycle_basis, bpivot)
        
        logger.info(f"H1 basis dimension: {dim_h1}")
        
        coord_mat = solve_z2_sequential(z_tilde, cycle_basis)
        annotation_dict = {}

        for idx, (vs, vd) in enumerate(residual_edges):
            h1_coeff = coord_mat[dim_h1:, idx]
            annotation_dict[(vs, vd)] = h1_coeff
        
        annotation_null_vector = np.zeros((dim_h1,), dtype=np.int8)
        return (annotation_dict, annotation_null_vector)
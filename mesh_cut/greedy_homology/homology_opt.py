from .linalg import get_Bopt_column
from .sp_tree import SpanningTree
from .graphbase import GraphBase
import numpy as np
import logging

logger = logging.getLogger(__name__)

class HomologyBasisOptimizer:
    def __init__(self, graphBase: GraphBase):
        self.graphBase = graphBase

    def compute_optimal_basis(self):
        cycles = []
        for v in range(0, self.graphBase.n_vertices):
            sp_sptree = SpanningTree(self.graphBase)

            sp_sptree.build_spt(v, True)

            # annotate all cycles & record cycle length
            residual_edges = sp_sptree.get_residual_edges()

            for idx, (vs, vd) in enumerate(residual_edges):
                path = sp_sptree.get_path(vs, vd)
                cycle_length = \
                        self.graphBase.get_path_length(path) + \
                        self.graphBase.get_path_length((vd, vs))

                annotation = \
                    sp_sptree.vertice_annotation[vs] ^ \
                    sp_sptree.vertice_annotation[vd] ^ \
                    self.graphBase.get_edge_annotation(vs, vd)

                cycles.append(
                    (cycle_length, path + [vs], annotation)
                )
        
        num_cycles = len(cycles)

        logger.info(f"Number of candidate cycles: {num_cycles}")
        assert(num_cycles != 0)

        cycles.sort(key= lambda x: x[0])

        anno_matrix = np.ndarray(
                        (cycles[0][2].shape[0], num_cycles),
                        dtype=np.int8
                    )
        for idx, cycleTuple in enumerate(cycles):
            anno_matrix[:, idx] = cycleTuple[2]
        
        pivots = get_Bopt_column(anno_matrix)
        assert(len(pivots) == cycles[0][2].shape[0])

        return [
            cycles[idx] for idx in pivots
        ]

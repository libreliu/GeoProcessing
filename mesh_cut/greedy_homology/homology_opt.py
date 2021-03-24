from .sp_tree import SpanningTree
from .graphbase import GraphBase
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
                    (cycle_length, path, annotation)
                )
        
        num_cycles = len(cycles)

        logger.info(f"Number of candidate cycles: {num_cycles}")

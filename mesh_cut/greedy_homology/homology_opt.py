from .sp_tree import SpanningTree
from .graphbase import GraphBase
import logging

logger = logging.getLogger(__name__)

class HomologyBasisOptimizer:
    def __init__(self, graphBase: GraphBase):
        self.graphBase = graphBase

    def compute_optimal_basis(self):
        for vs in range(0, self.graphBase.n_vertices):
            sp_sptree = SpanningTree(self.graphBase)

            sp_sptree.build_spt(vs)

            # annotate all cycles & record cycle length
            

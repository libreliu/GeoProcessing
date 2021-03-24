from .linalg import check_z2
from .graphbase import GraphBase
from .heapdict import heapdict
import numpy as np
import logging

logger = logging.getLogger(__name__)

class SpanningTree:
    """Represent a spanning tree on top of the graphBase"""
    def __init__(self, graphBase: GraphBase):
        self.graphBase = graphBase
    
        # child -> (parent, dist-between-child-and-parent)
        self.parent_tree = None
        self.root_id = None
        # (vs, vd) with vs <= vd
        self.edge_set = set()

        # edges that weren't chosen for the spanning tree
        self.residual_edges = None

        # Only available in SPT
        self.dists = None

        # vertice annotation
        self.vertice_annotation = None

    def get_residual_edges(self):
        if self.residual_edges is not None:
            return self.residual_edges

        assert(len(self.edge_set) > 0)
        self.residual_edges = list(self.graphBase.edge_set - self.edge_set)

        return self.residual_edges
    
    def get_path_to_root(self, node_id):
        assert(self.parent_tree is not None and self.root_id is not None)

        path = [node_id]
        if node_id == self.root_id:
            return path

        next_elem = self.parent_tree[node_id][0]
        while next_elem != self.root_id:
            path.append(next_elem)
            next_elem = self.parent_tree[next_elem][0]
        
        path.append(self.root_id)
        return path

    def get_path(self, start, end):
        """[start_idx, ..., end_idx]"""
        assert(self.parent_tree is not None)
        assert(start != end)

        # (start -> root_id)
        spath = self.get_path_to_root(start)
        # (end -> root_id)
        epath = self.get_path_to_root(end)

        # filter out the redundant by expanding from root to each vertices
        last_passage = -1

        while len(spath) > 0 and len(epath) > 0 \
            and spath[len(spath) - 1] == epath[len(epath) - 1]:
            
            last_passage = spath[len(spath) - 1]
            del spath[len(spath) - 1]
            del epath[len(epath) - 1]

        return spath + [last_passage] + epath[::-1]

    def build_mst(self, start: int = 0):
        """ Build MST using Prim """
        if self.parent_tree is not None:
            raise Exception("Tree already built.")
        
        self.parent_tree = {}
        self.root_id = start
        visited = set([start])

        # always tree -> non-tree, dist_heap contains tree dist to V - visited
        dist_heap = heapdict()
        for (vert_id, (_, dist)) in self.graphBase._v_pool.get(start).edges.items():
            # (dist, src_in_tree, target_outside_tree)
            # duplicate vert_id to ensure a stable sort
            dist_heap[vert_id] = (dist, start, vert_id)
        
        while len(dist_heap) > 0:
            vd, (dist, vs, _) = dist_heap.popitem()
            self.parent_tree[vd] = (vs, dist)
            self.edge_set.add(tuple(sorted((vs, vd))))

            # expand edge vs->vd
            visited.add(vd)

            # process vd
            for (vd_neigh, (_, neigh_dist)) in self.graphBase._v_pool.get(vd).edges.items():
                if vd_neigh not in visited:
                    if vd_neigh in dist_heap:
                        if dist_heap[vd_neigh][0] > neigh_dist:
                            dist_heap[vd_neigh] = (neigh_dist, vd, vd_neigh)
                    else:
                        dist_heap[vd_neigh] = (neigh_dist, vd, vd_neigh)
        
        if len(self.edge_set) != self.graphBase.n_vertices - 1:
            raise Exception("Mesh not connected.")

    def build_spt(self, start: int, annotate=True):
        """Build shortest path tree start from @start, uses Dijkstra"""
        if self.parent_tree is not None:
            raise Exception("Tree already built.")

        work_heap = heapdict()
        self.parent_tree = {}
        self.root_id = start
        
        if annotate:
            assert(self.graphBase.annotation_null_vector is not None)
            self.vertice_annotation = {}
            self.vertice_annotation[start] = self.graphBase.annotation_null_vector

        self.dists = [float('inf') for i in range(0, self.graphBase.n_vertices)]
        for i in range(0, self.graphBase.n_vertices):
            work_heap[i] = float('inf')

        self.dists[start] = 0

        while len(work_heap) > 0:
            vd, _ = work_heap.popitem()
            for (vd_neigh, (_, neigh_dist)) in self.graphBase._v_pool.get(vd).edges.items():
                alt = self.dists[vd] + neigh_dist
                if alt < self.dists[vd_neigh]:
                    self.dists[vd_neigh] = alt
                    self.parent_tree[vd_neigh] = (vd, neigh_dist)
                    if annotate:
                        self.vertice_annotation[vd_neigh] = \
                            (self.vertice_annotation[vd] ^ self.graphBase.get_edge_annotation(vd, vd_neigh))
                    
                    work_heap[vd_neigh] = alt
        
        self.edge_set = set(
            [tuple(sorted((k, v))) for k, (v, _) in self.parent_tree.items()]
        )

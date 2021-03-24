from .linalg import check_z2
import numpy as np
import logging
import openmesh as om

logger = logging.getLogger(__name__)

class GraphVertex:
    def __init__(self, idx: int, coord: np.ndarray):
        """
        coord: np.ndarray[np.float64]
        """
        self.id = idx
        self.coord = coord
        self.edges = {}
    
    def add_edge(self, vert: 'GraphVertex'):
        assert(vert.id not in self.edges)
        distance = np.sqrt(np.sum((vert.coord - self.coord) ** 2))
        self.edges[vert.id] = (vert, distance)
        vert.edges[self.id] = (self, distance)

    def add_edge_if_absent(self, vert: 'GraphVertex'):
        if vert.id not in self.edges:
            self.add_edge(vert)

class GraphPool:
    """Adjacency List"""
    def __init__(self):
        self.pool = {}
    
    def add(self, id, value):
        assert(id not in self.pool)
        self.pool[id] = value

    def get(self, id) -> GraphVertex:
        assert(id in self.pool)
        return self.pool[id]

class GraphBase:
    # --- Traversal mechanics ---
    def all_edges_iterator(self):
        """ returns (v0_id, v1_id, dist), v0 < v1 """
        for vert_id, vert in self._v_pool.pool.items():
            for neigh_id in vert.edges.keys():
                if neigh_id > vert_id:
                    yield (vert_id, neigh_id, vert.edges[neigh_id][1])
    
    def edge_list_from_vector(self, edge_vector: np.ndarray):
        assert(edge_vector.shape == (self.n_edges,))

        out_list = []
        for i in range(0, self.n_edges):
            if edge_vector[i] == 1:
                out_list.append(self.rev_edge_lookup[i])
        
        return out_list

    # Chain group mechanism
    def get_path_vector(self, path):
        """Column vector for path in a given chain group, works for Z_2
        Notice the path must not walk through same edge twice"""
        assert(len(path) >= 2)

        path_vec = np.zeros((self.n_edges), dtype=np.int8)
        for i in range(1, len(path)):
            e = tuple(sorted((path[i - 1], path[i])))
            e_idx = self.edge_lookup[e]
            path_vec[e_idx] += 1
        
        assert(check_z2(path_vec))
        return path_vec

    # ---------------------------

    def set_annotation(self, annotation, null_vector):
        assert(self.annotation is None and self.annotation_null_vector is None)
        self.annotation = annotation
        self.annotation_null_vector = null_vector
    
    def get_edge_annotation(self, vs, vd):
        assert(self.annotation is not None and self.annotation_null_vector is not None)

        edge_pair = tuple(sorted((vs, vd)))
        if edge_pair in self.annotation:
            return self.annotation[edge_pair]
        else:
            return self.annotation_null_vector

    def __init__(self, points: np.ndarray, fv_indices: np.ndarray):

        # -- Annotations --
        self.annotation = None
        self.annotation_null_vector = None
        # -----------------

        self._v_pool = GraphPool()
        self._fv_indices = fv_indices
        self._points = points

        # provide (vs, vd) -> edge index; vs < vd
        self.edge_lookup = {}
        self.rev_edge_lookup = {}
        self.edge_set = set()

        for idx, coord in enumerate(self._points):
            self._v_pool.add(
                idx,
                GraphVertex(idx, coord)
            )

        edge_idx = 0
        for face in self._fv_indices:
            v0 = self._v_pool.get(face[0])
            v1 = self._v_pool.get(face[1])
            v2 = self._v_pool.get(face[2])

            v0.add_edge_if_absent(v1)
            v0.add_edge_if_absent(v2)
            v1.add_edge_if_absent(v2)

            # must be sorted!
            e0 = tuple(sorted((face[0], face[1])))
            e1 = tuple(sorted((face[1], face[2])))
            e2 = tuple(sorted((face[0], face[2])))

            self.edge_set.add(e0)
            self.edge_set.add(e1)
            self.edge_set.add(e2)

            if e0 not in self.edge_lookup:
                self.edge_lookup[e0] = edge_idx
                self.rev_edge_lookup[edge_idx] = e0
                edge_idx += 1

            if e1 not in self.edge_lookup:
                self.edge_lookup[e1] = edge_idx
                self.rev_edge_lookup[edge_idx] = e1
                edge_idx += 1

            if e2 not in self.edge_lookup:
                self.edge_lookup[e2] = edge_idx
                self.rev_edge_lookup[edge_idx] = e2
                edge_idx += 1
        
        self.n_vertices = len(self._points)
        self.n_faces = len(self._fv_indices)
        self.n_edges = len(self.edge_set)

        assert(edge_idx == self.n_edges)
        
        # TODO: check if mesh is closed
        self.genus = 1 - (self.n_vertices - self.n_edges + self.n_faces) / 2
        assert(abs(round(self.genus) - self.genus) < 1e-3)
        self.genus = round(self.genus)

        logger.info(f"V={self.n_vertices}, E={self.n_edges}, F={self.n_faces}, genus={self.genus}")
    
    @staticmethod
    def from_openmesh(mesh: om.TriMesh, copy: bool = False):
        graphInst = GraphBase(mesh.points(), mesh.fv_indices())
        return graphInst
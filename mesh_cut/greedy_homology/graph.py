import openmesh as om
import numpy as np
from .heapdict import heapdict
from .tree import SpanningTree
import logging

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

class Graph:

    def __init__(self, points: np.ndarray, fv_indices: np.ndarray, coeff_field=np.int8):
        """Initialize a graph
        points: np.ndarray[np.float64]
        fv_indices: np.ndarray[np.int32]
        """
        self.min_tree = None
        self.residual_edges = None
        self.cycle_basis = None
        self.boundary_basis = None
        self.residual_edges = None
        self.dim_cycle = None

        self.edge_set = set()
        self.coeff_field = coeff_field

        self._points = points
        self.n_vertices = len(self._points)
        self._fv_indices = fv_indices
        self.n_faces = len(self._fv_indices)
        self._v_pool = GraphPool()

        # provide (vs, vd) -> edge index; vs < vd
        self.edge_lookup = {}

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
                edge_idx += 1

            if e1 not in self.edge_lookup:
                self.edge_lookup[e1] = edge_idx
                edge_idx += 1

            if e2 not in self.edge_lookup:
                self.edge_lookup[e2] = edge_idx
                edge_idx += 1
        
        self.n_edges = len(self.edge_set)
        assert(edge_idx == self.n_edges)
        
        # TODO: check if mesh is closed
        self.genus = 1 - (self.n_vertices - self.n_edges + self.n_faces) / 2
        assert(abs(round(self.genus) - self.genus) < 1e-3)
        self.genus = round(self.genus)

        logger.info(f"V={self.n_vertices}, E={self.n_edges}, F={self.n_faces}, genus={self.genus}")
    
    def all_edges_iterator(self):
        """ returns (v0_id, v1_id, dist), v0 < v1 """
        for vert_id, vert in self._v_pool.pool.items():
            for neigh_id in vert.edges.keys():
                if neigh_id > vert_id:
                    yield (vert_id, neigh_id, vert.edges[neigh_id][1])

    def build_mst(self, start: int = 0):
        """ Build MST using Prim """
        visited = set([start])
        self.min_tree = SpanningTree(start, self.n_vertices)

        # always tree -> non-tree, dist_heap contains tree dist to V - visited
        dist_heap = heapdict()
        for (vert_id, (_, dist)) in self._v_pool.get(start).edges.items():
            # (dist, src_in_tree, target_outside_tree)
            dist_heap[vert_id] = (dist, start, vert_id)
        
        while len(dist_heap) > 0:
            vd, (dist, vs, _) = dist_heap.popitem()
            self.min_tree.add_node(vd, vs, dist)

            # expand edge vs->vd
            visited.add(vd)

            # process vd
            for (vd_neigh, (_, neigh_dist)) in self._v_pool.get(vd).edges.items():
                if vd_neigh not in visited:
                    if vd_neigh in dist_heap:
                        if dist_heap[vd_neigh][0] > neigh_dist:
                            dist_heap[vd_neigh] = (neigh_dist, vd, vd_neigh)
                    else:
                        dist_heap[vd_neigh] = (neigh_dist, vd, vd_neigh)
        
        if len(self.min_tree.edge_set) != self.n_vertices - 1:
            raise Exception("Mesh not connected")

    def get_path_vector(self, path):
        """Column vector for path, works for Z_2 while ok for Z_n theoretically"""
        assert(len(path) >= 2)

        path_vec = np.zeros((self.n_edges), dtype=self.coeff_field)
        for i in range(1, len(path) - 1):
            e = tuple(sorted((path[i - 1], path[i])))
            e_idx = self.edge_lookup[e]
            path_vec[e_idx] += 1
        
        return path_vec

    def get_cycle_basis(self):
        """Returns a n_edges x d matrix"""
        if self.cycle_basis is not None:
            return self.cycle_basis

        if self.min_tree is None:
            self.build_mst()

        self.residual_edges = self.edge_set - self.min_tree.edge_set
        self.dim_cycle = len(self.residual_edges)
        logger.info(f"Residual edges: {self.dim_cycle}")

        # TODO: check if bool satisfy requirements
        self.cycle_basis = np.ndarray((self.n_edges, self.dim_cycle), dtype=self.coeff_field)
        for idx, (vs, vd) in enumerate(self.residual_edges):
            # TODO: write this
            path = self.min_tree.get_path(vs, vd)
            path_vector = self.get_path_vector(path)
            
            self.cycle_basis[..., idx] = path_vector

        return self.cycle_basis

    def get_boundary_basis(self):
        """Get boundary group basis (n_edges x n_faces)"""
        if self.boundary_basis is not None:
            return self.boundary_basis

        self.boundary_basis = np.zeros((self.n_edges, self.n_faces), dtype=self.coeff_field)
        for idx, face in enumerate(self._fv_indices):
            e0_idx = self.edge_lookup[tuple(sorted((face[0], face[1])))]
            e1_idx = self.edge_lookup[tuple(sorted((face[1], face[2])))]
            e2_idx = self.edge_lookup[tuple(sorted((face[0], face[2])))]

            # boundary(f_idx) w.r.t self.coeff_field
            self.boundary_basis[:, idx][e0_idx] = 1
            self.boundary_basis[:, idx][e1_idx] = 1
            self.boundary_basis[:, idx][e2_idx] = 1
        return self.boundary_basis

    def check_z2(self, A: np.ndarray):
        """Check if given matrix is in Z_2"""
        return ((A == 1) + (A == 0)).all()

    def get_Bopt_column(self, A_input: np.ndarray):
        """Get first rank(A) linear independent column vectors of A over Z_2"""
        A = A_input.copy()
        assert(self.check_z2(A))

        m, n = A.shape
        pivot_column = []

        # no need to care things upside working_row
        working_row = 0
        working_col = 0

        while working_row < m and working_col < n:
            # find element with first non-zero coeff this col
            for i in range(working_row, m):
                if A[i, working_col] != 0:
                    break
            else:
                # all elements in this column is zero
                working_col += 1
                continue

            if i > 0:
                # swap row working_row with row i (working row)
                A[[working_row, i], :] = A[[i, working_row], :]
            
            # (rows, cols) = (rows, cols) - row_vec * coeff_in_working_col
            for r in range(working_row + 1, m):

                A[r] -= A[working_row] * A[r, working_col]

            pivot_column.append(working_col)
            working_row += 1
            working_col += 1

        return pivot_column

    def get_h1_basis(self):
        cbasis = self.get_cycle_basis()
        bbasis = self.get_boundary_basis()
        base_mat = np.ndarray((self.n_edges, bbasis.shape[1] + cbasis.shape[1]), dtype=self.coeff_field)
        base_mat[:, 0:bbasis.shape[1]] = bbasis
        base_mat[:, bbasis.shape[1]:] = cbasis

        self.get_Bopt_column(bbasis)

    @staticmethod
    def from_objfile(filename: str):
        pass

    @staticmethod
    def from_openmesh(mesh: om.TriMesh, copy: bool = False):
        graphInst = Graph(mesh.points(), mesh.fv_indices())
        return graphInst
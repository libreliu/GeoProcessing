import openmesh as om
import numpy as np
from .heapdict import heapdict
from .tree import SpanningTree
import logging
from .linalg import get_Bopt_column, check_z2, solve_z2, solve_z2_sequential

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
        self.dim_cycle = None
        self.h1_basis = None
        self.hpivot = None
        self.z_tilde = None
        # (vs, vd) => annotation
        self.annotation_dict = None

        self.edge_set = set()
        self.coeff_field = coeff_field

        self._points = points
        self.n_vertices = len(self._points)
        self._fv_indices = fv_indices
        self.n_faces = len(self._fv_indices)
        self._v_pool = GraphPool()

        # provide (vs, vd) -> edge index; vs < vd
        self.edge_lookup = {}
        self.rev_edge_lookup = {}

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

    def edge_list_from_vector(self, edge_vector: np.ndarray):
        assert(edge_vector.shape == (self.n_edges,))

        out_list = []
        for i in range(0, self.n_edges):
            if edge_vector[i] == 1:
                out_list.append(self.rev_edge_lookup[i])
        
        return out_list


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
        for i in range(1, len(path)):
            e = tuple(sorted((path[i - 1], path[i])))
            e_idx = self.edge_lookup[e]
            path_vec[e_idx] += 1
        
        assert(check_z2(path_vec))
        return path_vec

    def get_cycle_basis(self):
        """Returns a n_edges x d matrix"""
        if self.cycle_basis is not None:
            return self.cycle_basis

        if self.min_tree is None:
            self.build_mst()

        self.residual_edges = list(self.edge_set - self.min_tree.edge_set)
        self.dim_cycle = len(self.residual_edges)
        logger.info(f"Residual edges: {self.dim_cycle}")

        # TODO: check if bool satisfy requirements
        self.cycle_basis = np.ndarray((self.n_edges, self.dim_cycle), dtype=self.coeff_field)
        for idx, (vs, vd) in enumerate(self.residual_edges):
            # TODO: write this
            path = self.min_tree.get_path(vs, vd)
            path_vector = self.get_path_vector(path)
            residual_edge_vector = self.get_path_vector([vs, vd])
            
            self.cycle_basis[..., idx] = path_vector + residual_edge_vector
            assert(check_z2(self.cycle_basis[..., idx]))

        return self.cycle_basis

    # TODO: reduce from here!
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

        self.pivot_boundary_basis = get_Bopt_column(self.boundary_basis)
        return (self.boundary_basis, self.pivot_boundary_basis)

    def get_h1_basis(self):
        if self.h1_basis != None and self.hpivot != None and self.z_tilde != None:
            return self.h1_basis, self.hpivot, self.z_tilde

        cbasis = self.get_cycle_basis()
        # DEBUG
        cpivot = get_Bopt_column(cbasis)

        bbasis, bpivot = self.get_boundary_basis()
        bcmatrix = np.ndarray((self.n_edges, bbasis.shape[1] + cbasis.shape[1]), dtype=self.coeff_field)
        bcmatrix[:, 0:bbasis.shape[1]] = bbasis
        bcmatrix[:, bbasis.shape[1]:] = cbasis

        bc_pivot = get_Bopt_column(bcmatrix)
        assert(bc_pivot[:len(bpivot)] == bpivot)
        self.z_tilde = bcmatrix[:, bc_pivot]

        self.hpivot = bc_pivot[len(bpivot):]
        self.h1_basis = bcmatrix[:, self.hpivot]

        return self.h1_basis, self.hpivot, self.z_tilde

    def get_cycle_annotations(self):
        """Get annotations for cycles"""
        if self.annotation_dict != None:
            return self.annotation_dict

        _, hpivot, z_tilde = self.get_h1_basis()
        dim_h1 = len(hpivot)
        if dim_h1 == 0:
            raise Exception("Well, no need to do trick on genus 0")
        dim_b1 = z_tilde.shape[1] - dim_h1
        cbasis = self.get_cycle_basis()

        coord_mat = solve_z2_sequential(z_tilde, cbasis)
        self.annotation_dict = {}

        # lookup residual edge according to cbasis
        assert(self.residual_edges != None)
        for idx, (vs, vd) in enumerate(self.residual_edges):
            h1_coeff = coord_mat[dim_h1:, idx]
            self.annotation_dict[(vs, vd)] = h1_coeff
        
        return self.annotation_dict

    def build_spt(self, start: int):
        """Build shortest path tree starting from start"""
        work_heap = heapdict()
        
        prevs = [None for i in range(0, self.n_vertices)]
        dists = [float('inf') for i in range(0, self.n_vertices)]
        for i in range(0, self.n_vertices):
            work_heap[i] = float('int')

        dists[start] = 0

        while len(work_heap) > 0:
            vd, _ = work_heap.popitem()
            for (vd_neigh, (_, neigh_dist)) in self._v_pool.get(vd).edges.items():
                alt = dists[vd] + neigh_dist
                if alt < dists[vd_neigh]:
                    dists[vd_neigh] = alt
                    prevs[vd_neigh] = vd
                    
                    work_heap[vd_neigh] = alt
        
        sptree = SpanningTree(start, self.n_vertices)
        # TODO: fix this
        sptree.parent_tree = prevs
        sptree.edge_set = set(
            [sorted((k, v)) for k, v in prevs.items()]
        )
        
        return dists, prevs, sptree

    def calculate_optimal_basis(self):
        pass

    @staticmethod
    def from_objfile(filename: str):
        pass

    @staticmethod
    def from_openmesh(mesh: om.TriMesh, copy: bool = False):
        graphInst = Graph(mesh.points(), mesh.fv_indices())
        return graphInst
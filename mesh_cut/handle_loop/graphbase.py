from .linalg import check_z2
import numpy as np
import logging
import openmesh as om
import meshpy, meshpy.tet, meshpy.geometry

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

    def get_path_length(self, path):
        """path: [vs, .., vd]"""
        dist = 0
        for i in range(1, len(path)):
            vs = self._v_pool.get(path[i - 1])
            vd = self._v_pool.get(path[i])

            dist += np.sqrt(np.sum((vs.coord - vd.coord) ** 2))

        return dist

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

    def __init__(self, points: np.ndarray, fv_indices: np.ndarray, volumetric: bool = False):

        # -- Annotations --
        self.annotation = None
        self.annotation_null_vector = None
        # -----------------

        self._v_pool = GraphPool()
        self._fv_indices = fv_indices
        self._points = points

        # provide (vs, vd) -> edge index; vs < vd
        self.edge_lookup = {}
        # edge_index -> (vs, vd)
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

        if not volumetric:
            if abs(self.genus - 0.5) < 1e-3:
                raise Exception("The mesh seems already planar!")
            assert(abs(round(self.genus) - self.genus) < 1e-3)

        self.genus = round(self.genus)

        logger.info(f"V={self.n_vertices}, E={self.n_edges}, F={self.n_faces}, genus={self.genus}")
    
    @staticmethod
    def from_openmesh(mesh: om.TriMesh, copy: bool = False):
        if copy:
            graphInst = GraphBase(np.copy(mesh.points()), np.copy(mesh.fv_indices()))
        else:
            graphInst = GraphBase(mesh.points(), mesh.fv_indices())
        return graphInst

    @staticmethod
    def volumetric_from_openmesh(mesh: om.TriMesh, copy: bool = False, box_margin: float = 0.5):
        if copy:
            points = np.copy(mesh.points())
            fv_indices = np.copy(mesh.fv_indices())
        else:
            points = mesh.points()
            fv_indices = mesh.fv_indices()

        num_points = points.shape[0]

        # find bounding box
        xmin = ymin = zmin = 100000.0
        xmax = ymax = zmax = -100000.0

        for point in points:
            # print(point)
            if point[0] < xmin:
                xmin = point[0]
            if point[1] < ymin:
                ymin = point[1]
            if point[2] < zmin:
                zmin = point[2]
            
            if point[0] > xmax:
                xmax = point[0]
            if point[1] > ymax:
                ymax = point[1]
            if point[2] > zmax:
                zmax = point[2]
        
        logger.info(f"AABBMin=({xmin},{ymin},{zmin}) AABBMax=({xmax},{ymax},{zmax})")

        boxPoints, boxFacets, _, boxFacetMarkers = meshpy.geometry.make_box(
            np.array([xmin - box_margin, ymin - box_margin, zmin - box_margin]), np.array([xmax + box_margin, ymax + box_margin, zmax + box_margin])
        )

        boxFacets = meshpy.geometry.offset_point_indices(boxFacets, len(points))

        meshInfo = meshpy.tet.MeshInfo()
        meshInfo.set_points(
            np.vstack((
                points,
                boxPoints
            ))
        )

        meshInfo.set_facets(
            [fv for fv in fv_indices.tolist()] + \
            [fv for fv in boxFacets]
        )

        # TODO: find a point inside the surface
        meshInfo.set_holes([(0.9, 0.9, 0.9)])

        mesh = meshpy.tet.build(meshInfo)

        mesh.write_vtk("tetgen_output.vtk")

        # tetra -> 4 faces, no orientation considered
        tet_points = np.asarray(mesh.points)
        assert(np.allclose(tet_points[0:len(points)], points))
        num_tets = len(mesh.elements)
        tet_fv_indices = np.ndarray((num_tets * 4, 3), dtype=np.int)
        for idx, tetra in enumerate(mesh.elements):
            tet_fv_indices[idx * 4] = [tetra[0], tetra[1], tetra[2]]
            tet_fv_indices[idx * 4 + 1] = [tetra[0], tetra[1], tetra[3]]
            tet_fv_indices[idx * 4 + 2] = [tetra[0], tetra[2], tetra[3]]
            tet_fv_indices[idx * 4 + 3] = [tetra[1], tetra[2], tetra[3]]

        graphInst = GraphBase(tet_points, tet_fv_indices, True)
        return graphInst

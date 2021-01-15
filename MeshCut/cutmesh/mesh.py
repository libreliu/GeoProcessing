import numpy as np
import pyvista as pv

class CutMesh:
    """An overlay to OpenMesh Mesh"""
    def __init__(self, om_mesh):

        # not moving, so a separate array is used to mark
        self.points = np.copy(om_mesh.points())
        self.n_points = om_mesh.n_vertices()
        #self.points_deleted = np.zeros(om_mesh.n_vertices(), dtype=np.uint8)

        self.edges = {k: set() for k in range(0, self.n_points)}

        for i in om_mesh.fv_indices():
            assert(i.shape == (3,))
            for m in range(0, 3):
                # duplicate is possible, and no need to warn
                self.edges[i[m]].add(i[(m + 1) % 3])
                self.edges[i[m]].add(i[(m + 2) % 3])

    def remove_edge(self, s, e):
        if s not in self.edges[e]:
            raise Exception("Removing a non-exist edge")
        else:
            assert(e in self.edges[s])
        
        self.edges[s].remove(e)
        self.edges[e].remove(s)

    def add_edge(self, s, e):
        pass

    def get_edge_iterator(self):
        """ An iterator that visits each edge once """

        def edge_iter():
            cmesh = self
            iterated_edge = set()

            for key, value_set in cmesh.edges.items():
                for value in value_set:
                    if (key, value) not in iterated_edge and (value, key) not in iterated_edge:
                        iterated_edge.add((key, value))
                        yield (key, value)
        
        return edge_iter()

    def adj_verts(self, s):
        """ Find vertices adjacent to s """
        return self.edges[s]

    def edge_adj_faces(self, s, e):
        """
        check and return if the edge is adjacent to faces
        returns a set, containing all triple nodes
        """
        left = self.edges[s]
        right = self.edges[e]

        return left & right

    def remove_face(self, v1, v2, v3):
        if (v2 not in self.edges[v1]) or (v3 not in self.edges[v1]):
            raise Exception("Removing a non-exist triangle")

        self.remove_edge(v1, v2)
        self.remove_edge(v1, v3)
        self.remove_edge(v2, v3)
    
    def get_num_edges(self):
        sum = 0
        for key, value in self.edges.items():
            sum += len(value)
        assert(sum % 2 == 0)
        return int(sum / 2)

    def to_vis_polydata(self):
        """
        Construct pyvista Polydata
        """
        poly = pv.PolyData()
        poly.points = np.copy(self.points)

        num_edges = self.get_num_edges()
        lines = np.full((num_edges, 3), 2, dtype=np.int_)
        idx = 0
        for v1, v2 in self.get_edge_iterator():
            lines[idx, 1] = v1
            lines[idx, 2] = v2
            idx += 1
        
        assert(num_edges == idx)
        
        poly.lines = lines
        return poly
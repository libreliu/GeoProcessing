#!/usr/bin/env python3

class SpanningTree:
    def __init__(self, root_id, n_vertices):
        # parent -> [childs]
        self.tree = {
            root_id: []
        }
        # child -> parent
        self.parent_tree = {}

        self.edge_set = set()
        self.root_id = root_id
        self.n_vertices = n_vertices
    
    def add_node(self, node_id, parent, dist):
        assert(parent in self.tree)
        assert(not any(
            map(lambda x: x[0] == node_id, self.tree[parent])
        ))

        self.tree[parent].append(
            (node_id, dist)
        )
        if node_id not in self.tree:
            self.tree[node_id] = []

        assert(node_id not in self.parent_tree)
        self.parent_tree[node_id] = parent

        # small number first
        vs, vd = sorted((node_id, parent))
        self.edge_set.add(
            (vs, vd)
        )

    def get_path_to_root(self, node_id):
        path = [node_id]
        if node_id == self.root_id:
            return path

        next_elem = self.parent_tree[node_id]
        while next_elem != self.root_id:
            path.append(next_elem)
            next_elem = self.parent_tree[next_elem]
        
        path.append(self.root_id)
        return path

    def get_path(self, start, end):
        """[start_idx, ..., end_idx]"""
        assert(start != end)
        # combine (start -> root_id) & (end -> root_id)
        spath = self.get_path_to_root(start)
        epath = self.get_path_to_root(end)

        return spath + epath[::-1][1:]

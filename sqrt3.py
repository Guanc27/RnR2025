# INTERPOLATORY 
# 1:3 division, flips edges, after two iterations \exists 9 triangles and edges are split to 3
# new vertices have deg 6 and deg of old ones don't change 
# Labsik et. al. also uses adaptive refinement (triangles which do not satisfy a flatness criterion is split)
# i.e. high local curv => more refinement 
# coord for new vertex is q^{k+1} = 32/81*(p1^k+p2^k+p3^k) - 1/81*(p4^k+p5^k+p6^k) - 1/81*(p7^k+...+p12^k)
# above computed based on x=5/3

# adv: 
    # allows computation of more refinement levels until prescribed mesh complexity is reached
    # doesn't hv the problem of holes, no edges are split 
    # sike, there's a bunch of holes, but seems to be boundary cases and where valence is high
# disadv: 
    # 

import numpy as np
import trimesh
import polyscope as ps
import polyscope.imgui as psim
import math

class Sqrt3Subdivision:
    def subdivide(self, V, F, iterations):
        """Perform interpolatory √3-subdivision on a triangular mesh.
        
        Suppose graph G = (V,E) where |V|=n, |E|=m
        - V: (n,3) numpy array of vertex positions.
        - F: (m,3) numpy array of face indices.
        """
        V_sub, F_sub = V.copy(), F.copy()
        for _ in range(iterations):
            V_sub, F_sub = self._subdivide_once(V_sub, F_sub)
        return V_sub, F_sub

    def _subdivide_once(self, V, F):
        # build edges_dict, which maps sorted edge (i, j) to a list of (face index, opposite vertex)
        edges_dict = {}
        for face_i, face in enumerate(F):
            for e in range(3):
                i0, i1, i2 = face[e], face[(e+1)%3], face[(e+2)%3]
                key = self.sort_edge(i0, i1)
                if key not in edges_dict:
                    edges_dict[key] = []
                edges_dict[key].append((face_i, i2))
        
        # build the adjacency dictionary (adj) mapping vertex index to its neighbors
        # adjacency array for a graph
        adj = {i: set() for i in range(len(V))}
        for face in F:
            v0, v1, v2 = face
            adj[v0].update([v1, v2])
            adj[v1].update([v0, v2])
            adj[v2].update([v0, v1])
        
        new_vertices = list(V) # contains original vertices
        
        # face_point_ids maps each face index to the index of its computed face point.
        face_point_ids = {}
        for fi, face in enumerate(F):
            fp = self.compute_face_point(V, face, edges_dict, adj)
            face_point_ids[fi] = len(new_vertices)
            new_vertices.append(fp)
        
        # build a mapping from each original vertex to a list of adjacent face point indices
        # adjacency list but for the new vertices
        vertex_face_map = {i: [] for i in range(len(V))}
        for fi, face in enumerate(F):
            for v in face:
                vertex_face_map[v].append(face_point_ids[fi])
        
        # For each original vertex, sort the associated face points in cyclic order.
        new_faces = []
        for v_idx, fp_list in vertex_face_map.items():
            pos = V[v_idx]
            angles = []
            
            for fp_idx in fp_list:
                fp = new_vertices[fp_idx]
                angle = math.atan2(fp[1] - pos[1], fp[0] - pos[0])
                angles.append((angle, fp_idx))
                
            # sort by angle to get cyclic order
            angles.sort(key=lambda x: x[0])
            sorted_fp = [fp_idx for (_, fp_idx) in angles]
            n = len(sorted_fp)
            
            if n < 2:
                continue
            # new faces by connecting the original vertex with each consecutive pair of face points
            # establishing new connections 
            for i in range(n):
                f1 = sorted_fp[i]
                f2 = sorted_fp[(i + 1) % n]
                new_faces.append([v_idx, f1, f2])
        
        return np.array(new_vertices, dtype=np.float32), np.array(new_faces, dtype=np.int32)

    def compute_face_point(self, V, face, edges_dict, adj):
        """
        Compute the new face point using a 12-point stencil.

        Let p1, p2, p3 be vertices of the triangle that identify a face 
        For each edge (p_i, p_j) of the triangle:
            - Let p_opposite be the vertex opposite that edge (from the neighboring face).
            (Used as one of the three points p4, p5, p6) - this should be same as butterfly
            - Let (p_extra1, p_extra2) be two extra neighbors for that edge,
            computed similarly to the extra neighbors in the Butterfly scheme.
            (Used as two stencil points per edge, total six points: p7,...,p12.)
        """
        p1 = V[face[0]]
        p2 = V[face[1]]
        p3 = V[face[2]]
        
        opposite_points = []  # will hold p4, p5, p6 (one per edge)
        extra_points = []     # will hold 6 points (two per edge - the longer edge)


        # For each edge of the triangle, get the opposite and extra neighbor points.
        vertices = [face[0], face[1], face[2]]
        for i in range(3):
            i0, i1 = vertices[i], vertices[(i+1)%3]
            edge_key = self.sort_edge(i0, i1)
            connected = edges_dict.get(edge_key, [])
            
            opp = None
            for (f_idx, opp_vertex) in connected:
                # assume that if the opposite vertex is not the third vertex of this face,
                # then it is the neighbor.
                if opp_vertex not in vertices:
                    opp = opp_vertex
                    break
            # For boundary edges, fallback to simple midpoint.
            if opp is None:
                opp = i0  # or use i1; here we simply reuse one endpoint

            opposite_points.append(V[opp])
            
            # For extra neighbors, mimic the butterfly extra neighbor search.
            # Here, we simply take two extra neighbors from the intersection of adjacencies.
            # First extra neighbor: common between i0 and opp (excluding i1)
            cand1 = adj[i0].intersection(adj[opp]) - {i1}
            p_extra1 = V[next(iter(cand1))] if cand1 else V[i0]
            # Second extra neighbor: common between i1 and opp (excluding i0)
            cand2 = adj[i1].intersection(adj[opp]) - {i0}
            p_extra2 = V[next(iter(cand2))] if cand2 else V[i1]
            extra_points.extend([p_extra1, p_extra2])
            
            # six extra neighbors, but found by three longer edges connecting the opposite vertices
        
        # apply weights as outlined in paper, a=32/81, b=1/81, c=2/81
        term1 = (32/81) * (p1 + p2 + p3)
        term2 = (1/81) * np.sum(np.array(opposite_points), axis=0)
        term3 = (2/81) * np.sum(np.array(extra_points), axis=0)
        q = term1 - term2 - term3
        #print(f"This is opposite points: {opposite_points}")
        #print(f"This is extra points: {extra_points}")
        return q
    
    def sort_edge(self, a, b):
            return (a, b) if a < b else (b, a)


#mesh = trimesh.load_mesh("./bunny/reconstruction/bun_zipper_res4.ply")
mesh = trimesh.load_mesh("./Aimshape2D/20_cow2.off")


sqrt3 = Sqrt3Subdivision()
V_sub, F_sub = sqrt3.subdivide(mesh.vertices, mesh.faces, iterations=3)
print("Original number of faces:", len(mesh.faces))
print("Subdivided number of faces:", len(F_sub))

ps.init()
ps.register_surface_mesh("Original Mesh", mesh.vertices, mesh.faces)
ps.register_surface_mesh("Subdivided Mesh", V_sub, F_sub)

def callback():
    mesh_enabled = True
    psim.PushItemWidth(100)
    clicked, new_state = psim.Checkbox("Show Original Mesh", mesh_enabled)
    if clicked:
        mesh_enabled = new_state
        ps.get_surface_mesh("Original Mesh").set_enabled(mesh_enabled)
        ps.get_surface_mesh("Subdivided Mesh").set_enabled(not mesh_enabled)
    psim.PopItemWidth()

ps.set_user_callback(callback)
ps.show()

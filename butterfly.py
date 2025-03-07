import thingi10k
import numpy as np
import polyscope as ps
import polyscope.imgui as psim 
import trimesh

class Butterfly_alg():    
    def subdivide_butterfly(self, V, F, iterations):
        """Butterfly subdivision on a triangular mesh."""
        V_sub, F_sub = V.copy(), F.copy() # or else V,F will be modified as well 
        
        for _ in range(iterations):
            V_sub, F_sub = self._subdivide_butterfly_once(V_sub, F_sub)
        return V_sub, F_sub

    def _subdivide_butterfly_once(self, V, F):
        edges_dict = {}
        
        def sort_edge(a, b) -> tuple[int, int]: 
            return (a, b) if a < b else (b, a)
        # sorts the edge in consistent order ()
        
    
        for face_i, face in enumerate(F):
            for e in range(3):
                i0, i1, i2 = face[e], face[(e+1)%3], face[(e+2)%3]
                # vertices are in matrix form
                # np.array([[0,0,0],
                #           [1,0,0],
                #           [0,1,0]])
                # faces are also in matrix form 
                # np.array([[v0,v1,v2],
                #           [v1,v2,v3]])
                edge_key = sort_edge(i0, i1)
                if edge_key not in edges_dict:
                    edges_dict[edge_key] = []
                edges_dict[edge_key].append((face_i, i2))
                # similar to adjacency matrix, but mapping each edge to the face that contains them 
        
        new_edge_verts, new_vertices, new_faces = {}, list(V), []
        # new_vertices will contain original v and new pts
        adj = {i: set() for i in range(len(V))}
        for face in F:
            v0, v1, v2 = face
            adj[v0].update([v1, v2])
            adj[v1].update([v0, v2])
            adj[v2].update([v0, v1])
            
            mid_idx = [] # the midpoint indices to preserve order 
            for e in range(3): # each face has three edges 
                i0, i1 = face[e], face[(e+1)%3] # gets the two vertices that form an edge
                edge_key = sort_edge(i0, i1)
                
                if edge_key not in new_edge_verts:
                    p0, p1 = V[i0], V[i1] # coord of the two vertices
                    connected = edges_dict[edge_key]
                    if len(connected) == 2:
                        oppA = connected[0][1]
                        oppB = connected[1][1]
                        
                        p2 = V[oppA]
                        p3 = V[oppB]
                        p4, p5, p6, p7 = self.find_extra_neighbors(i0, i1, oppA, oppB, adj)    
                        p4 = V[p4] if p4 is not None else V[i0]  
                        p5 = V[p5] if p5 is not None else V[i0]  
                        p6 = V[p6] if p6 is not None else V[i1]  
                        p7 = V[p7] if p7 is not None else V[i1]  
                        
                        w = 1/16
                        pm = 0.5 * (p0 + p1) + 2*w* (p2+p3) - w*(p4+p5+p6+p7)
                    else:
                        # if on the border, then just take the average
                        pm = 0.5 * (p0 + p1)
                    
                    new_edge_idx = len(new_vertices)
                    new_vertices.append(pm) # since index 0, contains coord
                    new_edge_verts[edge_key] = new_edge_idx # this maps edge to the new vertex index 
                
                else:
                    new_edge_idx = new_edge_verts[edge_key] # if alr computed then j use it
                mid_idx.append(new_edge_idx)

            v0, v1, v2, m01, m12, m20 = face[0], face[1], face[2], mid_idx[0], mid_idx[1], mid_idx[2]
            new_faces.extend([[v0, m01, m20], [m01, v1, m12], [m20, m12, v2], [m01, m12, m20]])
        return np.array(new_vertices, dtype=np.float32), np.array(new_faces, dtype=np.int32)

    def find_extra_neighbors(self, i0, i1, oppA, oppB, adj):
        """
        Find the four extra vertices needed for the butterfly subdivision scheme.
        
        Parameters:
        - i0, i1: indices of the two vertices that form the edge being subdivided
        - oppA, oppB: indices of the opposite vertices in the two triangles sharing the edge
        - adj: adjacency dictionary mapping vertex indices to their neighbors
        
        Returns:
        - p4, p5, p6, p7: coordinates of the four extra vertices
        """
        # Get the additional neighbors needed for the butterfly scheme
        # each other vertex is adjacent to at least two of the four given inputs 
        
        p4_cand = adj[i0].intersection(adj[oppA]) - {i1, oppB}
        p4_idx = next(iter(p4_cand)) if p4_cand else None
        
        p5_cand = adj[i0].intersection(adj[oppB]) - {i1, oppA}
        p5_idx = next(iter(p5_cand)) if p5_cand else None
        
        p6_cand = adj[i1].intersection(adj[oppB]) - {i0, oppA}
        p6_idx = next(iter(p6_cand)) if p6_cand else None

        p7_cand = adj[i1].intersection(adj[oppA]) - {i0, oppB}
        p7_idx = next(iter(p7_cand)) if p7_cand else None
        
        return p4_idx, p5_idx, p6_idx, p7_idx
        
#thingi10k.init()
#dataset = thingi10k.dataset()
#V, F = thingi10k.load_file(dataset[0]['file_path'])

mesh = trimesh.load_mesh("./bunny/reconstruction/bun_zipper_res4.ply")

butterfly = Butterfly_alg()
V_sub, F_sub = butterfly.subdivide_butterfly(mesh.vertices, mesh.faces, iterations=3) # # faces multiply by 4 each time
print(len(mesh.faces))
print(len(F_sub))

def compute_face_areas(V, F):
    """Compute the area of each face (triangle)."""
    v0 = V[F[:, 0]]
    v1 = V[F[:, 1]]
    v2 = V[F[:, 2]]
    face_areas = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    return face_areas

def compute_vertex_areas(V, F, face_areas):
    """
    Compute per-vertex area by distributing each face's area equally among its three vertices.
    """
    n = V.shape[0]
    A = np.zeros(n)
    for i, face in enumerate(F):
        for v in face:
            A[v] += face_areas[i] / 3.0
    return A

def compute_face_angles(V, F):
    """
    Compute the interior angles for each vertex in every face.
    Returns an array (M x 3) where each row contains the angles at the vertices of a face.
    """
    v0 = V[F[:, 0]]
    v1 = V[F[:, 1]]
    v2 = V[F[:, 2]]
    
    def angle(a, b, c):
        # Compute the angle at vertex 'a' in triangle (a, b, c).
        ab = b - a
        ac = c - a
        dot = np.einsum('ij,ij->i', ab, ac)
        norm_ab = np.linalg.norm(ab, axis=1)
        norm_ac = np.linalg.norm(ac, axis=1)
        cos_angle = dot / (norm_ab * norm_ac)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.arccos(cos_angle)
    
    angles0 = angle(v0, v1, v2)
    angles1 = angle(v1, v2, v0)
    angles2 = angle(v2, v0, v1)
    
    return np.stack([angles0, angles1, angles2], axis=1)

def compute_gaussian_curvature(V, F):
    """
    Compute per-vertex Gaussian curvature using the angle deficit method:
      K_i = (2*pi - sum(theta_i)) / A_i
    where A_i is the vertex area.
    """
    face_areas = compute_face_areas(V, F)
    vertex_areas = compute_vertex_areas(V, F, face_areas)
    face_angles = compute_face_angles(V, F)
    n = V.shape[0]
    angle_sum = np.zeros(n)
    for i, face in enumerate(F):
        angle_sum[face[0]] += face_angles[i, 0]
        angle_sum[face[1]] += face_angles[i, 1]
        angle_sum[face[2]] += face_angles[i, 2]
    
    gaussian_curvature = (2 * np.pi - angle_sum) / vertex_areas
    return gaussian_curvature

def compute_mean_curvature(V, F):
    """
    Compute per-vertex mean curvature.
    We approximate the Laplace–Beltrami operator using cotangent weights:
      delta v_i = (1/(2A_i)) * sum_{j in N(i)} (cot alpha + cot beta)(v_i - v_j)
    and then set:
      H_i = 0.5 * ||delta v_i||
    """
    n = V.shape[0]
    L = np.zeros((n, 3))  # Laplacian vector for each vertex
    W = {}  # Dictionary for cotangent weights keyed by an edge (i, j) with i < j

    def cot_angle(a, b, c):
        ba = b - a
        ca = c - a
        cos_val = np.dot(ba, ca)
        sin_val = np.linalg.norm(np.cross(ba, ca))
        return cos_val / sin_val if sin_val != 0 else 0

    # Accumulate cotangent weights for each edge.
    for face in F:
        i, j, k = face
        cot_i = cot_angle(V[i], V[j], V[k])
        cot_j = cot_angle(V[j], V[k], V[i])
        cot_k = cot_angle(V[k], V[i], V[j])
        
        key = tuple(sorted((j, k)))
        W[key] = W.get(key, 0) + cot_i

        key = tuple(sorted((i, k)))
        W[key] = W.get(key, 0) + cot_j

        key = tuple(sorted((i, j)))
        W[key] = W.get(key, 0) + cot_k

    # Build the Laplacian operator.
    for (i, j), w in W.items():
        L[i] += w * (V[i] - V[j])
        L[j] += w * (V[j] - V[i])
    
    face_areas = compute_face_areas(V, F)
    vertex_areas = compute_vertex_areas(V, F, face_areas)
    
    mean_curvature = np.zeros(n)
    for i in range(n):
        if vertex_areas[i] > 1e-10:
            delta_vi = L[i] / (2 * vertex_areas[i])
            mean_curvature[i] = 0.5 * np.linalg.norm(delta_vi)
        else:
            mean_curvature[i] = 0.0

    return mean_curvature

# Compute curvature on the subdivided mesh.
gaussian_curv_orig = compute_gaussian_curvature(mesh.vertices, mesh.faces)
mean_curv_orig     = compute_mean_curvature(mesh.vertices, mesh.faces)

gaussian_curv_sub = compute_gaussian_curvature(V_sub, F_sub)
mean_curv_sub     = compute_mean_curvature(V_sub, F_sub)

# === Visualization with Polyscope ===

ps.init()

mesh_orig = ps.register_surface_mesh("Original Mesh", mesh.vertices, mesh.faces)
mesh_orig.add_scalar_quantity("Gaussian Curvature", gaussian_curv_orig, defined_on='vertices')
mesh_orig.add_scalar_quantity("Mean Curvature", mean_curv_orig, defined_on='vertices')

mesh_sub = ps.register_surface_mesh("Subdivided Mesh", V_sub, F_sub)
mesh_sub.add_scalar_quantity("Gaussian Curvature", gaussian_curv_sub, defined_on='vertices')
mesh_sub.add_scalar_quantity("Mean Curvature", mean_curv_sub, defined_on='vertices')

# Register original mesh
#ps.register_surface_mesh("Bunny Mesh", mesh.vertices, mesh.faces, enabled=True)

# Register subdivided mesh
#ps.register_surface_mesh("Subdivided Bunny Mesh", V_sub, F_sub, enabled=True)


def callback():
    mesh_enabled = True
    psim.PushItemWidth(100)

    clicked, new_state = psim.Checkbox("Show Mesh", mesh_enabled)
    
    if clicked:
        mesh_enabled = new_state
        ps.get_surface_mesh("Bunny Mesh").set_enabled(mesh_enabled)
        ps.get_surface_mesh("Subdivided Bunny").set_enabled(not mesh_enabled)
            
    psim.PopItemWidth()
    
# simplify mesh (so can then use subdivision alg to smoothen it)
ps.set_user_callback(callback)
#simplify_mesh()
ps.show()

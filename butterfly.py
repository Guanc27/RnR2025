# INTERPOLATORY
# 1:4 division, most famous interpolatory subdivision scheme for triangular meshes
# smooth surfaces for vertices with valence 6 (deg 6) 
# coord for new vertex is q^{k+1} = 0.5*(p1^k+p2^k) + w/2*(p3^k+p4^k) - w*(p5^k+p6^k+p7^k+p8^k) where usually w=1/16

# adv: 
    # provides local rules to compute new vertex positions
# disadv: 
    # not optimal with vertices \neq 6
    # holes can also be produced when edge is split into triang but not neighboring one 
    # above is fixed by red-green triangulation
    
# Zorin et. al. have developed modified rules, but implement the most basic one here 
# Kobbelt generalize 4-point scheme for subdividing quadrilateral meshes with arbitrary topology


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
        valence=[]
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
            
        # i initially constructed adj as I computed the new neighbors, but could be problematic  
        for face in F:            
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
                        """valence_i0 = len(adj[i0])
                        factor = 1.0 if valence_i0 == 6 else (valence_i0 / 6.0)"""
                        
                        """if len(adj[i0]) < 6:
                            valence.append(len(adj[i0]))"""
                        
                        pm = 0.5 * (p0 + p1) + 2*w* (p2+p3) - (w)*(p4+p5+p6+p7)
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
        """
        # each other vertex is adjacent to at least two of the four given inputs 
        
        p4_cand = adj[i0].intersection(adj[oppA]) - {i1}

        p4_idx = next(iter(p4_cand)) if p4_cand else None
        
        p5_cand = adj[i0].intersection(adj[oppB]) - {i1}

        p5_idx = next(iter(p5_cand)) if p5_cand else None
        
        p6_cand = adj[i1].intersection(adj[oppB]) - {i0}
        p6_idx = next(iter(p6_cand)) if p6_cand else None

        p7_cand = adj[i1].intersection(adj[oppA]) - {i0}
        p7_idx = next(iter(p7_cand)) if p7_cand else None
        
        return p4_idx, p5_idx, p6_idx, p7_idx
        
#thingi10k.init()
#dataset = thingi10k.dataset()
#V, F = thingi10k.load_file(dataset[0]['file_path'])

mesh = trimesh.load_mesh("./Aimshape2D/20_cow2.off")
mesh1 = trimesh.load_mesh("./Aimshape2D/23-Egea/23_egea.off")
mesh2 = trimesh.load_mesh("./Aimshape2D/31-Horse/31_horse.off")

butterfly = Butterfly_alg()
V_sub, F_sub = butterfly.subdivide_butterfly(mesh.vertices, mesh.faces, iterations=1) # faces multiply by 4 each time
V_sub1, F_sub1 = butterfly.subdivide_butterfly(mesh1.vertices, mesh1.faces, iterations=1)
V_sub2, F_sub2 = butterfly.subdivide_butterfly(mesh2.vertices, mesh2.faces, iterations=1)
"""print(f"Vertices: {len(mesh.vertices)}")
print(f"Subdivided vertices: {len(V_sub)}")
print(f"Faces: {len(mesh.faces)}")
print(f"Subdivided faces: {len(F_sub)}")"""

ps.init()

offset = np.array([1, 0.0, 0.0])
mesh1.vertices = mesh1.vertices + offset 
V_sub1 = V_sub1 + offset

offset = np.array([0.0, 1, 0.0])
mesh2.vertices = mesh2.vertices + offset 
V_sub2 = V_sub2 + offset


mesh_orig = ps.register_surface_mesh("Original Mesh Cow", mesh.vertices, mesh.faces)
mesh_sub = ps.register_surface_mesh("Subdivided Mesh Cow", V_sub, F_sub)

mesh_orig1 = ps.register_surface_mesh("Original Mesh Head", mesh1.vertices, mesh1.faces)
mesh_sub1 = ps.register_surface_mesh("Subdivided Mesh Head", V_sub1, F_sub1)

mesh_orig1 = ps.register_surface_mesh("Original Mesh Horse", mesh2.vertices, mesh2.faces)
mesh_sub1 = ps.register_surface_mesh("Subdivided Mesh Horse", V_sub2, F_sub2)


#mesh_orig.add_scalar_quantity("Gaussian Curvature", gaussian_curv_orig, defined_on='vertices')
#mesh_orig.add_scalar_quantity("Mean Curvature", mean_curv_orig, defined_on='vertices')

#mesh_sub.add_scalar_quantity("Gaussian Curvature", gaussian_curv_sub, defined_on='vertices')
#mesh_sub.add_scalar_quantity("Mean Curvature", mean_curv_sub, defined_on='vertices')

def callback():
    mesh_enabled = True
    psim.PushItemWidth(100)

    clicked, new_state = psim.Checkbox("Show Mesh", mesh_enabled)
    
    if clicked:
        mesh_enabled = new_state
        ps.get_surface_mesh("Original Mesh Cow").set_enabled(mesh_enabled)
        ps.get_surface_mesh("Original Mesh Head").set_enabled(mesh_enabled)
        ps.get_surface_mesh("Original Mesh Horse").set_enabled(mesh_enabled)
        ps.get_surface_mesh("Loop Subdivided Mesh Cow").set_enabled(not mesh_enabled)
        ps.get_surface_mesh("Loop Subdivided Mesh Head").set_enabled(not mesh_enabled)
        ps.get_surface_mesh("Loop Subdivided Mesh Horse").set_enabled(not mesh_enabled)
            
    psim.PopItemWidth()
    
ps.set_user_callback(callback)
#simplify_mesh()
ps.show()

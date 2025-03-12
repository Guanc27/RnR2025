# APPROXIMATING
# 1:4, also rly famous by Charles Loop


import thingi10k
import numpy as np
import polyscope as ps
import polyscope.imgui as psim 
import trimesh

class Loop_alg():    
    def subdivide_loop(self, V, F, iterations):
        """Loop subdivision on a triangular mesh."""
        V_sub, F_sub = V.copy(), F.copy() # or else V,F will be modified as well 
        
        for _ in range(iterations):
            V_sub, F_sub = self._subdivide_loop_once(V_sub, F_sub)
        return V_sub, F_sub

    def _subdivide_loop_once(self, V, F):
        edges_dict = {}
        
        def sort_edge(a, b) -> tuple[int, int]: 
            return (a, b) if a < b else (b, a)
        
        # adjacency list info
        adj = {i: set() for i in range(len(V))}
        for face in F:
            v0, v1, v2 = face
            adj[v0].update([v1, v2])
            adj[v1].update([v0, v2])
            adj[v2].update([v0, v1])
        
        # edge to face mapping
        for face_i, face in enumerate(F):
            for e in range(3):
                i0, i1 = face[e], face[(e+1)%3]
                edge_key = sort_edge(i0, i1)
                if edge_key not in edges_dict:
                    edges_dict[edge_key] = []
                edges_dict[edge_key].append(face_i)
        
        # for existing vertices (new pos)
        new_vertex_positions = np.zeros_like(V)
        for i, neighbors in adj.items():
            n = len(neighbors)
            if n > 0:
                if n == 3:
                    beta = 3/16
                else:
                    beta = 3/(8*n) 
                
                new_pos = (1 - n*beta) * V[i]
                for neighbor in neighbors:
                    new_pos += beta * V[neighbor]
                
                new_vertex_positions[i] = new_pos
            else:
                new_vertex_positions[i] = V[i] # keep position if no neighbors
        
        # Create new vertices at the midpoints of edges
        new_edge_verts = {}
        new_vertices = list(new_vertex_positions)
        
        for edge_key, faces in edges_dict.items():
            i0, i1 = edge_key
            p0, p1 = V[i0], V[i1]
            
            if len(faces) == 2:
                oppA = None
                oppB = None
                for face_idx in faces:
                    face = F[face_idx]
                    for v in face:
                        if v != i0 and v != i1:
                            if oppA is None:
                                oppA = v
                            else:
                                oppB = v
                                break
                
                p2 = V[oppA]
                p3 = V[oppB]
                
                pm = 3/8 * (p0 + p1) + 1/8 * (p2 + p3)
            else:
                pm = 0.5 * (p0 + p1)
            
            new_edge_idx = len(new_vertices)
            new_vertices.append(pm)
            new_edge_verts[edge_key] = new_edge_idx
        
        new_faces = []
        for face in F:
            v0, v1, v2 = face
            
            e01 = new_edge_verts[sort_edge(v0, v1)]
            e12 = new_edge_verts[sort_edge(v1, v2)]
            e20 = new_edge_verts[sort_edge(v2, v0)]
            
            new_faces.extend([
                [v0, e01, e20],
                [e01, v1, e12],
                [e20, e12, v2],
                [e01, e12, e20]
            ])
        
        return np.array(new_vertices, dtype=np.float32), np.array(new_faces, dtype=np.int32)

mesh = trimesh.load_mesh("./Aimshape2D/20_cow2.off")
mesh1 = trimesh.load_mesh("./Aimshape2D/23-Egea/23_egea.off")
mesh2 = trimesh.load_mesh("./Aimshape2D/31-Horse/31_horse.off")

loop = Loop_alg()
V_sub, F_sub = loop.subdivide_loop(mesh.vertices, mesh.faces, iterations=1)
V_sub1, F_sub1 = loop.subdivide_loop(mesh1.vertices, mesh1.faces, iterations=1)
V_sub2, F_sub2 = loop.subdivide_loop(mesh2.vertices, mesh2.faces, iterations=1)
print(f"Vertices: {len(mesh.vertices)}")
print(f"Subdivided vertices: {len(V_sub)}")
print(f"Faces: {len(mesh.faces)}")
print(f"Subdivided faces: {len(F_sub)}")

ps.init()

offset = np.array([1, 0.0, 0.0])
mesh1.vertices = mesh1.vertices + offset 
V_sub1 = V_sub1 + offset

offset = np.array([0.0, 1, 0.0])
mesh2.vertices = mesh2.vertices + offset 
V_sub2 = V_sub2 + offset

mesh_orig = ps.register_surface_mesh("Original Mesh Cow", mesh.vertices, mesh.faces)
mesh_sub = ps.register_surface_mesh("Loop Subdivided Mesh Cow", V_sub, F_sub)

mesh_orig1 = ps.register_surface_mesh("Original Mesh Head", mesh1.vertices, mesh1.faces)
mesh_sub1 = ps.register_surface_mesh("Loop Subdivided Mesh Head", V_sub1, F_sub1)

mesh_orig2 = ps.register_surface_mesh("Original Mesh Horse", mesh2.vertices, mesh2.faces)
mesh_sub2 = ps.register_surface_mesh("Loop Subdivided Mesh Horse", V_sub2, F_sub2)

def callback():
    mesh_enabled = True
    psim.PushItemWidth(100)

    clicked, new_state = psim.Checkbox("Show Original Meshes", mesh_enabled)
    
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
ps.show()
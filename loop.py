# APPROXIMATING
# 



import numpy as np

def subdivide_mid_edge(V, F, iterations=1):
    V_sub, F_sub = V.copy(), F.copy()
    for _ in range(iterations):
        V_sub, F_sub = _subdivide_mid_edge_once(V_sub, F_sub)
    return V_sub, F_sub

def _subdivide_mid_edge_once(V, F):
    # Similar to butterfly but we only place midpoints, no fancy stencils
    edges_dict = {}
    def sort_edge(a,b):
        return (a,b) if a<b else (b,a)
    
    new_verts = list(V)
    new_edge_verts = {}
    new_faces = []
    
    for face in F:
        mid_idx = []
        for e in range(3):
            i0 = face[e]
            i1 = face[(e+1)%3]
            ek = sort_edge(i0, i1)
            if ek not in new_edge_verts:
                pm = 0.5*(V[i0] + V[i1])
                new_idx = len(new_verts)
                new_verts.append(pm)
                new_edge_verts[ek] = new_idx
            else:
                new_idx = new_edge_verts[ek]
            mid_idx.append(new_idx)
        
        v0, v1, v2 = face
        m01, m12, m20 = mid_idx
        # 4 new triangles
        new_faces.append([v0, m01, m20])
        new_faces.append([m01, v1, m12])
        new_faces.append([m20, m12, v2])
        new_faces.append([m01, m12, m20])
    
    return np.array(new_verts, dtype=np.float32), np.array(new_faces, dtype=np.int32)

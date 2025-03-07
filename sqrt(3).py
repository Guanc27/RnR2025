import numpy as np

def subdivide_sqrt3(V, F, iterations=1):
    """
    Pseudocode for Kobbelt's sqrt(3) subdivision on a triangular mesh.
    """
    V_sub, F_sub = V.copy(), F.copy()
    for _ in range(iterations):
        V_sub, F_sub = _subdivide_sqrt3_once(V_sub, F_sub)
    return V_sub, F_sub

def _subdivide_sqrt3_once(V, F):
    # 1) Insert face centroids
    centroids = (V[F[:,0]] + V[F[:,1]] + V[F[:,2]]) / 3.0
    base_idx = len(V)
    newV = np.vstack([V, centroids])  # (N + M, 3)
    
    # 2) Split each face into 3 smaller triangles
    newF = []
    for i, face in enumerate(F):
        c_idx = base_idx + i
        v0, v1, v2 = face
        newF.append([v0, v1, c_idx])
        newF.append([v1, v2, c_idx])
        newF.append([v2, v0, c_idx])
    
    newF = np.array(newF, dtype=np.int32)
    
    # 3) Edge flips:
    #    Typically, we flip the "old" edges (the ones in the original mesh).
    #    Implementation details can be tricky, requiring a halfedge structure to identify diagonals.
    #    For demonstration, we only show the conceptual step:
    #    
    #    for each old edge e in the original mesh:
    #        check the two new triangles that share e in newF
    #        flip e if it forms a diamond shape
    #    
    #    We'll skip explicit code here for brevity.
    
    return newV, newF

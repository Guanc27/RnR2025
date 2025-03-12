"""
Compare Butterfly, Loop, and √3 subdivision algorithms on:
  1) Computational time -> faster is better.
  2) Mesh complexity (#vertices, #edges, #faces) -> higher complexity may give more detail but costs processing time.
  3) Geometric differences (distance from original) -> preserving features vs. better smoothness.
  4) Qualitative properties (e.g. stationarity, primal vs. dual).
"""

import time
import numpy as np
import trimesh
from butterfly import Butterfly_alg
from loop import Loop_alg
from refinedSqrt3 import Sqrt3Subdivision  

butter = Butterfly_alg()
loop = Loop_alg()
sqrt = Sqrt3Subdivision()

def compute_num_edges(faces):
    """
    For a triangular mesh, count unique edges from the face array.
    This is to compute # of edges.
    """
    edges_set = set()
    for face in faces:
        e1 = tuple(sorted([face[0], face[1]]))
        e2 = tuple(sorted([face[1], face[2]]))
        e3 = tuple(sorted([face[2], face[0]]))
        edges_set.update([e1, e2, e3])
    return len(edges_set)

def compare_subdivision_algorithms(V, F, iterations, samples):
    mesh_original = trimesh.Trimesh(vertices=V, faces=F, process=False)
    
    t0 = time.time()
    Vb, Fb = butter.subdivide_butterfly(V, F, iterations)
    t_butterfly = time.time() - t0

    nbV = len(Vb)
    nbF = len(Fb)
    nbE = compute_num_edges(Fb)

    t1 = time.time()
    Vl, Fl = loop.subdivide_loop(V, F, iterations)
    t_loop = time.time() - t1

    nlV = len(Vl)
    nlF = len(Fl)
    nlE = compute_num_edges(Fl)

    t2 = time.time()
    Vs, Fs = sqrt.subdivide(V, F, iterations)
    t_sqrt = time.time() - t2

    nsV = len(Vs)
    nsF = len(Fs)
    nsE = compute_num_edges(Fs)
    mesh_butterfly = trimesh.Trimesh(vertices=Vb, faces=Fb, process=False)
    mesh_loop      = trimesh.Trimesh(vertices=Vl, faces=Fl, process=False)
    mesh_sqrt      = trimesh.Trimesh(vertices=Vs, faces=Fs, process=False)

    pts_butterfly = mesh_butterfly.sample(samples)
    pts_loop      = mesh_loop.sample(samples)
    pts_sqrt      = mesh_sqrt.sample(samples)

    dist_butterfly = mesh_original.nearest.on_surface(pts_butterfly)[1]
    dist_loop      = mesh_original.nearest.on_surface(pts_loop)[1]
    dist_sqrt      = mesh_original.nearest.on_surface(pts_sqrt)[1]

    avg_butterfly = np.mean(dist_butterfly)
    max_butterfly = np.max(dist_butterfly)
    avg_loop      = np.mean(dist_loop)
    max_loop      = np.max(dist_loop)
    avg_sqrt      = np.mean(dist_sqrt)
    max_sqrt      = np.max(dist_sqrt)

    print(f"\n--- Compare Subdivision (iterations={iterations}) ---")
    print("Butterfly Subdivision:")
    print(f"  Time (s):         {t_butterfly:.4f}")
    print(f"  #Vertices:        {nbV}")
    print(f"  #Edges:           {nbE}")
    print(f"  #Faces:           {nbF}")
    print(f"  Avg dist to orig: {avg_butterfly:.5f}")
    print(f"  Max dist to orig: {max_butterfly:.5f}")

    print("\nLoop Subdivision:")
    print(f"  Time (s):         {t_loop:.4f}")
    print(f"  #Vertices:        {nlV}")
    print(f"  #Edges:           {nlE}")
    print(f"  #Faces:           {nlF}")
    print(f"  Avg dist to orig: {avg_loop:.5f}")
    print(f"  Max dist to orig: {max_loop:.5f}")

    print("\n√3 Subdivision:")
    print(f"  Time (s):         {t_sqrt:.4f}")
    print(f"  #Vertices:        {nsV}")
    print(f"  #Edges:           {nsE}")
    print(f"  #Faces:           {nsF}")
    print(f"  Avg dist to orig: {avg_sqrt:.5f}")
    print(f"  Max dist to orig: {max_sqrt:.5f}")

    return {
        'Butterfly': {
            'time': t_butterfly,
            'V': nbV, 'E': nbE, 'F': nbF,
            'avg_dist': avg_butterfly, 'max_dist': max_butterfly
        },
        'Loop': {
            'time': t_loop,
            'V': nlV, 'E': nlE, 'F': nlF,
            'avg_dist': avg_loop, 'max_dist': max_loop
        },
        'Sqrt3': {
            'time': t_sqrt,
            'V': nsV, 'E': nsE, 'F': nsF,
            'avg_dist': avg_sqrt, 'max_dist': max_sqrt
        }
    }
    
def qual():
    print("\n--- Qualitative Comparison ---")
    print("Algorithm      Stationary  Interpol./Approx.  Primal/Dual  Smoothness")
    print("---------------------------------------------------------------------")
    print("Butterfly      yes         interpolatory       primal       C^1 in reg. areas")
    print("Loop           yes         approximating       primal       C^2 in reg. areas")
    print("√3             yes         interpolatory       primal       (varies, ideally C^1)")

if __name__ == "__main__":
    print("Script started")
    print("Loading mesh from './Aimshape2D/20_cow2.off'...")
    
    mesh = trimesh.load_mesh("./Aimshape2D/20_cow2.off")
    print("Mesh loaded: {} vertices, {} faces".format(len(mesh.vertices), len(mesh.faces)))
    
    # Set the sample count to 10 times the number of faces to thoroughly sample the surface.
    sample_count = len(mesh.faces) * 10
    print("Sampling {} points on the subdivided surface for geometry comparison...".format(sample_count))
    
    print("Comparing subdivision algorithms...")
    compare_results = compare_subdivision_algorithms(mesh.vertices, mesh.faces, iterations=2, samples=sample_count)
    
    print("Performing qualitative comparison...")
    qual()
    
    print("Comparison complete.")

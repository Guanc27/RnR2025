# Learned Info About CMD

**Caching:**  
  - Cache consists of locally stored copies, so calling `thingi10k.init()` won't re-download everything and allows downloads to resume.
- **Installation Tips:**  
  - If you can’t install with conda, try using pip and vice versa.
- **Handling Corrupted Files:**  
  - If you encounter issues, uninstall and then reinstall (e.g., for polyscope).

**General notes**
  - Tools/Libraries:
    - vedo, trimesh, matplotlib, or use `igl.writeobj` to convert to `.obj` file.
  - trimesh/polyscope registers TrackedArray (although v similar to array), the former enables caching, speeds up applications in processing vertices/faces
  
  
# Project PAIN Points

- **C/C++ Challenges:**  
  - Attempted using C/C++ initially, but had problems with dependencies (whether to use MingGW 64, MSVC for Microsoft VS, etc.)

- **Python/Conda Issues:**  
  - Switched to conda instead but after updating Python to 3.10 for thingi10k, ran into further dependency problems with VTK/vedo.
  - Encountered igl import errors and conflicts.
  - **Butterfly Subdivision:**  
    - Issues with matching new vertex indices with those of new edges.
    - Laptop can't handle more than 7 iterations. 

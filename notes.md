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
  

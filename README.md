## Project Outline

- Interested in using various subdivision algorithms to smoothen meshes for better processing and visuals.
- Curious as to how different existing algorithms compare with one another in terms of geometrical properties, runtime, memory space, etc.

## General notes

- Tools/Libraries:
  - vedo, trimesh, matplotlib, or use `igl.writeobj` to convert to `.obj` file.

- **Existing Databases:**
  - **RWTT:**  
    - Uses photo reconstruction methods.
    - Subdivision algorithms can exacerbate problems due to too many vertices.
  - **Digital Michelangelo Project:**  
    - 3D scanning of large figures.
    - Digitized 10 statues by Michelangelo.
  - **Stanford 3D Scanning Repository:**  
    - Provides raw 3D scans.
    - Their zippering and volumetric range image merging methods produce smooth manifolds and already reduce noise.
      - *Range Imaging:*  
        - Captures 3D shapes by measuring the distance from sensor to object.
        - Reconstructs a 3D model from partial scans.
        - The first step is to align each range image in the same coordinate system, but gaps may remain; hence, zippering them together is used, which is not always ideal for surface reconstruction.
      - *Additional Note:*  
        - Some 3D scanning techniques capture connectivity, so treating them as point clouds (i.e., just x, y, z coordinates) discards useful information.
  - Alternatively, use standard methodologies like quadric edge collapse decimation and voxel-based remeshing.
  - **Polytechnique**
    - Multiple datasets of triangular meshes (from aim@shape and Thingi10k repo)
      - genus 0, so surface closed, topologically equiv to spheres but triangulations

- **Testing Subdivision Algorithms:**
  - Either find a database with irregular meshes or decompose a mesh.
  - **Curvature Metrics:**  
    - **Gaussian Curvature:**  
      - Product of the two principal curvatures.
      - Gives insight into local geometric properties (e.g., ellipse, saddle-like shapes).
      - Example: Pre-div (-95870), post-div (-1.6e+8).
    - **Mean Curvature:**  
      - Average of the two principal curvatures.
      - Example: Pre-div (0.2), post-div (0).
  - To test smoothness, consider other metrics such as normal variation.
  - Subdivision either interpolation (includes original control pts like butterfly) or approximating (doesn't include the original points)

  # Project PAIN Points

- **C/C++ Challenges:**  
  - Attempted using C/C++ initially, but had problems with dependencies (whether to use MingGW 64, MSVC for Microsoft VS, etc.)

- **Python/Conda Issues:**  
  - Switched to conda instead but after updating Python to 3.10 for thingi10k, ran into further dependency problems with VTK/vedo.
  - Encountered igl import errors and conflicts.
  - **Butterfly Subdivision:**  
    - Issues with matching new vertex indices with those of new edges.
    - Laptop can't handle more than 7 iterations. 

## References

- [Stanford 3D Scanning Repository](https://graphics.stanford.edu/data/3Dscanrep/)
- (Unrelated but for surface construction) [Stanford CS468 Lecture Slides on Surface Reconstruction](https://graphics.stanford.edu/courses/cs468-12-spring/LectureSlides/03_Surface_Reconstruction.pdf)
- [Research Paper on Subdivision Algorithms](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=99ca8274377ee438fbb748438aa3057e7f6654a2)
- http://www.math.tau.ac.il/~niradyn/papers/butterfly.pdf
- https://people.eecs.berkeley.edu/~sequin/CS284/PAPERS/root3subdiv.pdf
- 
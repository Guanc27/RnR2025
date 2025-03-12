## Project Outline

- Interested in using various subdivision algorithms to smoothen meshes for better processing and visuals.
- Curious as to how different existing algorithms compare with one another in terms of geometrical properties, runtime, memory space, etc.
- Understanding different subdivision categories and pros/cons 

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
      - https://www.lix.polytechnique.fr/~amturing/Software/Datasets/Datasets.html
  - **Computer Graphics Group at MIT**
    - https://people.csail.mit.edu/sumner/research/deftransfer/data.html

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
  - Convergent analysis on the refinement meshes

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

- [Research Paper on Subdivision Algorithms](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=99ca8274377ee438fbb748438aa3057e7f6654a2)
- [Butterfly scheme] http://www.math.tau.ac.il/~niradyn/papers/butterfly.pdf
- [Interpolatory sqrt(3) subdivision] https://people.eecs.berkeley.edu/~sequin/CS284/PAPERS/root3subdiv.pdf
- [Loop subdivision] https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/thesis-10.pdf
- [Original sqrt(3) subdivision] https://dl.acm.org/doi/pdf/10.1145/344779.344835 
- https://cs.engr.uky.edu/~cheng/PUBL/Book_CCSS.pdf
- [General info] https://people.computing.clemson.edu/~dhouse/courses/405/notes/Catmull-Clark.pdf
- [Alg comparison] https://pubs.aip.org/aip/acp/article/1487/1/343/855266/A-comparison-of-surface-subdivision-algorithms-for?pdfCoverIconEvent=cite
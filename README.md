# Anatomy-Volume-Analysis

The MATLAB class **Anatomy Volume Analysis** provides a volume-oriented workflow for anatomical surface meshes.  
It detects boundary loops, patches holes, generates shrink-wrapped surfaces, reconstructs surfaces from point clouds, and fills regions with structured point sets for intersection/volume computations.  
The volume class is for geometry-processing tasks where robust, watertight models and volumetric point sets are required.

---

## Key Features
1. **Edge Detection:** Boundary vertex/edge extraction and automatic loop grouping  
2. **Hole Patching:** Fill holes with triangulated patches  
3. **Point Filling:** Fill reference and anatomical model with points (grid points)  
4. **Shrink-Wrap Surfaces:** Alpha shape hulls with UI or auto iteration to closed meshes  
5. **Reference Intersection Support:** Identification of reference model parts within alpha shape hulls and assignment of grid points for subsequent volume/intersection steps  
6. **Point Filtering:** Remove grid points outside anatomical model  
7. **Volume Reconstruction:** Refined alpha shape hull with grid points  

---

## How to Use
1. **Code structure:** The volume class can be integrated into a corresponding main script. The key parts of the code are extracted into individual functions for reuse. Provide the mesh as an STL file (triangulated).  
2. **Import:** Import the anatomical model data and any reference model in STL format.  
3. **Run Volume Calculations**  
4. **Results:** Use the grid points and closed hulls for voxel/point-based volume or intersection analysis subsequently.  

---

## Use Cases
- Generating watertight hulls from surface meshes  
- Preparing regions for voxel/point intersection and volume computation  
- Reconstructing surfaces from point clouds  
- Computing distance maps for quantitative mesh comparison  

---

## Benefits
- **Robust to real data:** Handles holes, open boundaries, and noisy points  
- **Interactive where needed:** Quick UI checks  
- **Automatable:** Workflow with auto-closure of alpha shapes  
- **Volume-ready outputs:** Inside/outside masks, alpha shapes, and dense grid points  
- **Modular:** Single functions can be dropped into custom workflows  

---

## Keywords
`anatomy`, `volume analysis`, `shrink wrap`, `alpha shape`, `mesh repair`, `hole filling`, `boundary detection`, `clustering`, `KDE`, `DBSCAN`, `point cloud`, `defect filling`, `watertight mesh`, `biomedical geometry`

---

## References
- N. Douillet, *‘‘Discrete contour mesh patch,’’* MATLAB Central File Exchange, 2025. [Online]. Available: [https://de.mathworks.com/matlabcentral/fileexchange/78901-discrete-contour-mesh-patch-2d-3d](https://de.mathworks.com/matlabcentral/fileexchange/78901-discrete-contour-mesh-patch-2d-3d)  
- Sven, *‘‘inpolyhedron - Are points inside a triangulated volume?’’* MATLAB Central File Exchange, 2024. [Online]. Available: [https://de.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume](https://de.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume)  
- S. de Wolski, *‘‘File Exchange Pick of the Week - INPOLYHEDRON,’’* 2013. [Online]. Available: [https://blogs.mathworks.com/pick/2013/09/06/inpolyhedron/](https://blogs.mathworks.com/pick/2013/09/06/inpolyhedron/)  

---

## Cite As
Carina Micheler (2025). [Anatomy Volume Analysis](https://www.mathworks.com/matlabcentral/fileexchange/181321-anatomy-volume-analysis), MATLAB Central File Exchange.

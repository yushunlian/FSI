# FSI
A C++ library for performing fluid-structure interaction using OpenFOAM's fluid solver and multi-degree of freedom(MDOF) system solver. 

* **MDOFBodyMotion** directory contains classes that hold the MDOF system properties like, mass, stifness, damping, motion state, etc. All the structural properties are stored in this clases including manipulating functions that operate on the data of the clases.  

* **MDOFBodyMotionSolver** holds a class for mesh motion solver. This class calculates and applies displacement to the grid in the nearfield of the structure. The class holds the patchs, interpolation distance and scale used to calculate the grid displacement.   

* **MDOFSolver** contains classes that are used to solve the MDOF systme given the forces and structural properties. This classes hold the MDOFBodyMotion object numerical scheme used to advace the solution of the structural systme like CentralDifference, CrackNicolsen, ...



## Plan












**Note**
- There are some modifications on the OpenFOAM-dev version for some of the clases. Future development need to acomodate these chages and fix if there is any backward compatability issue. The change includes class and inheritance structure for clases in **MDOFBodyMotionSolver** directory.


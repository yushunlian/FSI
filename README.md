# FSI
A C++ library for performing fluid-structure interaction using OpenFOAM's fluid solver and multi-degree of freedom(MDOF) system solver. 

* **mDoFSystemMotion** directory contains classes that hold the MDOF system properties like, mass, stifness, damping, motion state, etc. All the structural properties are stored in this clases including manipulating functions that operate on the data of the clases.  

* **mDoFSystemMotionSolver** holds a class for mesh motion solver. This class calculates and applies displacement to the grid in the nearfield of the structure. The class holds the patchs, interpolation distance and scale used to calculate the grid displacement.   

* **mDoFSolver** contains classes that are used to solve the MDOF systme given the forces and structural properties. This classes hold the MDOFBodyMotion object numerical scheme used to advace the solution of the structural systme like CentralDifference, CrackNicolsen, ...

* **mDoFSystemControl** Convergence control based on the mDoF motion state. The body linear and angular velocities are averaged over a specified time window and compared to specified velocity convergence criteria and the run stopped after writing the current time results, if the criteria are met. 

* **mDoFSystemState** Writes the mDoF motion state.

* **trush** Contains old files which have been deleted.

**Order of calling**
- *dynamicMesh.update()* -> *dynamicMotionSolverFvMesh.update()* -> *motionSolver.newPoints()* -> *displacementMotionSolver* -> *mDoFSystemMotionSolver.solve()* -> *mDoFSystemMotion.update(args)* -> *mDoFSolver.solve(args)*


**Note**
- [x] There are some modifications on the OpenFOAM-dev version for some of the clases showen above. Future development need to acomodate these chages and fix if there is any backward compatability issue. The changes include a new class inheritance structure for clases in **mDoFSystemMotionSolver** directory and others.


Task List
===========

- [ ] In the **mDoFSystemMotion** classes create properties and functions that will hold multi-degree-of-freedom data(structural data) for a building structure. The structural data is stored as an array each intries representing a value at a particular floor and values that are universal like damping, number of modes to include, etc. The story data includes floor elevation, center of mass, mass, torsional radius vector, mode shape matrix, diagonal matrix stifnness matrix, etc. 

- [ ] Implement a MDOF solver in **mDoFSolver** classes. The solver uses explicity/implicit technique to advance the lamped-mass system for the current time step. The load for the solver comes from the forces calculated on a patches that surounds each lamped mass. Each mass need to have it's own surounding patch look into the new update in [OpenFoam-8](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/src/rigidBodyDynamics/rigidBodyModel/rigidBodyModel.H). The forces on each patch is extracted using the 'forces' post processing tool of OpenFOAM look the implementation in [sixDoFRigidBodyMotionSolver.solve()](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/src/sixDoFRigidBodyMotion/sixDoFRigidBodyMotionSolver/sixDoFRigidBodyMotionSolver.C) function. Follow the following procedure to create the patches for each floor. 
  - First, create an "\*stl" surface deviding each floor for the building as a separate file and define a boundary condition for each in the corresponding velocity and pressure field files. 
  - List the patches in [dynamicMeshDict](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/tutorials/incompressible/pimpleFoam/RAS/wingMotion/wingMotion2D_pimpleFoam/constant/dynamicMeshDict) dictionary so that the motionSolver reads them. 
  - Create a separate forces function object for each patch in the [sixDoFRigidBodyMotionSolver.solve()](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/src/sixDoFRigidBodyMotion/sixDoFRigidBodyMotionSolver/sixDoFRigidBodyMotionSolver.C) file and feed to *mDoFSolver.solve(args)* function as 
  ```cpp
  
        dictionary forcesDicts(patches_.size());
        
        forAll(patches_, i)
        {
          forcesDict.add("type", functionObjects::forces::typeName);
          forcesDict.add("patches", patches_);
          forcesDict.add("rhoInf", rhoInf_);
          forcesDict.add("rho", rhoName_);
          forcesDict.add("CofR", centreOfRotation());

          functionObjects::forces f("forces", t, forcesDict);

          f.calcForcesMoment();
          ...
        }
  ```

- [ ] In the **mDoFSystemMotionSolver** classes first, implement a method that uses the displacement and rotation at each floor level and interpolate for intermediate elevations and create a fully diformed point field of the structural body(outer skin). Then, using is deformed configuraton defiend by the point field, solve the near field motion of the fluid grid.   

- [ ] Develope a **fsiPimpleFoam** solver that performes iterations in the fluid and MDOF solver. This solver needs to be implemented in $FOAM_APP folder as incompressible solver modiying the existing **pimpleFoam** solver. 






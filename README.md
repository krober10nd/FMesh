# FMesh
A 2D/3D, an-/isotropic triangular/tetrahedral mesh generator based on the DistMesh algorithm in modern Fortran 
<br /> 
<br /> 
<br /> 
(NOTE: this is a work in progress!).
=======================================


The software will need: 
---------------------------------------
1. GNU compiler (e.g., gfortran and gcc)
2. Qhull (static re-entrant libraries)
3. Python3+ (for visualization of mesh and for development of mesh size functions)

To install the software:
---------------------------------------
1. Enter into the work directory.
2. In the makefile, change the path to the include directory (CPP) and library directory (LDFLAGS). 
3. In make.inc, change the type of build configuration (i.e., dbg, optimised, etc.) 
4. Run make and the generated binary should be installed in work/bin.

To use the software:
---------------------------------------
It requires a Planar Straight Line graph (in 2D) or PSLG that describes the boundary of the domain. 
% The PSLG is a text file called "PSLG.txt" in the working directory in the following format: 

5  2      ! # of points and dimension <br /> 
-1.0  1.0 ! x and y coordinates of each point <br /> 
 1.0  1.0 <br /> 
 1.0 -1.0 <br /> 
-1.0 -1.0 <br /> 
-1.0  1.0 <br /> 
<br /> 
It requires the user define their own mesh size function in the module src/YourMeshSizeFunction.F90. The mesh size function takes a vector of coordinates and returns a desired element size and orientation (i.e., anisotropic) In the simplest case, you can simply write a function to return only a element size (isotropic). 

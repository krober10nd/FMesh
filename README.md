# FMesh
A 2D/3D, an-/isotropic triangular mesh generator based on the DistMesh algorithm in modern Fortran 
<br /> 
<br /> 
<br /> 
(NOTE: this is a work in progress!).
=======================================


The software will need: 
---------------------------------------
1. GNU compiler (e.g., gfortran and gcc)
2. Qhull (static re-entrant libraries)

To install the software:
---------------------------------------
1. Enter into the work directory.
2. In the makefile, change the path to the include directory (CPP) and library directory (LDFLAGS). 
3. In make.inc, change the type of build configuration (i.e., dbg, optimised, etc.) 
4. Run make and the generated binary should be installed in work/bin.

To use the software:
---------------------------------------
It requires a Planar Straight Line graph (in 2D) or PSLG that describes the boundary of the domain in either ccw or cw order. 
% The PSLG is a text file called "PSLG.txt" in the working directory in the following format: 

5  2      ! # of points and dimension <br /> 
-1.0  1.0 ! x and y coordinates of each point <br /> 
 1.0  1.0 <br /> 
 1.0 -1.0 <br /> 
-1.0 -1.0 <br /> 
-1.0  1.0 <br /> 
<br /> 

Currently the domain has to be a singly-connected polygon (but it likely easy to modify to support mutliply-connected polygons).

It requires the user define their own mesh size function. These are rasters/structured grids than span the entire meshing domain with a minimum grid spacing at least twice as small as the minimum element size. Three grid files that are required are 1) the size of the mesh in the major axis of the circum-ellipses, 2) the size of the mesh in the minor axis of the circum-ellipses, and 3) the angle the major axis of the circum-ellipse makes with the x-axis (in radians). Please see MakeMeshSizes.m for more details in creating these files. 

The program is run through the command line type (assuming you've named the files described above like below): 

./distmesh.x PSLG.txt MeshSize1.txt MeshSize2.txt MeshAngle.txt 

An example metric tensor is: 

x,y spans [-2 2]x[-2 2] <br /> 
MeshSize1(x,y) = 0.005 + 1.5*abs(1-(x.^2  + y.^2 ).^0.5);  
MeshSize2(x,y) = 0.1*(x.^2 + y.^2).^0.5 + 1.5*abs(1-(x.^2  + y.^2 ).^0.5); <br /> 
Angle(x,y) = atan(x(x,y)/y(x,y)) + 90  ; <br /> 



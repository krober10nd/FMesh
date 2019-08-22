# FMesh
A distance-based mesh generator (based on the DistMesh algorithm) written in modern Fortran (Note: this is a work in progress!).

The software will need: 
1. GNU compiler (e.g., gfortran and gcc)
2. Qhull (static re-entrant libraries)
3. Python3+ (for visualization of mesh)

To install the software:
1. Enter into the work directory.
2. In the makefile, change the path to the include directory (CPP) and library directory (LDFLAGS). 
3. In make.inc, change the type of build configuration (i.e., dbg, optimised, etc.) 
4. Run make and the generated binary should be installed in work/bin.



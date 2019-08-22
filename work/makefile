INCLUDE  := -I../include/libqhull_r
LDFLAGS  := -L../lib
LDLIBS   := -lqhullstatic_r
# fortran flags 
GFLAGS   := -O0 -fbacktrace -g -DREAL64 
# c flags
CFLAGS   := -O2 -g

distmesh.x: ctriangulate.o utils.o distmesh.o
	gfortran -o distmesh.x ctriangulate.o utils.o distmesh.o $(LDFLAGS) $(LDLIBS) $(INCLUDE) $(CFLAGS)  
utils.mod: utils.o utils.F90 ctriangulate.c
	gfortran -c $(GFLAGS) utils.F90
utils.o: utils.F90
	gfortran -c $(GFLAGS) utils.F90
distmesh.o: utils.mod distmesh.F90
	gfortran -c $(GFLAGS) distmesh.F90
ctriangulate.o: ctriangulate.c
	gcc -c $(CFLAGS) $(INCLUDE) ctriangulate.c
clean:
	rm utils.mod distmesh.o utils.o ctriangulate.o distmesh.x

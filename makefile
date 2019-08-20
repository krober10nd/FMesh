INCLUDE  := -I/Users/Keith/Desktop/include/libqhull_r
LDFLAGS  := -L/Users/Keith/Desktop/lib
LDLIBS   := -lqhullstatic_r
# fortran flags 
GFLAGS   := -O0 -g -fbacktrace -fbounds-check -cpp -DREAL64
# c flags
CFLAGS   := -O0 -g

testqhull.x: ctriangulate.o utils.o driver.o
	gfortran -o testqhull.x ctriangulate.o utils.o driver.o $(LDFLAGS) $(LDLIBS) $(INCLUDE) $(CFLAGS)  
utils.mod: utils.o utils.F90 ctriangulate.c
	gfortran -c $(GFLAGS) utils.F90
utils.o: utils.F90
	gfortran -c $(GFLAGS) utils.F90
driver.o: utils.mod driver.f90
	gfortran -c $(GFLAGS) driver.f90
ctriangulate.o: ctriangulate.c
	gcc -c $(CFLAGS) $(INCLUDE) ctriangulate.c
clean:
	rm utils.mod driver.o utils.o ctriangulate.o testqhull.x

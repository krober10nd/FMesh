# Start of the makefile
CFLAGS  := -I/Users/Keith/Desktop/include/libqhull_r
LDLIBS  := -lqhullstatic_r
LDFLAGS := -L/Users/Keith/Desktop/lib
GFLAGS   := -O0 -g -fbacktrace -fbounds-check -cpp

testqhull.x: ctriangulate.o utils.o driver.o
	gfortran -o testqhull.x $(GFLAGS) $(LDLIBS) $(LDFLAGS) ctriangulate.o utils.o driver.o
utils.mod: utils.o utils.f90 ctriangulate.c
	gfortran -c utils.f90
utils.o: utils.f90
	gfortran -c utils.f90
driver.o: utils.mod driver.f90
	gfortran -c driver.f90
ctriangulate.o: ctriangulate.c
	gcc -c $(CFLAGS) ctriangulate.c
clean:
	rm utils.mod driver.o utils.o ctriangulate.o testqhull.x
# End of the makefile

# Example make.inc file.  Not all options need to be set.

#-----
# Config info.

# Used to create a configuration-level and optimisation-level specific
# directory in which compiled objects are placed and to uniquely name binaries.

CONFIG = gnu
OPT = opt

#-----
# Compiler configuration.

# --- Preprocessor ---
# If CPP is defined, then this is used to do the preprocessing on Fortran
# source files and then the Fortran compiler is called.  If CPP is not defined,
# then CPPFLAGS is passed to the Fortran compiler and the Fortran source files
# are preprocessed and compiled in a single step.
# All C and C++ files are preprocessed and compiled in one step by the C and
# C++ compilers, with CPPFLAGS passed to the C and C++ compilers.
# NOTE: if you use cpp you will almost certainly have to provide the -C -P flags.
CPP = 

# Preprocessing flags
# Used for *.F, *.F90, *.c and *.cpp files.
CPPFLAGS = -I/usr/local/include/libqhull_r

# --- Fortran ---
# compiler
FC = gfortran
# compiler flags
FFLAGS = -O3 -fopenmp # -Wall -Wextra -pedantic
# flag to specify directory to put compiled .mod files (e.g. -module for ifort, -J for gfortran)
F90_MOD_FLAG = -J

# --- C ---
# compiler
CC = gcc
# compiler flags
CFLAGS = -O3 -fopenmp #-Wall -Wextra -pedantic
# Not all compilers (e.g. xlc) support -MM -MT that we require to produce the
# C dependencies.  Fortunately there is (usually) a version of (e.g.) gcc
# kicking around which does.  If CCD is defined, then use that to generate the
# C dependency files, otherwise CC is used.
CCD =

# --- C++ ---
# compiler
CXX = g++
# compiler flags
CXXFLAGS = -O3 -fopenmp
# Not all compilers (e.g. xlc++) support -MM -MT that we require to produce the
# C++ dependencies.  Fortunately there is (usually) a version of (e.g.) g++
# kicking around which does.  If CXXD is defined, then use that to generate the
# C++ dependency files, otherwise CXX is used.
CXXD =

# --- Linker ---
# linker
LD = gfortran
# linker flags
LDFLAGS = -L/usr/local/lib -fopenmp
# additional libraries to link to
LIBS = -lqhullstatic_r

# --- Archiver ---
# archive program (usually ar)
AR = ar
# archive flags (usually -rcs)
ARFLAGS = -rcs

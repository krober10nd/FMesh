!**************************************************
! serial DistMesh algorithm in modern FORTRAN
! to build mesh in 2d or 3d
! AUTHOR: kjr, universdade de sao paulo, 2019--
!**************************************************
PROGRAM DistMesh
USE utils 
IMPLICIT NONE

!! SPECIFY PROBELM SPECIFIC PARAMETERS HERE
DIM = 2
LMIN=0.25d0
!! 

! READ IN GEOMETRY 
CALL ReadPSLGtxt(PSLG,LMIN)

! STEP 1: Create initial triangulation to iterate on
CALL FormInitialTria2D(DIM,PSLG,LMIN,POINTS,NP,TRIAS,NF)

! STEP 2: Incrementally move points and retriangulate based 
! on "spring dynamics"




END PROGRAM 

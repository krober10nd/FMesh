!**************************************************
! serial DistMesh algorithm in modern FORTRAN
! to build mesh in 2d or 3d
! AUTHOR: kjr, universdade de sao paulo, 2019--
!**************************************************
PROGRAM DistMesh
USE utils 
IMPLICIT NONE

integer :: i 

!! SPECIFY PROBELM SPECIFIC PARAMETERS HERE
DIM = 2
LMIN=0.05d0
CALL ReadPSLGtxt(PSLG,LMIN) ! READ IN BOUNDARY DESCRIPTION 
!! 

! STEP 1-2: Create initial points to iterate on
CALL FormInitialPoints2D(DIM,PSLG,LMIN,POINTS,NP)

! STEP 3: Retriangulation by Delaunay algorithm 
CALL DelTriangulate(DIM,NP,POINTS,NF,TRIAS,IERR)

! DEBUG VISUALIZE INITIAL TRIANGULATION 
OPEN(UNIT=300,FILE="Points.txt",ACTION='WRITE')
DO i=1,NP
  WRITE(300,"(2F12.8)")POINTS(1,i),POINTS(2,i)
ENDDO
CLOSE(300)

OPEN(UNIT=301,FILE="Facets.txt",ACTION='WRITE')
DO i=1,NF
  WRITE(301,"(3I8)")TRIAS(1,i),TRIAS(2,i),TRIAS(3,i)
ENDDO
CLOSE(301)



END PROGRAM 

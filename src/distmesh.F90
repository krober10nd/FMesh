!**************************************************
! serial DistMesh algorithm in modern FORTRAN
! to build mesh in 2d or 3d
! AUTHOR: kjr, universdade de sao paulo, 2019--
!**************************************************
PROGRAM DistMesh
USE utils 
IMPLICIT NONE

INTEGER(idx_t) :: i

!! SPECIFY PROBELM SPECIFIC PARAMETERS HERE
DIM = 2
LMIN=0.25d0

BBOX(1) = -1.d0 
BBOX(2) = 1.d0 
BBOX(3) = -1.d0 
BBOX(4) = 1.d0
!!

CALL FormInitialPoints2D(DIM,BBOX,LMIN,POINTS,NP)

OPEN(UNIT=300,FILE="Points.txt",ACTION='WRITE')
DO i=1,NP
  WRITE(300,"(2F12.8)")POINTS(1,i),POINTS(2,i)
ENDDO
CLOSE(300)

CALL Triangulate(DIM,NP,POINTS,NF,TRIAS,IERR)

OPEN(UNIT=301,FILE="Facets.txt",ACTION='WRITE')
DO i=1,NF
  WRITE(301,"(3I8)")TRIAS(1,i),TRIAS(2,i),TRIAS(3,i)
ENDDO
CLOSE(301)



END PROGRAM 

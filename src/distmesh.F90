!**************************************************
! Serial DistMesh algorithm in FORTRAN
! to build a mesh in 2d or 3d.
! 
! AUTHOR: Keith J. Roberts, 
! PLACE: Universdade de Sao Paulo
! START: 2019-08-17 
!**************************************************
! NOTE: SWAP YOUR MESH SIZE IN THE MODULE YourMeshSize.F90
!
PROGRAM DistMesh
USE YourMeshSize, ONLY : MeshSize
USE utils 

IMPLICIT NONE

INTEGER I 

DIM = 2
LMIN=0.05

CALL ReadPSLGtxt(PSLG,LMIN)                                 ! Read in boundary description 

CALL FormInitialPoints2D(MeshSize,DIM,PSLG,LMIN,POINTS,NP)  ! Step 1-2: Create initial points to iterate on

WRITE(*,'(A)') "                                      "
WRITE(*,'(A)') "********************************************************"
WRITE(*,'(A)') "****************BEGIN ITERATING*************************"
WRITE(*,'(A)') "********************************************************"
WRITE(*,'(A)') "                                      "

! Main DistMesh loop
ALLOCATE(PointsOLD(NP,DIM)) 
PointsOld = -9999.0d0                                       ! For the first iteration 
DO

  IF(MAXVAL(SQRT(SUM((Points(1:NP,:)-PointsOld(1:NP,:))**2,2))/LMIN).GT.TTOL) THEN ! Any Large Movement?

    ! 3. Retriangulation by the Delaunay algorithm
    CALL DelTriaWElim(DIM,PSLG,NP,POINTS,NF,TRIAS,IERR)      
    WRITE(*,'(A,I8,A,I8,A)')"INFO: " & 
    //" THE MESH HAS ",NP," VERTICES " &
    //" AND ",NF," FACES "

    PointsOld=Points

    ! 4. Describe each bar by a unique pair of nodes
    CALL findUniqueBars(DIM,NF,TRIAS,NUMBARS,BARS)

    ! 5. Graphical output of the current mesh
    CALL WriteMesh(DIM,POINTS,NP,TRIAS,NF)

  ENDIF

  ! 6. Move mesh points based on bar lengths L and forces F
  CALL CalcForces(MeshSize,DIM,POINTS,NP,BARS,NUMBARS,FVEC)

  CALL ApplyForces(DIM,POINTS,NP,BARS,NUMBARS,FVEC) 
  EXIT 
  

ENDDO
!
!
!
END PROGRAM 

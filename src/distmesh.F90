!**************************************************
! serial DistMesh algorithm in modern FORTRAN
! to build mesh in 2d or 3d
! AUTHOR: kjr, universdade de sao paulo, 2019--
!**************************************************
PROGRAM DistMesh
USE utils 
USE YourMeshSize, ONLY : MeshSize

IMPLICIT NONE

DIM = 2
LMIN=0.05d0

CALL ReadPSLGtxt(PSLG,LMIN)                                 ! Read in boundary description 

CALL FormInitialPoints2D(MeshSize,DIM,PSLG,LMIN,POINTS,NP)  ! Step 1-2: Create initial points to iterate on

ALLOCATE(PointsOLD(NP,DIM)) 
PointsOld = -9999.0d0                                       ! For the first iteration 

DO
  IF(MAXVAL(SQRT(SUM((Points(1:NP,:)-PointsOld(1:NP,:))**2,2))/LMIN).GT.TTOL) THEN ! Any Large Movement?

    ! 3. Retriangulation by the Delaunay algorithm
    print *, "TRIA'ing NP ",NP," points "
    CALL DelTriaWElim(DIM,PSLG,NP,POINTS,NF,TRIAS,IERR)      
    print *, "PAST TRIA STEP "

    PointsOld=Points

    !% 4. Describe each bar by a unique pair of nodes
    CALL findUniqueBars(DIM,NF,TRIAS,NUMBARS,BARS)

    ! 5. Graphical output of the current mesh
    !CALL WriteMesh(DIM,POINTS,NP,TRIAS,NF)

    EXIT 

  ENDIF

  ! 6. Move mesh points based on bar lengths L and forces F
  !CALL CalculateEdgeLengths()

ENDDO
!
!
!
END PROGRAM 

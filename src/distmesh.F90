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
LMIN=0.0001d0

CALL ReadPSLGtxt(PSLG,LMIN)                                 ! Read in boundary description 

CALL FormInitialPoints2D(MeshSize,DIM,PSLG,LMIN,POINTS,NP)  ! Step 1-2: Create initial points to iterate on

STOP 

ALLOCATE(PointsOLD(NP,DIM)) 
PointsOld = -9999.0d0                                       ! For the first iteration 
DO
  ! 3. Retriangulation by the Delaunay algorithm
  IF(MAXVAL(SQRT(SUM((Points-PointsOld)**2,2))/LMIN).GT.TTOL) THEN ! Any Large Movement?

    CALL DelTriaWElim(DIM,PSLG,NP,POINTS,NF,TRIAS,IERR)      

    PointsOld=Points

    !% 4. Describe each bar by a unique pair of nodes
    !bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])]; % Interior bars duplicated
    !bars=unique(sort(bars,2),’rows’); %

    !SUBROUTINE WriteMeshData(DIM,POINTS,NP,FACETS,NF)
    ! 5. Graphical output of the current mesh
    CALL WriteMeshData(DIM,POINTS,NP,TRIAS,NF)
  ENDIF

  ! 6. Move mesh points based on bar lengths L and forces F
  CALL CalculateEdgeLengths()

ENDDO



END PROGRAM 

!**************************************************
! Serial DistMesh algorithm in FORTRAN
! to build a mesh in 2d or 3d.
! 
! AUTHOR: Keith J. Roberts, 
! PLACE: Universdade de Sao Paulo
! START: 2019-08-17 --
!**************************************************
! NOTE: SWAP THE MESH SIZE IN THE MODULE YourMeshSize.F90
!
PROGRAM DistMesh
USE YourMeshSize
USE utils 

IMPLICIT NONE

INTEGER I 
REAL(8) :: TS,TF
MaxIter =100 ! MAXIMUM NUMBER OF ITERATIONS

CALL ReadPSLGtxt(PSLG,LMIN)                                 ! Read in boundary description 

CALL FormInitialPoints2D(MeshSize,DIM,PSLG,LMIN,POINTS,NP)  ! Step 1-2: Create initial points to iterate on

CALL DelTriaWElim(DIM,PSLG,NP,POINTS,NF,TRIAS,IERR)         ! Compute Delaunay triangulation of point set with masking

CALL TriaToTria(NF,TRIAS,T2T,T2N)                               ! Calculate the triangle adj. matrix 

STOP 

WRITE(*,'(A)') "                                      "
WRITE(*,'(A)') "****************BEGIN ITERATING*************************"
WRITE(*,'(A)') "                                      "

ITER = 0 ! iteration counter 

DO ! distmesh loop
  CALL CPU_TIME(TS) 

  ! Perform edge flips to achieve Del. hood
  !CALL edgeFlips()  

  ! 4. Describe each bar by a unique pair of nodes
  CALL findUniqueBars(DIM,NF,TRIAS,NUMBARS,BARS)

  ! 5. Output of the current mesh
  IF(MOD(ITER,5).EQ.0) THEN
    CALL WriteMesh(DIM,POINTS,NP,TRIAS,NF,ITER)
  ENDIF

  ! 6. Calculate forces on bars
  CALL CalcForces(MeshSize,DIM,POINTS,NP,BARS,NUMBARS,FVEC)

  ! 7. Move points based on forces
  CALL ApplyForces(DIM,POINTS,NP,BARS,NUMBARS,FVEC) 

  ! 8. Bring outside points back to the boundary
  CALL ProjectPointsBack(DIM,PSLG,POINTS,NP)
  
  ITER = ITER + 1 
  WRITE(*,'(A,I4,A)') "INFO: ITERATION: ",ITER," COMPLETE"

  CALL CPU_TIME(TF) 
  WRITE(*,'(A,F12.8)') "INFO: ELAPSED TIME IS: ",TF-TS 
   
  ! 9. Termination criterion: reached max iterations. 
  IF(ITER.EQ.MaxIter) THEN
    WRITE(*,'(A)') "INFO: MAXIMUM NUMBER OF ITERATIONS REACHED..." 
    WRITE(*,'(A)') "********************************************************"
    EXIT 
  ENDIF

ENDDO

WRITE(*,'(A,I8,A,I8,A)')"INFO: " & 
//" THE MESH HAS ",NP," VERTICES " &
//" AND ",NF," FACES "
CALL WriteMesh(DIM,POINTS,NP,TRIAS,NF,ITER)

!
!
!
END PROGRAM 

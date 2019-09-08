!-----------------------------------------------------------------------
!> @brief Serial DistMesh algorithm in FORTRAN
!>        to build a mesh in 2d or 3d.
!-----------------------------------------------------------------------
PROGRAM DistMesh
!-----------------------------------------------------------------------
USE YourMeshSize ! contains your mesh size function queried during execution
USE utils 
USE vars 

IMPLICIT NONE
integer :: lastNumFlips
REAL(8) :: TS,TF

MaxIter = 100                                               ! Maximum number of iterations

CALL ReadPSLGtxt(PSLG,LMIN)                                 ! Read in boundary description 

CALL FormInitialPoints2D(MeshSize,DIM,PSLG,LMIN,POINTS,NP)  ! Create initial points to iterate on

CALL DelTriaWElim(DIM,PSLG,NP,POINTS,NF,TRIAS)              ! Compute Delaunay triangulation of point set with masking

CALL TriaToTria(NF,TRIAS,T2T,T2N)                           ! Calculate the triangle adj. matrices

WRITE(*,'(A)') "                                      "
WRITE(*,'(A)') "****************BEGIN ITERATING*************************"
WRITE(*,'(A)') "                                      "

ITER = 1 ! iteration counter 

DO 
  CALL CPU_TIME(TS) 

  ! Perform edge flips to achieve Del.-hood
  NUMFLIPS=99999;
  !LastNumFlips=99999
  DO WHILE(NUMFLIPS.GT.0) 
    CALL edgeFlipper(DIM,NP,POINTS,NF,TRIAS,T2N,T2T,NUMFLIPS) 
    IF(LastNumFlips.EQ.NumFlips) THEN 
      print *, "outta here"
      EXIT ! stuck lets get out of here! 
    ENDIF
    LastNumFlips=NumFlips 
  ENDDO
  
  ! Output of the current mesh
  IF(MOD(ITER,20).EQ.0.OR.ITER.EQ.1) THEN
    CALL WriteMesh(DIM,POINTS,NP,TRIAS,NF,ITER)
  ENDIF
  
  ! Describe each bar by a unique pair of nodes
  CALL findUniqueBars(DIM,NF,TRIAS,NUMBARS,BARS)

  ! Calculate forces on bars
  CALL CalcForces(MeshSize,DIM,POINTS,NP,BARS,NUMBARS,FVEC)

  ! Move points based on forces
  CALL ApplyForces(DIM,POINTS,NP,BARS,NUMBARS,FVEC) 

  ! Bring outside points back to the boundary
  CALL ProjectPointsBack(DIM,PSLG,POINTS,NP)
  
  CALL CPU_TIME(TF) 
  WRITE(*,'(A,F12.8)') "INFO: ELAPSED TIME IS: ",TF-TS 
   
  ! Termination criterion: reached max iterations. 
  IF(ITER.EQ.MaxIter) THEN
    WRITE(*,'(A)') "INFO: MAXIMUM NUMBER OF ITERATIONS REACHED..." 
    WRITE(*,'(A)') "********************************************************"
    EXIT 
  ENDIF
  
  ITER = ITER + 1 
  WRITE(*,'(A,I4,A)') "INFO: ITERATION: ",ITER," COMPLETE"

ENDDO

WRITE(*,'(A,I8,A,I8,A)')"INFO: " & 
//" FINAL MESH HAS ",NP," VERTICES " &
//" AND ",NF," FACES "

CALL WriteMesh(DIM,POINTS,NP,TRIAS,NF,ITER)

!-----------------------------------------------------------------------
END PROGRAM 
!-----------------------------------------------------------------------

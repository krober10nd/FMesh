!-----------------------------------------------------------------------
! PROGRAM DISTMESH 
!-----------------------------------------------------------------------
!> @brief Serial DistMesh algorithm to build a mesh in 2d or 3d.
!-----------------------------------------------------------------------
PROGRAM DISTMESH 
!-----------------------------------------------------------------------

USE msizes 
USE utils 
USE vars 

IMPLICIT NONE

SzFields=ParseInputs()                                       ! Parse Mesh.inp and load in all sizing fields 

CALL FormInitialPoints2D(SzFields,PSLG,POINTS,NP)            ! Create initial points to iterate on

CALL DelTriaWElim(PSLG,NP,POINTS,NF,TRIAS)                   ! Compute Delaunay triangulation of point set with masking

CALL TriaToTria(NF,TRIAS,T2T,T2N)                            ! Calculate the triangle adj. matrices

WRITE(*,'(A)') "                                      "
WRITE(*,'(A)') "****************BEGIN ITERATING*************************"
WRITE(*,'(A)') "                                      "

ITER    = 1                                                  ! Intialize iteration counter 

CALL CPU_TIME(TS) 
DO 
  CALL CPU_TIME(TSiter) 

  ! Perform edge flips to achieve Del.-hood
  NumFlips=99999; LastNumFlips=99999 ; ITFlips =1
  DO WHILE(NUMFLIPS.GT.0) 
    CALL edgeFlipper(SzFields,NP,POINTS,NF,TRIAS,T2N,T2T,NUMFLIPS) 
    IF(LastNumFlips.EQ.NumFlips) THEN 
      EXIT ! stuck..
    ENDIF
    ITFlips = ITFlips + 1
    IF(ITFlips.EQ.20) THEN
      EXIT ! stuck..
    ENDIF
    LastNumFlips=NumFlips 
    ! TODO check for any overlapped elements, will need to lower timestep?
  ENDDO
  
  ! Output of the current mesh
  IF(MOD(ITER,NSCREEN).EQ.0.OR.ITER.EQ.1) THEN
    CALL WriteMesh(POINTS,NP,TRIAS,NF,ITER)
  ENDIF
  
  ! Describe each bar by a unique pair of nodes
  CALL FindUniqueBars(NF,TRIAS,NUMBARS,BARS)

  ! Calculate forces on bars
  CALL CalcForces(SzFields,POINTS,BARS,NUMBARS,FVEC)

  ! Move points based on forces
  CALL ApplyForces(POINTS,NP,BARS,NUMBARS,FVEC) 

  ! Bring outside points back to the boundary
  CALL ProjectPointsBack(PSLG,POINTS,NP)
  
  CALL CPU_TIME(TFiter) 
  IF(MOD(ITER,NSCREEN).EQ.0.OR.ITER.EQ.1) THEN
    WRITE(*,'(A,F12.8)') "INFO: ELAPSED TIME IS: ",TFiter-TSiter 
  ENDIF
   
  ! Termination criterion: reached max iterations. 
  IF(ITER.EQ.MaxIter) THEN
    WRITE(*,'(A)') "INFO: MAXIMUM NUMBER OF ITERATIONS REACHED..." 
    WRITE(*,'(A)') "********************************************************"
    EXIT 
  ENDIF
  
  ITER = ITER + 1 
  WRITE(*,'(A,I4,A)') "INFO: ITERATION: ",ITER," COMPLETE"

ENDDO

CALL CPU_TIME(TF) 
WRITE(*,'(A,F12.8)') "INFO: ELAPSED TIME FOR MESH GENERATION IS: ",TF-TS 

WRITE(*,'(A,I8,A,I8,A)')"INFO: " & 
//" FINAL MESH HAS ",NP," VERTICES " &
//" AND ",NF," FACES "


!-----------------------------------------------------------------------
END PROGRAM DISTMESH
!-----------------------------------------------------------------------

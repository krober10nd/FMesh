!-----------------------------------------------------------------------
! PROGRAM DISTMESH 
!-----------------------------------------------------------------------
!> @brief Serial DistMesh algorithm to build a mesh in 2d or 3d.
!-----------------------------------------------------------------------
PROGRAM DISTMESH 
!-----------------------------------------------------------------------

USE YourMeshSize 
USE utils 
USE vars 

IMPLICIT NONE

!--------------------------------------------------
! THESE SHOULD BECOME INPUT FILE PARAMETERS
!--------------------------------------------------
MaxIter = 300 ! Maximum number of iterations
ITER    = 1 ! intialize iteration counter 
DELTAT  = 0.05d0 ! psuedo-timestep
NSCREEN = 5 ! number of times to write data to disk
!--------------------------------------------------

SzFields=ParseInputs()                                       ! Load in all sizing fields 

CALL FormInitialPoints2D(CalcMeshSize,SzFields,PSLG,POINTS,NP)! Create initial points to iterate on

stop

CALL DelTriaWElim(PSLG,NP,POINTS,NF,TRIAS)                   ! Compute Delaunay triangulation of point set with masking

CALL TriaToTria(NF,TRIAS,T2T,T2N)                            ! Calculate the triangle adj. matrices

WRITE(*,'(A)') "                                      "
WRITE(*,'(A)') "****************BEGIN ITERATING*************************"
WRITE(*,'(A)') "                                      "


DO 
  CALL CPU_TIME(TS) 

  ! Perform edge flips to achieve Del.-hood
  NumFlips=99999;
  LastNumFlips=99999
  DO WHILE(NUMFLIPS.GT.0) 
    CALL edgeFlipper(SzFields,NP,POINTS,NF,TRIAS,T2N,T2T,NUMFLIPS) 
    IF(LastNumFlips.EQ.NumFlips) THEN 
      EXIT ! stuck lets get out of here! 
    ENDIF
    LastNumFlips=NumFlips 
    ! TODO check for any overlapped elements, will need to lower timestep?
    ! TODO check if stuck repeating the same number of iterations 
  ENDDO
  
  ! Output of the current mesh
  IF(MOD(ITER,NSCREEN).EQ.0.OR.ITER.EQ.1) THEN
    CALL WriteMesh(POINTS,NP,TRIAS,NF,ITER)
  ENDIF
  
  ! Describe each bar by a unique pair of nodes
  CALL FindUniqueBars(NF,TRIAS,NUMBARS,BARS)

  ! Calculate forces on bars
  CALL CalcForces(SzFields,POINTS,NP,BARS,NUMBARS,FVEC)

  ! Move points based on forces
  CALL ApplyForces(POINTS,NP,BARS,NUMBARS,FVEC) 

  ! Bring outside points back to the boundary
  CALL ProjectPointsBack(PSLG,POINTS,NP)
  
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

CALL WriteMesh(POINTS,NP,TRIAS,NF,ITER)

!-----------------------------------------------------------------------
END PROGRAM DISTMESH
!-----------------------------------------------------------------------

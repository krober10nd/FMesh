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

! ARGUMENT PARSING GOES HERE !
! TODO MAKE THIS INTO A SUBROUTINE 
IF(COMMAND_ARGUMENT_COUNT().LT.1)THEN
  WRITE(*,*)'ERROR, COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
  STOP
ENDIF
CALL GET_COMMAND_ARGUMENT(1,pslgfname)   
CALL GET_COMMAND_ARGUMENT(2,sizefname)   
CALL GET_COMMAND_ARGUMENT(3,elongfname)                    
SzFx = LoadMeshSizes(sizefname)                              ! Load size function into memory
PSLG = ReadPSLGtxt(pslgfname,SzFx)                           ! Load in boundary description 
ElongFx = LoadElongation(elongfname)                         ! Load in elongation into memory
! ARGUMENT PARSING ENDS HERE !


CALL FormInitialPoints2D(MeshSize,SzFx,ElongFx,DIM,PSLG,POINTS,NP)   ! Create initial points to iterate on

IF(NP.LT.3) THEN 
  WRITE(*,'(A)') "FATAL: NOT ENOUGH INITIAL POINTS TO MESH" 
  STOP
ENDIF

CALL DelTriaWElim(DIM,PSLG,NP,POINTS,NF,TRIAS)               ! Compute Delaunay triangulation of point set with masking

CALL TriaToTria(NF,TRIAS,T2T,T2N)                            ! Calculate the triangle adj. matrices

WRITE(*,'(A)') "                                      "
WRITE(*,'(A)') "****************BEGIN ITERATING*************************"
WRITE(*,'(A)') "                                      "

MaxIter = 300 ! Maximum number of iterations
ITER    = 1 ! intialize iteration counter 
DELTAT  = 0.05d0 ! psuedo-timestep
NSCREEN = 5 ! number of times to write data to disk

DO 
  CALL CPU_TIME(TS) 

  ! Perform edge flips to achieve Del.-hood
  NUMFLIPS=99999;
  LastNumFlips=99999
  DO WHILE(NUMFLIPS.GT.0) 
    CALL edgeFlipper(DIM,ElongFx,NP,POINTS,NF,TRIAS,T2N,T2T,NUMFLIPS) 
    IF(LastNumFlips.EQ.NumFlips) THEN 
      EXIT ! stuck lets get out of here! 
    ENDIF
    LastNumFlips=NumFlips 
    ! check for any overlapped elements, will need to lower timestep?
    ! check if stuck repeating the same number of iterations 
  ENDDO
  
  ! Output of the current mesh
  IF(MOD(ITER,NSCREEN).EQ.0.OR.ITER.EQ.1) THEN
    CALL WriteMesh(DIM,POINTS,NP,TRIAS,NF,ITER)
  ENDIF
  
  ! Describe each bar by a unique pair of nodes
  CALL findUniqueBars(DIM,NF,TRIAS,NUMBARS,BARS)

  ! Calculate forces on bars
  CALL CalcForces(MeshSize,SzFx,ElongFx,DIM,POINTS,NP,BARS,NUMBARS,FVEC)

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
END PROGRAM DISTMESH
!-----------------------------------------------------------------------

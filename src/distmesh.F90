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
USE YourMeshSize, ONLY : MeshSize
USE utils 

IMPLICIT NONE

INTEGER I 
REAL(8) :: TS,TF,T1,T2

DIM = 2 ! 2 DIMENIONS 
LMIN=0.05 ! MINIMUM ELEMENT SIZE 
MaxIter = 50 ! MAXIMUM NUMBER OF ITERATIONS

CALL ReadPSLGtxt(PSLG,LMIN)                                 ! Read in boundary description 

CALL FormInitialPoints2D(MeshSize,DIM,PSLG,LMIN,POINTS,NP)  ! Step 1-2: Create initial points to iterate on

WRITE(*,'(A)') "                                      "
WRITE(*,'(A)') "****************BEGIN ITERATING*************************"
WRITE(*,'(A)') "                                      "

ITER = 0 ! iteration counter 

ALLOCATE(PointsOLD(NP,DIM)) 
PointsOld = -9999.0d0                                       ! For the first iteration 
DO ! distmesh loop

  CALL CPU_TIME(TS) 

  CALL CPU_TIME(T1)
  CALL DelTriaWElim(DIM,PSLG,NP,POINTS,NF,TRIAS,IERR)      
  CALL CPU_TIME(T2)

  print *, "DELTRAI ", T2-T1

  PointsOld(1:NP,:)=Points(1:NP,1:DIM)

  ! 4. Describe each bar by a unique pair of nodes
  CALL CPU_TIME(T1) 
  CALL findUniqueBars(DIM,NF,TRIAS,NUMBARS,BARS)
  CALL CPU_TIME(T2) 
  print *, "FINDBARS ",T2-T1

  ! 5. Output of the current mesh
  IF(MOD(ITER,10).EQ.0) THEN
    CALL WriteMesh(DIM,POINTS,NP,TRIAS,NF,ITER)
  ENDIF

  ! 6. Calculate forces on bars
  CALL CPU_TIME(T1)
  CALL CalcForces(MeshSize,DIM,POINTS,NP,BARS,NUMBARS,FVEC)
  CALL CPU_TIME(T2) 
  print *, "CalcForces ",T2-T1

  ! 7. Move points based on forces
  CALL CPU_TIME(T1)
  CALL ApplyForces(DIM,POINTS,NP,BARS,NUMBARS,FVEC) 
  CALL CPU_TIME(T2)
  print *, "ApplyForces ",T2-T1

  ! 8. Bring outside points back to the boundary
  CALL CPU_TIME(T1)
  CALL ProjectPointsBack(DIM,PSLG,POINTS,NP)
  CALL CPU_TIME(T2)
  print *, "ProjectPointsBack ",T2-T1
  
  WRITE(*,'(A,I4,A)') "INFO: ITERATION: ",ITER," COMPLETE!"
  ITER = ITER + 1 

  CALL CPU_TIME(TF) 
  print *, (TF-TS)
   
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

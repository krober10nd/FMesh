MODULE YourMeshSize
USE vars 
implicit none

! SIZING FUNCTION RELATED PARAMETERS GO HERE 
REAL(kind=real_t)   :: LMIN=0.001d0 ! minimum mesh size 
REAL(kind=real_t)   :: LMAX=0.55d0  ! maximum mesh size 
REAL(kind=real_t)   :: GRADE=0.35d0 ! mesh size gradation 
INTEGER(kind=idx_t) :: DIM=2        ! dimension of problem

contains

! DEFINE YOU MESH SIZE FUNCTION HERE
FUNCTION MeshSize(POINTS) 
real(kind=real_t):: MeshSize
real(kind=real_t),intent(in) :: points(2)
MeshSize = minval((/LMAX,(1.0d0-points(1))*GRADE + LMIN/))
END FUNCTION MeshSize

! YOUR ORIENTATION AND ELONGATION CALCULATIONS HERE
FUNCTION CalcMetricTensor(POINTS) RESULT(ME)
real(kind=real_t) :: ME(2,2) 
real(kind=real_t),intent(in) :: POINTS(2)
! just isotropic for now 
ME = 0.0d0 
ME(1,1) = 1.0d0 
ME(2,2) = 1.0d0 
END FUNCTION 

END MODULE YourMeshSize

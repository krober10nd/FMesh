MODULE YourMeshSize
USE vars 
implicit none

! SIZING FUNCTION RELATED PARAMETERS GO HERE 
REAL(kind=real_t)   :: LMIN=0.05d0 ! minimum mesh size 
REAL(kind=real_t)   :: LMAX=0.50d0  ! maximum mesh size 
REAL(kind=real_t)   :: GRADE=0.35d0 ! mesh size gradation 
INTEGER(kind=idx_t) :: DIM=2        ! dimension of problem

contains

! DEFINE YOU MESH SIZE FUNCTION HERE
FUNCTION MeshSize(POINTS) 
real(kind=real_t):: MeshSize
real(kind=real_t),intent(in) :: points(2)
MeshSize = minval((/LMAX,(5.0d0-points(1))*GRADE + LMIN/))
END FUNCTION MeshSize


! FX TO LOAD IN GRIDDED INTERP OF MESH SIZE FUNCTIONS 

! FX TO ACCEPT GRIDDED INTERP OF MESH SIZE FUNCTION AND A POINT AND RETURN A SIZE 

! YOUR ORIENTATION AND ELONGATION CALCULATIONS HERE
! IF ISOTROPIC SET TO IDENTITY MATRIX 
FUNCTION CalcMetricTensor(POINTS) RESULT(ME)
real(kind=real_t) :: ME(1:2,1:2) 
real(kind=real_t),intent(in) :: POINTS(2)
! just isotropic for now 
ME = 0.0d0 
ME(1,1) = 1.0d0
ME(2,2) = 1.0d0
! create a simple smooth variation in anisotropicness 
if(POINTS(1).GT.-2.0d0.AND.POINTS(1).LT.0.0d0) THEN 
  ME(2,2)= (2.0d0+points(1))*2.0d0+1.0d0
elseif(POINTS(1).GT.0.0d0.AND.POINTS(1).LT.2.0d0) THEN
  ME(2,2)= (2.0d0-points(1))*2.0d0+1.0d0
elseif(POINTS(1).LT.EPSILON(1.0d0).AND.POINTS(1).GT.-EPSILON(1.0d0)) THEN
  ME(2,2)=4.0d0
endif
ME(2,2)=MAXVAL((/ 1.0d0,ME(2,2)/))

!ME(2,2) = 10.d0
END FUNCTION 

END MODULE YourMeshSize

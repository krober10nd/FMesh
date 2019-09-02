MODULE YourMeshSize
USE UTILS 
implicit none

REAL(kind=real_t)   :: LMIN=0.001d0 ! minimum mesh size 
INTEGER(kind=idx_t) :: DIM=2 ! dimension of problem

contains

! DEFINE YOU MESH SIZE FUNCTION HERE
FUNCTION MeshSize(POINTS) 
real(kind=real_t):: MeshSize
real(kind=real_t),intent(in) :: points(2)


MeshSize = minval((/0.005d0,(1.0d0-points(1))*0.55d0 + LMIN/))
END FUNCTION MeshSize

end module YourMeshSize

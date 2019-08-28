MODULE YourMeshSize
USE UTILS 
implicit none

contains

! DEFINE YOU MESH SIZE FUNCTION HERE
FUNCTION MeshSize(POINTS) 
real(kind=real_t):: MeshSize
real(kind=real_t),intent(in) :: points(2)
MeshSize = (1.0d0-points(1))*0.15d0 + 0.05d0
END FUNCTION MeshSize

end module YourMeshSize

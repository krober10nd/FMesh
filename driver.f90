program testqhull
!
! dummy program to test calling qhull
!
use iso_c_binding, only : c_ptr, c_f_pointer
use utils
implicit none

integer(kind=idx_t) :: d=2   ! dimension of problem
integer(kind=idx_t) :: np=1000 ! number of points
integer(kind=idx_t) :: nf    ! number of faces
integer(kind=idx_t) :: ierr 

integer(kind=idx_t),allocatable :: trias(:,:) ! facet table [dim+1 x nf] array of point indices 
real(kind=real_t),allocatable   :: points(:,:) !  [dim x np] array of points

integer :: i,j

allocate(points(d,np))
points=0.d0
do i =1,np
  do j =1,d
    points(j,i)=RAND()
  enddo
enddo

!! this call creates the face table given the options 
ierr = triangulate(d,np,points,nf,trias)

print *, "NUM OF FACES IS ",nf 
! VERIFY OUTPUT
!!do i =1,nf
!!  do j =1,d
!!    print *, "FACE ", i, " HAS ",trias(j,i)
!!  enddo
!!enddo

end program testqhull

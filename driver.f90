program testqhull
use utils
implicit none

integer :: d=2
integer :: np=10
integer :: i

real(8),allocatable :: points(:,:)

integer :: ierr 

allocate(points(np,2))
points=0.d0
do i =1,np
  points(i,1)=RAND()
  points(i,2)=RAND()
enddo

ierr = triangulate(d,np,points)

end program testqhull

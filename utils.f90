module utils
use iso_c_binding, only: c_int32_t, c_int64_t, c_float, c_double, c_ptr
! calls qhull library from c
! PUT KIND PRECISION CODE HERE
implicit none
interface
function ctriangulate(DIM, NUMPOINTS, fpoints)bind(c,name='ctriangulate')
use iso_c_binding
implicit none
integer(c_int) :: ctriangulate
integer(c_int), intent(in) :: DIM
integer(c_int), intent(in) :: NUMPOINTS 
real(c_double), intent(in) :: fpoints(NUMPOINTS,DIM)
end function ctriangulate
end interface

contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function triangulate(DIM, NUMPOINTS, points)
use iso_c_binding
implicit none
integer(c_int), intent(in) :: DIM
integer(c_int), intent(in) :: NUMPOINTS 
real(c_double), intent(in) :: points(NUMPOINTS,DIM)
integer(c_int) :: ierr


print *, points(1,1)
print *, points(1,2) 
print *, "IN UTILS" 
ierr = ctriangulate(DIM,NUMPOINTS,points)

end function triangulate

end module utils

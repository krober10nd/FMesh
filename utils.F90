module utils
use iso_c_binding, only: c_int, c_int32_t, c_int64_t, c_float, c_double, c_ptr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calls qhull library from fortran language to produce facet table of delaunay triangulation 
!
! INPUTS: 
!
! OUTPUTS: 
!
! AUTHOR: Keith Jared Roberts, August/19/2019
!         USP, Brasil 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

!
! Width of elementary data types (should match those used in metis.h)
!
#ifdef INT64
    integer, parameter, public :: idx_t = c_int64_t ! <--- modify integer size here (c_int32_t or c_int64_t)
#else
    integer, parameter, public :: idx_t = c_int32_t 
#endif

#ifdef REAL64
    integer, parameter, public :: real_t = c_double  ! <--- modify real size here (c_float or c_double)
#else
    integer, parameter, public :: real_t = c_float
#endif

interface

function ctriangulate(DIM, NUMPOINTS, fpoints, ffacets)bind(c,name='ctriangulate')
import idx_t,real_t,c_ptr
implicit none
integer(kind=idx_t)             :: ctriangulate
integer(kind=idx_t), value, intent(in) :: DIM
integer(kind=idx_t), value, intent(in) :: NUMPOINTS 
real(kind=real_t), intent(in)   :: fpoints(DIM,NUMPOINTS)
type(c_ptr), intent(out) :: ffacets
end function ctriangulate

end interface

contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function triangulate(DIM, NUMPOINTS, points, facets)
implicit none
integer(kind=idx_t), intent(in) :: DIM
integer(kind=idx_t), intent(in) :: NUMPOINTS 
real(kind=real_t),   intent(in) :: points(DIM,NUMPOINTS)
type(c_ptr),intent(out) :: facets(*)
integer(kind=idx_t)             :: ierr

ierr = ctriangulate(DIM,NUMPOINTS,points,facets)

end function triangulate

end module utils

module utils
use iso_c_binding, only: c_int, c_int32_t, c_int64_t, c_float, c_double, c_ptr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROCEDURES FOR FortranMesh
! AUTHOR: Keith Jared Roberts, V1.0 August/19/2019
!         Universidade de São Paulo, Brasil 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

!
! Width of elementary data types change directive in makefile please 
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

function faces(DIM, NUMPOINTS, fpoints, NUMFACETS)bind(c,name='faces')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calls qhull library from fortran language to produce facet table of delaunay triangulation 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
import idx_t,real_t,c_ptr
implicit none
type(c_ptr) :: faces
integer(kind=idx_t), value, intent(in) :: DIM
integer(kind=idx_t), value, intent(in) :: NUMPOINTS 
integer(kind=idx_t), intent(inout) :: NUMFACETS
real(kind=real_t), intent(in)   :: fpoints(DIM,NUMPOINTS)
end function 

SUBROUTINE destroy_storage(p) BIND(C, NAME='destroy_storage')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! releases memory that was allocated in c 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
  IMPLICIT NONE
  TYPE(C_PTR), INTENT(IN), VALUE :: p
END SUBROUTINE destroy_storage


end interface

contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function triangulate(DIM, NUMPOINTS, points, NUMFACETS,facets)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calls qhull library from fortran language to produce facet table of delaunay triangulation 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use iso_c_binding, only : c_f_pointer,C_PTR
implicit none
integer(kind=idx_t), intent(in) :: DIM
integer(kind=idx_t), intent(in) :: NUMPOINTS 
integer(kind=idx_t), intent(inout):: NUMFACETS
real(kind=real_t),   intent(in) :: points(DIM,NUMPOINTS)
integer(kind=idx_t),intent(out),allocatable :: facets(:,:)

type(c_ptr) :: cfacetemp
integer(kind=idx_t),pointer :: ffacetemp(:)
integer(kind=idx_t):: ierr
integer :: i,j,k

! call qhull library 
cfacetemp = faces(DIM,NUMPOINTS,points,NUMFACETS)

CALL C_F_POINTER(cfacetemp, ffacetemp, [NUMFACETS*(DIM+1)])

! reshape vector ffacetemp to facets array
! TODO: check if it's possible to allocate to memory
allocate(facets(DIM+1,NUMFACETS))

facets=0 
k=0 
do i = 1,NUMFACETS
  do j = 1,DIM+1 
    k = k + 1
    facets(j,i)=ffacetemp(k)+1 ! ensure indexing starts at 1
  enddo
enddo

! IMPORTANT: RELEASE MEMORY MALLOC'ed in C 
CALL destroy_storage(cfacetemp)

! TODO: BETTER ERROR CHECKING (SEE ZOLTAN ERROR LEVELS FOR REFERENCE?)
ierr = 0 !! SUCCESS 

end function triangulate

end module utils

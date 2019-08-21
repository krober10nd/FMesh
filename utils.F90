module utils
use iso_c_binding, only: c_int, c_int32_t, c_int64_t, c_float, c_double, c_ptr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROCEDURES FOR FortranMesh
! AUTHOR: Keith Jared Roberts, V1.0 August/19/2019
!         Universidade de São Paulo, Brasil 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

! Width of elementary data types change directive in makefile please 
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

! Later on this should become a derived type 
integer(kind=idx_t),allocatable :: trias(:,:)  ! facet table [dim+1 x nf] array of point indices 
real(kind=real_t),allocatable   :: points(:,:) !  [dim x np] array of points
integer(kind=idx_t) :: NP    ! number of vertices/nodes/points in the mesh
integer(kind=idx_t) :: NF    ! number of facets in mesh 

INTEGER(idx_t) :: DIM
REAL(real_t)   :: LMIN
REAL(real_t)   :: BBOX(4) 

! parameters for distmesh
REAL(kind=real_t),PARAMETER :: DPTOL =0.001d0
REAL(kind=real_t),PARAMETER :: TTOL  =0.1d0
REAL(kind=real_t),PARAMETER :: FSCALE=1.2d0 
REAL(kind=real_t),PARAMETER :: DELTAT=0.2d0 
REAL(kind=real_t) :: GEPS 
REAL(kind=real_t) :: DEPS 
REAL(kind=real_t) :: EPS 

! error conditions 
integer(kind=idx_t) :: ierr 

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
integer(kind=idx_t), intent(inout)     :: NUMFACETS
real(kind=real_t), intent(in)          :: fpoints(DIM,NUMPOINTS)

end function 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!


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
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

end interface

contains 


subroutine formInitialPoints2D(DIM, BBOX, LMIN,IPTS,NP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form initial points to populate in the domain according to mesh size function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUTS arguments that will be passed later on when this becomes a subroutine 
INTEGER(idx_t),INTENT(IN) :: DIM
REAL(real_t),INTENT(IN) :: LMIN       
REAL(real_t),INTENT(IN) :: BBOX(*)   

! OUTPUTS 
REAL(real_t),INTENT(OUT),ALLOCATABLE :: IPTS(:,:)
INTEGER(idx_t),INTENT(OUT) :: NP 

! local to subroutine 
REAL(real_t),ALLOCATABLE :: xvec(:),yvec(:)
REAL(real_t),ALLOCATABLE :: xg(:,:)
REAL(real_t),ALLOCATABLE :: yg(:,:)
REAL(real_t) :: temp,junk
INTEGER(idx_t)           :: nx,ny,nz
INTEGER  :: i,j,k

EPS=EPSILON(1.0d0) ! machine precision on your computer 
GEPS=0.001d0*LMIN 
DEPS=SQRT(EPS)*LMIN 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1. Create initial distribution in bounding box (equilateral triangles)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! [x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
! x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
! p=[x(:),y(:)];                                       % List of node coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

temp=LMIN*(sqrt(3.0d0)/2.0d0)

nx=0
junk = bbox(1) ! x minimum 
DO
  if(junk.GT.BBOX(3)) exit
   nx = nx + 1 
   junk = junk + LMIN ! next grid point    
ENDDO 
print *, NX 
ny=0
junk = bbox(3) ! y minimum 
DO
  if(junk.GT.BBOX(4)) continue
   ny = ny + 1 
   junk = junk + temp ! next grid point 
ENDDO 


ALLOCATE(XG(NY,NX),YG(NY,NX))
XG = -99999.d0
YG = -99999.d0

ALLOCATE(XVEC(NX),YVEC(NY))

ALLOCATE(IPTS(2,NX*NY))
IPTS = -99999.d0

DO I = 1,NX
  XVEC(I) = LMIN*(I-1) + BBOX(1) ! plus x minimum 
ENDDO

DO J = 1,NY
  YVEC(J) = temp*(J-1) + BBOX(3)! plus y minimum 
ENDDO
print *, XVEC(NX)
print *, YVEC(NY) 

! now repmat xvec and yvec
DO I = 1,NY
  DO J = 1,NX
    XG(I,J) = XVEC(J)
    YG(I,J) = YVEC(I)
  ENDDO
ENDDO

!! SHIFT X-COORD TO FORM EQUILATERAL TRIAS
DO I = 2,NY,2
  DO J = 1,NX
    XG(I,J) = XG(I,J) + LMIN/2.d0
  ENDDO
ENDDO

! SAVE ALL INITIAL POINTS 
NP = 0 
DO I = 1,NX
  DO J = 1,NY 
    NP = NP + 1
    IPTS(1,NP) = XG(i,j)
    IPTS(2,NP) = YG(i,j) 
  ENDDO
ENDDO

! 2. Remove points outside the region, apply the rejection method
! THIS STEP REQUIRES 
! 1. A LINEAR INTERPOLANT TO DETERMINE THE PROBABILITY OF RETAINING THE POINT 
! 2. A FUNCTION TO GENERATE A "RANDOM" NUMBER
! 3. A WAY TO DETERMINE IF A POINT IS IN A POLYGON 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! p=p(feval(fd,p,varargin{:})<geps,:);                 % Keep only d<0 points.
! r0=1./feval(fh,p,varargin{:}).^2;                    % Probability to keep point.
! p=p(rand(size(p,1),1)<r0./max(r0),:);                % Rejection method.
! N=size(p,1);                                         % Number of initial points.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine FormInitialPoints2D
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!


SUBROUTINE Triangulate(DIM, NP, POINTS, NF, FACETS,IERR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calls qhull library from fortran language to output facet table of delaunay triangulation 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use iso_c_binding, only : c_f_pointer,C_PTR
implicit none

! INPUTS 
integer(kind=idx_t), intent(in) :: DIM
integer(kind=idx_t), intent(in) :: NP
real(kind=real_t),   intent(in) :: POINTS(DIM,NP)
integer(kind=idx_t), intent(inout):: NF

! OUTPUTS 
integer(kind=idx_t),intent(out),allocatable :: FACETS(:,:)
integer(kind=idx_t),INTENT(OUT) :: IERR

! LOCAL TO SUBROUTINE 
type(c_ptr) :: cfacetemp
integer(kind=idx_t),pointer :: ffacetemp(:)
integer :: i,j,k

! call qhull library 
cfacetemp = faces(DIM,NP,POINTS,NF)

CALL C_F_POINTER(cfacetemp, ffacetemp, [NF*(DIM+1)])

! reshape vector ffacetemp to facets array
! TODO: check if it's possible to allocate to memory with an error condition
allocate(facets(DIM+1,NF))

facets=0 
k=0 
do i = 1,NF
  do j = 1,DIM+1 
    k = k + 1
    facets(j,i)=ffacetemp(k)+1 ! ensure indexing starts at 1
  enddo
enddo

! IMPORTANT: RELEASE MEMORY MALLOC'ed in C 
CALL destroy_storage(cfacetemp)

! TODO: BETTER ERROR CHECKING (SEE ZOLTAN ERROR LEVELS FOR REFERENCE?)
ierr = 0 !! SUCCESS 

end subroutine Triangulate
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

end module utils

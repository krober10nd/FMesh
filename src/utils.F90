module utils
use iso_c_binding, only: c_int, c_int32_t, c_int64_t, c_float, c_double, c_ptr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROCEDURES FOR FMesh
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


! 2D domain boundary description 
TYPE BounDescrip
   INTEGER(kind=idx_t) :: NumVert
   INTEGER(kind=idx_t) :: DIM
   REAL(kind=real_t),ALLOCATABLE :: Vert(:,:)
ENDTYPE

TYPE(BounDescrip) :: PSLG 

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

subroutine ReadPSLGtxt(PSLG) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read in PSLG from a text file named into derived type BounDescrip
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUTS 

! OUTPUTS
TYPE(BounDescrip),INTENT(OUT) :: PSLG 

! LOCAL TO SUBROUTINE 
logical :: fileFound 
integer :: tempNP,tempDIM
integer :: i 

INQUIRE(FILE='PSLG.txt',EXIST=fileFound)

IF(fileFound) THEN
  OPEN(UNIT=1,FILE='PSLG.txt',ACTION='READ')
  READ(1,*) tempNP,tempDIM ! number of points in file  
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(A,I5,A,I5,A)') "INFO: Reading PSLG file called PSLG.txt with " &  
      //" ",tempNP," points in ",tempDIM," dimensions."
  WRITE(*,'(A)') "********************************************************"

  ALLOCATE(PSLG%Vert(tempDIM,tempNP))
  PSLG%NumVert = tempNP 
  PSLG%dim = tempDIM 
  DO I = 1,PSLG%NumVert
    READ(1,*) PSLG%Vert(1,I),PSLG%Vert(2,I)
  ENDDO
  CLOSE(1)
  WRITE(*,'(A)') "INFO: Read PSLG file succesfully!"
ELSE
  WRITE(*,'(A)') "FATAL: PSLG.txt file not found"
  STOP  
ENDIF
end subroutine ReadPSLGtxt 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

subroutine formInitialPoints2D(DIM, PSLG, LMIN,IPTS,NP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form initial points to populate in the domain according to mesh size function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUTS arguments that will be passed later on when this becomes a subroutine 
INTEGER(idx_t),INTENT(IN) :: DIM
TYPE(BounDescrip),INTENT(IN) :: PSLG 
REAL(real_t),INTENT(IN) :: LMIN       

! OUTPUTS 
REAL(real_t),INTENT(OUT),ALLOCATABLE :: IPTS(:,:)
INTEGER(idx_t),INTENT(OUT) :: NP 

! local to subroutine 
REAL(real_t),ALLOCATABLE :: xvec(:),yvec(:)
REAL(real_t),ALLOCATABLE :: xg(:,:)
REAL(real_t),ALLOCATABLE :: yg(:,:)
REAL(real_t) :: temp,junk
INTEGER(idx_t)           :: nx,ny,nz
INTEGER(idx_t)  :: i,j,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1. Create initial distribution in bounding box (equilateral triangles)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! [x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
! x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
! p=[x(:),y(:)];                                       % List of node coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

BBOX(1) = MINVAL(PSLG%Vert(1,:))
BBOX(2) = MAXVAL(PSLG%Vert(1,:))
BBOX(3) = MINVAL(PSLG%Vert(2,:))
BBOX(4) = MAXVAL(PSLG%Vert(2,:))

temp=LMIN*(sqrt(3.0d0)/2.0d0)
NX=FLOOR((BBOX(2) - BBOX(1))/LMIN)+1
NY=FLOOR((BBOX(4) - BBOX(3))/temp)+1

WRITE(*,'(A)') "********************************************************"
WRITE(*,'(A,F6.3,A,F6.3,A,F6.3,A,F6.3)') "INFO: Domain extents are: " &  
       //" XMIN: ",BBOX(1)," XMAX: ",BBOX(2), " " & 
       //" YMIN: ",BBOX(3)," YMAX: ",BBOX(4) 
WRITE(*,'(A)') "********************************************************"

ALLOCATE(XVEC(NX),YVEC(NY))

ALLOCATE(IPTS(2,NX*NY))
IPTS = -99999.d0

DO I = 1,NX
  XVEC(I) = LMIN*(I-1) + BBOX(1) ! plus x minimum 
ENDDO

DO J = 1,NY
  YVEC(J) = temp*(J-1) + BBOX(3)! plus y minimum 
ENDDO

ALLOCATE(XG(NY,NX))
ALLOCATE(YG(NY,NX))

CALL MeshGrid2D(XVEC,YVEC,XG,YG) 

!! SHIFT X-COORD TO FORM EQUILATERAL TRIAS
DO I = 2,NY,2
  DO J = 1,NX
    XG(I,J) = XG(I,J) + LMIN/2.d0
  ENDDO
ENDDO

! SAVE ALL INITIAL POINTS 
NP = 0 
DO I = 1,NY
  DO J = 1,NX 
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
CALL CalculateSignedDistance(PSLG) 
! Keep points with dist < geps 

! Calculate R0 by evaluating sizing function 

! Apply rejection method 


end subroutine FormInitialPoints2D
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE CalculateSignedDistance(PSLG)
implicit none 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Given a Planar straight line graph and a set of points, determine the nearest 
! Euclidean distance to the boundary using a KD-tree and sign the distance 
! negative if the point in question is inside the polygon defined by the PSLG. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE(BounDescrip),INTENT(IN) :: PSLG 

END SUBROUTINE CalculateSignedDistance
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

subroutine meshgrid2D(xgv, ygv, X, Y)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  real(kind=real_t),intent(in)   :: xgv(:), ygv(:)
  real(kind=real_t),intent(out)  :: X(:,:), Y(:,:)
  integer           :: sX, sY

  X(:,:) = spread( xgv, 1, size(ygv) )
  Y(:,:) = spread( ygv, 2, size(xgv) )

end subroutine
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

end module utils

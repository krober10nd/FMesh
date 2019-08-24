MODULE UTILS 
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

#ifdef REAL32
   integer, parameter, public :: real_t = c_float  ! <--- modify real size here (c_float or c_double)
#else
    integer, parameter, public :: real_t = c_double
#endif

! Later on this should become a derived type 
integer(kind=idx_t),allocatable :: trias(:,:)  ! facet table [dim+1 x nf] array of point indices 
real(kind=real_t),allocatable   :: points(:,:) !  [dim x np] array of points
integer(kind=idx_t) :: NP    ! number of vertices/nodes/points in the mesh
integer(kind=idx_t) :: NF    ! number of facets in mesh 

INTEGER(kind=idx_t) :: DIM
REAL(kind=real_t)   :: LMIN
REAL(kind=real_t)   :: BBOX(4) 

! parameters for distmesh
REAL(kind=real_t),PARAMETER :: DPTOL =0.001d0
REAL(kind=real_t),PARAMETER :: TTOL  =0.1d0
REAL(kind=real_t),PARAMETER :: FSCALE=1.2d0 
REAL(kind=real_t),PARAMETER :: DELTAT=0.2d0 
REAL(kind=real_t) :: GEPS 
REAL(kind=real_t) :: DEPS 
REAL(kind=real_t),PARAMETER :: EPS=EPSILON(1.d0) 

! 2D domain boundary description 
TYPE BounDescrip2D
   INTEGER(kind=idx_t) :: NumVert
   INTEGER(kind=idx_t) :: DIM
   REAL(kind=real_t),ALLOCATABLE :: Vert(:,:)
ENDTYPE

TYPE(BounDescrip2D) :: PSLG 

! error conditions 
integer(kind=idx_t) :: ierr 

INTERFACE

FUNCTION pnpoly(NUMVERT, VERTx, VERTy, TESTx, TESTy)bind(c,name='pnpoly')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALLS PNPOLY TO DETERMINE IF POINT IS IN A 2D multiply-connected POLYGON
! C CODE IS IN CFUNCTIONS.c
! https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html#The%20C%20Code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
import idx_t,real_t,c_ptr

implicit none

integer(kind=idx_t) :: pnpoly 

integer(kind=idx_t), value, intent(in) :: NUMVERT
real(kind=real_t),   value, intent(in) :: TESTx
real(kind=real_t),   value, intent(in) :: TESTy
real(kind=real_t),intent(in)           :: VERTx(NUMVERT)
real(kind=real_t),intent(in)           :: VERTy(NUMVERT)


END FUNCTION pnpoly
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

FUNCTION faces(DIM, NUMPOINTS, fpoints, NUMFACETS)bind(c,name='faces')
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

END FUNCTION faces
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


END INTERFACE



CONTAINS 



SUBROUTINE ReadPSLGtxt(PSLG,LMIN) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read in Planar Straight Line Graph (PSLG) from a text file named 
! 'PSLG.txt' into derived type BounDescrip2D. Interpolate points on boundary 
! so spacing between two points does nto exceed LMIN/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUTS 
REAL(real_t),INTENT(IN) :: LMIN       

! OUTPUTS
TYPE(BounDescrip2D),INTENT(OUT) :: PSLG 

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

CALL densify(LMIN/2.0d0,PSLG) 

END SUBROUTINE ReadPSLGtxt 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE densify(MAXDIFF,PSLG) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Edit vertex spacing on PSLG to ensure it is less than MAXDIFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUTS 
REAL(kind=real_t),INTENT(IN) :: MAXDIFF

! OUTPUTS
TYPE(BounDescrip2D), INTENT(INOUT) :: PSLG 

! LOCAL
REAL(kind=real_t),ALLOCATABLE :: temp(:,:)
REAL(kind=real_t) :: DX(PSLG%NumVert),DY(PSLG%NumVert)
REAL(kind=real_t),ALLOCATABLE :: ilat(:),ilon(:)

INTEGER(kind=idx_t) :: nx,ny,nout
INTEGER(kind=idx_t) :: NIN(PSLG%NumVert-1)
INTEGER(kind=idx_t) :: SUMIN
INTEGER(kind=idx_t) :: i,j,nstep,n,ni

nx = size(PSLG%Vert(1,:))
ny = size(PSLG%Vert(2,:))

DO i = 2,PSLG%NumVert
  DY(i)=ABS(PSLG%Vert(2,i)-PSLG%Vert(2,i-1));
  DX(i)=ABS(PSLG%Vert(1,i)-PSLG%Vert(1,i-1));
  NIN(i-1)=CEILING(MAXVAL((/DX(i),DY(i)/))/maxdiff)-1;
ENDDO

SUMIN=SUM(NIN)
IF(SUMIN.eq.0) then 
    WRITE(*,'(A)') "INFO: No insertion needed on PSLG."
    RETURN
ENDIF

NOUT=SUMIN+NX;

ALLOCATE(temp(2,NOUT)) 
temp = -99999.d0

n=1
do i=1,nx-1
    ni=nin(i)
    if(ni.EQ.0) THEN
        temp(2,n)=PSLG%Vert(2,i)
        temp(1,n)=PSLG%Vert(1,i)
        nstep=1;
    ELSE
        ilat=linspace(PSLG%Vert(2,i),PSLG%Vert(2,i+1),ni+2)
        ilon=linspace(PSLG%Vert(1,i),PSLG%Vert(1,i+1),ni+2)
        temp(2,n:n+ni)=ilat(1:ni+1)
        temp(1,n:n+ni)=ilon(1:ni+1)
        
        nstep=ni+1;
    endif
    n = n + nstep 
enddo
temp(2,nout)=PSLG%Vert(2,ny)
temp(1,nout)=PSLG%Vert(1,nx)

! debug 
OPEN(UNIT=303,FILE="debug.txt",ACTION='WRITE')
do i = 1,nout 
  WRITE(303,"(2F12.8)")temp(1,i),temp(2,i)
ENDDO
CLOSE(303)
DEALLOCATE(PSLG%Vert)
ALLOCATE(PSLG%Vert(PSLG%DIM,NOUT))
PSLG%Vert = temp 
PSLG%NumVert = NOUT
DEALLOCATE(temp)

end subroutine densify 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

subroutine formInitialTria2D(DIM,PSLG,LMIN,IPTS,NP,TRIAS,NF)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form initial Delaunay triangulation to fill in the domain according to mesh size function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUTS 
INTEGER(idx_t),INTENT(IN)      :: DIM
TYPE(BounDescrip2D),INTENT(IN) :: PSLG 
REAL(real_t),INTENT(IN)        :: LMIN       

! OUTPUTS 
REAL(real_t),INTENT(OUT),ALLOCATABLE :: IPTS(:,:)
INTEGER(kind=idx_t),INTENT(OUT),ALLOCATABLE:: TRIAS(:,:)
INTEGER(idx_t),INTENT(OUT) :: NP 
INTEGER(kind=idx_t), INTENT(INOUT):: NF

! local to subroutine 
REAL(real_t),ALLOCATABLE :: xvec(:),yvec(:)
REAL(real_t),ALLOCATABLE :: xg(:,:)
REAL(real_t),ALLOCATABLE :: yg(:,:)
REAL(real_t),ALLOCATABLE :: Dist(:)
REAL(real_t)             :: temp,junk
INTEGER(idx_t)           :: nx,ny,nz
INTEGER(idx_t)           :: i,j,k

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
DO J = 1,NX
  DO I = 2,NY,2
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
CALL CalculateSignedDistance(IPTS,PSLG,Dist) 

! Keep points with dist < geps 

! Calculate R0 by evaluating sizing function 

! Apply rejection method 


CALL DelTriangulate(DIM,NP,POINTS,NF,TRIAS,IERR)


! DEBUG VISUALIZE INITIAL TRIANGULATION 
OPEN(UNIT=300,FILE="Points.txt",ACTION='WRITE')
DO i=1,NP
  WRITE(300,"(2F12.8)")POINTS(1,i),POINTS(2,i)
ENDDO
CLOSE(300)

OPEN(UNIT=301,FILE="Facets.txt",ACTION='WRITE')
DO i=1,NF
  WRITE(301,"(3I8)")TRIAS(1,i),TRIAS(2,i),TRIAS(3,i)
ENDDO
CLOSE(301)

END SUBROUTINE FormInitialTria2D
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE CalculateSignedDistance(IPTS,PSLG,SignedDistance)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Given a PSLG and a set of points, determine the nearest 
! Euclidean distance to the PSLG using a KD-tree and sign the distance 
! negative if the point in question is inside the polygon defined by the PSLG. 
! and vice-versa otherwise.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use kdtree2_module
implicit none 

! INPUTS 
TYPE(BounDescrip2D),INTENT(IN) :: PSLG 
REAL(real_t),INTENT(IN),ALLOCATABLE :: IPTS(:,:)

! OUTPUTS 
REAL(real_t),INTENT(OUT),ALLOCATABLE :: SignedDistance(:)

! LOCAL 
TYPE(kdtree2), POINTER :: TREE=>null()
TYPE(kdtree2_result), ALLOCATABLE :: KDRESULTS(:)
INTEGER(kind=idx_t) :: tempSZ
INTEGER(kind=idx_t) :: INoOUT 
INTEGER(kind=idx_t) :: I

! BUILD KD-TREE with PSLG vertices
tree => kdtree2_create(PSLG%Vert,rearrange=.true.,sort=.true.)

! Loop over all the initial points 
tempSZ = SIZE(IPTS,2) 
ALLOCATE(SignedDistance(tempSZ))

ALLOCATE(KDRESULTS(1))

INoOUT = -999
DO I =1,tempSZ
  call kdtree2_n_nearest(tp=tree,qv=IPTS(:,I),nn=1, & 
                               results=KDRESULTS)

  SignedDistance(I) = SQRT(KDRESULTS(1)%DIS)
  !print *, I
  !print *, PSLG%Vert(1,I), PSLG%Vert(2,I) 
  !print *, IPTS(1,I),IPTS(2,I)
  !print *, "****************************"

  ! Determine if point is in the PSLG defined polygon
  CALL FPNPOLY(PSLG,IPTS(:,I),INoOUT)
  
  IF(INoOUT.EQ.1) THEN 
    SignedDistance(I)= -SignedDistance(I)
  ENDIF

ENDDO

call kdtree2_destroy(tp=tree)

DEALLOCATE(KDRESULTS)

print *, tempSZ
! DEBUG VISUALIZE INITIAL TRIANGULATION 
OPEN(UNIT=305,FILE="SignedDistance.txt",ACTION='WRITE')
DO i=1,tempSZ
  WRITE(305,"(F12.8)") SignedDistance(I)
ENDDO
CLOSE(305)

END SUBROUTINE CalculateSignedDistance
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE FPNPOLY(PSLG,TEST,INoOUT) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fortran subroutine to call inpoly PNPOLY written in cfunctions.c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use iso_c_binding, only : c_f_pointer,C_PTR
implicit none 

! INPUTS
TYPE(BounDescrip2D), INTENT(IN) :: PSLG 
REAL(kind=real_t), INTENT(IN) :: TEST(2)

! OUTPUTS 
INTEGER(kind=idx_t),INTENT(OUT) :: INoOUT 

!FUNCTION pnpoly(NUMVERT, VERTx, VERTy, TESTx, TESTy)bind(c,name='pnpoly')
INoOUT = pnpoly(PSLG%NumVert, PSLG%Vert(1,:), PSLG%Vert(2,:), &
                TEST(1),TEST(2))

END SUBROUTINE 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!


SUBROUTINE DelTriangulate(DIM, NP, POINTS, NF, FACETS,IERR)
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

END SUBROUTINE DelTriangulate
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE meshgrid2D(xgv, ygv, X, Y)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create a struture grid of points using the vectors xgv and ygv
! this mimics what the MATLAB command does 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
  ! INPUTS 
  real(kind=real_t),intent(in)   :: xgv(:), ygv(:)
  ! OUTPUTS
  real(kind=real_t),intent(out)  :: X(:,:), Y(:,:)

  X(:,:) = spread( xgv, 1, size(ygv) )
  Y(:,:) = spread( ygv, 2, size(xgv) )

END SUBROUTINE
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

FUNCTION linspace(a, b, n) result(s)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creates a vector s that spans the interval (a,b) with n points. 
! this mimics what the MATLAB command does 
! REFERENCE: https://github.com/certik/fortran-utils
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
! INPUTS 
real(kind=real_t), intent(in) :: a, b
integer(kind=idx_t), intent(in) :: n
! OUTPUTS 
real(kind=real_t) :: s(n)

s = meshexp(a, b, 1.0d0, n-1)
END FUNCTION
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

FUNCTION meshexp(rmin, rmax, a, N) result(mesh)
! REFERENCE: https://github.com/certik/fortran-utils
! Generates exponential mesh of N elements on [rmin, rmax]
!
! Arguments
! ---------
!
! The domain [rmin, rmax], the mesh will contain both enkind=real_toints:
real(kind=real_t), intent(in) :: rmin, rmax
!
! The ratio of the rightmost to leftmost element lengths in the mesh (for a > 1
! this means the "largest/smallest"); The only requirement is a > 0. For a == 1
! a uniform mesh will be returned:
real(kind=real_t), intent(in) :: a
!
! The number of elements in the mesh:
integer(kind=idx_t), intent(in) :: N
!
! Returns
! -------
!
! The generated mesh:
real(kind=real_t) :: mesh(N+1)
!
! Note: Every exponential mesh is fully determined by the set of parameters
! (rmin, rmax, a, N). Use the get_meshexp_pars() subroutine to obtain them
! from the given mesh.
!
! Example
! -------
!
! real(kind=real_t) :: r(11)
! r = meshexp(0._kind=real_t, 50._kind=real_t, 1e9_kind=real_t, 10)

integer :: i
real(kind=real_t) :: alpha, beta
if (a < 0.0d0) then
    WRITE(*,'(A)')"meshexp: a > 0 required"
else if (abs(a - 1) < tiny(1.0d0)) then
    alpha = (rmax - rmin) / N
    do i = 1, N+1
        mesh(i) = alpha * (i-1.0d0) + rmin
    end do
else
    if (N > 1) then
        beta = log(a) / (N-1)
        alpha = (rmax - rmin) / (exp(beta*N) - 1)
        do i = 1, N+1
            mesh(i) = alpha * (exp(beta*(i-1)) - 1) + rmin
        end do
    else if (N == 1) then
        mesh(1) = rmin
        mesh(2) = rmax
    else
        WRITE(*,'(A)') "meshexp: N >= 1 required"
    end if
end if
END FUNCTION
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

END MODULE UTILS

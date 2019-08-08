!**************************************************
! serial DistMesh algorithm in modern FORTRAN
! INPUTS: 
! OUTPUTS: 
! AUTHOR: kjr, universdade de sao paulo, 2019--
!**************************************************
PROGRAM DistMesh

use iso_c_binding, only: c_int32_t, c_int64_t, c_float, c_double, c_ptr

IMPLICIT NONE

! Width of elementary data types
#ifdef INT64
    integer, parameter :: idx_t = c_int64_t ! <--- modify integer size here (c_int32_t or c_int64_t)
#else
    integer, parameter :: idx_t = c_int32_t ! <--- modify integer size here (c_int32_t or c_int64_t)
#endif

#ifdef REAL64
    integer, parameter :: real_t = c_double  ! <--- modify real size here (c_float or c_double)
#else
    integer, parameter :: real_t = c_float
#endif

! parameters for distmesh to function
REAL(real_t),PARAMETER :: DPTOL =0.001d0
REAL(real_t),PARAMETER :: TTOL  =0.1d0
REAL(real_t),PARAMETER :: FSCALE=1.2d0 
REAL(real_t),PARAMETER :: DELTAT=0.2d0 
REAL(real_t) :: GEPS 
REAL(real_t) :: DEPS 
REAL(real_t) :: EPS 

! to form initial points
REAL(real_t),ALLOCATABLE :: XG(:,:)
REAL(real_t),ALLOCATABLE :: YG(:,:)
REAL(real_t),ALLOCATABLE :: IPTS(:,:)
INTEGER(idx_t)           :: NX,NY 

! arguments that will be passed later on when this becomes a subroutine 
REAL(real_t) :: LMIN = 0.50d0 ! minimum element size 
REAL(real_t) :: BBOX(4)    ! bounding box

INTEGER  :: I,J,K ! all counters go here

EPS=EPSILON(1.0d0) ! machine precision on your computer 

! this will be passed later on the development process 
BBOX = (/-1.0d0,1.0d0,-1.0d0,1.0d0/)  ! TL TR BL BR 

GEPS=0.001d0*LMIN 
DEPS=SQRT(EPS)*LMIN 

! STEP ONE (THIS WILL BECOME A PRE_PROCESSING STEP LATER ON)
! 1. Create initial distribution in bounding box (equilateral triangles)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! [x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
! x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
! p=[x(:),y(:)];                                       % List of node coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

NX = CEILING((BBOX(2) - BBOX(1))/LMIN) 
NY = CEILING((BBOX(4) - BBOX(3))/LMIN) 

ALLOCATE(XG(NX,NY),YG(NX,NY))
XG = -99999.d0
YG = -99999.d0

ALLOCATE(IPTS(NX*NY,2))
IPTS = -99999.d0

DO I = 1,NX
  DO J = 1,NY
    XG(I,J) = LMIN*(I-1)
  ENDDO
ENDDO

DO J = 1,NX
  DO I = 1,NY
    YG(I,J) = LMIN*(J-1) 
  ENDDO
ENDDO

! SHIFT X-COORD TO FORM EQUILATERAL TRIAS
DO I = 2,NX,2
  DO J = 1,NY
    XG(I,J) = XG(I,J) + LMIN/2.d0
  ENDDO
ENDDO

! SAVE ALL INITIAL POINTS 
DO I = 1,NX
  DO J = 1,NY 
    IPTS(I,1) = XG(I,J)
    IPTS(I,2) = YG(I,J) 
  ENDDO
ENDDO

! DEBUG 
!OPEN(UNIT=300,FILE="XG.txt",ACTION='WRITE')
!OPEN(UNIT=301,FILE="YG.txt",ACTION='WRITE')
!
!DO I=1,NX 
!  DO J=1,NY 
!    WRITE(300,"(F12.8)",ADVANCE='NO') XG(I,J)
!    WRITE(301,"(F12.8)",ADVANCE='NO') YG(I,J)
!  ENDDO
!  WRITE(300,"(/)")
!  WRITE(301,"(/)")
!ENDDO
!
!CLOSE(300)
!CLOSE(301)


! STEP TWO (THIS WILL BECOME A PRE_PROCESSING STEP LATER ON)
! THIS STEP REQUIRES 
! 1. A LINEAR INTERPOLANT TO DETERMINE THE PROBABILITY OF RETAINING THE POINT 
! 2. A FUNCTION TO GENERATE A "RANDOM" NUMBER
! 3. A WAY TO DETERMINE IF A POINT IS IN A POLYGON 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2. Remove points outside the region, apply the rejection method
! p=p(feval(fd,p,varargin{:})<geps,:);                 % Keep only d<0 points.
! r0=1./feval(fh,p,varargin{:}).^2;                    % Probability to keep point.
! p=p(rand(size(p,1),1)<r0./max(r0),:);                % Rejection method.
! N=size(p,1);                                         % Number of initial points.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END PROGRAM 

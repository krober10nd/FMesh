MODULE UTILS 
use iso_c_binding, only: c_int, c_int32_t, c_int64_t, c_float, c_double, c_ptr

!**************************************************
! PROCEDURES FOR FMesh
! 
! AUTHOR: Keith J. Roberts, 
! PLACE: Universdade de Sao Paulo
! START: 2019-08-17 
!**************************************************

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
integer(kind=idx_t),allocatable :: trias(:,:)  ! facet table [nf x dim+1] array of point indices 
integer(kind=idx_t),allocatable :: t2t(:,:)  ! triangle-to-triangle neighbors [nf x dim+1] of connected trias
integer(kind=idx_t),allocatable :: t2n(:,:)  ! triangle-to-vertex neighbors [nf x dim+1] of connected trias (see TriaToTria)
integer(kind=idx_t),allocatable :: bars(:,:)   ! unique bars of the triangulation
real(kind=real_t),allocatable   :: points(:,:) !  [dim x np] array of points
real(kind=real_t),allocatable   :: pointsOld(:,:) ! [dim x np] array of points from previous iteration
real(kind=real_t),allocatable   :: fvec(:,:) ! [numbars x 2] array of spring forces on bars
integer(kind=idx_t) :: NP    ! number of vertices/nodes/points in the mesh
integer(kind=idx_t) :: NF    ! number of facets in mesh 
integer(kind=idx_t) :: NumBars ! number of edges in the mesh 

REAL(kind=real_t)   :: BBOX(4) ! bounding box coords. Top left and bot. right. 

! some parameters for distmesh
REAL(kind=real_t),PARAMETER :: DPTOL =0.001d0
REAL(kind=real_t),PARAMETER :: TTOL  =0.1d0
REAL(kind=real_t),PARAMETER :: FSCALE=1.2d0 
REAL(kind=real_t),PARAMETER :: DELTAT=0.1d0
REAL(kind=real_t) :: GEPS 
REAL(kind=real_t) :: DEPS 
REAL(kind=real_t),PARAMETER :: EPS=EPSILON(1.d0) 
INTEGER(kind=idx_t) :: iter,MaxIter

! 2D domain boundary description 
TYPE BounDescrip2D
   INTEGER(kind=idx_t) :: NumVert
   INTEGER(kind=idx_t) :: NumVert0
   INTEGER(kind=idx_t) :: DIM
   REAL(kind=real_t),ALLOCATABLE :: Vert(:,:)  ! densified poly
   REAL(kind=real_t),ALLOCATABLE :: Vert0(:,:) ! org. poly
ENDTYPE

TYPE(BounDescrip2D) :: PSLG 

! error conditions 
integer(kind=idx_t) :: ierr 

! AN ISOTROPIC MESH SIZE FUNCTION 
! this is defined by the user in YourMeshSize.F90 module
ABSTRACT INTERFACE
     function IsoSZ(P)
        import real_t
        real(kind=real_t) :: sz
        real(kind=real_t), intent (in) :: p(2)
     end function IsoSZ
END INTERFACE


INTERFACE

FUNCTION pnpoly(NUMVERT, VERTx, VERTy, TESTx, TESTy)bind(c,name='pnpoly')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALLS PNPOLY TO DETERMINE IF POINT IS IN A 2D multiply-connected POLYGON
! C CODE IS IN CFUNCTIONS.c
! https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html#The%20C%20Code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
import idx_t,real_t

implicit none

integer(kind=idx_t) :: pnpoly 

integer(kind=idx_t), value, intent(in) :: NUMVERT
real(kind=real_t),   value, intent(in) :: TESTx
real(kind=real_t),   value, intent(in) :: TESTy
real(kind=real_t),intent(in)           :: VERTx(NUMVERT)
real(kind=real_t),intent(in)           :: VERTy(NUMVERT)

END FUNCTION pnpoly
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

END INTERFACE

INTERFACE 

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

SUBROUTINE ProjectPointsBack(DIM,PSLG,POINTS,NP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update vertex locations according to calc'ed forces. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE 

! INPUTS 
INTEGER(idx_t),INTENT(IN) :: DIM , NP 
TYPE(BounDescrip2D),INTENT(IN) :: PSLG 

! OUTPUTS
REAL(real_t),INTENT(INOUT),ALLOCATABLE :: POINTS(:,:)

! LOCAL 
INTEGER(idx_t) :: I
INTEGER(idx_t) :: NOUT
REAL(real_t),ALLOCATABLE :: Dist(:)
REAL(real_t),ALLOCATABLE :: dgradx(:),dgrady(:),dgrad2(:)
REAL(real_t),ALLOCATABLE :: tempDistx(:),tempDisty(:)
REAL(real_t),ALLOCATABLE :: tempP1(:,:),tempP2(:,:)
REAL(real_t),ALLOCATABLE :: PointsOut(:,:),DistOut(:)
INTEGER(idx_t),ALLOCATABLE :: INoOUT(:) 

DEPS = SQRT(EPSILON(1.d0)) 

!  %7. Bring outside points back to the boundary
!  d = feval(obj.fd,p,obj,[],1); ix = d > 0;                  % Find points outside (d>0)
ALLOCATE(INoOUT(NP)) 
INoOUT = 0
DO I = 1,NP 
  CALL FPNPOLY(PSLG,POINTS(I,:),INoOUT(I))
ENDDO

NOUT=0
ALLOCATE(PointsOut(NP,DIM),DistOut(NP))
PointsOut=-99999.d0 
DistOut=-9999.d0
DO I=1,NP 
  IF(INoOUT(I).EQ.0) THEN 
    NOUT=NOUT+1
    PointsOut(NOUT,:) = Points(I,:)
  ENDIF
ENDDO

!  if sum(ix) > 0
IF(NOUT.GT.0) THEN  

  CALL CalcSignedDistance(PointsOut,PSLG,DistOut) 

  ALLOCATE(dgradx(NOUT),dgrady(NOUT),dgrad2(NOUT))
  ALLOCATE(tempDistx(NOUT),tempDisty(NOUT))
  ALLOCATE(tempP1(NOUT,DIM),tempP2(NOUT,DIM))

!      dgradx = (feval(obj.fd,[p(ix,1)+deps,p(ix,2)],obj,[])...%,1)...
!                -d(ix))/deps; % Numerical
  tempP1 = PointsOUT 
  tempP1(:,1) = tempP1(:,1) + DEPS 
  CALL CalcSignedDistance(tempP1,PSLG,tempDistx) 
  dgradx= (tempDistx - DistOut(1:NOUT))/DEPS 

!      dgrady = (feval(obj.fd,[p(ix,1),p(ix,2)+deps],obj,[])...%,1)...
!                -d(ix))/deps; % gradient
  tempP2 = PointsOUT 
  tempP2(:,2) = tempP1(:,2) + DEPS 
  CALL CalcSignedDistance(tempP2,PSLG,tempDisty) 
  dgrady= (tempDisty - DistOut(1:NOUT))/DEPS 

!      dgrad2 = dgradx.^+2 + dgrady.^+2;
!      p(ix,:) = p(ix,:)-[d(ix).*dgradx./dgrad2,...
!                         d(ix).*dgrady./dgrad2];
  dgrad2 = dgradx**2 + dgrady**2 
  NOUT=0
  DO I =1,NP 
    IF(INoOUT(I).EQ.0) THEN
      NOUT=NOUT+1
      POINTS(I,1) = POINTS(I,1) - DistOut(NOUT)*dgradx(NOUT)/dgrad2(NOUT)
      POINTS(I,2) = POINTS(I,2) - DistOut(NOUT)*dgrady(NOUT)/dgrad2(NOUT) 
    ENDIF
  ENDDO

ENDIF

! RELEASE ALLOCATABLES 
DEALLOCATE(DistOut,PointsOut)
DEALLOCATE(dgradx,dgrady,dgrad2)
DEALLOCATE(tempDistx,tempDisty) 
DEALLOCATE(tempP1,tempP2)

END SUBROUTINE 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE ApplyForces(DIM,POINTS,NP,BARS,NUMBARS,FVEC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update vertex locations according to calc'ed forces. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE 

! INPUTS 
INTEGER(idx_t),INTENT(IN) :: DIM,NP,NUMBARS 
INTEGER(idx_t),INTENT(IN),ALLOCATABLE :: BARS(:,:)
REAL(real_t),INTENT(IN),ALLOCATABLE :: FVEC(:,:)

! OUTPUTS 
REAL(real_t),INTENT(INOUT),ALLOCATABLE :: POINTS(:,:) 

! LOCAL 
INTEGER(idx_t) :: I 

DO I = 1,NUMBARS
  POINTS(BARS(I,1),1) = POINTS(BARS(I,1),1) + DELTAT * FVEC(I,1)
  POINTS(BARS(I,1),2) = POINTS(BARS(I,1),2) + DELTAT * FVEC(I,2)

  POINTS(BARS(I,2),1) = POINTS(BARS(I,2),1) - DELTAT * FVEC(I,1)
  POINTS(BARS(I,2),2) = POINTS(BARS(I,2),2) - DELTAT * FVEC(I,2)
ENDDO

END SUBROUTINE 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE CalcForces(HFX,DIM,POINTS,NP,BARS,NUMBARS,FVEC) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate forcing terms for each nodes based on mesh size function 
! and actual length of the edge 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUTS
REAL(real_t) :: HFX  ! Mesh Size function 
INTEGER(idx_t),INTENT(IN) :: DIM,NP,NUMBARS
INTEGER(idx_t),INTENT(IN),ALLOCATABLE :: BARS(:,:)
REAL(real_t),INTENT(IN),ALLOCATABLE :: POINTS(:,:)

! OUTPUTS
REAL(real_t),INTENT(OUT),ALLOCATABLE :: FVEC(:,:)

! LOCAL 
INTEGER(idx_t) :: i 
REAL(real_t)   :: barvec(NUMBARS,DIM)
REAL(real_t)   :: L(NUMBARS),L0(NUMBARS),LN(NUMBARS)
REAL(real_t)   :: FORCES(NUMBARS)
REAL(real_t)   :: HBARS(NUMBARS)
REAL(real_t)   :: midpt(2),a,b,SCALE_FACTOR

!barvec=p(bars(:,1),:)-p(bars(:,2),:); % List of bar vectors
!L=sqrt(sum(barvec.^2,2)); % L = Bar lengths
!hbars=feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2,varargin{:});
!L0=hbars*Fscale*median(L)/median(HBARS); ! From Roberts et al. 2019
!LN = L./L0;  
!F    = (1-LN.^4).*exp(-LN.^4)./LN; % Bossens-Heckbert edge force
!Fvec = F*[1,1].*barvec;

ALLOCATE(FVEC(NUMBARS,2))

FORCES=0.0d0
DO I = 1,NUMBARS
  BARVEC(I,:)=POINTS(bars(I,1),:) - POINTS(bars(I,2),:)
  L(I)=SQRT(BARVEC(I,1)**2.0d0 + BARVEC(I,2)**2.0d0)
  midpt=(POINTS(BARS(I,1),:) + POINTS(BARS(I,2),:))/2.0d0
  HBARS(I)=HFX(midpt)
ENDDO 

a=MEDIAN(L,1,NUMBARS) ; b = MEDIAN(HBARS,1,NUMBARS)
SCALE_FACTOR = a/b 

DO I = 1,NUMBARS
  L0(I)=HBARS(I)*FSCALE*SCALE_FACTOR
  LN(I)=L(I)/L0(I)
  ! Bossens Heckbert 
  FORCES(I)=(1-LN(I)**4)*EXP(-LN(I)**4)/LN(I)
  ! Linear spring (Hooke's Law)
  !FORCES(I)=MAXVAL( (/1.0d0-LN(I),0.0d0/))
  FVEC(I,1)=FORCES(I)*BARVEC(I,1)
  FVEC(I,2)=FORCES(I)*BARVEC(I,2)
ENDDO

END SUBROUTINE CalcForces
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!


SUBROUTINE sortRows(NumRows,ROWS,IDX) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Given an array of int indices (i.e., a ROW), sort the ROW in ascending 
! lexiographic  order and return the mapping to achieve the sorted ROWS array.  
! assumes the ROWS are the first two columns of the array ROWS (i.e., ROWS(:,1:2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUT 
INTEGER(idx_t),INTENT(IN),ALLOCATABLE :: ROWS(:,:)
INTEGER(idx_t),INTENT(IN) :: NumRows

! OUTPUT 
INTEGER(idx_t),INTENT(OUT),ALLOCATABLE :: IDX(:)

! LOCAL 
INTEGER(idx_t) :: i,k,temp
INTEGER(idx_t) :: junkShort(2) 
INTEGER(8)     :: junkLong 
INTEGER(8),ALLOCATABLE :: tmpB1(:),tmpB2(:)
INTEGER(idx_t),ALLOCATABLE :: tROWS(:,:)

ALLOCATE(tROWS(SIZE(ROWS,1),SIZE(ROWS,2)))
tROWS = ROWS

ALLOCATE(tmpB1(NumROWS)) 

! 2. Convert to unique integer number
DO I = 1,NumRows
  tmpB1(I) = TRANSFER((/tROWS(I,1),tROWS(I,2)/),junkLONG)  
ENDDO

! 3. Apply sorting method to sort in ascending order 
! need the array idx which maps the sort 
allocate(idx(numrows))
! original indexing 
do i =1,numRows
  idx(i) = i
enddo
CALL quickSortINTLONG(tmpB1,1,NumROWS,idx) 

END SUBROUTINE sortRows
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!


SUBROUTINE findUniqueBars(DIM,NF,FACES,NumBars,BARS) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Find the unique edges of the triangulation 
! bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
! bars=unique(sort(bars,2),'rows');      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUT 
INTEGER(idx_t),INTENT(IN) :: DIM,NF 
INTEGER(idx_T),INTENT(IN),ALLOCATABLE :: FACES(:,:)

! OUTPUT 
INTEGER(idx_t),INTENT(OUT),ALLOCATABLE :: BARS(:,:)
INTEGER(idx_t),INTENT(OUT) :: NumBars

INTEGER(idx_t) :: i,k,temp
INTEGER(idx_t) :: junkShort(2) 
INTEGER(8)     :: junkLong 
INTEGER(8),ALLOCATABLE :: tmpB1(:),tmpB2(:)
INTEGER(idx_t),ALLOCATABLE :: idx(:)

! 1. Determine non unique edges of mesh 
ALLOCATE(BARS(3*NF,2))
BARS = 0
K = 0 
DO I = 1,NF 
  K = K + 1
  BARS(K,1) = FACES(I,1) 
  BARS(K,2) = FACES(I,2)
ENDDO
DO I = 1,NF 
  K = K + 1
  BARS(K,1) = FACES(I,1) 
  BARS(K,2) = FACES(I,3)
ENDDO
DO I = 1,NF 
  K = K + 1
  BARS(K,1) = FACES(I,2) 
  BARS(K,2) = FACES(I,3)
ENDDO

NumBars = K 
! 2. Ensure first column has a smaller value than 2nd column 
DO I = 1,NumBars
  IF(BARS(I,1).GT.BARS(I,2)) THEN 
    temp = BARS(I,1) 
    BARS(I,1) = BARS(I,2) 
    BARS(I,2) = temp 
  ENDIF  
ENDDO

ALLOCATE(tmpB1(NumBars)) 

! 3. Convert to unique integer number
DO I = 1,NumBars
  tmpB1(I) = TRANSFER((/BARS(I,1),BARS(I,2)/),junkLONG)  
ENDDO

! 4. Apply sorting method to sort in ascending order 
allocate(idx(numbars))
CALL quickSortINTLONG(tmpB1,1,NumBars,idx) 
deallocate(idx)

ALLOCATE(tmpB2(NumBars))
tmpB2=0
! 5. Only keep unique entries 
tmpB2(1) = tmpB1(1) 
K=1 
DO I = 2,NumBars-1
  IF(tmpB1(I).NE.tmpB2(K)) THEN 
    K = K + 1 
    tmpB2(K) = tmpB1(I) 
  ENDIF
ENDDO
NumBars = K 

! 6. Convert back to the bar format
BARS = 0 
DO I = 1,NumBars 
  BARS(I,:) = TRANSFER(tmpB2(I),junkShort)
ENDDO 

!OPEN(UNIT=300,FILE="BARS.txt",ACTION='WRITE')
!DO i=1,NumBars
!  WRITE(300,"(2I8)")BARS(I,1),BARS(I,2)
!ENDDO
!CLOSE(300)

END SUBROUTINE findUniqueBars 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE WriteMesh(DIM,POINTS,NP,FACETS,NF,ITER)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write the mesh to disk 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUT 
REAL(real_t),INTENT(IN),ALLOCATABLE :: POINTS(:,:) 
INTEGER(idx_t),INTENT(IN),ALLOCATABLE :: FACETS(:,:) 
INTEGER(idx_t),INTENT(IN) :: NF,NP,DIM,ITER

! OUTPUT 
! WRITES SEVERAL TEXT FILES WITH POINTS AND FACETS 

! LOCAL 
INTEGER(idx_t) :: i 
CHARACTER(14) :: filename

IF(DIM.EQ.2) THEN 
  WRITE(filename,'(a,i4.4,a)') "POINTS",ITER,".TXT"
  ! 2D VISUALIZATION GOES HERE 
  OPEN(UNIT=300,FILE=filename,ACTION='WRITE')
  DO i=1,NP
    WRITE(300,"(2F16.8)")POINTS(i,1),POINTS(i,2)
  ENDDO
  CLOSE(300)
  
  WRITE(filename,'(a,i4.4,a)') "FACETS",ITER,".TXT"
  OPEN(UNIT=301,FILE=filename,ACTION='WRITE')
  DO i=1,NF
    WRITE(301,"(3I8)")FACETS(i,1),FACETS(i,2),FACETS(i,3)
  ENDDO
  CLOSE(301)
ELSE 

    ! 3D VISUALIZATION GOES HERE 
ENDIF

END SUBROUTINE WriteMesh
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

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

  ALLOCATE(PSLG%Vert(tempNP,tempDIM))
  PSLG%NumVert = tempNP 
  PSLG%dim = tempDIM 
  DO I = 1,PSLG%NumVert
    READ(1,*) PSLG%Vert(I,1),PSLG%Vert(I,2)
  ENDDO
  CLOSE(1)
  WRITE(*,'(A)') "INFO: Read PSLG file succesfully!"
  WRITE(*,'(A)') "********************************************************"
ELSE
  WRITE(*,'(A)') "FATAL: PSLG.txt file not found"
  WRITE(*,'(A)') "********************************************************"
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

nx = size(PSLG%Vert(:,1))
ny = size(PSLG%Vert(:,2))

DO i = 2,PSLG%NumVert
  DY(i)=ABS(PSLG%Vert(i,2)-PSLG%Vert(i-1,2));
  DX(i)=ABS(PSLG%Vert(i,1)-PSLG%Vert(i-1,1));
  NIN(i-1)=CEILING(MAXVAL((/DX(i),DY(i)/))/maxdiff)-1;
ENDDO

SUMIN=SUM(NIN)
IF(SUMIN.eq.0) then 
    WRITE(*,'(A)') "                                      "
    WRITE(*,'(A)') "INFO: No insertion needed on PSLG."
    WRITE(*,'(A)') "                                      "
    RETURN
ENDIF

NOUT=SUMIN+NX;

ALLOCATE(temp(NOUT,2)) 
temp = -99999.d0

n=1
do i=1,nx-1
    ni=nin(i)
    if(ni.EQ.0) THEN
        temp(n,2)=PSLG%Vert(i,2)
        temp(n,1)=PSLG%Vert(i,1)
        nstep=1;
    ELSE
        ilat=linspace(PSLG%Vert(i,2),PSLG%Vert(i+1,2),ni+2)
        ilon=linspace(PSLG%Vert(i,1),PSLG%Vert(i+1,1),ni+2)
        temp(n:n+ni,2)=ilat(1:ni+1)
        temp(n:n+ni,1)=ilon(1:ni+1)
        
        nstep=ni+1;
    endif
    n = n + nstep 
enddo
temp(nout,2)=PSLG%Vert(ny,2)
temp(nout,1)=PSLG%Vert(nx,1)

!! NEED TO TEST MULTIPLY CONNECTED POLYGONS
!!! debug 
!OPEN(UNIT=303,FILE="DensifiedPSLG.txt",ACTION='WRITE')
!do i = 1,nout 
!  WRITE(303,"(2F12.8)")temp(i,1),temp(i,2)
!ENDDO
!CLOSE(303)
ALLOCATE(PSLG%Vert0(PSLG%NumVert,PSLG%DIM))
PSLG%Vert0 = PSLG%Vert
PSLG%NumVert0 = PSLG%NumVert 
DEALLOCATE(PSLG%Vert)
ALLOCATE(PSLG%Vert(NOUT,PSLG%DIM))
PSLG%Vert = temp 
PSLG%NumVert = NOUT
DEALLOCATE(temp)

end subroutine densify 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

subroutine FormInitialPoints2D(HFX,DIM,PSLG,LMIN,IPTS,NP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form initial points to fill in the domain according to mesh size function
! This corresponds with steps 1 to 2 of the distmesh algorithm. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUTS 
REAL(real_t)                   :: HFX  ! Mesh Size function 
INTEGER(idx_t),INTENT(IN)      :: DIM
TYPE(BounDescrip2D),INTENT(IN) :: PSLG 
REAL(real_t),INTENT(IN)        :: LMIN       

! OUTPUTS 
REAL(real_t),INTENT(OUT),ALLOCATABLE :: IPTS(:,:)
INTEGER(idx_t),INTENT(OUT) :: NP 

! local to subroutine 
REAL(real_t),ALLOCATABLE :: xvec(:),yvec(:)
REAL(real_t),ALLOCATABLE :: xg(:,:)
REAL(real_t),ALLOCATABLE :: yg(:,:)
!REAL(real_t),ALLOCATABLE :: Dist(:)
INTEGER(idx_t),ALLOCATABLE:: INoOUT(:) 
REAL(real_t),ALLOCATABLE :: R0(:)
REAL(real_t)             :: H
REAL(real_t)             :: a,b
REAL(real_t)             :: temp
INTEGER(idx_t)           :: nx,ny,nz
INTEGER(idx_t)           :: i,j,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1. Create initial distribution in bounding box (equilateral triangles)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! [x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
! x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
! p=[x(:),y(:)];                                       % List of node coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
BBOX(1) = MINVAL(PSLG%Vert(:,1))
BBOX(2) = MAXVAL(PSLG%Vert(:,1))
BBOX(3) = MINVAL(PSLG%Vert(:,2))
BBOX(4) = MAXVAL(PSLG%Vert(:,2))

temp=LMIN*(sqrt(3.0d0)/2.0d0)
NX=FLOOR((BBOX(2) - BBOX(1))/LMIN)+1
NY=FLOOR((BBOX(4) - BBOX(3))/temp)+1

WRITE(*,'(A)') "                                      "
WRITE(*,'(A)') "********************************************************"
WRITE(*,'(A,F6.3,A,F6.3,A,F6.3,A,F6.3)') "INFO: Domain extents are: " &  
       //" XMIN: ",BBOX(1)," XMAX: ",BBOX(2), " " & 
       //" YMIN: ",BBOX(3)," YMAX: ",BBOX(4) 
WRITE(*,'(A)') "********************************************************"
WRITE(*,'(A)') "                                      "

ALLOCATE(XVEC(NX),YVEC(NY))

ALLOCATE(IPTS(NX*NY,2))
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
    IPTS(NP,1) = XG(i,j)
    IPTS(NP,2) = YG(i,j) 
  ENDDO
ENDDO


! 2. Remove points outside the region, apply the rejection method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! p=p(feval(fd,p,varargin{:})<geps,:);                 % Keep only d<0 points.
! r0=1./feval(fh,p,varargin{:}).^2;                    % Probability to keep point.
! p=p(rand(size(p,1),1)<r0./max(r0),:);                % Rejection method.
! N=size(p,1);                                         % Number of initial points.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Can we just call inpoly instead?
!CALL CalcSignedDistance(IPTS,PSLG,Dist) 
ALLOCATE(INoOUT(NP)) 
INoOUT = 0
DO I = 1,NP 
  CALL FPNPOLY(PSLG,IPTS(I,:),INoOUT(I))
ENDDO

!GEPS = 0.01d0 * LMIN  ! this is how it was in DistMesh

! Remove points with dist > geps (mark them)
DO I = 1,NP
!  IF(Dist(I).GT.GEPS) THEN 
  IF(INoOUT(I).EQ.0) THEN 
    IPTS(I,1) = -9999.d0
    IPTS(I,2) = -9999.d0
  ENDIF
ENDDO

!DEALLOCATE(DIST,XG,YG,XVEC,YVEC) 
DEALLOCATE(INoOUT,XG,YG,XVEC,YVEC) 

! push them to the back of the array
CALL PushZerosToBackREAL(IPTS,NP)
!
! Apply rejection method 
! p=p(rand(size(p,1),1)<r0./max(r0),:);  
ALLOCATE(r0(NP)) 
DO I = 1,NP
  H=HFX(IPTS(I,:))
  r0(i)  = 1.0d0/(H**2.0d0)
ENDDO
a = maxval(r0) 
DO I = 1,NP 
  b = r0(i)/a
  if(rand().gt.b) then
    ipts(I,1) = -9999.d0 
    ipts(I,2) = -9999.d0
  endif
ENDDO

CALL PushZerosToBackREAL(IPTS,NP)

print *, NP 

END SUBROUTINE FormInitialPoints2D
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE PushZerosToBackREAL(ARR,NP) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Push array entries with zero back of array and return updated NP size 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUTS 
REAL(kind=real_t),ALLOCATABLE,INTENT(INOUT) :: ARR(:,:) 

! OUTPUTS 
INTEGER(kind=idx_t),INTENT(INOUT) :: NP 

! LOCAL 
INTEGER(kind=idx_t) :: I,J,temp, k

K=0 
DO I = 1,NP 
  IF(ARR(I,1).GT.-9000.0d0) THEN 
    K = K + 1 
    DO J =1,SIZE(ARR,2)
      ARR(K,J)=ARR(I,J) 
    ENDDO
  ENDIF
ENDDO

temp = K ! temporary NP size  

DO WHILE (K < NP ) 
  K = K + 1 
  DO J = 1,SIZE(ARR,2) 
    ARR(K,J) = -9999.d0 
  ENDDO
ENDDO  

NP = temp 

END SUBROUTINE PushZerosToBackREAL
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE PushZerosToBackINT(ARR,NP) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Push array entries with zero back of array and return updated NP size 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUTS 
INTEGER(kind=idx_t),ALLOCATABLE,INTENT(INOUT) :: ARR(:,:) 

! OUTPUTS 
INTEGER(kind=idx_t),INTENT(INOUT) :: NP 

! LOCAL 
INTEGER(kind=idx_t) :: I,J,temp, k

K=0 
DO I = 1,NP 
  IF(ARR(I,1).GT.-9000) THEN 
    K = K + 1 
    DO J =1,SIZE(ARR,2)
      ARR(K,J)=ARR(I,J) 
    ENDDO
  ENDIF
ENDDO

temp = K ! temporary NP size  

DO WHILE (K < NP ) 
  K = K + 1 
  DO J = 1,SIZE(ARR,2) 
    ARR(K,J) = -9999
  ENDDO
ENDDO  

NP = temp 

END SUBROUTINE PushZerosToBackInt
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE CalcSignedDistance(IPTS,PSLG,SignedDistance)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Given a PSLG and a set of points, determine the nearest 
! Euclidean distance to the PSLG using a KD-tree and sign the distance 
! negative if the point in question is inside the polygon defined by the PSLG. 
! and vice-versa otherwise.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use kdtree2_module
USE OMP_LIB
implicit none 

! INPUTS 
TYPE(BounDescrip2D),INTENT(IN) :: PSLG 
REAL(real_t),INTENT(IN),ALLOCATABLE :: IPTS(:,:)

! OUTPUTS 
REAL(real_t),INTENT(OUT),ALLOCATABLE :: SignedDistance(:)

! LOCAL 
TYPE(kdtree2), POINTER :: TREE=>null()
TYPE(kdtree2_result),ALLOCATABLE :: KDRESULTS(:)
INTEGER(kind=idx_t) :: tempSZ
INTEGER(kind=idx_t) :: INoOUT
INTEGER(kind=idx_t) :: I
integer :: nthread,myid  
REAL(real_t) :: t1,t2

! For every point we compute the distance to every point in the PSLG 

! BUILD KD-TREE with PSLG vertices
tree => kdtree2_create(TRANSPOSE(PSLG%Vert),rearrange=.true.,sort=.true.)

! Loop over all the initial points 
tempSZ = SIZE(IPTS,1) 
ALLOCATE(SignedDistance(tempSZ))
ALLOCATE(KDRESULTS(1)) 
INoOUT = -999

DO I =1,tempSZ
  call kdtree2_n_nearest(tp=tree,qv=IPTS(I,:),nn=1,results=KDRESULTS)
  SignedDistance(I) = SQRT(KDRESULTS(1)%DIS)
  !! Determine if point is in the PSLG defined polygon
  CALL FPNPOLY(PSLG,IPTS(I,:),INoOUT)
  IF(INoOUT.EQ.1) THEN 
    SignedDistance(I)= -SignedDistance(I)
  ENDIF
ENDDO

call kdtree2_destroy(tp=tree)

DEALLOCATE(KDRESULTS)

!! DEBUG VISUALIZE SDF
!print *, "WRITING" 
!OPEN(UNIT=305,FILE="SignedDistance.txt",ACTION='WRITE')
!DO i=1,tempSZ
!  WRITE(305,"(3F12.8)") IPTS(i,1),IPTS(i,2),SignedDistance(I)
!ENDDO
!CLOSE(305)
!stop

END SUBROUTINE CalcSignedDistance
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
INoOUT = pnpoly(PSLG%NumVert0, PSLG%Vert0(:,1), PSLG%Vert0(:,2), &
                TEST(1),TEST(2))

END SUBROUTINE 
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE DelTriaWElim(DIM,PSLG, NP, POINTS, NF, FACETS,IERR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calls qhull library from fortran language to output facet table of delaunay triangulation 
! removes triangles with centroids outside of the PSLG 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use iso_c_binding, only : c_f_pointer,C_PTR
implicit none

! INPUTS 
integer(kind=idx_t), intent(in) :: DIM
integer(kind=idx_t), intent(in) :: NP
real(kind=real_t),   intent(in),ALLOCATABLE :: POINTS(:,:)
integer(kind=idx_t), intent(inout):: NF
TYPE(BounDescrip2D) :: PSLG 

! OUTPUTS 
integer(kind=idx_t),intent(out),allocatable :: FACETS(:,:)
integer(kind=idx_t),INTENT(OUT) :: IERR

! LOCAL TO SUBROUTINE 
INTEGER(kind=idx_t) :: INoOUT 
type(c_ptr) :: cfacetemp
integer(kind=idx_t),pointer :: ffacetemp(:)
real(real_t),allocatable :: CENTROIDS(:,:)
real(8) :: blah1,blah2

integer :: i,j,k

! call qhull library 
CALL CPU_TIME(blah1) 
cfacetemp = faces(DIM,NP,TRANSPOSE(POINTS),NF)
CALL CPU_TIME(blah2) 
!print *, blah2 - blah1 

CALL C_F_POINTER(cfacetemp, ffacetemp, [NF*(DIM+1)])

! reshape vector ffacetemp to facets array
! TODO: check if it's possible to allocate to memory with an error condition
allocate(facets(NF,DIM+1))

facets=0 
k=0 
do i = 1,NF
  do j = 1,DIM+1 
    k = k + 1
    facets(i,j)=ffacetemp(k)+1 ! ensure indexing starts at 1
  enddo
enddo

! IMPORTANT: RELEASE MEMORY MALLOC'ed in C 
CALL destroy_storage(cfacetemp)

ALLOCATE(CENTROIDS(NF,DIM)) 
CALL CalcBaryCenter(DIM,POINTS,NP,FACETS,NF,CENTROIDS)  

! Remove triangles whose centroid is out of the domain 
DO I = 1,NF
  CALL FPNPOLY(PSLG,CENTROIDS(I,:),INoOUT)
  IF(INoOUT.EQ.0) THEN 
    FACETS(I,:)=-9000
  ENDIF
ENDDO

CALL PushZerosToBackINT(FACETS,NF) 

! TODO: BETTER ERROR CHECKING (SEE ZOLTAN ERROR LEVELS FOR REFERENCE?)
ierr = 0 !! SUCCESS 

END SUBROUTINE DelTriaWElim
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE CalcBaryCenter(DIM,POINTS,NP,FACES,NF,CENTROIDS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the centroid of the triangles described in the arrays POINTS FACES
! For triangles: 
! centroids = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none 

! INPUTS 
REAL(real_t),INTENT(IN),ALLOCATABLE   :: POINTS(:,:)
INTEGER(idx_t),INTENT(IN),ALLOCATABLE :: FACES(:,:)
INTEGER(idx_t),INTENT(IN)           :: NP,NF,DIM

! OUTPUTS
REAL(real_t),INTENT(OUT),ALLOCATABLE :: CENTROIDS(:,:)

! LOCAL
INTEGER(idx_t) :: tempF(DIM+1) 
REAL(real_t)   :: temp(DIM,DIM+1) 
INTEGER(idx_t) :: i,j,k

ALLOCATE(CENTROIDS(NF,DIM)) 

DO I = 1,NF
  tempF = FACES(I,:)
  DO J = 1,DIM+1 ! for each face
    temp(:,J) = POINTS(tempF(J),:)
  ENDDO
  DO K = 1,DIM
    CENTROIDS(I,K) = SUM(temp(K,:))/DBLE(DIM+1)
  ENDDO
ENDDO

!OPEN(UNIT=303,FILE="debug.txt",ACTION='WRITE')
!do i = 1,nf
!  WRITE(303,"(2F12.8)")CENTROIDS(1,i),CENTROIDS(2,i)
!ENDDO
!CLOSE(303)

END SUBROUTINE CalcBaryCenter
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

SUBROUTINE meshgrid2D(xgv, ygv, X, Y)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create a struture grid of points using the vectors xgv and ygv
! this mimics what the MATLAB command does 
! [xg,yg]=meshgrid(xvec,yvec)
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
! vec = linspace(a,b,n) 
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

SUBROUTINE TriaToTria(triangle_num,triangle_node, triangle_neighbor,triangle_notsharedvertex )

!*****************************************************************************80
!
!! TriaToTria determines triangle neighbors in a mesh (triangle_neighbor)
!  and the vertex index (in the neighbor triangle) not on the shared edge
!  between a given triangle and its neighbor (triangle_notsharedvertex)
!
! THIS CODE HAS BEEN ADAPTED FROM:
! https://people.sc.fsu.edu/~jburkardt/f_src/triangulation_triangle_neighbors/triangulation_triangle_neighbors.html
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Last Modified:
!
!    1 September 2019
!
!  Author:
!
!    John Burkardt
!  
!  Adapated/Modified by: 
! 
!    Keith Roberts 
  implicit none

  ! INPUTS 
  integer(idx_t),intent(in)::triangle_num
  integer(idx_t),intent(in),ALLOCATABLE::triangle_node(:,:)

  ! OUTPUTS 
  integer(idx_t),intent(out),ALLOCATABLE::triangle_neighbor(:,:)
  integer(idx_t),intent(out),ALLOCATABLE::triangle_notsharedvertex(:,:)

  ! LOCALS
  integer(idx_t) ::triangle_order=3
  integer(idx_t),allocatable :: idx(:)
  integer(idx_t),allocatable :: col(:,:),tcol(:,:)
  integer(idx_t) :: Nei(3), Cur(3)
  integer(idx_t) i
  integer(idx_t) icol
  integer(idx_t) j
  integer(idx_t) k
  integer(idx_t) p,pp
  integer(idx_t) side1
  integer(idx_t) side2
  integer(idx_t) tri
  integer(idx_t) tri1
  integer(idx_t) tri2
  logical :: found 

!  Step 1.
!  From the list of nodes for triangle T, of the form: (I,J,K)
!  construct the three neighbor relations:
!
!    (I,J,3,T) or (J,I,3,T),
!    (J,K,1,T) or (K,J,1,T),
!    (K,I,2,T) or (I,K,2,T)
!
!  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
!
  allocate(col(3*triangle_num,4))
  allocate(tcol(3*triangle_num,4))

  do tri = 1, triangle_num

    i = triangle_node(tri,1)
    j = triangle_node(tri,2)
    k = triangle_node(tri,3)

    if ( i < j ) then
      col(3*(tri-1)+1,1:4) = (/ i, j, 3, tri /)
    else
      col(3*(tri-1)+1,1:4) = (/ j, i, 3, tri /)
    end if

    if ( j < k ) then
      col(3*(tri-1)+2,1:4) = (/ j, k, 1, tri /)
    else
      col(3*(tri-1)+2,1:4) = (/ k, j, 1, tri /)
    end if

    if ( k < i ) then
      col(3*(tri-1)+3,1:4) = (/ k, i, 2, tri /)
    else
      col(3*(tri-1)+3,1:4) = (/ i, k, 2, tri /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1 and 2; the routine we call here
!  sorts on rows 1 through 4 but that won't hurt us.
!
!  What we need is to find cases where two triangles share an edge.
!  Say they share an edge defined by the nodes I and J.  Then there are
!  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
!  we make sure that these two columns occur consecutively.  That will
!  make it easy to notice that the triangles are neighbors.
!
! Sort the COL array 
  call sortRows(3*Triangle_num,col,IDX) 

  ! Sort the COL array 
  tCOL=0
  DO I = 1, 3*Triangle_num
    tCOL(I,1:4) = COL(IDX(I),1:4) 
  ENDDO
  COL = tCOL 
  deallocate(tCOL)

!
!  Step 3. Neighboring triangles show up as consecutive columns with
!  identical first two entries.  Whenever you spot this happening,
!  make the appropriate entries in TRIANGLE_NEIGHBOR.
!
  allocate(triangle_neighbor(triangle_num,triangle_order))

  triangle_neighbor = -1 

  icol = 1

  do

    if ( 3 * triangle_num <= icol ) then
      exit
    end if

    if ( col(icol,1) /= col(icol+1,1) .or. col(icol,2) /= col(icol+1,2) ) then
      icol = icol + 1
      cycle
    end if

    side1 = col(icol,3)
    tri1 = col(icol,4)
    side2 = col(icol+1,3)
    tri2 = col(icol+1,4)

    triangle_neighbor(tri1,side1) = tri2
    triangle_neighbor(tri2,side2) = tri1

    icol = icol + 2

  end do

  ! compute t2n 
  !find(~ismember(t(t2t(j,k),:), t(j,:)))    % Should be = t2n(j,k)
  allocate(triangle_NotSharedVertex(triangle_num,triangle_order)) 
  
  triangle_NotSharedVertex = -1 

  do tri = 1,triangle_num
    Cur(:) = triangle_node(tri,1:3) 
    do k = 1,3 ! adj
      if(triangle_neighbor(tri,k).eq.-1) cycle 
      Nei(:) = triangle_node(triangle_neighbor(tri,k),1:3)

      do p = 1,3
        found=.false. 
        do pp = 1,3 
          if(Nei(p).eq.Cur(pp)) then
            ! Nei and Curr share this vertex 
            found=.true.
          endif
        enddo
        ! if we haven't found this vertex yet, it must be *not shared*
        if(found.eqv..false.) then 
          triangle_notSharedVertex(tri,k) = p
        endif
      enddo

    enddo
  enddo

END SUBROUTINE TriaToTria
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

recursive subroutine quicksortINT(a, first, last)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  implicit none
  integer(4)  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksortINT(a, first, i-1)
  if (j+1 < last)  call quicksortINT(a, j+1, last)
end subroutine quicksortINT
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

recursive subroutine quicksortINTLONG(a, first, last,idx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  implicit none
  integer(8)  a(*), x, t
  integer first, last
  integer(idx_t),intent(inout),allocatable :: idx(:)
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     t = idx(i);idx(i) = idx(j); idx(j) = t ! save the sort map
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksortINTLONG(a, first, i-1,idx)
  if (j+1 < last)  call quicksortINTLONG(a, j+1, last,idx)
end subroutine quicksortINTLONG
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

recursive subroutine quicksortREAL(a, first, last)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! quicksort.f -*-f90-*-
! Author: t-nissie
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  implicit none
  real(8)  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksortREAL(a, first, i-1)
  if (j+1 < last)  call quicksortREAL(a, j+1, last)
end subroutine quicksortREAL
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

function median(a,first,last)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the median of a vector of reals 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE

! INPUTS 
REAL(real_t),intent(in) :: a(*)
INTEGER(idx_t),intent(in) :: first, last

! OUTPUTS
REAL(real_t) :: median 

! LOCALS 
INTEGER(idx_t) :: lng
REAL(real_t),ALLOCATABLE :: at(:)

lng = last - first + 1

allocate(at(lng))
at=a(first:last)

! 1. sort the numbers
CALL quicksortREAL(at,first,last)

! 2. If length is odd 
if(MOD(LNG,2).EQ.1) then
  ! then odd length vec
  median = at(floor(lng/2.0d0)+1)
else
  ! then even length vec
  median = ( at((lng/2)+1) + a((lng/2)-1) )/2.0d0
endif

deallocate(at)

end function median
!*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*!

END MODULE UTILS


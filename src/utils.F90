!-----------------------------------------------------------------------
!  MODULE UTILS
!-----------------------------------------------------------------------
!> @author Keith J. Roberts, Universidade de Sao Paulo, krober@usp.br 
!>
!> @brief The program called DistMesh uses the procedures contained herein.
!>        First public functions are declared followed by private functions.
!-----------------------------------------------------------------------
MODULE UTILS 
!-----------------------------------------------------------------------
use YourMeshSize
use vars

implicit none

!----------------------------begin of c-interface-----------------------
INTERFACE
!-----------------------------------------------------------------------
!>@brief call the inpolygon algorithm written in c from fortran.
!-----------------------------------------------------------------------
FUNCTION pnpoly(NUMVERT, VERTx, VERTy, TESTx, TESTy)bind(c,name='pnpoly')
import idx_t,real_t
implicit none
integer(kind=idx_t) :: pnpoly 
integer(kind=idx_t), value, intent(in) :: NUMVERT
real(kind=real_t),   value, intent(in) :: TESTx
real(kind=real_t),   value, intent(in) :: TESTy
real(kind=real_t),intent(in)           :: VERTx(NUMVERT)
real(kind=real_t),intent(in)           :: VERTy(NUMVERT)
END FUNCTION pnpoly
!-----------------------------------------------------------------------
END INTERFACE
!-----------------------------------------------------------------------
!> @brief calls qhull library from fortran language to produce facet table
!>        of delaunay triangulation 
!-----------------------------------------------------------------------
INTERFACE 
FUNCTION faces(DIM, NUMPOINTS, fpoints, NUMFACETS)bind(c,name='faces')
import idx_t,real_t,c_ptr
implicit none
type(c_ptr) :: faces
integer(kind=idx_t), value, intent(in) :: DIM
integer(kind=idx_t), value, intent(in) :: NUMPOINTS 
integer(kind=idx_t), intent(inout)     :: NUMFACETS
real(kind=real_t), intent(in)          :: fpoints(DIM,NUMPOINTS)
END FUNCTION faces
!-----------------------------------------------------------------------
END INTERFACE 
!-----------------------------------------------------------------------
!> @brief releases memory that was allocated in c in the faces function call 
!-----------------------------------------------------------------------
INTERFACE
SUBROUTINE destroy_storage(p) BIND(C, NAME='destroy_storage')
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
IMPLICIT NONE
TYPE(C_PTR), INTENT(IN), VALUE :: p
END SUBROUTINE destroy_storage
!-----------------------------------------------------------------------
END INTERFACE
!---------------------end of c-interface------------------------------------


PUBLIC edgeFlipper,ProjectPointsBack,ReadPSLGtxt,FormInitialPoints2D
PUBLIC DelTriaWElim,TriaToTria,findUniqueBars,CalcForces,ApplyForces
PUBLIC WriteMesh


PRIVATE sortRows,quickSortINTLONG,densify,MeshGrid2D,PushZerosToBackREAL
PRIVATE PushZerosToBackINT,FPNPOLY,CalcSignedDistance,CalcBaryCenter
PRIVATE linspace,meshexp,quicksortINT,quickSortREAL,median

!---------------------end of data declarations--------------------------------

CONTAINS 



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! P U B L I C  F U N C T I O N S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
FUNCTION ParseInputs() RESULT(SizingFields)
!-----------------------------------------------------------------------
IMPLICIT NONE 

TYPE(MetricTensor) :: SizingFields 
LOGICAL :: FileFound=.FALSE.

! check if it exists 
INQUIRE(FILE="Mesh.inp",EXIST=fileFound)

IF(fileFound) THEN
  OPEN(UNIT=1,FILE="Mesh.inp",ACTION='READ')
  READ(1,*) MaxIter
  READ(1,*) DeltaT 
  READ(1,*) NSCREEN  
ELSE
  WRITE(*,'(3A)') "FATAL: Mesh.inp control file not found."
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(A)') "                                                        "
  STOP  
ENDIF
CLOSE(1)


IF(COMMAND_ARGUMENT_COUNT().LT.1)THEN
  WRITE(*,*)'ERROR, COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
  STOP
ENDIF


CALL GET_COMMAND_ARGUMENT(2,sizefname)   
CALL GET_COMMAND_ARGUMENT(3,elongfname)                    
CALL GET_COMMAND_ARGUMENT(4,anglefname)                    

SzFx    = LoadMeshSizes(sizefname)                           ! Load size function into memory
ElongFx = LoadMeshElong(elongfname)                          ! Load in elongation factors into memory
AngleFx = LoadMeshAngle(anglefname)                          ! Load in angles into memory  

CALL GET_COMMAND_ARGUMENT(1,pslgfname)   
PSLG = ReadPSLGtxt(pslgfname,SzFx)                           ! Load in boundary description 

SizingFields%Iso = SzFx 
SizingFields%Elong = ElongFx
SizingFields%Angle = AngleFx

!-----------------------------------------------------------------------
END FUNCTION ParseInputs
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Flips edges to achieve Delaunay triangulation based on metric tensor.
!> @param[in]   MetricTensor contains the spatially-variable elongation factor
!> @param[in]     NP         number of points in triangulation
!> @param[in]     POINTS     vertices of the triangulation
!> @param[in]     NF         number of faces in the triangulation
!> @param[inout]  FACETS     table vertex indices describing faces 
!> @param[inout]  T2N        table of vertex index not shared in nei. tria.
!> @param[inout]  T2T        table of tria.-to-tria. adj.
!> @param[out]    NUMFLIPS   number of edge flips that were performed 
!-----------------------------------------------------------------------
SUBROUTINE edgeFlipper(ElongFx,NP,POINTS,NF,FACETS,T2N,T2T,NUMFLIPS) 
!-----------------------------------------------------------------------
INTEGER(idx_t),INTENT(IN)                :: NP,NF 
TYPE(MetricTensor),INTENT(IN)            :: ElongFx
REAL(real_t),  INTENT(IN),ALLOCATABLE    :: POINTS(:,:)
INTEGER(idx_t),INTENT(INOUT),ALLOCATABLE :: FACETS(:,:)
INTEGER(idx_t),INTENT(INOUT),ALLOCATABLE :: T2N(:,:)
INTEGER(idx_t),INTENT(INOUT),ALLOCATABLE :: T2T(:,:)
INTEGER(idx_t),INTENT(INOUT)             :: numflips 

INTEGER(idx_t) :: i
INTEGER(idx_t) :: t1
INTEGER(idx_t) :: t2
INTEGER(idx_t) :: n1,n2
INTEGER(idx_t) :: tix11,tix12,tix13,tix21,tix22,tix23
INTEGER(idx_t) :: newt(1:2,1:3) 
INTEGER(idx_t) :: nbt,nbn
INTEGER(idx_t) :: mod3x1(1:3),mod3x2(1:3),mod3x3(1:3)
REAL(real_t)   :: cp1,cp2
REAL(real_t)   :: temp1(1:1,1:2),temp2(1:1,1:1),temp3(1:1,1:2),temp4(1:1,1:1)
REAL(real_t)   :: edgeAV(1:2,1:1),edgeBV(1:2,1:1),edgeCV(1:2,1:1),edgeDV(1:2,1:1)
REAL(real_t)   :: edgeCV_t(1:1,1:2),edgeAV_t(1:1,1:2)
REAL(real_t)   :: ME(1:2,1:2)
REAL(real_t)   :: pmid(2)
REAL(real_t)   :: test(1:1,1:1)
LOGICAL        :: FLIP

mod3x1(1:3)=(/2,3,1/)
mod3x2(1:3)=(/3,1,2/)
mod3x3(1:3)=(/1,2,3/) 

numflips = 0 

DO T1 = 1,NF 
  DO N1 = 1,3
    FLIP=.FALSE. 
    T2=T2T(T1,N1) 
    IF(T2.GE.1) THEN 
      n2    = T2N(t1,n1) 

      tix11 = mod3x1(n1) 
      tix12 = mod3x2(n1) 
      tix13 = mod3x3(n1)

      tix21 = mod3x1(n2) 
      tix22 = mod3x2(n2) 
      tix23 = mod3x3(n2)
      
      newt(1,1:3) = (/FACETS(t1,1),FACETS(t1,2),FACETS(t1,3) /)
      newt(2,1:3) = (/FACETS(t2,1),FACETS(t2,2),FACETS(t2,3) /)

      ! compute midpoint of quadilateral of connected triangles 
      pmid(1) = POINTS(NEWT(1,1),1)+POINTS(NEWT(1,2),1)+POINTS(NEWT(1,3),1)
      pmid(1) = pmid(1)+POINTS(NEWT(2,1),1)+POINTS(NEWT(2,2),1)+POINTS(NEWT(2,3),1)
      pmid(1) = pmid(1)/6.0d0
      
      pmid(2) = POINTS(NEWT(1,1),2)+POINTS(NEWT(1,2),2)+POINTS(NEWT(1,3),2)
      pmid(2) = pmid(2)+POINTS(NEWT(2,1),2)+POINTS(NEWT(2,2),2)+POINTS(NEWT(2,3),2)
      pmid(2) = pmid(2)/6.0d0

      ! eval metric tensor at mdpt 
      ME = CalcMetricTensor(pmid,ElongFx) 

      ! vectors of shape 2x1 after transpose from 1x2 shape 
      temp1 = POINTS(newt(1,tix13):newt(1,tix13),1:2) - POINTS(newt(1,tix11):newt(1,tix11),1:2)
      edgeBV = TRANSPOSE(temp1) 
      temp1  = POINTS(newt(2,tix23):newt(2,tix23),1:2) - POINTS(newt(2,tix21):newt(2,tix21),1:2)
      edgeDV = TRANSPOSE(temp1)
      temp1 = POINTS(newt(1,tix13):newt(1,tix13),1:2) - POINTS(newt(2,tix21):newt(2,tix21),1:2)
      edgeCV = TRANSPOSE(temp1)
      temp1  = POINTS(newt(2,tix23):newt(2,tix23),1:2) - POINTS(newt(1,tix11):newt(1,tix11),1:2)
      edgeAV = TRANSPOSE(temp1)

      ! cross-->(v1.X*v2.Y) - (v1.Y*v2.X);
      cp1 = (edgeAV(1,1)*edgeBV(2,1)) - (edgeAV(2,1)*edgeBV(1,1))  
      cp2 = (edgeCV(1,1)*edgeDV(2,1)) - (edgeCV(2,1)*edgeDV(1,1))  

      EDGECV_t = TRANSPOSE(EDGECV) ! 1 x 2
      EDGEAV_t = TRANSPOSE(EDGEAV)

      ! Del. criterion
      temp1 = MATMUL(EDGECV_t,ME) !1x2 x 2x2 result is 1x2 
      temp2 = MATMUL(temp1,edgeDV)!1x2 x 2x1 result is 1x1

      temp3 = MATMUL(EDGEAV_t,ME) !1x2 x 2x2 result is 1x2
      temp4 = MATMUL(temp3,EDGEBV)!1x2 x 2x1 result is 1x1

      test = (cp1*temp2) + (cp2*temp4) 
      
      IF(test(1,1).GT.0.d0) THEN 
        FLIP=.TRUE.
      ENDIF
       
      IF(FLIP.eqv..true.) THEN 
        numflips = numflips + 1
        ! swap edge 
        newt(1,tix12) = newt(2,n2) 
        newt(2,tix22) = newt(1,n1)

        ! Insert new triangles (with flipped edges)
        FACETS(T1,:) = NEWT(1,:) 
        FACETS(T2,:) = NEWT(2,:) 
        
        ! Update t2t and t2n
        NBT = T2T(t2,tix21)
        NBN = T2N(t2,tix21) 

        T2T(T1,N1)=NBT 
        T2N(T1,N1)=NBN 
        
        IF(NBT.GE.1) THEN 
          T2T(nbt,nbn)=t1
          T2N(nbt,nbn)=n1
        ENDIF

        nbt=t2t(t1,tix11) 
        nbn=t2n(t1,tix11) 
        T2T(t2,n2)=nbt
        T2N(t2,n2)=nbn

        if(nbt.GE.1) THEN
          T2T(nbt,nbn)=t2
          T2N(nbt,nbn)=n2
        endif
        
        T2T(t1,tix11)=t2
        T2N(t1,tix11)=tix21
        T2T(t2,tix21)=t1
        T2N(t2,tix21)=tix11

      ENDIF
    ENDIF
  ENDDO
ENDDO

!  print *, "NUMFLIPS : ", numflips

!-----------------------------------------------------------------------
END SUBROUTINE edgeFlipper 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Projects vertices that have exited the domain back into the domain
!>        using the path of steepest descent in the signed distance function.
!> @param[in]          NP     number of points in the mesh
!> @param[in]        PSLG     planar-straight line graph 
!> @param[inout]   POINTS     coordinates of points in the mesh.
!-----------------------------------------------------------------------
SUBROUTINE ProjectPointsBack(PSLG,POINTS,NP)
!-----------------------------------------------------------------------
IMPLICIT NONE 

INTEGER(idx_t),INTENT(IN) :: NP 
TYPE(BounDescrip2D),INTENT(IN) :: PSLG 
REAL(real_t),INTENT(INOUT),ALLOCATABLE :: POINTS(:,:)

INTEGER(idx_t) :: I
INTEGER(idx_t) :: NOUT
REAL(real_t),ALLOCATABLE :: Dist(:)
REAL(real_t),ALLOCATABLE :: dgradx(:),dgrady(:),dgrad2(:)
REAL(real_t),ALLOCATABLE :: tempDistx(:),tempDisty(:)
REAL(real_t),ALLOCATABLE :: tempP1(:,:),tempP2(:,:)
REAL(real_t),ALLOCATABLE :: PointsOut(:,:),DistOut(:)
INTEGER(idx_t),ALLOCATABLE :: INoOUT(:) 
INTEGER(idx_t) :: DIM =2

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
  dgradx= (tempDistx(1:NOUT) - DistOut(1:NOUT))/DEPS 

!      dgrady = (feval(obj.fd,[p(ix,1),p(ix,2)+deps],obj,[])...%,1)...
!                -d(ix))/deps; % gradient
  tempP2 = PointsOUT 
  tempP2(:,2) = tempP1(:,2) + DEPS 
  CALL CalcSignedDistance(tempP2,PSLG,tempDisty) 
  dgrady= (tempDisty(1:NOUT) - DistOut(1:NOUT))/DEPS 

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

DEALLOCATE(dgradx,dgrady,dgrad2)
DEALLOCATE(tempDistx,tempDisty) 
DEALLOCATE(tempP1,tempP2)

ENDIF

! RELEASE ALLOCATABLES 
DEALLOCATE(InoOut,DistOut,PointsOut)

!-----------------------------------------------------------------------
END SUBROUTINE ProjectPointsBack 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Update vertex locations according to calculated forces. 
!-----------------------------------------------------------------------
SUBROUTINE ApplyForces(POINTS,NP,BARS,NUMBARS,FVEC)
!-----------------------------------------------------------------------
IMPLICIT NONE 

! INPUTS 
INTEGER(idx_t),INTENT(IN) :: NP,NUMBARS 
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

!-----------------------------------------------------------------------
END SUBROUTINE ApplyForces
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Calculate forcing terms for each vertex based on mesh size function 
!>        and actual length of the edge. 
!>        You can comment out the default Bossens-Heckbert potential function
!>        for the spring-based force function. 
!-----------------------------------------------------------------------
SUBROUTINE CalcForces(MeshSizes,POINTS,NP,BARS,NUMBARS,FVEC) 
!-----------------------------------------------------------------------
implicit none 

! INPUTS
TYPE(MetricTensor),INTENT(IN) :: MeshSizes
INTEGER(idx_t),INTENT(IN) :: NP,NUMBARS
INTEGER(idx_t),INTENT(IN),ALLOCATABLE :: BARS(:,:)
REAL(real_t),INTENT(IN),ALLOCATABLE :: POINTS(:,:)

! OUTPUTS
REAL(real_t),INTENT(OUT),ALLOCATABLE :: FVEC(:,:)

! LOCAL 
INTEGER(idx_t) :: i 
REAL(real_t)   :: barvec(NUMBARS,2)
REAL(real_t)   :: L(1:NUMBARS,1:1),L0(1:NUMBARS,1:1),LN(1:NUMBARS,1:1)
REAL(real_t)   :: FORCES(1:NUMBARS)
REAL(real_t)   :: HBARS(1:NUMBARS,1:1)
REAL(real_t)   :: temp1(1:1,1:2),temp2(1:1,1:1),VEC1(1:2,1:1),VEC1_t(1:1,1:2)
REAL(real_t)   :: midpt(1:2,1:1),a,b,SCALE_FACTOR
REAL(real_t)   :: ME(1:2,1:2)
INTEGER(idx_t) :: DIM=2 

DIM = 2
ALLOCATE(FVEC(NUMBARS,2))
! L(jj,1) = sqrt(rij'*M*rij);
FORCES=0.0d0
DO I = 1,NUMBARS
  ! To calculate edgelength, assume edge extrues along ideal length 
  MIDPT(1:2,1)=(POINTS(BARS(I,1),:) + POINTS(BARS(I,2),:))/2.0d0
  ME = CalcMetricTensor(MIDPT(1:2,1),MeshSizes) ! query metric tensor at midpoint
  HBARS(I,1)=CalcMeshSize(MIDPT(1:2,1),MeshSizes) ! query the ideal isotropic element size

  ! assume ideal edge extrudes at ideal angle 
  A = CalcMeshAngle(MIDPT(1:2,1),MeshSizes) 
  A = COS(A) ! in radians 
  A = A*(180.d0/3.14d0) ! in degrees  
  VEC1_t(1,1) = MIDPT(1,1) - (MIDPT(1,1)-A*HBARS(I,1)) ! cos(a)

  A = CalcMeshAngle(MIDPT(1:2,1),MeshSizes) 
  A = SIN(A) ! in radians 
  A = A*(180.d0/3.14d0) ! in degrees  
  VEC1_t(1,2) = MIDPT(2,1) - (MIDPT(2,1)-A*HBARS(I,1)) ! sin(a)

  VEC1 = transpose(VEC1_t) 
  temp1(1:1,1:2)=MATMUL(VEC1_t(1:1,1:2),ME(1:2,1:2))
  temp2(1:1,1:1)=MATMUL(temp1(1:1,1:2),VEC1(1:2,1:1))
  HBARS(I:I,1:1)=SQRT(temp2) 

  ! compute the actual length in metric space 
  BARVEC(I,:)=POINTS(BARS(I,1),:) - POINTS(BARS(I,2),:)
  VEC1_t(1:1,1:2)=BARVEC(I:I,1:2)
  VEC1=TRANSPOSE(VEC1_t) 
  temp1(1:1,1:2)=MATMUL(VEC1_t(1:1,1:2),ME(1:2,1:2))
  temp2(1:1,1:1)=MATMUL(temp1(1:1,1:2),VEC1(1:2,1:1))
  L(I:I,1:1)=SQRT(temp2) 
ENDDO 

a=MEDIAN(L,1,NUMBARS) ; b = MEDIAN(HBARS,1,NUMBARS)
SCALE_FACTOR = a/b 

DO I = 1,NUMBARS
  L0(I,1)=HBARS(I,1)*FSCALE*SCALE_FACTOR
  LN(I,1)=L(I,1)/L0(I,1)
  ! Bossens Heckbert 
  FORCES(I)=(1-LN(I,1)**4)*EXP(-LN(I,1)**4)/LN(I,1)
  ! Linear spring (Hooke's Law)
  !FORCES(I)=MAXVAL( (/1.0d0-LN(I,1),0.0d0/))
  FVEC(I,1)=FORCES(I)*BARVEC(I,1)
  FVEC(I,2)=FORCES(I)*BARVEC(I,2)
ENDDO


!-----------------------------------------------------------------------
END SUBROUTINE CalcForces
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief TriaToTria determines triangle neighbors in a mesh (triangle_neighbor)
!>         and the vertex index (in the neighbor triangle) not on the shared edge
!>         between a given triangle and its neighbor (triangle_notsharedvertex)
!>         THIS CODE HAS BEEN ADAPTED FROM:
!> https://people.sc.fsu.edu/~jburkardt/f_src/triangulation_triangle_neighbors/triangulation_triangle_neighbors.html
!-----------------------------------------------------------------------
SUBROUTINE TriaToTria(triangle_num,triangle_node, triangle_neighbor,triangle_notsharedvertex )
!-----------------------------------------------------------------------
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

  triangle_neighbor = 0

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
  
  triangle_NotSharedVertex = 0

  do tri = 1,triangle_num
    Cur(:) = triangle_node(tri,1:3) 
    do k = 1,3 ! adj
      if(triangle_neighbor(tri,k).eq.0) cycle 
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

!-----------------------------------------------------------------------
END SUBROUTINE TriaToTria
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Computes the unique bars of the mesh 
!>        Replicates the following MATLAB 
!>        bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         
!>        bars=unique(sort(bars,2),'rows');      
!-----------------------------------------------------------------------
SUBROUTINE FindUniqueBars(NF,FACES,NumBars,BARS) 
!-----------------------------------------------------------------------
implicit none 

! INPUT 
INTEGER(idx_t),INTENT(IN) :: NF 
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

!-----------------------------------------------------------------------
END SUBROUTINE findUniqueBars 
!-----------------------------------------------------------------------


!> @brief Write the mesh to disk.
!-----------------------------------------------------------------------
SUBROUTINE WriteMesh(POINTS,NP,FACETS,NF,ITER)
!-----------------------------------------------------------------------
implicit none 

! INPUT 
REAL(real_t),INTENT(IN),ALLOCATABLE :: POINTS(:,:) 
INTEGER(idx_t),INTENT(IN),ALLOCATABLE :: FACETS(:,:) 
INTEGER(idx_t),INTENT(IN) :: NF,NP,ITER

! OUTPUT 
! WRITES SEVERAL TEXT FILES WITH POINTS AND FACETS 

! LOCAL 
INTEGER(idx_t) :: i 
CHARACTER(14) :: filename
INTEGER(idx_t) :: DIM=2

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

!-----------------------------------------------------------------------
END SUBROUTINE WriteMesh
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Read in Planar Straight Line Graph (PSLG) from a text file named 
!>        'PSLG.txt' into derived type BounDescrip2D. Interpolate points on boundary 
!>         so spacing between two points does nto exceed LMIN/2
!-----------------------------------------------------------------------
FUNCTION ReadPSLGtxt(fname,SzFx) RESULT(PSLG) 
!-----------------------------------------------------------------------
implicit none 

! INPUTS 
CHARACTER(*),intent(in) :: fname
TYPE(GridData),intent(in) :: SzFx       

! OUTPUTS
TYPE(BounDescrip2D) :: PSLG 

! LOCAL TO SUBROUTINE 
INTEGER(idx_t) :: LenFN
logical        :: fileFound 
integer        :: tempNP,tempDIM
integer        :: i 

! get length of the fname containing the pslg data 
LenFN = LengthString(fname)
IF(LenFN .lt. 1) THEN 
  write(*,'(A)') "FATAL: FNAME of PSLG is invalid." 
  stop
ENDIF !error out

! check if it exists 
INQUIRE(FILE=fname(1:LenFN),EXIST=fileFound)

IF(fileFound) THEN
  OPEN(UNIT=1,FILE=fname(1:LenFN),ACTION='READ')
  READ(1,*) tempNP,tempDIM ! number of points in file  
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(A,A,A,I5,A,I5,A)') "INFO: Reading PSLG file called ",fname(1:LenFN)," with " &  
      //" ",tempNP," points in ",tempDIM," dimensions."

  ALLOCATE(PSLG%Vert(tempNP,tempDIM))
  PSLG%NumVert = tempNP 
  PSLG%dim = tempDIM 
  DO I = 1,PSLG%NumVert
    READ(1,*) PSLG%Vert(I,1),PSLG%Vert(I,2)
  ENDDO
  CLOSE(1)
  WRITE(*,'(3A)') "INFO: Read PSLG file called ",fname(1:LenFN)," succesfully!"
  WRITE(*,'(A)') "********************************************************"
ELSE
  WRITE(*,'(3A)') "FATAL: ",fname(1:LenFN) ," file not found."
  WRITE(*,'(A)') "********************************************************"
  STOP  
ENDIF

CALL densify(SzFx%LMIN/2.0d0,PSLG) 

!-----------------------------------------------------------------------
END FUNCTION ReadPSLGtxt 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Form initial points to fill in the domain according to mesh size function
!>         This corresponds with steps 1 to 2 of the distmesh algorithm. 
!-----------------------------------------------------------------------
subroutine FormInitialPoints2D(SzFields,PSLG,IPTS,NP)
implicit none 

! INPUTS 
TYPE(MetricTensor),INTENT(IN)  :: SzFields 
TYPE(BounDescrip2D),INTENT(IN) :: PSLG 

! OUTPUTS 
REAL(real_t),INTENT(OUT),ALLOCATABLE :: IPTS(:,:)
INTEGER(idx_t),INTENT(OUT) :: NP 

! local to subroutine 
REAL(real_t)   :: temp1(1:1,1:2),temp2(1:1,1:1),VEC1(1:2,1:1),VEC1_t(1:1,1:2)
REAL(real_t)   :: ME(1:2,1:2)
REAL(real_t)   :: LMIN
REAL(real_t),ALLOCATABLE :: xvec(:),yvec(:)
REAL(real_t),ALLOCATABLE :: xg(:,:)
REAL(real_t),ALLOCATABLE :: yg(:,:)
INTEGER(idx_t),ALLOCATABLE:: INoOUT(:) 
REAL(real_t),ALLOCATABLE :: R0(:,:)
REAL(real_t)             :: H(1:1,1:1)
REAL(real_t)             :: a,b
REAL(real_t)             :: temp
INTEGER(idx_t)           :: nx,ny,nz
INTEGER(idx_t)           :: i,j,k
INTEGER(idx_t)           :: DIM

LMIN = SzFields%Iso%LMIN
DIM  = PSLG%DIM

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
ALLOCATE(INoOUT(NP)) 
INoOUT = 0
DO I = 1,NP 
  CALL FPNPOLY(PSLG,IPTS(I,:),INoOUT(I))
ENDDO

! Remove points with dist > geps (mark them)
DO I = 1,NP
  IF(INoOUT(I).EQ.0) THEN 
    IPTS(I,1) = -9999.d0
    IPTS(I,2) = -9999.d0
  ENDIF
ENDDO

DEALLOCATE(INoOUT,XG,YG,XVEC,YVEC) 

! push them to the back of the array
CALL PushZerosToBackREAL(IPTS,NP)

! Apply rejection method 
! p=p(rand(size(p,1),1)<r0./max(r0),:);  
ALLOCATE(r0(1:NP,1:1)) 

DO I = 1,NP

  H=CalcMeshSize(IPTS(I,:),SzFields)

  ME=CalcMetricTensor(IPTS(I,:),SzFields)
  
  ! assume ideal edge extrudes at ideal angle 
  A = CalcMeshAngle(IPTS(I,:),SzFields) 
  A = COS(A) ! in radians 
  A = A*(180.d0/3.14d0) ! in degrees  
  VEC1_t(1,1) = IPTS(1,1) - (IPTS(1,1)-A*H(1,1)) ! cos(a)
  A = SIN(A) ! in radians 
  A = A*(180.d0/3.14d0) ! in degrees  
  VEC1_t(1,2) = IPTS(1,2) - (IPTS(1,2)-A*H(1,1)) ! sin(a)
  VEC1 = transpose(VEC1_t) 
  temp1(1:1,1:2)=MATMUL(VEC1_t(1:1,1:2),ME(1:2,1:2))
  temp2(1:1,1:1)=MATMUL(temp1(1:1,1:2),VEC1(1:2,1:1))
  H=SQRT(temp2) 
  r0(i,1)  = 1.0d0/(H(1,1)**2.0d0)
ENDDO
a = maxval(r0) 
DO I = 1,NP 
  b = r0(i,1)/a
  if(rand().gt.b) then
    ipts(I,1) = -9999.d0 
    ipts(I,2) = -9999.d0
  endif
ENDDO

CALL PushZerosToBackREAL(IPTS,NP)

IF(NP.LT.3) THEN 
  WRITE(*,'(A)') "FATAL: NOT ENOUGH INITIAL POINTS TO MESH" 
  STOP
ENDIF

WRITE(*,'(A,I8,A)') "INFO: There are ",NP," initial points in the domain."
WRITE(*,'(A)') "********************************************************"
WRITE(*,'(A)') "                                      "


!-----------------------------------------------------------------------
END SUBROUTINE FormInitialPoints2D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief calls qhull library from fortran language to output facet table of delaunay triangulation 
!>        removes triangles with centroids outside of the PSLG 
!>        Organizes the vertices in ccw order.
!-----------------------------------------------------------------------
SUBROUTINE DelTriaWElim(PSLG, NP, POINTS, NF, FACETS)
!-----------------------------------------------------------------------
use iso_c_binding, only : c_f_pointer,C_PTR
implicit none

! INPUTS 
integer(kind=idx_t), intent(in) :: NP
real(kind=real_t),   intent(in),ALLOCATABLE :: POINTS(:,:)
integer(kind=idx_t), intent(inout):: NF
TYPE(BounDescrip2D) :: PSLG 

! OUTPUTS 
integer(kind=idx_t),intent(out),allocatable :: FACETS(:,:)

! LOCAL TO SUBROUTINE 
INTEGER(kind=idx_t) :: INoOUT 
type(c_ptr) :: cfacetemp
integer(kind=idx_t),pointer :: ffacetemp(:)
real(real_t),allocatable :: CENTROIDS(:,:)

integer :: i,j,k
integer :: dim 

DIM = PSLG%DIM 

! call qhull library 
cfacetemp = faces(DIM,NP,TRANSPOSE(POINTS),NF)

CALL C_F_POINTER(cfacetemp, ffacetemp, [NF*(DIM+1)])

! reshape vector ffacetemp to facets array
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

! organize vertices in ccw order
CALL WindingOrder2D(FACETS,NF,Points) 

!-----------------------------------------------------------------------
END SUBROUTINE DelTriaWElim
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! P R I V A T E  F U N C T I O N S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------



!> @brief Given an array of triangular facets, organize them in a 
!>        CCW order by calculating the cross-product and reversing 
!>        the vertex order in the facet if necessary.
!-----------------------------------------------------------------------
SUBROUTINE WindingOrder2D(FACETS,NF,POINTS) 
!-----------------------------------------------------------------------
implicit none 
INTEGER(idx_t),INTENT(INOUT),ALLOCATABLE :: FACETS(:,:) 
INTEGER(idx_t),INTENT(IN) :: NF
REAL(real_t),INTENT(IN),ALLOCATABLE :: POINTS(:,:)

REAL(real_t) :: AtoB(2),BtoC(2)
REAL(real_t) :: crossz 
integer(idx_t):: i
integer(idx_t):: nm1,nm2,nm3
integer(idx_t):: temp

do i = 1,nf
  nm1=facets(i,1) ! a
  nm2=facets(i,2) ! b
  nm3=facets(i,3) ! c
  AtoB = points(nm2,:) - points(nm1,:) ! b - a
  BtoC = points(nm3,:) - points(nm2,:) ! c - b
  crossz = AtoB(1)*BtoC(2) - AtoB(2)*BtoC(1) 
  ! then cw,reverse   
  if(crossz.lt.0.d0) then
    temp=facets(i,1)  
    facets(i,1)=facets(i,3) 
    facets(i,3)=temp
    !print *, "reversed"
  endif
enddo
!-----------------------------------------------------------------------
END SUBROUTINE WindingOrder2D
!-----------------------------------------------------------------------


!> @brief Given an array of int indices (i.e., a ROW), sort the ROW in ascending 
!>        lexiographic  order and return the mapping to achieve the sorted ROWS array.  
!>        assumes the ROWS are the first two columns of the array ROWS (i.e., ROWS(:,1:2)
!-----------------------------------------------------------------------
SUBROUTINE sortRows(NumRows,ROWS,IDX) 
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
END SUBROUTINE sortRows
!-----------------------------------------------------------------------


!> @brief Edit vertex spacing on PSLG to ensure it is less than MAXDIFF
!-----------------------------------------------------------------------
SUBROUTINE densify(MAXDIFF,PSLG) 
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
end subroutine densify 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Push array entries with zero back of array and return updated NP size 
!-----------------------------------------------------------------------
SUBROUTINE PushZerosToBackREAL(ARR,NP) 
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

!-----------------------------------------------------------------------
END SUBROUTINE PushZerosToBackREAL
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Push array entries with zero back of array and return updated NP size 
!-----------------------------------------------------------------------
SUBROUTINE PushZerosToBackINT(ARR,NP) 
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
END SUBROUTINE PushZerosToBackInt
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Given a PSLG and a set of points, determine the nearest 
!         Euclidean distance to the PSLG using a KD-tree and sign the distance 
!         negative if the point in question is inside the polygon defined by the PSLG. 
!         and vice-versa otherwise.
!-----------------------------------------------------------------------
SUBROUTINE CalcSignedDistance(IPTS,PSLG,SignedDistance)
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
END SUBROUTINE CalcSignedDistance
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Calls the c code pnpoly from cfunctions.c to determine if a point is in a 2D mulitply-
!>        connected polygon. 
!>        https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html#The%20C%20Code
!-----------------------------------------------------------------------
SUBROUTINE FPNPOLY(PSLG,TEST,INoOUT) 
use iso_c_binding, only : c_f_pointer,C_PTR
implicit none 

TYPE(BounDescrip2D), INTENT(IN) :: PSLG 
REAL(kind=real_t), INTENT(IN) :: TEST(2)
INTEGER(kind=idx_t),INTENT(OUT) :: INoOUT 

INoOUT = pnpoly(PSLG%NumVert0, PSLG%Vert0(:,1), PSLG%Vert0(:,2), &
                TEST(1),TEST(2))

!-----------------------------------------------------------------------
END SUBROUTINE 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief  Calculate the centroid of the triangles described in the arrays POINTS FACES
!>         For triangles: 
!>         centroids = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
!-----------------------------------------------------------------------
SUBROUTINE CalcBaryCenter(DIM,POINTS,NP,FACES,NF,CENTROIDS)
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

!-----------------------------------------------------------------------
END SUBROUTINE CalcBaryCenter
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Create a strutured grid of points using the vectors xgv and ygv
!>        this mimics what the MATLAB command does 
!>        [xg,yg]=meshgrid(xvec,yvec)
!-----------------------------------------------------------------------
SUBROUTINE meshgrid2D(xgv, ygv, X, Y)
!-----------------------------------------------------------------------
  implicit none
  ! INPUTS 
  real(kind=real_t),intent(in)   :: xgv(:), ygv(:)
  ! OUTPUTS
  real(kind=real_t),intent(out)  :: X(:,:), Y(:,:)

  X(:,:) = spread( xgv, 1, size(ygv) )
  Y(:,:) = spread( ygv, 2, size(xgv) )

!-----------------------------------------------------------------------
END SUBROUTINE
!-----------------------------------------------------------------------


!> @brief Creates a vector s that spans the interval (a,b) with n points. 
!> this mimics what the MATLAB command does 
!> REFERENCE: https://github.com/certik/fortran-utils
!> vec = linspace(a,b,n) 
!-----------------------------------------------------------------------
FUNCTION linspace(a, b, n) result(s)
!-----------------------------------------------------------------------
implicit none
! INPUTS 
real(kind=real_t), intent(in) :: a, b
integer(kind=idx_t), intent(in) :: n
! OUTPUTS 
real(kind=real_t) :: s(n)

s = meshexp(a, b, 1.0d0, n-1)

!-----------------------------------------------------------------------
END FUNCTION
!-----------------------------------------------------------------------


!> @brief REFERENCE: https://github.com/certik/fortran-utils
!> Generates exponential mesh of N elements on [rmin, rmax]
!>
!> Arguments
!> ---------
!>
!>The domain [rmin, rmax], the mesh will contain both enkind=real_toints:
!-----------------------------------------------------------------------
FUNCTION meshexp(rmin, rmax, a, N) result(mesh)
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
END FUNCTION
!-----------------------------------------------------------------------


!> @brief  quicksort.f -*-f90-*-
!> Author: t-nissie
!> License: GPLv3
!> Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!-----------------------------------------------------------------------
recursive subroutine quicksortINT(a, first, last)
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
end subroutine quicksortINT
!-----------------------------------------------------------------------


!> @brief quicksort.f -*-f90-*-
!> Author: t-nissie
!> License: GPLv3
!> Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!-----------------------------------------------------------------------
recursive subroutine quicksortINTLONG(a, first, last,idx)
!-----------------------------------------------------------------------
  implicit none
  integer(8)  a(*), x, t
  integer first, last,tt
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
     tt = idx(i);idx(i) = idx(j); idx(j) = tt ! save the sort map
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksortINTLONG(a, first, i-1,idx)
  if (j+1 < last)  call quicksortINTLONG(a, j+1, last,idx)
!-----------------------------------------------------------------------
end subroutine quicksortINTLONG
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief quicksort.f -*-f90-*-
!> Author: t-nissie
!> License: GPLv3
!> Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!-----------------------------------------------------------------------
recursive subroutine quicksortREAL(a, first, last)
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
end subroutine quicksortREAL
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief the median of a vector of reals 
!-----------------------------------------------------------------------
function median(a,first,last)
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
end function median
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
END MODULE UTILS
!-----------------------------------------------------------------------


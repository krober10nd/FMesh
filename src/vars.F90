!-----------------------------------------------------------------------
!  MODULE VARS 
!-----------------------------------------------------------------------
!> @brief variables and types for fmesh
!-----------------------------------------------------------------------
MODULE VARS 
!-----------------------------------------------------------------------
use iso_c_binding, only: c_int, c_int32_t, c_int64_t, c_float, c_double, c_ptr

implicit none

integer, parameter, public :: idx_t = c_int32_t
integer, parameter, public :: real_t = c_double

REAL(real_t) :: TSiter,TFiter ! < variables for timing each iteration of distmesh
REAL(real_t) :: TS,TF ! < variables for timing all iterations of distmesh

CHARACTER(100) :: pslgfname !< filename of PSLG data
CHARACTER(100) :: sizefname !< filename of mesh sizes data.
CHARACTER(100) :: elongfname1 !< filename of elongation data.
CHARACTER(100) :: elongfname2 !< filename of elongation data.
CHARACTER(100) :: anglefname !< filename of angle data.

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
integer(kind=idx_t) :: NumFlips ! counts the number of edge flips per iteration
integer(kind=idx_t) :: lastNumFlips ! in the last iteration, the number of edge flips
integer(kind=idx_t) :: ITFlips ! number of iterations spent flipping

REAL(kind=real_t)   :: BBOX(4) ! bounding box coords. Top left and bot. right. 

! some parameters for distmesh
REAL(kind=real_t) :: DELTAT ! important
REAL(kind=real_t),PARAMETER :: DPTOL =0.001d0
REAL(kind=real_t),PARAMETER :: TTOL  =0.1d0
REAL(kind=real_t),PARAMETER :: FSCALE=1.2d0 
REAL(kind=real_t) :: GEPS 
REAL(kind=real_t) :: DEPS 
REAL(kind=real_t),PARAMETER :: EPS=EPSILON(1.d0) 
INTEGER(kind=idx_t) :: iter,MaxIter,nscreen

!-----------------------------------------------------------------------
TYPE BounDescrip2D !< 2D boundary description as a winding polygon
   INTEGER(kind=idx_t) :: NumVert !< number of vertices along polygon
   INTEGER(kind=idx_t) :: NumVert0 !< original number of vertices along polygon
   INTEGER(kind=idx_t) :: DIM !< dimension of problem
   REAL(kind=real_t),ALLOCATABLE :: Vert(:,:)  !< densified poly
   REAL(kind=real_t),ALLOCATABLE :: Vert0(:,:) !< org. poly
ENDTYPE
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
TYPE GridData !< container for raster fields
    REAL(real_t),ALLOCATABLE :: vals(:,:)  !< discrete values of size function 
    REAL(real_t) ::   delta !< grid spacing in 
    INTEGER(idx_t) :: ni,nj  !< dimension of raster 
    REAL(real_t) :: x0y0(2) !< bottom left corner coordinate
    REAL(real_t) :: LMIN !< minimum element size in domain
END TYPE 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
TYPE MetricTensor !< container for all the raster fields
   TYPE(GridData) :: Iso !< isotropic sizes 
   TYPE(GridData) :: Angle !< angle of elongation 
   TYPE(GridData) :: Elong1 !< factor of elongation along principal
   TYPE(GridData) :: Elong2 !< factor of elongation along transverse
END TYPE 
!-----------------------------------------------------------------------

TYPE(GridData) :: SzFx !< for the size of elements in the domain 

TYPE(GridData) :: ElongFx1 !< for the elongation of elements in the domain 

TYPE(GridData) :: ElongFx2 !< for the elongation of elements in the domain 

TYPE(GridData) :: AngleFx !< for the orientation of elements in the domain 

TYPE(MetricTensor) :: SzFields !< container for all atributes of the metric tensor

TYPE(BounDescrip2D) :: PSLG !< read in from disk


PUBLIC LengthString 

!---------------------end of data declarations--------------------------------

CONTAINS 


!-----------------------------------------------------------------------
INTEGER(idx_t) FUNCTION LengthString(String) Result(StringLength) 
!-----------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(*),INTENT(IN) :: String
INTEGER(idx_t) :: I
StringLength = -1 !default is -1 just because.
DO I=LEN(String), 1, -1
  IF(String(I:I) .ne. ' ') THEN
     StringLength = I
     EXIT
  ENDIF
ENDDO
!-----------------------------------------------------------------------
END FUNCTION LengthString 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
END MODULE VARS 
!-----------------------------------------------------------------------

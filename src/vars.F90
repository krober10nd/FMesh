MODULE VARS 
use iso_c_binding, only: c_int, c_int32_t, c_int64_t, c_float, c_double, c_ptr

!**************************************************
! VARIABLES AND TYPES FOR FMesh
! 
! AUTHOR: Keith J. Roberts, 
! PLACE: Universdade de Sao Paulo
! START: 2019-08-17 
!**************************************************
IMPLICIT NONE 

!!!!!!
! only support 32-bit integers and 64-bit reals 
!!!!!!
integer, parameter, public :: idx_t = c_int32_t

integer, parameter, public :: real_t = c_double

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

integer(kind=idx_t) :: NUMFLIPS 

REAL(kind=real_t)   :: BBOX(4) ! bounding box coords. Top left and bot. right. 

! some parameters for distmesh
REAL(kind=real_t),PARAMETER :: DPTOL =0.001d0
REAL(kind=real_t),PARAMETER :: TTOL  =0.1d0
REAL(kind=real_t),PARAMETER :: FSCALE=1.2d0 
REAL(kind=real_t),PARAMETER :: DELTAT=0.10d0
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


END MODULE VARS 

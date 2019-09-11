!-----------------------------------------------------------------------
!  MODULE YourMeshSize
!-----------------------------------------------------------------------
!> @brief Sizing functions for meshes. 
!>        These sizing functions are bilinear gridded interpolants. 
!-----------------------------------------------------------------------
MODULE YourMeshSize
!-----------------------------------------------------------------------
USE vars 
implicit none

REAL(kind=real_t)   :: LMIN=0.10d0  !< minimum mesh size 
REAL(kind=real_t)   :: LMAX=0.10d0  !< maximum mesh size  
REAL(kind=real_t)   :: GRADE=0.35d0 !< mesh size gradation (in decimal percent) 
INTEGER(kind=idx_t) :: DIM=2        !< dimension of problem (2 or 3)

PUBLIC LoadMeshSizes,MeshSize,CalcMetricTensor
PRIVATE LinearInterp2D

!---------------------end of data declarations--------------------------------

contains


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! P U B L I C  F U N C T I O N S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Load the gridded data into a GridData type  
!-----------------------------------------------------------------------
FUNCTION LoadMeshSizes(fname) Result(SzFx) 
!-----------------------------------------------------------------------
implicit none 
type(GridData) :: SzFx
CHARACTER(80),intent(in) :: fname

integer(idx_t) :: LenFN 
logical        :: fileFound
integer(idx_t) :: nx,ny 
real(real_t)   :: dx 
integer(idx_t) :: i,j

! get length of the fname containing the pslg data 
LenFN = LengthString(fname)
IF(LenFN.lt. 1) THEN 
  write(*,'(A)') "FATAL: FNAME of file with mesh sizes is invalid." 
  stop
ENDIF !error out

! check if it exists 
INQUIRE(FILE=fname(1:LenFN),EXIST=fileFound)

IF(fileFound) THEN
  OPEN(UNIT=1,FILE=fname(1:LenFN),ACTION='READ')
  ! READ HEADER 
  READ(1,*) SzFx%Ni,SzFx%Nj,SzFx%delta,SzFx%x0y0(1),SzFx%x0y0(2)
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(3A,I5,A)') "INFO: Reading mesh sizes file called ",fname(1:LenFN)," with " &  
      //" ",SzFx%Ni*SzFx%Nj," points."

  WRITE(*,'(A)') "********************************************************"

  ALLOCATE(SzFx%Vals(SzFx%Ni,SzFx%Nj))
  SzFx%Vals=-9999.d0 
  DO J=1,SzFx%Nj
    DO I=1,SzFx%Ni 
      READ(1,*) SzFx%Vals(i,j)
    ENDDO
  ENDDO
ELSE
  WRITE(*,'(3A)') "FATAL: ",fname(1:LenFN) ," file not found."
  WRITE(*,'(A)') "********************************************************"
  STOP  
ENDIF

!-----------------------------------------------------------------------
END FUNCTION LoadMeshSizes 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Calls the bilinear interpolant to determine mesh sizes
!-----------------------------------------------------------------------
FUNCTION MeshSize(POINTS,SzFx) 
!-----------------------------------------------------------------------
real(real_t)            :: MeshSize
type(GridData),intent(in)    :: SzFx
real(real_t),intent(in) :: points(2)

MeshSize = LinearInterp2D(points,SzFx)

!-----------------------------------------------------------------------
END FUNCTION MeshSize
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief YOUR ORIENTATION AND ELONGATION CALCULATIONS HERE
!         IF ISOTROPIC SET TO IDENTITY MATRIX 
!-----------------------------------------------------------------------
FUNCTION CalcMetricTensor(POINTS) RESULT(ME)
!-----------------------------------------------------------------------
real(kind=real_t) :: ME(1:2,1:2) 
real(kind=real_t),intent(in) :: POINTS(2)
! just isotropic for now 
ME = 0.0d0 
ME(1,1) = 1.0d0
ME(2,2) = 1.0d0
! create a simple smooth variation in anisotropicness 
if(POINTS(1).GT.-2.0d0.AND.POINTS(1).LT.0.0d0) THEN 
  ME(2,2)= (2.0d0+points(1))*50.0d0+1.0d0
elseif(POINTS(1).GT.0.0d0.AND.POINTS(1).LT.2.0d0) THEN
  ME(2,2)= (2.0d0-points(1))*50.0d0+1.0d0
elseif(POINTS(1).LT.EPSILON(1.0d0).AND.POINTS(1).GT.-EPSILON(1.0d0)) THEN
  ME(2,2)=100.0d0
endif
ME(2,2)=MAXVAL((/ 1.0d0,ME(2,2)/))

!ME(2,2) = 10.d0
!-----------------------------------------------------------------------
END FUNCTION CalcMetricTensor 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! P R I V A T E  F U N C T I O N S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Given a gridded dataset defined on a structured grid 
!>        performs a bilinear interpolation.
!-----------------------------------------------------------------------
FUNCTION LinearInterp2D(Qp,MeshSizes) RESULT(InterpVal)
!-----------------------------------------------------------------------
IMPLICIT NONE 

REAL(real_t)  :: InterpVal
TYPE(GridData),INTENT(IN)  :: MeshSizes
REAL(real_t),INTENT(IN)  :: QP(2)

REAL(real_t)   :: interp_lat, interp_lon
REAL(real_t)   :: ri,rj,dlat,dlon,dataint,trueval
REAL(real_t)   :: wx,wy
INTEGER(idx_t) :: ileft,iright,jbot,jtop
INTEGER(idx_t) :: i,j

rj=(Qp(2)-MeshSizes%x0y0(2))/(MeshSizes%Ni+1)
ri=(Qp(1)-MeshSizes%x0y0(1))/(MeshSizes%Nj+1)

ileft=int(ri)
iright=ileft+1

jbot=int(rj)
jtop=jbot+1

wx=ri-int(ri)
wy=rj-int(rj)

! perform bilinear interpolation
InterpVal=(1.-wy)*(1.-wx)*MeshSizes%vals(ileft,jbot) + &
        (1.-wy)*wx*MeshSizes%vals(iright,jbot) + & 
        wy*(1.-wx)*MeshSizes%vals(ileft,jtop) + & 
        wy*wx*MeshSizes%vals(iright,jtop) 

!-----------------------------------------------------------------------
END FUNCTION LinearInterp2D
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
END MODULE YourMeshSize
!-----------------------------------------------------------------------

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

REAL(kind=real_t)   :: LMIN !< minimum mesh size (set when reading mesh size function from file) 


PUBLIC LoadMeshSizes,LoadMeshElong,LoadMeshAngle
PUBLIC CalcMeshSize,CalcMeshElong,CalcMeshAngle
PUBLIC CalcMetricTensor


PRIVATE LinearInterp2D,MatInv2

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
  WRITE(*,'(A)') "                                                        "
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(3A,I5,A,F12.8)') "INFO: Reading mesh sizes file called ",fname(1:LenFN)," with " &  
      //" ",SzFx%Ni*SzFx%Nj," points and grid spacing ",SzFx%delta

  ALLOCATE(SzFx%Vals(SzFx%Ni,SzFx%Nj))
  SzFx%Vals=-9999.d0 
  DO I=1,SzFx%Ni
    READ(1,*) SzFx%Vals(I,:)
  ENDDO
  
  SzFx%LMIN = MINVAL(SzFx%Vals) 
  WRITE(*,'(A,F12.8)') "INFO: The minimum element size is ",SzFx%LMIN 
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(A)') "                                                        "
ELSE
  WRITE(*,'(3A)') "FATAL: ",fname(1:LenFN) ," file not found."
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(A)') "                                                        "
  STOP  
ENDIF

!-----------------------------------------------------------------------
END FUNCTION LoadMeshSizes 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Load the gridded data into a GridData type  
!-----------------------------------------------------------------------
FUNCTION LoadMeshElong(fname) Result(SzFx) 
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
  write(*,'(A)') "FATAL: FNAME of file with mesh elongations is invalid." 
  stop
ENDIF !error out

! check if it exists 
INQUIRE(FILE=fname(1:LenFN),EXIST=fileFound)

IF(fileFound) THEN
  OPEN(UNIT=1,FILE=fname(1:LenFN),ACTION='READ')
  ! READ HEADER 
  READ(1,*) SzFx%Ni,SzFx%Nj,SzFx%delta,SzFx%x0y0(1),SzFx%x0y0(2)
  WRITE(*,'(A)') "                                                        "
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(3A,I5,A,F12.8)') "INFO: Reading elongation file called ",fname(1:LenFN)," with " &  
      //" ",SzFx%Ni*SzFx%Nj," points and grid spacing ",SzFx%delta

  ALLOCATE(SzFx%Vals(SzFx%Ni,SzFx%Nj))
  SzFx%Vals=-9999.d0 
  DO I=1,SzFx%Ni
    READ(1,*) SzFx%Vals(I,:)
  ENDDO
  
  SzFx%LMIN = MINVAL(SzFx%Vals) 
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(A)') "                                                        "
ELSE
  WRITE(*,'(3A)') "FATAL: ",fname(1:LenFN) ," file not found."
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(A)') "                                                        "
  STOP  
ENDIF

!-----------------------------------------------------------------------
END FUNCTION LoadMeshElong
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Load the gridded data into a GridData type  
!-----------------------------------------------------------------------
FUNCTION LoadMeshAngle(fname) Result(SzFx) 
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
  write(*,'(A)') "FATAL: FNAME of file with mesh angles is invalid." 
  stop
ENDIF !error out

! check if it exists 
INQUIRE(FILE=fname(1:LenFN),EXIST=fileFound)

IF(fileFound) THEN
  OPEN(UNIT=1,FILE=fname(1:LenFN),ACTION='READ')
  ! READ HEADER 
  READ(1,*) SzFx%Ni,SzFx%Nj,SzFx%delta,SzFx%x0y0(1),SzFx%x0y0(2)
  WRITE(*,'(A)') "                                                        "
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(3A,I5,A,F12.8)') "INFO: Reading angle file called ",fname(1:LenFN)," with " &  
      //" ",SzFx%Ni*SzFx%Nj," points and grid spacing ",SzFx%delta

  ALLOCATE(SzFx%Vals(SzFx%Ni,SzFx%Nj))
  SzFx%Vals=-9999.d0 
  DO I=1,SzFx%Ni
    READ(1,*) SzFx%Vals(I,:)
  ENDDO
  
  SzFx%LMIN = MINVAL(SzFx%Vals) 
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(A)') "                                                        "
ELSE
  WRITE(*,'(3A)') "FATAL: ",fname(1:LenFN) ," file not found."
  WRITE(*,'(A)') "********************************************************"
  WRITE(*,'(A)') "                                                        "
  STOP  
ENDIF

!-----------------------------------------------------------------------
END FUNCTION LoadMeshAngle
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Calls the bilinear interpolant to determine mesh sizes
!-----------------------------------------------------------------------
FUNCTION CalcMeshSize(POINTS,SzFx) Result(MeshSize)
!-----------------------------------------------------------------------
real(real_t):: MeshSize
type(MetricTensor),intent(in):: SzFx
real(real_t),intent(in) :: points(2)

MeshSize = LinearInterp2D(points,SzFx%Iso)

!-----------------------------------------------------------------------
END FUNCTION CalcMeshSize
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Calls the bilinear interpolant to determine mesh elongation
!-----------------------------------------------------------------------
FUNCTION CalcMeshElongation(POINTS,SzFx) Result(MeshElong)
!-----------------------------------------------------------------------
real(real_t):: MeshElong
type(MetricTensor),intent(in):: SzFx
real(real_t),intent(in) :: points(2)

MeshElong = LinearInterp2D(points,SzFx%Elong)

!-----------------------------------------------------------------------
END FUNCTION CalcMeshElongation
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Calls the bilinear interpolant to determine mesh angle
!-----------------------------------------------------------------------
FUNCTION CalcMeshAngle(POINTS,SzFx) Result(MeshAngle)
!-----------------------------------------------------------------------
real(real_t):: MeshAngle
type(MetricTensor),intent(in):: SzFx
real(real_t),intent(in) :: points(2)

MeshAngle = LinearInterp2D(points,SzFx%Angle)

!-----------------------------------------------------------------------
END FUNCTION CalcMeshAngle
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief YOUR ORIENTATION AND ELONGATION CALCULATIONS HERE
!         IF ISOTROPIC SET TO IDENTITY MATRIX 
!-----------------------------------------------------------------------
FUNCTION CalcMetricTensor(POINTS,SzFx) RESULT(ME)
!-----------------------------------------------------------------------
real(kind=real_t) :: ME(1:2,1:2) 
real(kind=real_t),intent(in) :: points(2)
type(MetricTensor),intent(in) :: SzFx

real(real_t) :: MeshAngle,MeshElong
real(real_t) :: rot(1:2,1:2),rot_inv(1:2,1:2),elong(1:2,1:2)
real(real_t) :: temp1(1:2,1:2)

MeshAngle = CalcMeshAngle(points,SzFx) 
MeshElong = CalcMeshElong(points,SzFx) 

A = COS(MeshAngle)*(180.0d0/3.14d0)
B = SIN(MeshAngle)*(180.0d0/3.14d0

rot(1,1) =  A
rot(2,2) =  A
rot(1,2) =  B
rot(2,1) = -B
 
rot_inv = MatInv(Rot) ; 

elong(1:2,1:2) = 0.d0  ! off-diag
elong(1,1) = MeshElong ! only principal direc
elong(2,2) = 1.0d0 

temp1 = MatMul(rot,elong)
ME    = MatMul(temp1,rot_inv) 

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

rj=(Qp(2)-MeshSizes%x0y0(2))/MeshSizes%delta+1
ri=(Qp(1)-MeshSizes%x0y0(1))/MeshSizes%delta+1

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
!> @brief Performs a direct calculation of the inverse of a 2×2 matrix.
!-----------------------------------------------------------------------
  pure function matinv2(A) result(B)
!-----------------------------------------------------------------------
    complex(wp), intent(in) :: A(2,2)   !! Matrix
    complex(wp)             :: B(2,2)   !! Inverse matrix
    complex(wp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
!-----------------------------------------------------------------------
  end function matinv2
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
END MODULE YourMeshSize
!-----------------------------------------------------------------------

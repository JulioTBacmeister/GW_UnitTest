module ncread_mod

!
! This module 
!
use shr_kind_mod, only: r8 => shr_kind_r8
use physconst, only: rearth
use netcdf


implicit none
private


! Public interface
public :: ncread_topo
public :: ncread_topo_latlon
public :: ncread_era5_se_ic
public :: ncread_camsnap
public :: ncread_xy_data


!==========================================================================
contains
!==========================================================================

subroutine ncread_topo(ncfile, mxdis, angll, aniso, anixy, hwdth, clngt, gbxar, isovar, isowgt, sgh )

  character(len=*), intent(in) :: ncfile
  real(r8) , allocatable, intent(out) :: mxdis(:,:), angll(:,:), aniso(:,:), anixy(:,:), hwdth(:,:), clngt(:,:)
  real(r8) , allocatable, intent(out) :: gbxar(:), isowgt(:), isovar(:), sgh(:)
  
  integer ncid,status, dimid,varid  ! for netCDF data file

  integer :: ncol, nrdgs
  logical :: llatlon_x

  status = nf90_open( trim(ncfile) , 0, ncid)
  write(*,*) "STATUS" ,status,NF90_NOERR
  IF (STATUS .NE. NF90_NOERR) STOP "1" ! CALL HANDLE_ERR(STATUS)

  status = NF90_INQ_DIMID(ncid, 'ncol', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=ncol)
  IF (status .NE. NF90_NOERR) STOP "2" ! CALL HANDLE_ERR(status)

  status = NF90_INQ_DIMID(ncid, 'nrdg', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=nrdgs)
  IF (status .NE. NF90_NOERR) STOP "3" ! CALL HANDLE_ERR(status)

  allocate(  mxdis(ncol,nrdgs),  angll(ncol,nrdgs),  aniso(ncol,nrdgs),  anixy(ncol,nrdgs), &
       hwdth(ncol,nrdgs),  clngt(ncol,nrdgs)  )
  allocate(  gbxar(ncol),   isowgt(ncol),   isovar(ncol),  sgh(ncol)  )
  
  status = NF90_INQ_VARID(ncid, 'MXDIS', varid)
  status = NF90_GET_VAR(ncid, varid, mxdis )
  write(*,*) " MXDIS ", minval(mxdis ),maxval(mxdis )
  
  status = NF90_INQ_VARID(ncid, 'ANGLL', varid)
  status = NF90_GET_VAR(ncid, varid, angll )
  
  status = NF90_INQ_VARID(ncid, 'ANISO', varid)
  status = NF90_GET_VAR(ncid, varid, aniso )
  
  status = NF90_INQ_VARID(ncid, 'ANIXY', varid)
  status = NF90_GET_VAR(ncid, varid, anixy )
  
  status = NF90_INQ_VARID(ncid, 'HWDTH', varid)
  status = NF90_GET_VAR(ncid, varid, hwdth )

  status = NF90_INQ_VARID(ncid, 'CLNGT', varid)
  status = NF90_GET_VAR(ncid, varid, clngt )

  status = NF90_INQ_VARID(ncid, 'ISOWGT', varid)
  status = NF90_GET_VAR(ncid, varid, isowgt )

  status = NF90_INQ_VARID(ncid, 'ISOVAR', varid)
  status = NF90_GET_VAR(ncid, varid, isovar )

  status = NF90_INQ_VARID(ncid, 'GBXAR', varid)
  status = NF90_GET_VAR(ncid, varid, gbxar )

  status = NF90_INQ_VARID(ncid, 'SGH', varid)
  status = NF90_GET_VAR(ncid, varid, sgh )
  write(*,*) " SGH   ", minval(sgh ),maxval(sgh )

  ! Scale gbxar as it is done in CAM. This needs to checked out. (3/1/25)
  gbxar = gbxar * (rearth/1000._r8)*(rearth/1000._r8) ! transform to km^2


end subroutine ncread_topo

!====================================================================
subroutine ncread_topo_latlon(ncfile, mxdis, angll, aniso, anixy, hwdth, clngt, gbxar, isovar, isowgt, sgh )

  character(len=*), intent(in) :: ncfile
  real(r8) , allocatable, intent(out) :: mxdis(:,:), angll(:,:), aniso(:,:), anixy(:,:), hwdth(:,:), clngt(:,:)
  real(r8) , allocatable, intent(out) :: gbxar(:), isowgt(:), isovar(:), sgh(:)

  ! local varaibles
  real(r8) , allocatable :: mxdisR(:,:,:), angllR(:,:,:), anisoR(:,:,:)
  real(r8) , allocatable :: anixyR(:,:,:), hwdthR(:,:,:), clngtR(:,:,:)
  real(r8) , allocatable :: gbxarR(:,:), isowgtR(:,:), isovarR(:,:), sghR(:,:)
  
  integer ncid,status, dimid,varid  ! for netCDF data file

  integer :: ncol, nrdgs, nx, ny
  logical :: llatlon_x

  status = nf90_open( trim(ncfile) , 0, ncid)
  write(*,*) "STATUS" ,status,NF90_NOERR
  IF (STATUS .NE. NF90_NOERR) STOP "1" ! CALL HANDLE_ERR(STATUS)

  status = NF90_INQ_DIMID(ncid, 'lon', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=nx)
  status = NF90_INQ_DIMID(ncid, 'lat', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=ny)
  ncol=nx*ny
  write(*,*) "nx,ny,ncol",nx,ny,ncol
  IF (status .NE. NF90_NOERR) STOP "2" ! CALL HANDLE_ERR(status)

  status = NF90_INQ_DIMID(ncid, 'nrdg', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=nrdgs)
  IF (status .NE. NF90_NOERR) STOP "3" ! CALL HANDLE_ERR(status)

  allocate(  mxdis(ncol,nrdgs),  angll(ncol,nrdgs),  aniso(ncol,nrdgs),  anixy(ncol,nrdgs), &
       hwdth(ncol,nrdgs),  clngt(ncol,nrdgs)  )
  allocate(  gbxar(ncol),   isowgt(ncol),   isovar(ncol),  sgh(ncol)  )

  allocate(  mxdisR(nx,ny,nrdgs),  angllR(nx,ny,nrdgs),  anisoR(nx,ny,nrdgs),  anixyR(nx,ny,nrdgs), &
       hwdthR(nx,ny,nrdgs),  clngtR(nx,ny,nrdgs)  )
  allocate(  gbxarR(nx,ny),   isowgtR(nx,ny),   isovarR(nx,ny),  sghR(nx,ny)  )
  
  status = NF90_INQ_VARID(ncid, 'MXDIS', varid)
  status = NF90_GET_VAR(ncid, varid, mxdisR )
  mxdis =reshape( mxdisR , [ncol,nrdgs] )
  write(*,*) " MXDIS(R) ", minval(mxdisR ),maxval(mxdisR )


  
  
  status = NF90_INQ_VARID(ncid, 'ANGLL', varid)
  status = NF90_GET_VAR(ncid, varid, angllR )
  mxdis =reshape( mxdisR , [ncol,nrdgs] )
  
  status = NF90_INQ_VARID(ncid, 'ANISO', varid)
  status = NF90_GET_VAR(ncid, varid, anisoR )
  aniso =reshape( anisoR , [ncol,nrdgs] )
  
  status = NF90_INQ_VARID(ncid, 'ANIXY', varid)
  status = NF90_GET_VAR(ncid, varid, anixyR )
  anixy =reshape( anixyR , [ncol,nrdgs] )
  
  status = NF90_INQ_VARID(ncid, 'HWDTH', varid)
  status = NF90_GET_VAR(ncid, varid, hwdthR )
  hwdth =reshape( hwdthR , [ncol,nrdgs] )

  status = NF90_INQ_VARID(ncid, 'CLNGT', varid)
  status = NF90_GET_VAR(ncid, varid, clngtR )
  clngt =reshape( clngtR , [ncol,nrdgs] )

  status = NF90_INQ_VARID(ncid, 'ISOWGT', varid)
  status = NF90_GET_VAR(ncid, varid, isowgtR )
  isowgt =reshape( isowgtR , [ncol] )

  status = NF90_INQ_VARID(ncid, 'ISOVAR', varid)
  status = NF90_GET_VAR(ncid, varid, isovarR )
  isovar =reshape( isovarR , [ncol] )

  status = NF90_INQ_VARID(ncid, 'GBXAR', varid)
  status = NF90_GET_VAR(ncid, varid, gbxarR )
  gbxar =reshape( gbxarR , [ncol] )

  status = NF90_INQ_VARID(ncid, 'SGH', varid)
  status = NF90_GET_VAR(ncid, varid, sghR )
  sgh =reshape( sghR , [ncol] )
  
  write(*,*) " SGH(col)   ", minval(sgh ),maxval(sgh )
  write(*,*) " MXDIS(col) ", minval(mxdis ),maxval(mxdis )

  ! Scale gbxar as it is done in CAM. This needs to checked out. (3/1/25)
  gbxar = gbxar * (rearth/1000._r8)*(rearth/1000._r8) ! transform to km^2


end subroutine ncread_topo_latlon

!====================================================================
subroutine ncread_xy_data(ncfile , ncol, pver, ntim, &
     hyai , hybi , hyam , hybm , lon , lat , &
     PS, U , V , T , Q , ZETA ) !! , Q )

  use sphere_ops, only: sphere_curl2_xy

  character(len=*), intent(in) :: ncfile
  
  integer, intent(out) :: ncol,pver,ntim
  
  real(r8) , allocatable, intent(out) :: hyai(:), hybi(:), hyam(:), hybm(:)
  real(r8) , allocatable, intent(out) :: lon(:), lat(:)
  real(r8) , allocatable, intent(out) :: U(:,:,:), V(:,:,:), T(:,:,:), PS(:,:)
  real(r8) , allocatable, intent(out) :: Q(:,:,:), ZETA(:,:,:)
  
  !----- local variables --------------------------
  real(r8) , allocatable :: lon_R(:), lat_R(:)
  real(r8) , allocatable :: lon_xy(:,:), lat_xy(:,:)
  real(r8) , allocatable :: U_R(:,:,:,:), V_R(:,:,:,:), T_R(:,:,:,:), PS_R(:,:,:)
  real(r8) , allocatable :: ZETA_tk(:,:)

  integer :: i,j,k,icol,itim,nx,ny

  integer :: ncid,status, dimid,varid  ! for netCDF data file
  
  status = nf90_open( trim(ncfile) , 0, ncid)
  IF (STATUS .NE. NF90_NOERR) STOP "cant open "// trim(ncfile) 
  
  status = NF90_INQ_DIMID(ncid, 'lev', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=pver)
  
  status = NF90_INQ_DIMID(ncid, 'lon', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=nx)

  status = NF90_INQ_DIMID(ncid, 'lat', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=ny)
  
  status = NF90_INQ_DIMID(ncid, 'time', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=ntim)

  write(*,*) ncfile
  
  write(*,*) nx,ny,pver,ntim
  
  allocate(  U_R(nx,ny,pver,ntim) , V_R(nx,ny,pver,ntim) , T_R(nx,ny,pver,ntim)  , PS_R(nx,ny,ntim) )
  allocate(  lat_R(ny) , lon_R(nx)  )
  allocate(  lat_xy(nx,ny) , lon_xy(nx,ny), ZETA_tk(nx,ny)  )
  allocate(  hyai(pver+1) , hybi(pver+1) , hyam(pver) , hybm(pver)  )

  ncol=nx*ny
  allocate(  U(ncol,pver,ntim) , V(ncol,pver,ntim) , T(ncol,pver,ntim)  , PS(ncol,ntim) )
  allocate(  Q(ncol,pver,ntim), ZETA(ncol,pver,ntim) )
  allocate(  lat(ncol) , lon(ncol)  )

  status = NF90_INQ_VARID(ncid, 'lat', varid)
  status = NF90_GET_VAR(ncid, varid, lat_R )
  write(*,*) "Got ","lat ",minval(lat_R), maxval(lat_R)
  
  status = NF90_INQ_VARID(ncid, 'lon', varid)
  status = NF90_GET_VAR(ncid, varid, lon_R )
  
  status = NF90_INQ_VARID(ncid, 'hyai', varid)
  status = NF90_GET_VAR(ncid, varid, hyai )
  
  status = NF90_INQ_VARID(ncid, 'hybi', varid)
  status = NF90_GET_VAR(ncid, varid, hybi )
  
  status = NF90_INQ_VARID(ncid, 'hyam', varid)
  status = NF90_GET_VAR(ncid, varid, hyam )
  
  status = NF90_INQ_VARID(ncid, 'hybm', varid)
  status = NF90_GET_VAR(ncid, varid, hybm )
  
  status = NF90_INQ_VARID(ncid, 'PS', varid)
  status = NF90_GET_VAR(ncid, varid, PS_R )
  PS =reshape( PS_R , [ncol,ntim] )
  write(*,*) "Got ","PS  ",minval(PS), maxval(PS)
  
  status = NF90_INQ_VARID(ncid, 'U', varid)
  status = NF90_GET_VAR(ncid, varid, U_R )
  U =reshape( U_R , [ncol,pver,ntim] )
  write(*,*) "Got ","U   ",minval(U), maxval(U)
  
  status = NF90_INQ_VARID(ncid, 'V', varid)
  status = NF90_GET_VAR(ncid, varid, V_R )
  V =reshape( V_R , [ncol,pver,ntim] )
  write(*,*) "Got ","V   ",minval(V), maxval(V)
  
  status = NF90_INQ_VARID(ncid, 'T', varid)
  status = NF90_GET_VAR(ncid, varid, T_R )
  T =reshape( T_R , [ncol,pver,ntim] )
  write(*,*) "Got ","T   ",minval(T), maxval(T)
  
  !status = NF90_INQ_VARID(ncid, 'Q', varid)
  !status = NF90_GET_VAR(ncid, varid, Q )
  Q(:,:,:) = 0._r8
  write(*,*) "Got ","Q   ",minval(Q), maxval(Q)
  
  do j=1,ny
     do i=1,nx
        lat_xy(i,j)=lat_R(j)
        lon_xy(i,j)=lon_R(i)
     end do
  end do
  icol=0
  do j=1,ny
     do i=1,nx
        icol=icol+1
        lat(icol) = lat_xy(i,j)
        lon(icol) = lon_xy(i,j)
     end do
  end do

#if 1
  do itim=1,ntim
     do k=1,pver
        call sphere_curl2_xy(U_R(:,:,k,itim) , V_R(:,:,k,itim) , lat_R, lon_R, ZETA_tk, .TRUE.)
        ZETA(:,k,itim) = reshape( ZETA_tk , [ncol] )
     end do
  end do
  write(*,*) "Got ","ZETA ",minval(ZETA), maxval(ZETA)
#endif  
  
end subroutine ncread_xy_data

!==========================================================================

subroutine ncread_era5_se_ic(ncfile , ncol, pver, ntim, &
     hyai , hybi , hyam , hybm , lon , lat , &
     PS, U , V , T  , Q )

  character(len=*), intent(in) :: ncfile
  
  integer, intent(out) :: ncol,pver,ntim
  
  real(r8) , allocatable, intent(out) :: hyai(:), hybi(:), hyam(:), hybm(:)
  real(r8) , allocatable, intent(out) :: lon(:), lat(:)
  real(r8) , allocatable, intent(out) :: U(:,:,:), V(:,:,:), T(:,:,:), Q(:,:,:), PS(:,:)

  integer ncid,status, dimid,varid  ! for netCDF data file
  
  status = nf90_open( trim(ncfile) , 0, ncid)
  IF (STATUS .NE. NF90_NOERR) STOP "cant open "// trim(ncfile) 
  
  status = NF90_INQ_DIMID(ncid, 'lev', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=pver)
  
  status = NF90_INQ_DIMID(ncid, 'ncol', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=ncol)
  
  status = NF90_INQ_DIMID(ncid, 'time', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=ntim)

  write(*,*) ncfile
  
  write(*,*) ntim,ncol,pver

  
  allocate(  U(ncol,pver,ntim) , V(ncol,pver,ntim) , T(ncol,pver,ntim) , Q(ncol,pver,ntim) , PS(ncol,ntim) )
  allocate(  lat(ncol) , lon(ncol)  )
  allocate(  hyai(pver+1) , hybi(pver+1) , hyam(pver) , hybm(pver)  )

  status = NF90_INQ_VARID(ncid, 'lat', varid)
  status = NF90_GET_VAR(ncid, varid, lat )
  write(*,*) "Got ","lat ",minval(lat), maxval(lat)
  
  status = NF90_INQ_VARID(ncid, 'lon', varid)
  status = NF90_GET_VAR(ncid, varid, lon )
  
  status = NF90_INQ_VARID(ncid, 'hyai', varid)
  status = NF90_GET_VAR(ncid, varid, hyai )
  
  status = NF90_INQ_VARID(ncid, 'hybi', varid)
  status = NF90_GET_VAR(ncid, varid, hybi )
  
  status = NF90_INQ_VARID(ncid, 'hyam', varid)
  status = NF90_GET_VAR(ncid, varid, hyam )
  
  status = NF90_INQ_VARID(ncid, 'hybm', varid)
  status = NF90_GET_VAR(ncid, varid, hybm )
  
  status = NF90_INQ_VARID(ncid, 'PS', varid)
  status = NF90_GET_VAR(ncid, varid, PS )
  write(*,*) "Got ","PS  ",minval(PS), maxval(PS)
  
  status = NF90_INQ_VARID(ncid, 'U', varid)
  status = NF90_GET_VAR(ncid, varid, U )
  write(*,*) "Got ","U   ",minval(U), maxval(U)
  
  status = NF90_INQ_VARID(ncid, 'V', varid)
  status = NF90_GET_VAR(ncid, varid, V )
  write(*,*) "Got ","V   ",minval(V), maxval(V)
  
  status = NF90_INQ_VARID(ncid, 'T', varid)
  status = NF90_GET_VAR(ncid, varid, T )
  write(*,*) "Got ","T   ",minval(T), maxval(T)
  
  status = NF90_INQ_VARID(ncid, 'Q', varid)
  status = NF90_GET_VAR(ncid, varid, Q )
  write(*,*) "Got ","Q   ",minval(Q), maxval(Q)
  
  

  
end subroutine ncread_era5_se_ic

!==========================================================================

subroutine ncread_camsnap(ncfile , ncol, pver, ntim, &
     hyai , hybi , hyam , hybm , lon , lat , &
     PS, U , V , T  , Q , ZETA , &
     zm, zi, nm, ni, rhoi, pint, piln )

  character(len=*), intent(in) :: ncfile
  
  integer, intent(out) :: ncol,pver,ntim
  
  real(r8) , allocatable, intent(out) :: hyai(:), hybi(:), hyam(:), hybm(:)
  real(r8) , allocatable, intent(out) :: lon(:), lat(:)
  real(r8) , allocatable, intent(out) :: U(:,:,:), V(:,:,:), T(:,:,:), Q(:,:,:), PS(:,:), ZETA(:,:,:)
  real(r8) , allocatable, intent(out) :: zm(:,:,:), zi(:,:,:), nm(:,:,:), ni(:,:,:), rhoi(:,:,:)
  real(r8) , allocatable, intent(out) :: pint(:,:,:), piln(:,:,:)

  integer ncid,status, dimid,varid  ! for netCDF data file
  
  status = nf90_open( trim(ncfile) , 0, ncid)
  IF (STATUS .NE. NF90_NOERR) STOP "cant open "// trim(ncfile) 
  
  status = NF90_INQ_DIMID(ncid, 'lev', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=pver)
  
  status = NF90_INQ_DIMID(ncid, 'ncol', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=ncol)
  
  status = NF90_INQ_DIMID(ncid, 'time', dimid)
  status = NF90_INQUIRE_DIMENSION(ncid, dimid, len=ntim)

  write(*,*) ncfile
  
  write(*,*) ntim,ncol,pver

  
  allocate(  U(ncol,pver,ntim) , V(ncol,pver,ntim) , T(ncol,pver,ntim) , Q(ncol,pver,ntim) , PS(ncol,ntim) )
  allocate(  ZETA(ncol,pver,ntim)  )
  allocate(  zm(ncol,pver,ntim) , zi(ncol,pver+1,ntim) , nm(ncol,pver,ntim) , ni(ncol,pver+1,ntim)  )
  allocate(  pint(ncol,pver+1,ntim) , piln(ncol,pver+1,ntim)    )
  allocate(  rhoi(ncol,pver+1,ntim)   )
  allocate(  lat(ncol) , lon(ncol)  )
  allocate(  hyai(pver+1) , hybi(pver+1) , hyam(pver) , hybm(pver)  )

  status = NF90_INQ_VARID(ncid, 'lat', varid)
  status = NF90_GET_VAR(ncid, varid, lat )
  write(*,*) "Got ","lat ",minval(lat), maxval(lat)
  
  status = NF90_INQ_VARID(ncid, 'lon', varid)
  status = NF90_GET_VAR(ncid, varid, lon )
  
  status = NF90_INQ_VARID(ncid, 'hyai', varid)
  status = NF90_GET_VAR(ncid, varid, hyai )
  
  status = NF90_INQ_VARID(ncid, 'hybi', varid)
  status = NF90_GET_VAR(ncid, varid, hybi )
  
  status = NF90_INQ_VARID(ncid, 'hyam', varid)
  status = NF90_GET_VAR(ncid, varid, hyam )
  
  status = NF90_INQ_VARID(ncid, 'hybm', varid)
  status = NF90_GET_VAR(ncid, varid, hybm )
  
  status = NF90_INQ_VARID(ncid, 'PS', varid)
  status = NF90_GET_VAR(ncid, varid, PS )
  write(*,*) "Got ","PS  ",minval(PS), maxval(PS)
  
  status = NF90_INQ_VARID(ncid, 'UEGW', varid)
  status = NF90_GET_VAR(ncid, varid, U )
  write(*,*) "Got ","U   ",minval(U), maxval(U)
  
  status = NF90_INQ_VARID(ncid, 'VEGW', varid)
  status = NF90_GET_VAR(ncid, varid, V )
  write(*,*) "Got ","V   ",minval(V), maxval(V)
  
  status = NF90_INQ_VARID(ncid, 'TEGW', varid)
  status = NF90_GET_VAR(ncid, varid, T )
  write(*,*) "Got ","T   ",minval(T), maxval(T)

  status = NF90_INQ_VARID(ncid, 'QEGW', varid)
  status = NF90_GET_VAR(ncid, varid, Q )
  write(*,*) "Got ","Q   ",minval(Q), maxval(Q)

  status = NF90_INQ_VARID(ncid, 'ZIEGW', varid)
  status = NF90_GET_VAR(ncid, varid, zi )
  write(*,*) "Got ","zi   ",minval(zi), maxval(zi)

  status = NF90_INQ_VARID(ncid, 'ZMEGW', varid)
  status = NF90_GET_VAR(ncid, varid, zm )
  write(*,*) "Got ","zm   ",minval(zm), maxval(zm)
  
  status = NF90_INQ_VARID(ncid, 'NMEGW', varid)
  status = NF90_GET_VAR(ncid, varid, nm )
  write(*,*) "Got ","nm   ",minval(nm), maxval(nm)
  
  status = NF90_INQ_VARID(ncid, 'NIEGW', varid)
  status = NF90_GET_VAR(ncid, varid, ni )
  write(*,*) "Got ","ni   ",minval(ni), maxval(ni)

  status = NF90_INQ_VARID(ncid, 'RHOIEGW', varid)
  status = NF90_GET_VAR(ncid, varid, rhoi )
  write(*,*) "Got ","rhoi   ",minval(rhoi), maxval(rhoi)

  status = NF90_INQ_VARID(ncid, 'PINTEGW', varid)
  status = NF90_GET_VAR(ncid, varid, pint )
  write(*,*) "Got ","pint   ",minval(pint), maxval(pint)

  status = NF90_INQ_VARID(ncid, 'PILNEGW', varid)
  status = NF90_GET_VAR(ncid, varid, piln )
  write(*,*) "Got ","piln   ",minval(piln), maxval(piln)

  status = NF90_INQ_VARID(ncid, 'VORT4GW', varid)
  IF (STATUS .NE. NF90_NOERR) then
     write(*,*) "cant find VORT in "// trim(ncfile)
     zeta = U*0.
  else
     status = NF90_GET_VAR(ncid, varid, zeta  )
  end IF
  write(*,*) "Range of ","ZETA   ",minval(zeta), maxval(zeta)
   
end subroutine ncread_camsnap

end module ncread_mod

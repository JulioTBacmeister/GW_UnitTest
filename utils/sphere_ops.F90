module sphere_ops
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: Re => rearth

  implicit none
  private
  public :: sphere_curl2_xy
  public :: sphere_fronto_wrapper

  !!integer, parameter :: dp = selected_real_kind(15, 300)
contains

  !--------------------------------------------------------------------
  !> Convenience wrapper: inputs in (nx,ny) order → internally transpose
  !!   fx_xy(nx,ny), fy_xy(nx,ny)
  !!   lat(ny), lon(nx)
  !!   curlz_xy(nx,ny)  (result transposed back)
  !--------------------------------------------------------------------
  subroutine sphere_curl2_xy(fx_xy, fy_xy, lat, lon, curlz_xy, wrap)
    real(r8), intent(in)  :: fx_xy(:,:), fy_xy(:,:)          ! (nx,ny)
    real(r8), intent(in)  :: lat(:), lon(:)                  ! lat(ny), lon(nx)
    real(r8), intent(out) :: curlz_xy(size(fx_xy,1), size(fx_xy,2))  ! (nx,ny)
    !logical,  intent(in), optional :: wrap
    logical,  intent(in) :: wrap

    integer :: nx, ny
    real(r8), allocatable :: fx_yx(:,:), fy_yx(:,:), curlz_yx(:,:)

    nx = size(fx_xy,1)
    ny = size(fx_xy,2)

    if (size(fy_xy,1) /= nx .or. size(fy_xy,2) /= ny) stop "sphere_curl2_xy: size mismatch"
    if (size(lat) /= ny) stop "sphere_curl2_xy: lat size mismatch"
    if (size(lon) /= nx) stop "sphere_curl2_xy: lon size mismatch"

    allocate(fx_yx(ny,nx), fy_yx(ny,nx), curlz_yx(ny,nx))

    ! (nx,ny) → (ny,nx)
    fx_yx = transpose(fx_xy)
    fy_yx = transpose(fy_xy)

    ! Call core that expects (ny,nx)
    call sphere_curl2(fx_yx, fy_yx, lat, lon, curlz_yx, wrap)

    ! Back to (nx,ny)
    curlz_xy = transpose(curlz_yx)

  end subroutine sphere_curl2_xy

  !--------------------------------------------------------------------
  !> Curl (vertical component) of a horizontal vector field on a sphere
  !! Inputs:
  !!   fx(ny,nx), fy(ny,nx)  : vector components on lat-lon grid
  !!   lat(ny), lon(nx)      : degrees
  !!   Re                    : Earth radius (same units as desired output denom)
  !!   wrap (optional)       : if .true. (default), periodic in longitude
  !! Output:
  !!   curlz(ny,nx)
  !!
  !! Discretization matches the Python:
  !!   ζ ≈ (1/(Re*cosφ)) [ d(fy)/dlon - ( (cosφ_{j+1} fx_{j+1} - cosφ_{j-1} fx_{j-1})/dlat ) ]
  !! using centered differences and cyclic shifts.
  !--------------------------------------------------------------------
  subroutine sphere_curl2(fx, fy, lat, lon, curlz, wrap)
    real(r8), intent(in)  :: fx(:,:), fy(:,:)
    real(r8), intent(in)  :: lat(:), lon(:)
    real(r8), intent(out) :: curlz(size(fx,1), size(fx,2))
    logical,  intent(in), optional :: wrap
    integer :: j, i
    real(r8) :: dfdx, termy


    integer :: ny, nx
    real(r8), allocatable :: rlat(:), rlon(:), coslat(:)
    real(r8) :: dlat, dlon, deg2rad
    real(r8), allocatable :: coslat_jn(:), coslat_js(:)
    real(r8), allocatable :: coslat_2d(:,:), coslat_jn_2d(:), coslat_js_2d(:)
    logical :: lwrap

    ny = size(lat)
    nx = size(lon)
    lwrap = .true.; if (present(wrap)) lwrap = wrap

    deg2rad = acos(-1.0_r8) / 180.0_r8

    allocate(rlat(ny), rlon(nx), coslat(ny), coslat_jn(ny), coslat_js(ny))
    rlat = deg2rad * lat
    rlon = deg2rad * lon

    ! Uniform-step assumption, matching your Python (uses [3]-[1])
    if (ny >= 3) then
      dlat = rlat(3) - rlat(1)
    else
      dlat = 0.0_r8
    end if
    if (nx >= 3) then
      dlon = rlon(3) - rlon(1)
    else
      dlon = 0.0_r8
    end if

    coslat = cos(rlat)

    ! North/South shifted cos(lat): j+1 and j-1 (cyclic for simplicity)
    coslat_jn = cshift(coslat, shift = +1)   ! j+1
    coslat_js = cshift(coslat, shift = -1)   ! j-1

    ! Build 2D spreads for broadcasting-like ops
    allocate(coslat_2d(ny, nx))
    coslat_2d = spread(coslat, dim=2, ncopies=nx)
    ! (For the terms that multiply fx at j±1 we just spread the 1D factors on the fly.)

    ! Guard against zero spacings (degenerate grids)
    if (dlat == 0.0_r8 .or. dlon == 0.0_r8) then
      curlz = 0.0_r8
      ! Apply the same boundary zeroing as in Python
      if (ny >= 1) curlz(1, :) = 0.0_r8
      if (ny >= 2) curlz(ny-1, :) = 0.0_r8
      return
    end if

    if (lwrap) then
      ! Periodic in longitude: use cyclic shifts (like np.roll)
      ! d(fy)/dx ~ (fy(i+1) - fy(i-1)) / dlon
      ! the “x” direction is longitude → second dimension
      !curlz = (  ( cshift(fy, shift=-1, dim=2) - cshift(fy, shift=+1, dim=2) ) / dlon  &
      !         - ( spread(coslat_jn, dim=2, ncopies=nx) * cshift(fx, shift=-1, dim=1)  &
      !           - spread(coslat_js, dim=2, ncopies=nx) * cshift(fx, shift=+1, dim=1) ) / dlat ) &
      !        / ( Re * coslat_2d )
      curlz = (  ( cshift(fy, shift=+1, dim=2) - cshift(fy, shift=-1, dim=2) ) / dlon  &
               - ( spread(coslat_jn, dim=2, ncopies=nx) * cshift(fx, shift=+1, dim=1)  &
                 - spread(coslat_js, dim=2, ncopies=nx) * cshift(fx, shift=-1, dim=1) ) / dlat ) &
              / ( Re * coslat_2d )

    else
      ! Non-periodic in longitude:
      ! Use centered diffs internally; for i=1 and i=nx, fall back to one-sided (or copy)
      curlz = 0.0_r8
      do j = 1, ny
        do i = 1, nx
          ! d(fy)/dlon
          if (i == 1) then
            dfdx = (fy(j,2) - fy(j,1)) / dlon
          else if (i == nx) then
            dfdx = (fy(j,nx) - fy(j,nx-1)) / dlon
          else
            dfdx = (fy(j,i+1) - fy(j,i-1)) / (2.0_r8*dlon)
          end if

          ! (cosφ_{j+1} fx_{j+1} - cosφ_{j-1} fx_{j-1}) / dlat
          if (j == 1) then
            termy = ( coslat(2) * fx(2,i) - coslat(1) * fx(1,i) ) / dlat
          else if (j == ny) then
            termy = ( coslat(ny) * fx(ny,i) - coslat(ny-1) * fx(ny-1,i) ) / dlat
          else
            termy = ( coslat(j+1) * fx(j+1,i) - coslat(j-1) * fx(j-1,i) ) / (2.0_r8 * (dlat/2.0_r8))
            ! The above simplifies back to (cos_{j+1} fx_{j+1} - cos_{j-1} fx_{j-1})/dlat
          end if

          curlz(j,i) = ( dfdx - termy ) / ( Re * coslat(j) )
        end do
      end do
    end if

    ! Match your Python boundary handling:
    !   curlf_z[0:1,:] = 0.0    → j=1 in Fortran
    !   curlf_z[ny-2:ny-1,:] = 0.0 → j=ny-1 in Fortran (penultimate row)
    if (ny >= 1) curlz(1,   :) = 0.0_r8
    if (ny >= 2) curlz(ny-1,:) = 0.0_r8

  end subroutine sphere_curl2

  subroutine sphere_grad2( f, lat, lon, fx, fy )
    real(r8), intent(in)  :: f(:,:)
    real(r8), intent(in)  :: lat(:), lon(:)
    real(r8), intent(out) :: fx(size(f,1), size(f,2))
    real(r8), intent(out) :: fy(size(f,1), size(f,2))
    integer :: ny, nx, i, j
    real(r8), allocatable :: rlat(:), rlon(:), coslat(:)
    real(r8) :: dlat, dlon, deg2rad
    logical :: lwrap
    
    ny = size(lat)
    nx = size(lon)

    deg2rad = acos(-1.0_r8) / 180.0_r8

    allocate(rlat(ny), rlon(nx), coslat(ny) )
    rlat = deg2rad * lat
    rlon = deg2rad * lon

    ! Uniform-step assumption, matching your Python (uses [3]-[1])
    if (ny >= 3) then
      dlat = rlat(3) - rlat(1)
    else
      dlat = 0.0_r8
    end if
    if (nx >= 3) then
      dlon = rlon(3) - rlon(1)
    else
      dlon = 0.0_r8
    end if

    coslat = cos(rlat)
    
    do j=1,ny
       do i=2,nx-1
          fx(i,j) = ( f(i+1,j) - f(i-1,j) ) /( Re * dlon * coslat(j) )
       end do
    end do
    do j=1,ny
       i=1
       fx(i,j) = ( f(i+1,j) - f(nx,j) ) /( Re * dlon * coslat(j) )
       i=nx
       fx(i,j) = ( f(1,j) - f(i-1,j) ) /( Re * dlon * coslat(j) )
    end do

    do j=2,ny-1
       do i=1,nx
          fy(i,j) = ( f(i,j+1) - f(i,j-1) ) /( Re * dlat )
       end do
    end do
    do i=1,nx
       j=1
       fy(i,j) = ( f(i,j+1) - 0._r8 ) /( Re * 0.5 * dlat )
       j=ny
       fy(i,j) = ( 0._r8  - f(i,j-1) ) /( Re * 0.5 * dlat )
    end do

    deallocate(rlat , rlon, coslat )
    
  end subroutine sphere_grad2
  
  subroutine sphere_fronto( th, u, v, lat, lon, frontgf, frontga )
    real(r8), intent(in)  :: th(:,:,:), u(:,:,:), v(:,:,:)
    real(r8), intent(in)  :: lat(:), lon(:)
    real(r8), intent(out) :: frontgf(size(th,1), size(th,2), size(th,3))
    real(r8), intent(out) :: frontga(size(th,1), size(th,2), size(th,3))

    !--- local 
    integer :: ny, nx, nz, i, j, k
    real(r8), allocatable :: rlat(:), rlon(:), coslat(:), tanlat(:)
    real(r8) :: dlat, dlon, deg2rad
    logical  :: lwrap
    real(r8) :: thx(size(th,1), size(th,2), size(th,3))
    real(r8) :: ux(size(th,1), size(th,2), size(th,3))
    real(r8) :: vx(size(th,1), size(th,2), size(th,3))
    real(r8) :: thy(size(th,1), size(th,2), size(th,3))
    real(r8) :: uy(size(th,1), size(th,2), size(th,3))
    real(r8) :: vy(size(th,1), size(th,2), size(th,3))
    real(r8) :: grx(size(th,1), size(th,2) )
    real(r8) :: gry(size(th,1), size(th,2) )

    nx = size(th,1)
    ny = size(th,2)
    nz = size(th,3)
    

    do k=1,nz
       call sphere_grad2( th(:,:,k) , lat, lon, grx, gry )
       thx(:,:,k)=grx
       thy(:,:,k)=gry
    end do
    do k=1,nz
       call sphere_grad2( u(:,:,k) , lat, lon, grx, gry )
       ux(:,:,k)=grx
       uy(:,:,k)=gry
    end do
    do k=1,nz
       call sphere_grad2( v(:,:,k) , lat, lon, grx, gry )
       vx(:,:,k)=grx
       vy(:,:,k)=gry
    end do

    deg2rad = acos(-1.0_r8) / 180.0_r8
    allocate(rlat(ny), rlon(nx), tanlat(ny), coslat(ny) )
    rlat = deg2rad * lat
    rlon = deg2rad * lon

    tanlat=tan(rlat)

    
    do k=1,nz
       do j=1,ny
          frontgf(:,j,k) = -(thx(:,j,k) **2 *( ux(:,j,k) - v(:,j,k)*tanlat(j)/Re ) ) &
          - (thy(:,j,k)**2 * vy(:,j,k)) &
          -(thx(:,j,k)*thy(:,j,k) * &
          ( vx(:,j,k) + uy(:,j,k) + u(:,j,k)*tanlat(j)/Re ) )

       end do
    end do
    
    frontga=frontgf
    
    deallocate(rlat, rlon, tanlat, coslat )

    
  end subroutine sphere_fronto

  subroutine sphere_fronto_wrapper( th, u, v, lat, lon, frontgf, frontga )
    real(r8), intent(in)  :: th(:,:,:), u(:,:,:), v(:,:,:)
    real(r8), intent(in)  :: lat(:), lon(:)
    real(r8), intent(out) :: frontgf(size(th,1), size(th,2), size(th,3))
    real(r8), intent(out) :: frontga(size(th,1), size(th,2), size(th,3))

    !----- local ---------

    real(r8), allocatable  :: thR(:,:,:,:), uR(:,:,:,:), vR(:,:,:,:)
    real(r8), allocatable  :: frontgfR(:,:,:,:), frontgaR(:,:,:,:)
    integer :: nt,nz,ny,nx,ncol,i,j,k,n

    nt   = size(th,3)
    ncol = size(th,1)
    nz   = size(th,2)
    nx   = size(lon)
    ny   = size(lat) 
    

    
    allocate( thR(nx,ny,nz,nt), uR(nx,ny,nz,nt), vR(nx,ny,nz,nt) )
    allocate( frontgfR(nx,ny,nz,nt), frontgaR(nx,ny,nz,nt)  )

    
    uR =reshape( u , [nx,ny,nz,nt] )
    vR =reshape( v , [nx,ny,nz,nt] )
    thR =reshape( th , [nx,ny,nz,nt] )

    
    write(*,*) "In frontogen ","TH      ",minval(TH),    maxval(TH),    shape(TH)

    do n=1,nt
       call sphere_fronto( thR(:,:,:,n), uR(:,:,:,n), vR(:,:,:,n), &
            lat, lon, frontgfR(:,:,:,n), frontgaR(:,:,:,n) )
    end do

    frontgf =reshape( frontgfR , [ncol,nz,nt] )
    frontga =reshape( frontgaR , [ncol,nz,nt] )
    
    deallocate( thR, uR, vR )
    deallocate( frontgfR , frontgaR  )

    
  end subroutine sphere_fronto_wrapper
  
end module sphere_ops

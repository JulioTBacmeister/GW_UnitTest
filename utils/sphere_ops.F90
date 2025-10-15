module sphere_ops
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: Re => rearth

  implicit none
  private
  public :: sphere_curl2_xy

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

end module sphere_ops

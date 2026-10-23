module sphere_se_ops

  !-----------------------------------------------------------------------
  ! Horizontal differential operators for unstructured ("SE") column grids
  ! such as ne30pg3, where the columns carry only lat(ncol)/lon(ncol) and
  ! there is no (i,j) structure to difference along.
  !
  ! This is the unstructured counterpart of sphere_ops, which only works on
  ! regular lat-lon grids. The public routines mirror the sphere_* ones:
  !
  !    se_grad2            <->  sphere_grad2
  !    se_curl2            <->  sphere_curl2
  !    se_fronto           <->  sphere_fronto
  !    se_fronto_wrapper   <->  sphere_fronto_wrapper
  !
  ! Method
  ! ------
  ! Moving least squares on a nearest-neighbour stencil. For each column i
  ! we take its K nearest neighbours, project them into the local tangent
  ! plane at column i (east/north coordinates, metres), and fit
  !
  !     f(j) - f(i)  ~  c1*x + c2*y + c3*x^2 + c4*x*y + c5*y^2
  !
  ! by inverse-distance-weighted least squares. The fit is forced through
  ! the centre value (no constant term), so the gradient at i is just
  !
  !     df/dx = sum_j cgx(i,j) * ( f(nbr(i,j)) - f(i) )
  !
  ! The stencil indices and the coefficients cgx/cgy depend only on the
  ! grid, so se_ops_init builds them once and every later call is a cheap
  ! sparse stencil apply. The fit is done in coordinates normalised by the
  ! local stencil radius to keep the 5x5 normal-equation matrix well
  ! conditioned; if it is singular anyway (degenerate stencil) the column
  ! falls back to a plain linear fit.
  !
  ! Accuracy is second order for the gradient on a quasi-uniform grid. The
  ! cubed-sphere panel edges and corners of a pgN grid are handled without
  ! special-casing, because the fit uses the actual column positions.
  !-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst,    only: Re => rearth
  use cam_abortutils, only: endrun

  implicit none
  private

  public :: se_ops_init
  public :: se_ops_finalize
  public :: se_ops_ready
  public :: se_grad2
  public :: se_curl2
  public :: se_fronto
  public :: se_fronto_wrapper

  interface se_grad2
     module procedure se_grad2_1d
     module procedure se_grad2_2d
  end interface se_grad2

  ! Number of terms in the quadratic fit (no constant term).
  integer, parameter :: nterm = 5

  ! Default stencil size. Needs to exceed nterm comfortably for a stable
  ! least-squares fit; 12 is about two rings on a pgN grid.
  integer, parameter :: kstencil_default = 12

  ! Grid-dependent state, built by se_ops_init.
  integer              :: ncol_s   = 0
  integer              :: kstencil = 0
  integer, allocatable :: nbr(:,:)        ! (ncol,K) neighbour column indices
  real(r8),allocatable :: cgx(:,:)        ! (ncol,K) d/dx stencil weights
  real(r8),allocatable :: cgy(:,:)        ! (ncol,K) d/dy stencil weights
  real(r8),allocatable :: tanlat_s(:)     ! (ncol) tan(latitude), for metric terms
  logical              :: ready = .false.

contains

!==========================================================================

  logical function se_ops_ready()
    se_ops_ready = ready
  end function se_ops_ready

!==========================================================================

  subroutine se_ops_init(lat, lon, kstencil_in, verbose)

    !--------------------------------------------------------------------
    ! Build the neighbour stencils and least-squares gradient weights.
    ! lat/lon are column centres in DEGREES, as carried by the driver.
    ! Call once per grid; calling again rebuilds from scratch.
    !--------------------------------------------------------------------

    real(r8), intent(in) :: lat(:), lon(:)
    integer,  intent(in), optional :: kstencil_in
    logical,  intent(in), optional :: verbose

    real(r8), allocatable :: rlat(:), rlon(:)
    real(r8), allocatable :: xyz(:,:)          ! (3,ncol) unit vectors
    integer,  allocatable :: bin_start(:), bin_list(:), bin_count(:)
    integer,  allocatable :: cand(:)
    real(r8), allocatable :: best_d(:)
    integer,  allocatable :: best_i(:)

    real(r8) :: pi, deg2rad, dmean, dbin, radius
    real(r8) :: ev(3), nv(3), dv(3), pxy(2)
    real(r8) :: h, rho, w, dx, dy, d2
    integer  :: nlat_bin, nlon_bin, nbin
    integer  :: i, j, k, m, n, jb, ib, jj, ii, ibin, ncand, nbest, istat
    integer  :: jb0, jb1, nlon_search
    logical  :: lverbose
    logical  :: ok

    real(r8), allocatable :: Aloc(:,:), wloc(:)
    real(r8) :: Nmat(nterm,nterm), Ninv(nterm,nterm)
    real(r8) :: M2(2,2), M2inv(2,2), det2

    !--------------------------------------------------------------------
    lverbose = .true. ; if (present(verbose)) lverbose = verbose
    kstencil = kstencil_default
    if (present(kstencil_in)) kstencil = max(nterm+2, kstencil_in)

    ncol_s = size(lat)
    if (size(lon) /= ncol_s) call endrun('se_ops_init: lat/lon size mismatch')
    if (ncol_s < kstencil+1) call endrun('se_ops_init: too few columns for stencil')

    call se_ops_finalize()

    pi = acos(-1.0_r8)
    deg2rad = pi/180.0_r8

    allocate(rlat(ncol_s), rlon(ncol_s), xyz(3,ncol_s), stat=istat)
    if (istat /= 0) call endrun('se_ops_init: allocation failure (coords)')

    rlat = deg2rad*lat
    rlon = deg2rad*lon

    do i = 1, ncol_s
       xyz(1,i) = cos(rlat(i))*cos(rlon(i))
       xyz(2,i) = cos(rlat(i))*sin(rlon(i))
       xyz(3,i) = sin(rlat(i))
    end do

    allocate(nbr(ncol_s,kstencil), cgx(ncol_s,kstencil), cgy(ncol_s,kstencil), &
             tanlat_s(ncol_s), stat=istat)
    if (istat /= 0) call endrun('se_ops_init: allocation failure (stencil)')
    nbr = 0 ; cgx = 0._r8 ; cgy = 0._r8
    tanlat_s = tan(rlat)

    !--------------------------------------------------------------------
    ! Bucket the columns into lat-lon bins so the neighbour search is
    ! local rather than O(ncol^2).
    !--------------------------------------------------------------------
    dmean = sqrt(4.0_r8*pi/real(ncol_s,r8))      ! mean angular spacing (rad)
    dbin  = 3.0_r8*dmean
    nlat_bin = max(1, int(pi/dbin))
    nlon_bin = max(1, int(2.0_r8*pi/dbin))
    nbin = nlat_bin*nlon_bin

    allocate(bin_count(nbin), bin_start(nbin+1), bin_list(ncol_s), stat=istat)
    if (istat /= 0) call endrun('se_ops_init: allocation failure (bins)')

    bin_count = 0
    do i = 1, ncol_s
       call bin_of(rlat(i), rlon(i), nlat_bin, nlon_bin, pi, jb, ib)
       ibin = (jb-1)*nlon_bin + ib
       bin_count(ibin) = bin_count(ibin) + 1
    end do

    bin_start(1) = 1
    do n = 1, nbin
       bin_start(n+1) = bin_start(n) + bin_count(n)
    end do

    bin_count = 0
    do i = 1, ncol_s
       call bin_of(rlat(i), rlon(i), nlat_bin, nlon_bin, pi, jb, ib)
       ibin = (jb-1)*nlon_bin + ib
       bin_list(bin_start(ibin)+bin_count(ibin)) = i
       bin_count(ibin) = bin_count(ibin) + 1
    end do

    !--------------------------------------------------------------------
    ! Per-column stencil and weights.
    !--------------------------------------------------------------------
    allocate(cand(ncol_s), best_d(kstencil), best_i(kstencil), stat=istat)
    if (istat /= 0) call endrun('se_ops_init: allocation failure (search)')
    allocate(Aloc(kstencil,nterm), wloc(kstencil), stat=istat)
    if (istat /= 0) call endrun('se_ops_init: allocation failure (lsq)')

    do i = 1, ncol_s

       ! --- gather candidates within an expanding search radius ---------
       radius = 3.0_r8*dmean
       nbest = 0
       do n = 1, 8
          ncand = 0
          call bin_of(rlat(i), rlon(i), nlat_bin, nlon_bin, pi, jb, ib)
          jb0 = max(1, jb - (int(radius/dbin)+1))
          jb1 = min(nlat_bin, jb + (int(radius/dbin)+1))
          ! longitude span widens towards the poles
          if (abs(cos(rlat(i))) < 1.0e-6_r8) then
             nlon_search = nlon_bin
          else
             nlon_search = min(nlon_bin, int(radius/(dbin*abs(cos(rlat(i)))))+2)
          end if
          do jj = jb0, jb1
             do m = -nlon_search, nlon_search
                ii = modulo(ib-1+m, nlon_bin) + 1
                ibin = (jj-1)*nlon_bin + ii
                do k = bin_start(ibin), bin_start(ibin+1)-1
                   j = bin_list(k)
                   if (j == i) cycle
                   ncand = ncand + 1
                   cand(ncand) = j
                end do
             end do
             if (nlon_search >= nlon_bin) exit   ! already took every longitude
          end do
          if (ncand >= kstencil) exit
          radius = 2.0_r8*radius
       end do

       if (ncand < kstencil) &
            call endrun('se_ops_init: could not find enough neighbours')

       ! --- select the K nearest (chord distance on the unit sphere) ----
       nbest = 0
       do n = 1, ncand
          j = cand(n)
          dv = xyz(:,j) - xyz(:,i)
          d2 = dv(1)*dv(1) + dv(2)*dv(2) + dv(3)*dv(3)
          call insert_best(d2, j, best_d, best_i, nbest, kstencil)
       end do

       nbr(i,:) = best_i(:)

       ! --- local tangent-plane coordinates -----------------------------
       ! east and north unit vectors at column i
       ev(1) = -sin(rlon(i)) ; ev(2) =  cos(rlon(i)) ; ev(3) = 0._r8
       nv(1) = -sin(rlat(i))*cos(rlon(i))
       nv(2) = -sin(rlat(i))*sin(rlon(i))
       nv(3) =  cos(rlat(i))

       h = 0._r8
       do k = 1, kstencil
          j = nbr(i,k)
          dv = xyz(:,j) - xyz(:,i)
          dx = Re*( dv(1)*ev(1) + dv(2)*ev(2) + dv(3)*ev(3) )
          dy = Re*( dv(1)*nv(1) + dv(2)*nv(2) + dv(3)*nv(3) )
          Aloc(k,1) = dx
          Aloc(k,2) = dy
          h = h + sqrt(dx*dx + dy*dy)
       end do
       h = h/real(kstencil,r8)
       if (h <= 0._r8) call endrun('se_ops_init: degenerate stencil (zero radius)')

       ! normalise, build the quadratic design matrix and the weights
       do k = 1, kstencil
          dx = Aloc(k,1)/h
          dy = Aloc(k,2)/h
          rho = sqrt(dx*dx + dy*dy)
          if (rho <= 0._r8) call endrun('se_ops_init: duplicate column position')
          wloc(k)   = 1._r8/rho              ! inverse-distance weighting
          Aloc(k,1) = dx
          Aloc(k,2) = dy
          Aloc(k,3) = dx*dx
          Aloc(k,4) = dx*dy
          Aloc(k,5) = dy*dy
       end do

       ! --- weighted normal equations M = A^T W^2 A ---------------------
       do m = 1, nterm
          do n = 1, nterm
             Nmat(m,n) = 0._r8
             do k = 1, kstencil
                Nmat(m,n) = Nmat(m,n) + wloc(k)*wloc(k)*Aloc(k,m)*Aloc(k,n)
             end do
          end do
       end do

       call invert_sym(Nmat, Ninv, nterm, ok)

       if (ok) then
          ! rows 1 and 2 of Ninv A^T W^2 give d/dx and d/dy in normalised
          ! coordinates; divide by h to get them per metre.
          do k = 1, kstencil
             cgx(i,k) = 0._r8
             cgy(i,k) = 0._r8
             do m = 1, nterm
                cgx(i,k) = cgx(i,k) + Ninv(1,m)*Aloc(k,m)
                cgy(i,k) = cgy(i,k) + Ninv(2,m)*Aloc(k,m)
             end do
             cgx(i,k) = cgx(i,k)*wloc(k)*wloc(k)/h
             cgy(i,k) = cgy(i,k)*wloc(k)*wloc(k)/h
          end do
       else
          ! Fall back to a plain linear fit on the same stencil.
          do m = 1, 2
             do n = 1, 2
                M2(m,n) = 0._r8
                do k = 1, kstencil
                   M2(m,n) = M2(m,n) + wloc(k)*wloc(k)*Aloc(k,m)*Aloc(k,n)
                end do
             end do
          end do
          det2 = M2(1,1)*M2(2,2) - M2(1,2)*M2(2,1)
          if (abs(det2) <= 0._r8) call endrun('se_ops_init: singular stencil')
          M2inv(1,1) =  M2(2,2)/det2 ; M2inv(1,2) = -M2(1,2)/det2
          M2inv(2,1) = -M2(2,1)/det2 ; M2inv(2,2) =  M2(1,1)/det2
          do k = 1, kstencil
             cgx(i,k) = ( M2inv(1,1)*Aloc(k,1) + M2inv(1,2)*Aloc(k,2) )*wloc(k)*wloc(k)/h
             cgy(i,k) = ( M2inv(2,1)*Aloc(k,1) + M2inv(2,2)*Aloc(k,2) )*wloc(k)*wloc(k)/h
          end do
       end if

    end do

    deallocate(rlat, rlon, xyz, bin_count, bin_start, bin_list, cand, &
               best_d, best_i, Aloc, wloc)

    ready = .true.

    if (lverbose) then
       write(*,*) 'se_ops_init: unstructured gradient operators built'
       write(*,*) '   ncol              = ', ncol_s
       write(*,*) '   stencil size K    = ', kstencil
       write(*,*) '   mean spacing (km) = ', dmean*Re/1000._r8
       write(*,*) '   lat/lon bins      = ', nlat_bin, nlon_bin
    end if

  end subroutine se_ops_init

!==========================================================================

  subroutine se_ops_finalize()
    if (allocated(nbr))      deallocate(nbr)
    if (allocated(cgx))      deallocate(cgx)
    if (allocated(cgy))      deallocate(cgy)
    if (allocated(tanlat_s)) deallocate(tanlat_s)
    ready = .false.
  end subroutine se_ops_finalize

!==========================================================================

  subroutine se_grad2_1d(f, fx, fy)

    !--- gradient of a single level: f(ncol) -> fx(ncol), fy(ncol) [per m]
    real(r8), intent(in)  :: f(:)
    real(r8), intent(out) :: fx(:), fy(:)

    integer :: i, k, j
    real(r8) :: gx, gy, df

    call check_ready(size(f), 'se_grad2')

    do i = 1, ncol_s
       gx = 0._r8
       gy = 0._r8
       do k = 1, kstencil
          j  = nbr(i,k)
          df = f(j) - f(i)
          gx = gx + cgx(i,k)*df
          gy = gy + cgy(i,k)*df
       end do
       fx(i) = gx
       fy(i) = gy
    end do

  end subroutine se_grad2_1d

!==========================================================================

  subroutine se_grad2_2d(f, fx, fy)

    !--- gradient of every level: f(ncol,nz) -> fx(ncol,nz), fy(ncol,nz)
    real(r8), intent(in)  :: f(:,:)
    real(r8), intent(out) :: fx(:,:), fy(:,:)

    integer :: k

    call check_ready(size(f,1), 'se_grad2')

    do k = 1, size(f,2)
       call se_grad2_1d(f(:,k), fx(:,k), fy(:,k))
    end do

  end subroutine se_grad2_2d

!==========================================================================

  subroutine se_curl2(fx, fy, curlz)

    !--------------------------------------------------------------------
    ! Vertical component of the curl of a horizontal vector field whose
    ! components are given in the local east (fx) / north (fy) directions:
    !
    !   zeta = [ d(fy)/dlam - d(cos(phi) fx)/dphi ] / ( Re cos(phi) )
    !        =   dfy/dx - dfx/dy + fx*tan(phi)/Re
    !
    ! which is the same convention as sphere_ops::sphere_curl2. Input may
    ! be a single level (ncol) or all levels (ncol,nz).
    !--------------------------------------------------------------------
    real(r8), intent(in)  :: fx(:,:), fy(:,:)
    real(r8), intent(out) :: curlz(:,:)

    real(r8), allocatable :: fx_x(:,:), fx_y(:,:), fy_x(:,:), fy_y(:,:)
    integer :: nz, k

    call check_ready(size(fx,1), 'se_curl2')
    nz = size(fx,2)

    allocate(fx_x(ncol_s,nz), fx_y(ncol_s,nz), fy_x(ncol_s,nz), fy_y(ncol_s,nz))

    call se_grad2_2d(fx, fx_x, fx_y)
    call se_grad2_2d(fy, fy_x, fy_y)

    do k = 1, nz
       curlz(:,k) = fy_x(:,k) - fx_y(:,k) + fx(:,k)*tanlat_s(:)/Re
    end do

    deallocate(fx_x, fx_y, fy_x, fy_y)

  end subroutine se_curl2

!==========================================================================

  subroutine se_fronto(th, u, v, frontgf, frontga)

    !--------------------------------------------------------------------
    ! Frontogenesis function on an unstructured grid. Same expression as
    ! sphere_ops::sphere_fronto, with the lat-lon finite differences
    ! replaced by the least-squares gradients.
    !
    !   F = -( dth/dx^2 * ( du/dx - v tan(phi)/Re ) )
    !       -( dth/dy^2 * dv/dy )
    !       -( dth/dx dth/dy * ( dv/dx + du/dy + u tan(phi)/Re ) )
    !
    ! th, u, v are (ncol,nz).
    !
    ! NOTE: frontga is returned as a copy of frontgf, mirroring what
    ! sphere_fronto does. It is a placeholder there too -- the true
    ! frontogenesis ANGLE is not computed by either routine. Only frontgf
    ! is consumed by gw_front_calc; frontga is diagnostic output only.
    !--------------------------------------------------------------------
    real(r8), intent(in)  :: th(:,:), u(:,:), v(:,:)
    real(r8), intent(out) :: frontgf(:,:), frontga(:,:)

    real(r8), allocatable :: thx(:,:), thy(:,:)
    real(r8), allocatable :: ux(:,:),  uy(:,:)
    real(r8), allocatable :: vx(:,:),  vy(:,:)
    integer :: nz, k

    call check_ready(size(th,1), 'se_fronto')
    nz = size(th,2)

    allocate(thx(ncol_s,nz), thy(ncol_s,nz), ux(ncol_s,nz), uy(ncol_s,nz), &
             vx(ncol_s,nz),  vy(ncol_s,nz))

    call se_grad2_2d(th, thx, thy)
    call se_grad2_2d(u,  ux,  uy )
    call se_grad2_2d(v,  vx,  vy )

    do k = 1, nz
       frontgf(:,k) = -( thx(:,k)**2 * ( ux(:,k) - v(:,k)*tanlat_s(:)/Re ) ) &
                      -( thy(:,k)**2 * vy(:,k) )                             &
                      -( thx(:,k)*thy(:,k) *                                 &
                         ( vx(:,k) + uy(:,k) + u(:,k)*tanlat_s(:)/Re ) )
    end do

    frontga = frontgf

    deallocate(thx, thy, ux, uy, vx, vy)

  end subroutine se_fronto

!==========================================================================

  subroutine se_fronto_wrapper(th, u, v, frontgf, frontga)

    !--------------------------------------------------------------------
    ! (ncol,nz,nt) driver-facing wrapper, matching the shape convention of
    ! sphere_ops::sphere_fronto_wrapper. se_ops_init must have been called
    ! for this grid first.
    !--------------------------------------------------------------------
    real(r8), intent(in)  :: th(:,:,:), u(:,:,:), v(:,:,:)
    real(r8), intent(out) :: frontgf(:,:,:), frontga(:,:,:)

    integer :: n, nt

    call check_ready(size(th,1), 'se_fronto_wrapper')
    nt = size(th,3)

    write(*,*) "In se frontogen ","TH      ",minval(th), maxval(th), shape(th)

    do n = 1, nt
       call se_fronto( th(:,:,n), u(:,:,n), v(:,:,n), &
                       frontgf(:,:,n), frontga(:,:,n) )
    end do

  end subroutine se_fronto_wrapper

!==========================================================================
! Private helpers
!==========================================================================

  subroutine check_ready(ncol_in, caller)
    integer,          intent(in) :: ncol_in
    character(len=*), intent(in) :: caller
    if (.not. ready) call endrun(trim(caller)//': se_ops_init has not been called')
    if (ncol_in /= ncol_s) then
       write(*,*) trim(caller)//': ncol = ', ncol_in, ' but operators built for ', ncol_s
       call endrun(trim(caller)//': ncol does not match the initialised grid')
    end if
  end subroutine check_ready

!--------------------------------------------------------------------------

  subroutine bin_of(rlat, rlon, nlat_bin, nlon_bin, pi, jb, ib)
    !--- lat/lon bucket indices for a point (radians) ---
    real(r8), intent(in)  :: rlat, rlon, pi
    integer,  intent(in)  :: nlat_bin, nlon_bin
    integer,  intent(out) :: jb, ib

    jb = int( (rlat + 0.5_r8*pi)/pi * real(nlat_bin,r8) ) + 1
    jb = max(1, min(nlat_bin, jb))
    ib = int( modulo(rlon, 2.0_r8*pi)/(2.0_r8*pi) * real(nlon_bin,r8) ) + 1
    ib = max(1, min(nlon_bin, ib))
  end subroutine bin_of

!--------------------------------------------------------------------------

  subroutine insert_best(d2, j, best_d, best_i, nbest, kmax)
    !--- keep the kmax smallest distances seen so far (insertion sort) ---
    real(r8), intent(in)    :: d2
    integer,  intent(in)    :: j, kmax
    real(r8), intent(inout) :: best_d(kmax)
    integer,  intent(inout) :: best_i(kmax)
    integer,  intent(inout) :: nbest

    integer :: p, q

    if (nbest < kmax) then
       nbest = nbest + 1
       best_d(nbest) = d2
       best_i(nbest) = j
       ! bubble the new entry up into place
       do p = nbest, 2, -1
          if (best_d(p) < best_d(p-1)) then
             call swap_entry(best_d, best_i, p, p-1)
          else
             exit
          end if
       end do
    else if (d2 < best_d(kmax)) then
       best_d(kmax) = d2
       best_i(kmax) = j
       do p = kmax, 2, -1
          if (best_d(p) < best_d(p-1)) then
             call swap_entry(best_d, best_i, p, p-1)
          else
             exit
          end if
       end do
    end if

    q = 0   ! (silence unused warnings on some compilers)

  end subroutine insert_best

!--------------------------------------------------------------------------

  subroutine swap_entry(d, idx, a, b)
    real(r8), intent(inout) :: d(:)
    integer,  intent(inout) :: idx(:)
    integer,  intent(in)    :: a, b
    real(r8) :: dt
    integer  :: it
    dt = d(a) ; d(a) = d(b) ; d(b) = dt
    it = idx(a) ; idx(a) = idx(b) ; idx(b) = it
  end subroutine swap_entry

!--------------------------------------------------------------------------

  subroutine invert_sym(M, Minv, n, ok)
    !--- Gauss-Jordan inverse with partial pivoting of a small matrix ---
    integer,  intent(in)  :: n
    real(r8), intent(in)  :: M(n,n)
    real(r8), intent(out) :: Minv(n,n)
    logical,  intent(out) :: ok

    real(r8) :: A(n,2*n), piv, fac, amax
    integer  :: i, j, k, ip

    ok = .true.
    A = 0._r8
    A(:,1:n) = M
    do i = 1, n
       A(i,n+i) = 1._r8
    end do

    do k = 1, n
       ! partial pivot
       ip = k
       amax = abs(A(k,k))
       do i = k+1, n
          if (abs(A(i,k)) > amax) then
             amax = abs(A(i,k))
             ip = i
          end if
       end do
       if (amax <= 1.0e-12_r8) then
          ok = .false.
          Minv = 0._r8
          return
       end if
       if (ip /= k) then
          do j = 1, 2*n
             fac = A(k,j) ; A(k,j) = A(ip,j) ; A(ip,j) = fac
          end do
       end if
       piv = A(k,k)
       A(k,:) = A(k,:)/piv
       do i = 1, n
          if (i /= k) then
             fac = A(i,k)
             if (fac /= 0._r8) A(i,:) = A(i,:) - fac*A(k,:)
          end if
       end do
    end do

    Minv = A(:,n+1:2*n)

  end subroutine invert_sym

!==========================================================================

end module sphere_se_ops

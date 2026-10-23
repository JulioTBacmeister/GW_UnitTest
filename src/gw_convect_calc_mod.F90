#define UNITTEST
!++ jtb this just mostly shuts down outfld calls
module gw_convect_calc_mod

!-------------------------------------------------------------------------
! Offline driver-side wrapper for the Beres (2004) convective gravity wave
! source, gw_convect::gw_beres_src.
!
! This is the analog of gw_movmtn_calc_mod / gw_front_calc_mod: it packages
! the source call, the drag-profile solve, the momentum/energy fixers and the
! diagnostic output that live inline in gw_drag::gw_tend (the
! "use_gw_convect_dp" block) into something the unit-test driver can call.
!
! Unlike the moving-mountain source, Beres cannot be run with a dummy lookup
! table -- the source spectrum IS the table. So gw_init_beres here is a real
! reader for gw_drag_file, written against plain netCDF rather than the PIO
! calls used by the CAM original.
!-------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8, cl=>shr_kind_cl
use gw_common, only: pver , GWBand
use gw_common, only: west, east, south, north
use physics_types,  only: physics_ptend
use ppgrid, only: pcnst, pcols
use cam_abortutils, only: endrun
use gw_convect,  only: BeresSourceDesc
use nc_flexout_mod

implicit none
private

public :: gw_convect_calc
public :: set_band_convect
public :: set_vramp
public :: report_from_within
public :: gw_init_beres

! A mid-scale "band". Same construction as the frontal band -- Beres uses
! band_mid in CAM, exactly as the C&M frontal source does.
type(GWBand) :: band_mid

real(r8) , pointer :: vramp(:) => null()

! Beres (deep convection) settings and source-spectrum table. The table is
! read once, on the first call, and kept for the life of the run.
type(BeresSourceDesc) :: beres_dp_desc
logical :: beres_dp_ready = .false.

logical, parameter :: masterproc =.TRUE.

contains
!==========================================================================

subroutine gw_convect_calc( &
   ncol, lchnk, dt, pref_edge, &
   u, v, t, netdt, &
   p, piln, zm, zi, &
   nm, ni, rhoi, kvtt, q, dse, effgw_beres_dp, &
   gw_drag_file, &
   !++jtb  ptend defined in massively hacked physics_types
   ptend, flx_heat, lat)

   use coords_1d,  only: Coords1D
   use gw_convect, only: gw_beres_src
   use gw_common,  only: gw_drag_prof, energy_change
   use gw_common,  only: calc_taucd, momentum_flux, momentum_fixer
   use namelist_mod, only: gw_apply_tndmax

   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: lchnk        ! chunk identifier
   real(r8),         intent(in) :: dt           ! Time step.

   real(r8),         intent(in) :: pref_edge(pver+1)  ! reference pressure at edges (Pa)


   real(r8),         intent(in) :: u(ncol,pver)    ! Midpoint zonal winds. ( m s-1)
   real(r8),         intent(in) :: v(ncol,pver)    ! Midpoint meridional winds. ( m s-1)
   real(r8),         intent(in) :: t(ncol,pver)    ! Midpoint temperatures. (K)
   ! Convective heating rate (K s-1). This is TTEND_DP off of pbuf in CAM.
   real(r8),         intent(in) :: netdt(ncol,pver)
   type(Coords1D),   intent(in) :: p               ! Pressure coordinates.
   real(r8),         intent(in) :: piln(ncol,pver+1)  ! Log of interface pressures.
   real(r8),         intent(in) :: zm(ncol,pver)   ! Midpoint altitudes above ground (m).
   real(r8),         intent(in) :: zi(ncol,pver+1) ! Interface altitudes above ground (m).
   real(r8),         intent(in) :: nm(ncol,pver)   ! Midpoint Brunt-Vaisalla frequencies (s-1).
   real(r8),         intent(in) :: ni(ncol,pver+1) ! Interface Brunt-Vaisalla frequencies (s-1).
   real(r8),         intent(in) :: rhoi(ncol,pver+1) ! Interface density (kg m-3).
   real(r8),         intent(in) :: kvtt(ncol,pver+1) ! Molecular thermal diffusivity.
   real(r8),         intent(in) :: q(:,:,:)        ! Constituent array.
   real(r8),         intent(in) :: dse(ncol,pver)  ! Dry static energy.


   real(r8),         intent(in) :: effgw_beres_dp  ! Tendency efficiency.

   ! Lookup table of Beres source spectra (CAM's gw_drag_file).
   character(len=*), intent(in) :: gw_drag_file

   !++jtb Using massively hacked physics_types
   type(physics_ptend), intent(inout):: ptend   ! Parameterization net tendencies.

   real(r8),        intent(out) :: flx_heat(pcols)

   ! Latitude in DEGREES. Optional. When present it is used only to zero the
   ! efficiency at the poles, as gw_tend does. Note CAM works in radians here;
   ! the offline driver carries degrees, so the test is written for degrees.
   real(r8),         intent(in), optional :: lat(ncol)

   !---------------------------Local storage-------------------------------

   integer :: k, m, istat, itime

   real(r8), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(r8), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(r8), allocatable :: phase_speeds(:,:)

   ! Efficiency for a gravity wave source.
   real(r8) :: effgw(ncol)

   ! Indices of top gravity wave source level and lowest level where wind
   ! tendencies are allowed.
   integer :: src_level(ncol)
   integer :: tend_level(ncol)

   ! Projection of wind at midpoints and interfaces.
   real(r8) :: ubm(ncol,pver)
   real(r8) :: ubi(ncol,pver+1)

   ! Unit vectors of source wind (zonal and meridional components).
   real(r8) :: xv(ncol)
   real(r8) :: yv(ncol)

   ! Convective source heating depth (m) and max heating rate (K day-1).
   real(r8) :: hdepth(ncol)
   real(r8) :: maxq0(ncol)

   real(r8) :: utgw(ncol,pver)       ! zonal wind tendency
   real(r8) :: vtgw(ncol,pver)       ! meridional wind tendency
   real(r8) :: ttgw(ncol,pver)       ! temperature tendency
   real(r8) :: qtgw(ncol,pver,pcnst) ! constituents tendencies

   ! Effective gravity wave diffusivity at interfaces.
   real(r8) :: egwdffi(ncol,pver+1)
   real(r8) :: egwdffi_tot(ncol,pver+1)

   ! Temperature tendencies from diffusion and kinetic energy.
   real(r8) :: dttdf(ncol,pver)
   real(r8) :: dttke(ncol,pver)

   ! Provisional absolute wave stress from gw_drag_prof
   real(r8) :: tau_diag(ncol,pver+1)

   ! Reynolds stress for waves propagating in each cardinal direction.
   real(r8) :: taucd(ncol,pver+1,4)
   ! Momentum fluxes used by fixer.
   real(r8) :: um_flux(ncol), vm_flux(ncol)

   ! Energy change used by fixer.
   real(r8) :: de(ncol)

   ! Net stress summed over the spectrum, for diagnostic output.
   real(r8) :: tausum(ncol,pver+1)

   !----------------------------------------------------------------------------
   itime = 1

   !------------------------------------------------------------------
   ! One-time setup: source level, heating depth floor, lookup table.
   ! gw_init (use_gw_convect_dp block) does this at CAM startup.
   !------------------------------------------------------------------
   if (.not. beres_dp_ready) then

      ! Deep convection uses storm-relative winds.
      beres_dp_desc%storm_shift = .true.

      ! Set the deep source level to 700 hPa.
      do k = 0, pver
         if (pref_edge(k+1) < 70000._r8) beres_dp_desc%k = k+1
      end do

      ! Don't use deep convection heating depths below this limit.
      beres_dp_desc%min_hdepth = 1000._r8

      call gw_init_beres(gw_drag_file, band_mid, beres_dp_desc)

      if (masterproc) then
         write (*,*) 'Beres deep level =', beres_dp_desc%k
         write (*,*) ' '
      end if

      beres_dp_ready = .true.

   end if

   ! Allocate wavenumber fields.
   allocate(tau(ncol,-band_mid%ngwv:band_mid%ngwv,pver+1),stat=istat)
   call alloc_err(istat,'gw_convect_calc','tau',ncol*(band_mid%ngwv**2+1)*(pver+1))
   allocate(gwut(ncol,pver,-band_mid%ngwv:band_mid%ngwv),stat=istat)
   call alloc_err(istat,'gw_convect_calc','gwut',ncol*pver*(band_mid%ngwv**2+1))
   allocate(phase_speeds(ncol,-band_mid%ngwv:band_mid%ngwv),stat=istat)
   call alloc_err(istat,'gw_convect_calc','phase_speeds',ncol*(band_mid%ngwv**2+1))

   write(*,*) " ooooohhh  ... in gw_convect_calc"
   write(*,*) " in convect_calc shape( tau ) ",shape( tau )
   write(*,*) "   ... ngwv ", band_mid%ngwv
   write(*,*) " range netdt        " , minval( netdt ) , maxval( netdt )

   ! initialize accumulated momentum fluxes and tendencies
   tau_diag = -9999._r8
   egwdffi_tot = 0._r8

   ! Efficiency of gravity wave momentum transfer.
   ! The pole test is really only to remove the pole points.
   if (present(lat)) then
      where (90._r8 - abs(lat) >= 4*epsilon(1._r8))
         effgw = effgw_beres_dp
      elsewhere
         effgw = 0._r8
      end where
   else
      effgw = effgw_beres_dp
   end if
   write(*,*) " range effgw      " , minval( effgw )  , maxval( effgw )

   ! Determine wave sources for Beres deep scheme
   call gw_beres_src(ncol, band_mid, beres_dp_desc, &
        u, v, netdt, zm, src_level, tend_level, tau, &
        ubm, ubi, xv, yv, phase_speeds, hdepth, maxq0)

   !-------------------------------------------------------------------
   ! gw_beres_src leaves ubi(:,pver+1) undefined -- in CAM ubi is a shared
   ! work array in gw_tend and retains whatever the previous scheme put
   ! there. Nothing downstream reads it (src_level/tend_level from Beres
   ! are always <= pver-1), but define it so this module never carries an
   ! uninitialized value.
   !-------------------------------------------------------------------
   ubi(:,pver+1) = ubm(:,pver)

   write(*,*) " Now back in gw_convect_calc        "
   write(*,*) " range src_level    " , minval( src_level ) , maxval( src_level )
   write(*,*) " range tend_level   " , minval( tend_level ) , maxval( tend_level )
   write(*,*) " range phase_speeds " , minval( phase_speeds ) , maxval( phase_speeds )
   write(*,*) " range tau          " , minval( tau )  , maxval(tau )
   write(*,*) " range hdepth       " , minval( hdepth ) , maxval( hdepth )
   write(*,*) " band_mid%effkwv    " , band_mid%effkwv

   tausum = sum( tau, dim=2 )
   call ncfile_put_col3d('TAU_SRC_BERES',  tausum , itime, 'N m-2', 'spectrum-summed stress after source' )

   ! Solve for the drag profile with Beres source spectrum.
   call gw_drag_prof(ncol, band_mid, p, src_level, tend_level, dt, &
        t, vramp,    &
        piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
        effgw,   phase_speeds,       kvtt, q,  dse,  tau,  utgw,  vtgw, &
        ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
        lapply_effgw_in=gw_apply_tndmax , tau_diag=tau_diag )
   write(*,*) " in gw_convect_calc after gw_drag_prof "

   ! Project stress into directional components.
   taucd = calc_taucd(ncol, band_mid%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)
   write(*,*) " in gw_convect_calc after calc_taucd "

   !  add the diffusion coefficients
   do k = 1, pver+1
      egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
   end do

   ! Store constituents tendencies
   do m = 1, pcnst
      do k = 1, pver
         ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
      end do
   end do

   ! Find momentum flux, and use it to fix the wind tendencies below
   ! the gravity wave region.
   call momentum_flux(tend_level, taucd, um_flux, vm_flux)
   call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

   ! Add the momentum tendencies to the output tendency arrays.
   do k = 1, pver
      ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
      ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
   end do

   !-------------------------------------------------------------------
   ! gw_tend runs energy_change/energy_fixer here, using the flx_heat
   ! carried in from the orographic scheme. flx_heat is intent(out) per
   ! call in this driver, so -- as in gw_front_calc -- the fixer is left
   ! out and flx_heat is just reported back.
   !
   !call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
   !     ptend%v(:ncol,:), ptend%s(:ncol,:)+ttgw, de)
   !call energy_fixer(tend_level, p, de-flx_heat(:ncol), ttgw)
   !-------------------------------------------------------------------

   do k = 1, pver
      ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
   end do

   ! Calculate energy change for output to CAM's energy checker.
   call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:), de)
   flx_heat(:ncol) = de

   !-----------------------------------------------------------------
   ! Diagnostics. HDEPTH/MAXQ0/NETDT mirror gw_tend's outfld calls.
   !-----------------------------------------------------------------
   tausum = sum( tau, dim=2 )

   call ncfile_put_col3d('NETDT',  netdt , itime, 'K s-1', 'heating rate driving Beres source' )
   call ncfile_put_col2d('HDEPTH', hdepth/1000._r8 , itime, 'km', 'heating depth' )
   call ncfile_put_col2d('MAXQ0',  maxq0 , itime, 'K day-1', 'max column heating rate' )
   call ncfile_put_col2d('K_SRC_BERES', 1.0_r8 * src_level , itime, '1', 'source level idx for Beres GW' )

   call ncfile_put_col3d('UBM_BERES', ubm , itime, 'ms-1', 'on wave wind at midlayer' )
   call ncfile_put_col3d('TAU_BERES', tausum , itime, 'N m-2', 'spectrum-summed stress profile - Beres' )
   call ncfile_put_col3d('TAU_DIAG_BERES',  tau_diag , itime, 'N m-2', 'pre-pixie stress profile - Beres'  )
   call ncfile_put_col3d('UTGW_BERES', utgw , itime, 'ms-2', 'zonal wind tendency - Beres' )
   call ncfile_put_col3d('VTGW_BERES', vtgw , itime, 'ms-2', 'meridional wind tendency - Beres' )

   ! Tau in each direction.
   call ncfile_put_col3d('TAUE_BERES', taucd(:,:,east)  , itime, 'N m-2', 'E stress profile - Beres' )
   call ncfile_put_col3d('TAUW_BERES', taucd(:,:,west)  , itime, 'N m-2', 'W stress profile - Beres' )
   call ncfile_put_col3d('TAUN_BERES', taucd(:,:,north) , itime, 'N m-2', 'N stress profile - Beres' )
   call ncfile_put_col3d('TAUS_BERES', taucd(:,:,south) , itime, 'N m-2', 'S stress profile - Beres' )

   write(*,*) " wheeeewww   ... survived gw_convect_calc"

   deallocate(tau, gwut, phase_speeds)

end subroutine gw_convect_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_band_convect(band_mid_in)
  type(GWBand) , intent(in) :: band_mid_in
! Just to bring band_mid into module
band_mid=band_mid_in
write(*,*) " Band Convect (Beres) ",band_mid%ngwv
end subroutine set_band_convect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_vramp()

! The driver calls this once per time step; only allocate the first time.
if (.not. associated(vramp)) allocate(vramp(pver))
vramp(:) = 1._r8

end subroutine set_vramp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine report_from_within()

write(*,*) " Inside gw_convect_calc_mod: Band Mid ",band_mid%ngwv
write(*,*) " Inside gw_convect_calc_mod: pcols ",pcols
end subroutine report_from_within


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alloc_err( istat, routine, name, nelem )

  !-----------------------------------------------------------------------
  ! Purpose:
  ! Issue error message after non-zero return from an allocate statement.
  !
  ! Author: B. Eaton
  !-----------------------------------------------------------------------

  integer, intent(in) ::&
       istat           ! status from allocate statement
  character(len=*), intent(in) ::&
       routine,       &! routine that called allocate
       name            ! name of array
  integer, intent(in) ::&
       nelem           ! number of elements attempted to allocate
  !-----------------------------------------------------------------------

  integer :: iulog
  iulog=6

  if ( istat .ne. 0 ) then
     write(iulog,*)'ERROR trying to allocate memory in routine: ' &
          //trim(routine)
     write(iulog,*)'  Variable name: '//trim(name)
     write(iulog,*)'  Number of elements: ',nelem
     call endrun ('ALLOC_ERR')
  end if
end subroutine alloc_err


!==============================================================
subroutine gw_init_beres(file_name, band, desc)

  !-------------------------------------------------------
  ! Read the Beres source-spectrum lookup table.
  !
  ! Same job as the routine of this name in {orig}/gw_drag.F90, but using
  ! plain netCDF instead of PIO, in the style of utils/ncread_mod.F90.
  !
  ! The file (e.g. newmfspectra40_dc25.nc) holds mfcc(PS,MW,HD) in CDL
  ! order, which is mfcc(HD,MW,PS) in Fortran -- matching the shape of
  ! desc%mfcc. The phase-speed dimension is clipped to the model's ngwv,
  ! centered on c=0, which is what the offset start index below does.
  !-------------------------------------------------------
  use netcdf

  character(len=*), intent(in) :: file_name
  type(GWBand), intent(in) :: band
  type(BeresSourceDesc), intent(inout) :: desc

  ! netCDF ids and error code.
  integer :: ncid, status, dimid, varid

  ! Number of wavenumbers in the input file.
  integer :: ngwv_file

  integer :: stat
  character(len=cl) :: msg

  !----------------------------------------------------------------------
  ! read in look-up table for source spectra
  !-----------------------------------------------------------------------

  status = nf90_open( trim(file_name), nf90_nowrite, ncid )
  if (status /= nf90_noerr) then
     write(*,*) 'gw_init_beres: cannot open '//trim(file_name)
     call endrun('gw_init_beres: cannot open gw_drag_file')
  end if

  ! Get HD (heating depth) dimension.
  status = nf90_inq_dimid(ncid, 'HD', dimid)
  call check_nc(status, 'inq_dimid HD in '//trim(file_name))
  status = nf90_inquire_dimension(ncid, dimid, len=desc%maxh)
  call check_nc(status, 'inquire_dimension HD in '//trim(file_name))

  ! Get MW (mean wind) dimension.
  status = nf90_inq_dimid(ncid, 'MW', dimid)
  call check_nc(status, 'inq_dimid MW in '//trim(file_name))
  status = nf90_inquire_dimension(ncid, dimid, len=desc%maxuh)
  call check_nc(status, 'inquire_dimension MW in '//trim(file_name))

  ! Get PS (phase speed) dimension.
  status = nf90_inq_dimid(ncid, 'PS', dimid)
  call check_nc(status, 'inq_dimid PS in '//trim(file_name))
  status = nf90_inquire_dimension(ncid, dimid, len=ngwv_file)
  call check_nc(status, 'inquire_dimension PS in '//trim(file_name))

  ! Number in each direction is half of total (and minus phase speed of 0).
  desc%maxuh = (desc%maxuh-1)/2
  ngwv_file = (ngwv_file-1)/2

  if (ngwv_file < band%ngwv) then
     write(*,*) 'gw_init_beres: table ngwv = ', ngwv_file, ' model ngwv = ', band%ngwv
     call endrun('gw_init_beres: PhaseSpeed in lookup table file does not &
          &cover the whole spectrum implied by the model''s ngwv.')
  end if

  ! Allocate hd and get data.
  allocate(desc%hd(desc%maxh), stat=stat, errmsg=msg)
  call alloc_err(stat, 'gw_init_beres', 'hd', desc%maxh)

  status = nf90_inq_varid(ncid, 'HD', varid)
  call check_nc(status, 'Error finding HD in '//trim(file_name))
  status = nf90_get_var(ncid, varid, desc%hd, start=[1], count=[desc%maxh])
  call check_nc(status, 'Error reading HD from '//trim(file_name))

  ! While not currently documented in the file, it uses kilometers. Convert
  ! to meters.
  desc%hd = desc%hd*1000._r8

  ! Allocate mfcc. "desc%maxh" and "desc%maxuh" are from the file, but the
  ! model determines wavenumber dimension.
  allocate(desc%mfcc(desc%maxh,-desc%maxuh:desc%maxuh,&
       -band%ngwv:band%ngwv), stat=stat, errmsg=msg)
  call alloc_err(stat, 'gw_init_beres', 'mfcc', &
       desc%maxh*(2*desc%maxuh+1)*(2*band%ngwv+1))

  ! Get mfcc data.
  status = nf90_inq_varid(ncid, 'mfcc', varid)
  call check_nc(status, 'Error finding mfcc in '//trim(file_name))
  status = nf90_get_var(ncid, varid, desc%mfcc, &
       start=[1,1,ngwv_file-band%ngwv+1], count=shape(desc%mfcc))
  call check_nc(status, 'Error reading mfcc from '//trim(file_name))

  status = nf90_close(ncid)

  if (masterproc) then
     write(*,*) "Read in source spectra from file."
     write(*,*) "   file  = ", trim(file_name)
     write(*,*) "   maxh, maxuh, ngwv_file = ", desc%maxh, desc%maxuh, ngwv_file
     write(*,*) "   hd   max, min = ", maxval(desc%hd),", ",minval(desc%hd)
     write(*,*) "   mfcc max, min = ", maxval(desc%mfcc),", ",minval(desc%mfcc)
  endif

contains

  subroutine check_nc(stat_in, context)
    integer, intent(in) :: stat_in
    character(len=*), intent(in) :: context
    if (stat_in /= nf90_noerr) then
       write(*,*) 'gw_init_beres: '//trim(context)
       write(*,*) '   '//trim(nf90_strerror(stat_in))
       call endrun('gw_init_beres: netCDF error')
    end if
  end subroutine check_nc

end subroutine gw_init_beres
!==========================================================================


end module gw_convect_calc_mod

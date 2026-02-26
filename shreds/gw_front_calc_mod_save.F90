#define UNITTEST
!++ jtb this just mostly shuts down outfld calls
module gw_front_calc_mod

use shr_kind_mod,   only: r8=>shr_kind_r8, cl=>shr_kind_cl
use gw_common, only: pver , GWBand
use physics_types,  only: physics_ptend
use ppgrid, only: pcnst, pcols
use cam_abortutils, only: endrun
use gw_front,  only: CMSourceDesc,gaussian_cm_desc
use nc_flexout_mod

implicit none
private

public :: gw_front_calc
public :: set_band_front
public :: set_vramp
public :: report_from_within

! A mid-scale "band".
type(GWBand) :: band_front

real(r8) , pointer :: vramp(:)

logical, parameter :: masterproc =.TRUE.


contains
!==========================================================================

subroutine gw_front_calc( &
   ncol, lchnk, dt, pref_edge, &
   u, v, t, frontgf, &
   p, piln, zm, zi, &
   nm, ni, rhoi, kvtt, q, dse, effgw_front, &
   frontgfc, taubgnd, front_gaussian_width, &
   !++jtb  ptend defined in massively hacked physics_types
   ptend, flx_heat) 

   use coords_1d,  only: Coords1D
   use gw_front,  only: gw_cm_src, gaussian_cm_desc
   use gw_common,  only: gw_drag_prof, energy_change
   use gw_common,  only: calc_taucd, momentum_flux, momentum_fixer
   use namelist_mod, only: gw_apply_tndmax

   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: lchnk        ! chunk identifier
   real(r8),         intent(in) :: dt           ! Time step.
   ! Some namelist driven params
   real(r8),         intent(in) :: frontgfc, taubgnd, front_gaussian_width

   real(r8),         intent(in) :: pref_edge(pver+1)  ! reference pressure at edges (Pa)

   
   real(r8),         intent(in) :: u(ncol,pver)    ! Midpoint zonal winds. ( m s-1)
   real(r8),         intent(in) :: v(ncol,pver)    ! Midpoint meridional winds. ( m s-1)
   real(r8),         intent(in) :: t(ncol,pver)    ! Midpoint temperatures. (K)
   real(r8),         intent(in) :: frontgf(ncol,pver) ! 3D frontogenesis (units??)
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


   real(r8),         intent(in) :: effgw_front     ! Tendency efficiency.

   !++jtb Using massively hacked physics_types
   type(physics_ptend), intent(inout):: ptend   ! Parameterization net tendencies. 

   real(r8),        intent(out) :: flx_heat(pcols)

   !---------------------------Local storage-------------------------------

   ! Frontogenesis wave settings.
   type(CMSourceDesc) :: cm_desc

   integer :: k, m, nn, istat, itime

   real(r8), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(r8), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(r8), allocatable :: phase_speeds(:,:)

   ! Efficiency for a gravity wave source.
   real(r8) :: effgw(ncol), effgw_cm

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

   ! Averages over source region.
   real(r8) :: ubmsrc(ncol) ! On-ridge wind.
   real(r8) :: usrc(ncol)   ! Zonal wind.
   real(r8) :: vsrc(ncol)   ! Meridional wind.
   real(r8) :: nsrc(ncol)   ! B-V frequency.
   real(r8) :: rsrc(ncol)   ! Density.

   ! normalized wavenumber
   real(r8) :: m2src(ncol)

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

   integer :: kfront,kbot_front
   
   character(len=1) :: cn
   character(len=9) :: fname(4)
   !----------------------------------------------------------------------------
   itime = 1

   
   do k = 0, pver
      ! Check frontogenesis at 600 hPa.
      if (pref_edge(k+1) < 60000._r8) kfront = k+1
   end do

   ! Source waves from 500 hPa.
   kbot_front = maxloc(pref_edge, 1, (pref_edge < 50000._r8)) - 1

   if (masterproc) then
      write (*,*) 'KFRONT      =',kfront
      write (*,*) 'KBOT_FRONT  =',kbot_front
      write(*,*) ' '
   end if
   

   cm_desc = gaussian_cm_desc(band_front, kbot_front, kfront, frontgfc, &
        taubgnd, front_gaussian_width)

   

   ! Allocate wavenumber fields.
   allocate(tau(ncol,band_front%ngwv:band_front%ngwv,pver+1),stat=istat)
   call alloc_err(istat,'gw_rdg_calc','tau',ncol*(band_front%ngwv**2+1)*(pver+1))
   allocate(gwut(ncol,pver,band_front%ngwv:band_front%ngwv),stat=istat)
   call alloc_err(istat,'rdg_calc','gwut',ncol*pver*(band_front%ngwv**2+1))
   allocate(phase_speeds(ncol,band_front%ngwv:band_front%ngwv),stat=istat)
   call alloc_err(istat,'rdg_calc','phase_speeds',ncol*(band_front%ngwv**2+1))


   
   write(*,*) " ooooohhh  ... in gw_front_calc"
   ! return
   
   ! initialize accumulated momentum fluxes and tendencies
   tau_diag = -9999._r8


     ! Efficiency of gravity wave momentum transfer.
     effgw = effgw_cm
     ! Frontogenesis is too high at the poles (at least for the FV
     ! dycore), so introduce a polar taper.
     !if (gw_polar_taper) effgw = effgw * cos(state1%lat(:ncol))

     ! Determine the wave source for C&M background spectrum
     call gw_cm_src(ncol, band_front, cm_desc, u, v, frontgf(:ncol,:), &
          src_level, tend_level, tau, ubm, ubi, xv, yv, phase_speeds)

     ! Solve for the drag profile with C&M source spectrum.
     call gw_drag_prof(ncol, band_front, p, src_level, tend_level,  dt, &
          t, vramp,   &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,   phase_speeds,       kvtt, q,  dse,  tau,  utgw,  vtgw, &
          ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
          lapply_effgw_in=gw_apply_tndmax)

     ! Project stress into directional components.
     taucd = calc_taucd(ncol, band_front%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

     !  add the diffusion coefficients
     do k = 1, pver+1
        egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
     end do

     !Add the constituent tendencies
     do m=1, pcnst
        do k = 1, pver
           ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
        end do
     end do

     ! Find momentum flux, and use it to fix the wind tendencies below
     ! the gravity wave region.
     call momentum_flux(tend_level, taucd, um_flux, vm_flux)
     call momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

     ! add the momentum tendencies to the output tendency arrays
     do k = 1, pver
        ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
        ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
     end do

     ! Find energy change in the current state, and use fixer to apply
     ! the difference in lower levels.
     !call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
     !     ptend%v(:ncol,:), ptend%s(:ncol,:)+ttgw, de)
     !call energy_fixer(tend_level, p, de-flx_heat(:ncol), ttgw)

     do k = 1, pver
        ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
     end do

   ! Calculate energy change for output to CAM's energy checker.
   call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:), de)
   flx_heat(:ncol) = de

   write(*,*) " wheeeewww   ... survived gw_movtmn_calc"


   deallocate(tau, gwut, phase_speeds)

end subroutine gw_front_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_band_front(band_front_in)
  type(GWBand) , intent(in) :: band_front_in
! Just to bring band_front into module
band_front=band_front_in
write(*,*) " Band MovMTN ",band_front%ngwv
end subroutine set_band_front

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_vramp()

allocate(vramp(pver))
vramp(:) = 1._r8

end subroutine set_vramp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine report_from_within()

write(*,*) " Inside gw_rdg_calc_mod: Band MovMtn ",band_front%ngwv
write(*,*) " Inside gw_rdg_calc_mod: pcols ",pcols
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
subroutine gw_init_front(file_name, band, desc)

  !-------------------------------------------------------
  ! Highly abbridged version of this subr that resides in
  ! {orig}/gw_drag.F90
  !-------------------------------------------------------
  character(len=*), intent(in) :: file_name
  type(GWBand), intent(in) :: band
  type(CMSourceDesc), intent(inout) :: desc
  

end subroutine gw_init_front
!==========================================================================


end module gw_front_calc_mod

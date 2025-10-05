#define UNITTEST
!++ jtb this just mostly shuts down outfld calls
module gw_movmtn_calc_mod

use shr_kind_mod,   only: r8=>shr_kind_r8, cl=>shr_kind_cl
use gw_common, only: pver , GWBand
use physics_types,  only: physics_ptend
use ppgrid, only: pcnst, pcols
use cam_abortutils, only: endrun
use gw_movmtn,  only: MovMtnSourceDesc

implicit none
private

public :: gw_movmtn_calc
public :: set_band_movmtn
public :: set_vramp
public :: report_from_within

! A mid-scale "band" with only one flow dependent phase speed (l = 0).
type(GWBand) :: band_movmtn

real(r8) , pointer :: vramp(:)

contains
!==========================================================================

subroutine gw_movmtn_calc( &
   ncol, lchnk, dt, pref_edge, &
   u, v, t, vort4gw, &
   p, piln, zm, zi, &
   nm, ni, rhoi, kvtt, q, dse, effgw_movmtn_pbl, &
   !++jtb  ptend defined in massively hacked physics_types
   ptend, flx_heat) 

   use coords_1d,  only: Coords1D
   use gw_movmtn,  only: gw_movmtn_src
   use gw_common,  only: gw_drag_prof, energy_change
   use namelist_mod, only: alpha_gw_movmtn, movmtn_plaunch, movmtn_psteer, movmtn_source 
   use namelist_mod, only: gw_apply_tndmax

   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: lchnk        ! chunk identifier
   real(r8),         intent(in) :: dt           ! Time step.

   real(r8),         intent(in) :: pref_edge(pver+1)  ! reference pressure at edges (Pa)

   
   real(r8),         intent(in) :: u(ncol,pver)    ! Midpoint zonal winds. ( m s-1)
   real(r8),         intent(in) :: v(ncol,pver)    ! Midpoint meridional winds. ( m s-1)
   real(r8),         intent(in) :: t(ncol,pver)    ! Midpoint temperatures. (K)
   real(r8),         intent(in) :: vort4gw(ncol,pver) ! 3D vorticity field (s-1)
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


   real(r8),         intent(in) :: effgw_movmtn_pbl  ! Tendency efficiency.

   !++jtb Using massively hacked physics_types
   type(physics_ptend), intent(inout):: ptend   ! Parameterization net tendencies. 

   real(r8),        intent(out) :: flx_heat(pcols)

   !---------------------------Local storage-------------------------------

   ! Moving mountain settings and table.
   type(MovMtnSourceDesc) :: movmtn_desc

   integer :: k, m, nn, istat

   real(r8), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(r8), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(r8), allocatable :: phase_speeds(:,:)

   ! Isotropic source flag [anisotropic orography].
   integer  :: isoflag(ncol)

   ! Efficiency for a gravity wave source.
   real(r8) :: effgw(ncol)

   ! Indices of top gravity wave source level and lowest level where wind
   ! tendencies are allowed.
   integer :: src_level(ncol)
   integer :: tend_level(ncol)
   integer :: bwv_level(ncol)
   integer :: tlb_level(ncol)

   
   ! Indices of source(launch) level and steering level
   integer ::  movmtn_ksteer, movmtn_klaunch

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

   ! Wave stress in zonal/meridional direction
   real(r8) :: taummx(ncol,pver+1)
   real(r8) :: taummy(ncol,pver+1)
   ! Provisional absolute wave stress from gw_drag_prof
   real(r8) :: tau_diag(ncol,pver+1)

   ! Energy change used by fixer.
   real(r8) :: de(ncol)
   ! Dummy placeholders for source variables
   real(r8) :: upwp_clubb_gw(ncol,pver), vpwp_clubb_gw(ncol,pver), xpwp_clubb(ncol,pver) 
   real(r8) :: ttend_clubb(ncol,pver), ttend_dp(ncol,pver) , hdepth(ncol)

   character(len=1) :: cn
   character(len=9) :: fname(4)
   !----------------------------------------------------------------------------
   
   ! Allocate wavenumber fields.
   allocate(tau(ncol,band_movmtn%ngwv:band_movmtn%ngwv,pver+1),stat=istat)
   call alloc_err(istat,'gw_rdg_calc','tau',ncol*(band_movmtn%ngwv**2+1)*(pver+1))
   allocate(gwut(ncol,pver,band_movmtn%ngwv:band_movmtn%ngwv),stat=istat)
   call alloc_err(istat,'rdg_calc','gwut',ncol*pver*(band_movmtn%ngwv**2+1))
   allocate(phase_speeds(ncol,band_movmtn%ngwv:band_movmtn%ngwv),stat=istat)
   call alloc_err(istat,'rdg_calc','phase_speeds',ncol*(band_movmtn%ngwv**2+1))


   !=============================================================
   ! None of this needed for the gw_movmtn_pbl source but here for
   ! future devel0pment
   call gw_init_movmtn("DummyName.nc", band_movmtn, movmtn_desc)
   do k = 0, pver
      ! 950 hPa index
      if (pref_edge(k+1) < 95000._r8) movmtn_desc%k = k+1
   end do

   ! Don't use deep convection heating depths below this limit.
   movmtn_desc%min_hdepth = 1._r8
   !=============================================================
   
   ! Find indices of sterring and launch levels
   do k = 1, pver
      ! Find steering level
      if ( (pref_edge(k+1) >= movmtn_psteer).and.(pref_edge(k) < movmtn_psteer) ) then
         movmtn_ksteer = k
      end if
   end do
   do k = 1, pver
      ! Find launch level
      if ( (pref_edge(k+1) >= movmtn_plaunch).and.(pref_edge(k) < movmtn_plaunch ) ) then
         movmtn_klaunch = k
      end if
   end do

   
   write(*,*) " ooooohhh  ... in gw_movtmn_calc"
   ! return
   
   ! initialize accumulated momentum fluxes and tendencies
   taummx = 0._r8
   taummy = 0._r8
   tau_diag = -9999._r8

   upwp_clubb_gw = 0.
   vpwp_clubb_gw = 0.
   ttend_clubb = 0.
   ttend_dp = 0.

   xpwp_clubb(:ncol,:) = sqrt( upwp_clubb_gw(:ncol,:)**2 + vpwp_clubb_gw(:ncol,:)**2 )

   call gw_movmtn_src(ncol, lchnk, band_movmtn , movmtn_desc, &
        u, v, ttend_dp(:ncol,:), ttend_clubb(:ncol,:), xpwp_clubb(:ncol,:), vort4gw(:ncol,:), &
        zm, alpha_gw_movmtn, movmtn_source, movmtn_ksteer, movmtn_klaunch, src_level, tend_level, &
        tau, ubm, ubi, xv, yv, &
        phase_speeds, hdepth)
   !-------------------------------------------------------------
   ! gw_movmtn_src returns wave-relative wind profiles ubm,ubi
   ! and unit vector components describing direction of wavevector
   ! and application of wave-drag force. I believe correct setting
   ! for c is c=0, since it is incorporated in ubm and (xv,yv)
   !--------------------------------------------------------------
   write(*,*) " Now back in gw_movmtn_calc  topi       " 
   write(*,*) " range src_level    " , minval( src_level ) , maxval( src_level )
   write(*,*) " range tend_level   " , minval( tend_level ) , maxval( tend_level )
   write(*,*) " range phase_speeds " , minval( phase_speeds ) , maxval( phase_speeds )
   write(*,*) " range tau          " , minval( tau )  , maxval(tau )
   write(*,*) " band_movmtn%effkwv " , band_movmtn%effkwv


#if 0
   call outfld('SRC_LEVEL_MOVMTN', real(src_level,r8), ncol, lchnk)
   call outfld('TND_LEVEL_MOVMTN', real(tend_level,r8), ncol, lchnk)
   call outfld('UBI_MOVMTN', ubi, ncol, lchnk)
   call outfld('UBM_MOVMTN', ubm, ncol, lchnk)
#endif
   
   effgw = effgw_movmtn_pbl
   write(*,*) " range effgw      " , minval( effgw )  , maxval( effgw )
   call gw_drag_prof(ncol, band_movmtn, p, src_level, tend_level, dt, &
        t, vramp,    &
        piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
        effgw,   phase_speeds,       kvtt, q,  dse,  tau,  utgw,  vtgw, &
        ttgw, qtgw, egwdffi,  gwut, dttdf, dttke,            &
        lapply_effgw_in=gw_apply_tndmax )

   ! Project stress into directional components.
   !! taucd = calc_taucd(ncol, band_movmtn%ngwv, tend_level, tau, phase_speeds, xv, yv, ubi)

   !  add the diffusion coefficients
   do k = 1, pver+1
      egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
   end do

   ! Add the tendencies from isotropic residual to the totals.
   do k = 1, pver
      !++jtb: Using massively hacked physics_types ...
      ! physics tendencies
      ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
      ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
      ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
   end do

   do m = 1, pcnst
      do k = 1, pver
         ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
      end do
   end do

   do k = 1, pver+1
      taummx(:,k) =  tau(:,0,k)*xv
      taummy(:,k) =  tau(:,0,k)*yv
   end do

   write(20) ubm
   write(20) vort4gw
   write(20) tau(:,0,:)

   ! Calculate energy change for output to CAM's energy checker.
   call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:), de)
   flx_heat(:ncol) = de

   write(*,*) " wheeeewww   ... survived gw_movtmn_calc"


   deallocate(tau, gwut, phase_speeds)

end subroutine gw_movmtn_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_band_movmtn(band_movmtn_in)
  type(GWBand) , intent(in) :: band_movmtn_in
! Just to bring band_movmtn into module
band_movmtn=band_movmtn_in
write(*,*) " Band MovMTN ",band_movmtn%ngwv
end subroutine set_band_movmtn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_vramp()

allocate(vramp(pver))
vramp(:) = 1._r8

end subroutine set_vramp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine report_from_within()

write(*,*) " Inside gw_rdg_calc_mod: Band MovMtn ",band_movmtn%ngwv
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
subroutine gw_init_movmtn(file_name, band, desc)

  !-------------------------------------------------------
  ! Highly abbridged version of this subr that resides in
  ! {orig}/gw_drag.F90
  !-------------------------------------------------------
  character(len=*), intent(in) :: file_name
  type(GWBand), intent(in) :: band
  type(MovMtnSourceDesc), intent(inout) :: desc

  ! PIO variable ids and error code.
  integer :: mfccid, uhid, hdid, stat

  ! Number of wavenumbers in the input file.
  integer :: ngwv_file

  ! Full path to gw_drag_file.
  character(len=cl) :: file_path

  character(len=cl) :: msg

  !----------------------------------------------------------------------
  ! read in look-up table for source spectra
  !-----------------------------------------------------------------------

  ! HD (heating depth) dimension.
  desc%maxh = 15 !get_pio_dimlen(gw_file_desc, "HD", file_path)

  ! MW (mean wind) dimension.
  desc%maxuh = 241 ! get_pio_dimlen(gw_file_desc, "MW", file_path)

  ! Get PS (phase speed) dimension.
  ngwv_file = 0 !get_pio_dimlen(gw_file_desc, "PS", file_path)

  ! Number in each direction is half of total (and minus phase speed of 0).
  desc%maxuh = (desc%maxuh-1)/2
  ngwv_file = (ngwv_file-1)/2

  ! Allocate hd and get data.

  allocate(desc%hd(desc%maxh), stat=stat, errmsg=msg)

  ! While not currently documented in the file, it uses kilometers. Convert
  ! to meters.
  desc%hd = 1000._r8

 ! Allocate wind.
  allocate(desc%uh(desc%maxuh), stat=stat, errmsg=msg)
  desc%uh = 0._r8

  ! Allocate mfcc. "desc%maxh" and "desc%maxuh".
  allocate(desc%mfcc(desc%maxh,-desc%maxuh:desc%maxuh,&
       -band%ngwv:band%ngwv), stat=stat, errmsg=msg)
  desc%mfcc = -999._r8
  

end subroutine gw_init_movmtn
!==========================================================================


end module gw_movmtn_calc_mod

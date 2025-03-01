#define UNITTEST
!++ jtb this just mostly shuts down outfld calls
module gw_rdg_calc_mod

use shr_kind_mod,   only: r8=>shr_kind_r8, cl=>shr_kind_cl
use gw_common, only: pver , GWBand
use physics_types,  only: physics_ptend
use ppgrid, only: pcnst, pcols
use cam_abortutils, only: endrun

implicit none
private

public :: gw_rdg_calc
public :: set_band_rdg
public :: set_vramp
public :: report_from_within

! anisotropic ridge fields
integer, parameter :: prdg = 16


! A mid-scale "band" with only stationary waves (l = 0).
type(GWBand) :: band_oro

real(r8) , pointer :: vramp(:)

contains
!==========================================================================

subroutine gw_rdg_calc( &
   type, ncol, lchnk, n_rdg, dt, &
   u, v, t, p, piln, zm, zi, &
   nm, ni, rhoi, kvtt, q, dse, &
   effgw_rdg, effgw_rdg_max, &
   effgw_rdg_resid, luse_gw_rdg_resid, &
   hwdth, clngt, gbxar, &
   mxdis, angll, anixy, &
   isovar, isowgt, &
   rdg_cd_llb, trpd_leewv, &
   !++jtb  ptend defined in massively hacked physics_types
   ptend, flx_heat) 

   use coords_1d,  only: Coords1D
   use gw_rdg,     only: gw_rdg_src, gw_rdg_resid_src, gw_rdg_belowpeak, gw_rdg_break_trap, gw_rdg_do_vdiff
   use gw_common,  only: gw_drag_prof, energy_change

   character(len=5), intent(in) :: type         ! BETA or GAMMA
   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: lchnk        ! chunk identifier
   integer,          intent(in) :: n_rdg
   real(r8),         intent(in) :: dt           ! Time step.

   real(r8),         intent(in) :: u(ncol,pver)    ! Midpoint zonal winds. ( m s-1)
   real(r8),         intent(in) :: v(ncol,pver)    ! Midpoint meridional winds. ( m s-1)
   real(r8),         intent(in) :: t(ncol,pver)    ! Midpoint temperatures. (K)
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


   real(r8),         intent(in) :: effgw_rdg       ! Tendency efficiency.
   real(r8),         intent(in) :: effgw_rdg_max
   real(r8),         intent(in) :: effgw_rdg_resid  ! Tendency efficiency.
   logical,          intent(in) :: luse_gw_rdg_resid ! On-Off switch 
   real(r8),         intent(in) :: hwdth(ncol,prdg) ! width of ridges.
   real(r8),         intent(in) :: clngt(ncol,prdg) ! length of ridges.
   real(r8),         intent(in) :: gbxar(ncol)      ! gridbox area

   real(r8),         intent(in) :: mxdis(ncol,prdg) ! Height estimate for ridge (m).
   real(r8),         intent(in) :: angll(ncol,prdg) ! orientation of ridges.
   real(r8),         intent(in) :: anixy(ncol,prdg) ! Anisotropy parameter.

   real(r8),         intent(in) :: isovar(ncol)     ! sqrt of residual variance
   real(r8),         intent(in) :: isowgt(ncol)     ! area frac of residual variance

   real(r8),         intent(in) :: rdg_cd_llb      ! Drag coefficient for low-level flow
   logical,          intent(in) :: trpd_leewv

   !++jtb Using massively hacked physics_types
   type(physics_ptend), intent(inout):: ptend   ! Parameterization net tendencies. 

   real(r8),        intent(out) :: flx_heat(pcols)

   !---------------------------Local storage-------------------------------

   integer :: k, m, nn, istat

   real(r8), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(r8), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(r8), allocatable :: phase_speeds(:,:)

   ! Isotropic source flag [anisotropic orography].
   integer  :: isoflag(ncol)

   ! horiz wavenumber [anisotropic orography].
   real(r8) :: kwvrdg(ncol)

   ! Efficiency for a gravity wave source.
   real(r8) :: effgw(ncol)

   ! Indices of top gravity wave source level and lowest level where wind
   ! tendencies are allowed.
   integer :: src_level(ncol)
   integer :: tend_level(ncol)
   integer :: bwv_level(ncol)
   integer :: tlb_level(ncol)

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

   ! Top of low-level flow layer.
   real(r8) :: tlb(ncol)

   ! Bottom of linear wave region.
   real(r8) :: bwv(ncol)

   ! Froude numbers for flow/drag regimes
   real(r8) :: Fr1(ncol)
   real(r8) :: Fr2(ncol)
   real(r8) :: Frx(ncol)

   ! Wave Reynolds stresses at source level
   real(r8) :: tauoro(ncol)
   real(r8) :: taudsw(ncol)

   ! Surface streamline displacement height for linear waves.
   real(r8) :: hdspwv(ncol)

   ! Surface streamline displacement height for downslope wind regime.
   real(r8) :: hdspdw(ncol)

   ! Wave breaking level
   real(r8) :: wbr(ncol)

   real(r8) :: utgw(ncol,pver)       ! zonal wind tendency
   real(r8) :: vtgw(ncol,pver)       ! meridional wind tendency
   real(r8) :: ttgw(ncol,pver)       ! temperature tendency
   real(r8) :: qtgw(ncol,pver,pcnst) ! constituents tendencies

   ! Effective gravity wave diffusivity at interfaces.
   real(r8) :: egwdffi(ncol,pver+1)

   ! Temperature tendencies from diffusion and kinetic energy.
   real(r8) :: dttdf(ncol,pver)
   real(r8) :: dttke(ncol,pver)

   ! Wave stress in zonal/meridional direction
   real(r8) :: taurx(ncol,pver+1)
   real(r8) :: taurx0(ncol,pver+1)
   real(r8) :: taury(ncol,pver+1)
   real(r8) :: taury0(ncol,pver+1)
   ! Provisional absolute wave stress from gw_drag_prof
   real(r8) :: tau_diag(ncol,pver+1)

   ! U,V tendency accumulators
   real(r8) :: utrdg(ncol,pver)
   real(r8) :: vtrdg(ncol,pver)
   real(r8) :: ttrdg(ncol,pver)

   ! Energy change used by fixer.
   real(r8) :: de(ncol)

   character(len=1) :: cn
   character(len=9) :: fname(4)
   !----------------------------------------------------------------------------
   
   ! Allocate wavenumber fields.
   allocate(tau(ncol,band_oro%ngwv:band_oro%ngwv,pver+1),stat=istat)
   call alloc_err(istat,'gw_rdg_calc','tau',ncol*(band_oro%ngwv**2+1)*(pver+1))
   allocate(gwut(ncol,pver,band_oro%ngwv:band_oro%ngwv),stat=istat)
   call alloc_err(istat,'rdg_calc','gwut',ncol*pver*(band_oro%ngwv**2+1))
   allocate(phase_speeds(ncol,band_oro%ngwv:band_oro%ngwv),stat=istat)
   call alloc_err(istat,'rdg_calc','phase_speeds',ncol*(band_oro%ngwv**2+1))

   ! initialize accumulated momentum fluxes and tendencies
   taurx = 0._r8
   taury = 0._r8
   ttrdg = 0._r8
   utrdg = 0._r8
   vtrdg = 0._r8
   tau_diag = -9999._r8

   do nn = 1, n_rdg
      kwvrdg  = 0.001_r8 / ( hwdth(:,nn) + 0.001_r8 ) ! this cant be done every time step !!!
      isoflag = 0
      effgw   = effgw_rdg * ( hwdth(1:ncol,nn)* clngt(1:ncol,nn) ) / gbxar(1:ncol)
      effgw   = min( effgw_rdg_max , effgw )

      call gw_rdg_src(ncol, band_oro, p, &
         u, v, t, mxdis(:,nn), angll(:,nn), anixy(:,nn), kwvrdg, isoflag, zi, nm, &
         src_level, tend_level, bwv_level, tlb_level, tau, ubm, ubi, xv, yv,  &
         ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, phase_speeds)

      call gw_rdg_belowpeak(ncol, band_oro, rdg_cd_llb, &
         t, mxdis(:,nn), anixy(:,nn), kwvrdg, &
         zi, nm, ni, rhoi, &
         src_level, tau, &
         ubmsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, &
         tauoro, taudsw, hdspwv, hdspdw)

#ifdef UNITTEST
      if (nn == 1) then
         write(20) tau(:,0,:)
      end if
#endif

      call gw_rdg_break_trap(ncol, band_oro, &
         zi, nm, ni, ubm, ubi, rhoi, kwvrdg , bwv, tlb, wbr, &
         src_level, tlb_level, hdspwv, hdspdw,  mxdis(:,nn), &
         tauoro, taudsw, tau, &
         ldo_trapped_waves=trpd_leewv)


#ifdef UNITTEST
      if (nn == 1) then
         write(20) tau(:,0,:)
      end if
#endif
      
      call gw_drag_prof(ncol, band_oro, p, src_level, tend_level, dt, &
         t, vramp,    &
         piln, rhoi, nm, ni, ubm, ubi, xv, yv,   &
         effgw, phase_speeds, kvtt, q, dse, tau, utgw, vtgw, &
         ttgw, qtgw, egwdffi,   gwut, dttdf, dttke, &
         kwvrdg=kwvrdg, &
         satfac_in = 1._r8, lapply_vdiff=gw_rdg_do_vdiff , tau_diag=tau_diag )

      ! Add the tendencies from each ridge to the totals.
      do k = 1, pver
         ! diagnostics
         utrdg(:,k) = utrdg(:,k) + utgw(:,k)
         vtrdg(:,k) = vtrdg(:,k) + vtgw(:,k)
         ttrdg(:,k) = ttrdg(:,k) + ttgw(:,k)
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
         taurx0(:,k) =  tau(:,0,k)*xv
         taury0(:,k) =  tau(:,0,k)*yv
         taurx(:,k)  =  taurx(:,k) + taurx0(:,k)
         taury(:,k)  =  taury(:,k) + taury0(:,k)
      end do

#ifdef UNITTEST
      if (nn == 1) then
         write(20) bwv
         write(20) tlb
         write(20) wbr
         write(20) ubmsrc
         write(20) nsrc
         write(20) tauoro
         write(20) taudsw
         write(20) ubm
         write(20) tau(:,0,:)
         write(20) tau_diag
      end if
#else
      if (nn == 1) then
         call outfld('BWV_HT1', bwv,     ncol, lchnk)
         call outfld('TLB_HT1', tlb,     ncol, lchnk)
         call outfld('WBR_HT1', wbr,     ncol, lchnk)
         call outfld('TAUDSW1', taudsw,  ncol, lchnk)
         call outfld('TAUORO1', tauoro,  ncol, lchnk)
         call outfld('UBMSRC1', ubmsrc,  ncol, lchnk)
         call outfld('USRC1',   usrc,    ncol, lchnk)
         call outfld('VSRC1',   vsrc,    ncol, lchnk)
         call outfld('NSRC1'  , nsrc,    ncol, lchnk)
         ! Froude numbers
         call outfld('Fr1_DIAG' , Fr1,    ncol, lchnk)
         call outfld('Fr2_DIAG' , Fr2,    ncol, lchnk)
         call outfld('Frx_DIAG' , Frx,    ncol, lchnk)
         ! Ridge quantities - don't change.  Written for convenience
         call outfld('MXDIS1' , mxdis(:,nn) ,  ncol, lchnk)
         call outfld('ANGLL1' , angll(:,nn) ,  ncol, lchnk)
         call outfld('ANIXY1' , anixy(:,nn) ,  ncol, lchnk)
         call outfld('HWDTH1' , hwdth(:,nn) ,  ncol, lchnk)
         call outfld('CLNGT1' , clngt(:,nn) ,  ncol, lchnk)
         call outfld('GBXAR1' , gbxar ,        ncol, lchnk)
         call outfld('TAUM1_DIAG' , tau_diag ,  ncol, lchnk)
         call outfld('TAU1RDG'//trim(type)//'M', tau(:,0,:),  ncol, lchnk)
         call outfld('UBM1'//trim(type),         ubm,         ncol, lchnk)
         call outfld('UBT1RDG'//trim(type),      gwut,        ncol, lchnk)
      end if

      if (nn <= 6) then
         write(cn, '(i1)') nn
         call outfld('TAU'//cn//'RDG'//trim(type)//'X', taurx0,  ncol, lchnk)
         call outfld('TAU'//cn//'RDG'//trim(type)//'Y', taury0,  ncol, lchnk)
         call outfld('UT'//cn//'RDG'//trim(type),       utgw,    ncol, lchnk)
         call outfld('VT'//cn//'RDG'//trim(type),       vtgw,    ncol, lchnk)
      end if
#endif
      
   end do ! end of loop over multiple ridges

#ifndef UNITTEST
   call outfld('TAUARDG'//trim(type)//'X', taurx,  ncol, lchnk)
   call outfld('TAUARDG'//trim(type)//'Y', taury,  ncol, lchnk)
#endif
   
   !if (luse_gw_rdg_resid == .true.) then ! is this line a possible problem??
   if (luse_gw_rdg_resid .eqv. .true.) then ! gfortran wants this!
   ! Add additional GW from residual variance. Assumed isotropic
      kwvrdg  = 0.001_r8 / ( 100._r8 )
      effgw   = effgw_rdg_resid * isowgt    !1.0_r8 * isowgt
      tauoro = 0._r8

      call gw_rdg_resid_src(ncol, band_oro, p, &
         u, v, t, isovar, kwvrdg, zi, nm, &
         src_level, tend_level, tau, ubm, ubi, xv, yv,  &
         ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, phase_speeds, tauoro )

      call gw_drag_prof(ncol, band_oro, p, src_level, tend_level, dt, &
         t, vramp,    &
         piln, rhoi, nm, ni, ubm, ubi, xv, yv,   &
         effgw, phase_speeds, kvtt, q, dse, tau, utgw, vtgw, &
         ttgw, qtgw, egwdffi,   gwut, dttdf, dttke, &
         kwvrdg=kwvrdg, &
         satfac_in = 1._r8, lapply_vdiff=gw_rdg_do_vdiff , tau_diag=tau_diag )

      ! Add the tendencies from isotropic residual to the totals.
      do k = 1, pver
         ! diagnostics
         utrdg(:,k) = utrdg(:,k) + utgw(:,k)
         vtrdg(:,k) = vtrdg(:,k) + vtgw(:,k)
         ttrdg(:,k) = ttrdg(:,k) + ttgw(:,k)
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
         taurx0(:,k) =  tau(:,0,k)*xv
         taury0(:,k) =  tau(:,0,k)*yv
         taurx(:,k)  =  taurx(:,k) + taurx0(:,k)
         taury(:,k)  =  taury(:,k) + taury0(:,k)
      end do

#ifndef UNITTEST
      call outfld('UTGWORO_RESID', utgw,  ncol, lchnk)
      call outfld('VTGWORO_RESID', vtgw,  ncol, lchnk)

      call outfld('TAUDIAG_RESID', tau_diag,  ncol, lchnk)
      call outfld('TAUORO_RESID', tauoro ,  ncol, lchnk)
      call outfld('TAURESID'//trim(type)//'M', tau(:,0,:),  ncol, lchnk)
      call outfld('TAURESID'//trim(type)//'X', taurx,  ncol, lchnk)
      call outfld('TAURESID'//trim(type)//'Y', taury,  ncol, lchnk)

      call outfld('UBMRESID'//trim(type),      ubm,         ncol, lchnk)
      call outfld('UBIRESID'//trim(type),      ubi,         ncol, lchnk)
      call outfld('SRC_LEVEL_RESID'//trim(type),      1._r8*src_level ,         ncol, lchnk)
#endif
      ! end of residual variance calc
   end if

   ! Calculate energy change for output to CAM's energy checker.
   call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
          ptend%v(:ncol,:), ptend%s(:ncol,:), de)
   flx_heat(:ncol) = de

   if (trim(type) == 'BETA') then
      fname(1) = 'TAUGWX'
      fname(2) = 'TAUGWY'
      fname(3) = 'UTGWORO'
      fname(4) = 'VTGWORO'
   else if (trim(type) == 'GAMMA') then
      fname(1) = 'TAURDGGMX'
      fname(2) = 'TAURDGGMY'
      fname(3) = 'UTRDGGM'
      fname(4) = 'VTRDGGM'
   else
      call endrun('gw_rdg_calc: FATAL: type must be either BETA or GAMMA'&
                  //' type= '//type)
   end if

#ifndef UNITTEST
   call outfld(fname(1), taurx(:,pver+1), ncol, lchnk)
   call outfld(fname(2), taury(:,pver+1), ncol, lchnk)
   call outfld(fname(3), utrdg,  ncol, lchnk)
   call outfld(fname(4), vtrdg,  ncol, lchnk)
   call outfld('TTGWORO', ttrdg / cpair,  ncol, lchnk)
#endif
   

   deallocate(tau, gwut, phase_speeds)

end subroutine gw_rdg_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_band_rdg(band_oro_in)
  type(GWBand) , intent(in) :: band_oro_in
! Just to bring band_oro into module
band_oro=band_oro_in
write(*,*) " Band ORO ",band_oro%ngwv
end subroutine set_band_rdg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_vramp()

allocate(vramp(pver))
vramp(:) = 1._r8

end subroutine set_vramp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine report_from_within()

write(*,*) " Inside gw_rdg_calc_mod: Band ORO ",band_oro%ngwv
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
    
end module gw_rdg_calc_mod

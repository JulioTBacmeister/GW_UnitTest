module namelist_mod

use shr_kind_mod, only : r8 => shr_kind_r8

use namelist_utils,  only: find_group_name
use units,           only: getunit, freeunit
use cam_abortutils,  only: endrun

implicit none

! Declare all variables

! === cam_initfiles_nl ===
character(len=256) :: bnd_topo, ncdata_root, gw_drag_file, gw_drag_file_mm, ncdata_type
character(len=256) :: calculation_type,casename,ncout_root
real(8) :: scale_dry_air_mass

character(len=256) :: ncdata

! === gw_drag_nl ===
real(8) :: alpha_gw_movmtn, effgw_beres_dp, effgw_cm, effgw_movmtn_pbl
real(8) :: effgw_rdg_beta, effgw_rdg_beta_max, effgw_rdg_resid
real(8) :: front_gaussian_width, frontgfc, gw_dc, gw_dc_long
real(8) :: gw_oro_south_fac, gw_prndl, gw_qbo_hdepth_scaling
real(8) :: movmtn_plaunch, movmtn_psteer, rdg_beta_cd_llb, taubgnd
integer :: movmtn_source, n_rdg_beta, pgwv, pgwv_long
integer :: start_year, start_month, start_day
integer :: start_hour, nsteps
integer :: dt_step
logical :: use_topo_file, gw_apply_tndmax, gw_limit_tau_without_eff
logical :: gw_lndscl_sgh, gw_polar_taper, gw_top_taper
logical :: tau_0_ubc, trpd_leewv_rdg_beta, use_gw_rdg_beta
logical :: use_gw_rdg_gamma, use_gw_rdg_resid

! === phys_ctl_nl ===
character(len=256) :: cam_chempkg, cam_physpkg, deep_scheme, eddy_scheme
character(len=256) :: macrop_scheme, microp_scheme, radiation_scheme
character(len=256) :: shallow_scheme, waccmx_opt
integer :: cam_snapshot_after_num, cam_snapshot_before_num, cld_macmic_num_steps
integer :: srf_flux_avg
logical :: convproc_do_aer, do_clubb_sgs, do_hb_above_clubb
logical :: history_aero_optics, history_aerosol, history_amwg
logical :: history_budget, history_chemistry, history_chemspecies_srf
logical :: history_clubb, history_dust, history_eddy, history_vdiag
logical :: history_waccm, history_waccmx, use_gw_convect_dp
logical :: use_gw_convect_sh, use_gw_front, use_gw_front_igw
logical :: use_gw_movmtn_pbl, use_gw_oro, use_hemco
logical :: use_hetfrz_classnuc, use_simple_phys, use_subcol_microp

public :: drv_readnl
public :: atm_readnl

! Single public line for all variables
public :: bnd_topo, ncdata, scale_dry_air_mass, use_topo_file,ncdata_type, &
     alpha_gw_movmtn, effgw_beres_dp, effgw_cm, effgw_movmtn_pbl, &
     effgw_rdg_beta, effgw_rdg_beta_max, effgw_rdg_resid, front_gaussian_width, &
     frontgfc, gw_apply_tndmax, gw_dc, gw_dc_long, gw_drag_file, &
     gw_drag_file_mm, gw_limit_tau_without_eff, gw_lndscl_sgh, gw_oro_south_fac, &
     gw_polar_taper, gw_prndl, gw_qbo_hdepth_scaling, gw_top_taper, &
     movmtn_plaunch, movmtn_psteer, movmtn_source, n_rdg_beta, pgwv, &
     pgwv_long, rdg_beta_cd_llb, tau_0_ubc, taubgnd, trpd_leewv_rdg_beta, &
     use_gw_rdg_beta, use_gw_rdg_gamma, use_gw_rdg_resid, &
     cam_chempkg, cam_physpkg, cam_snapshot_after_num, cam_snapshot_before_num, &
     cld_macmic_num_steps, convproc_do_aer, deep_scheme, do_clubb_sgs, &
     do_hb_above_clubb, eddy_scheme, history_aero_optics, history_aerosol, &
     history_amwg, history_budget, history_chemistry, history_chemspecies_srf, &
     history_clubb, history_dust, history_eddy, history_vdiag, history_waccm, &
     history_waccmx, macrop_scheme, microp_scheme, radiation_scheme, &
     shallow_scheme, srf_flux_avg, use_gw_convect_dp, use_gw_convect_sh, &
     use_gw_front, use_gw_front_igw, use_gw_movmtn_pbl, use_gw_oro, &
     use_hemco, use_hetfrz_classnuc, use_simple_phys, use_subcol_microp, &
     waccmx_opt, calculation_type, ncdata_root, ncout_root, casename

!==================================================
contains
!==================================================

subroutine atm_readnl(nlfile ) !, bnd_topo, ncdata )

  ! File containing namelist input.
  character(len=*), intent(in) :: nlfile

  ! File names for code:
  ! character(len=*), intent(out) :: bnd_topo, ncdata

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'gen_readnl'

  logical :: use_topo_file
  real(r8) ::  scale_dry_air_mass
  
  namelist /cam_initfiles_nl_camsnap/ bnd_topo, ncdata_root, scale_dry_air_mass, use_topo_file, ncdata_type
  namelist /cam_initfiles_nl_ERA5/ bnd_topo, ncdata_root, scale_dry_air_mass, use_topo_file, ncdata_type
  namelist /cam_initfiles_nl_xy/ bnd_topo, ncdata_root, scale_dry_air_mass, use_topo_file, ncdata_type
  namelist /cam_initfiles_nl_xympas/ bnd_topo, ncdata_root, scale_dry_air_mass, use_topo_file, ncdata_type

  namelist /gw_drag_nl/ alpha_gw_movmtn, effgw_beres_dp, effgw_cm, effgw_movmtn_pbl, &
     effgw_rdg_beta, effgw_rdg_beta_max, effgw_rdg_resid, front_gaussian_width, &
     frontgfc, gw_apply_tndmax, gw_dc, gw_dc_long, gw_drag_file, &
     gw_drag_file_mm, gw_limit_tau_without_eff, gw_lndscl_sgh, gw_oro_south_fac, &
     gw_polar_taper, gw_prndl, gw_qbo_hdepth_scaling, gw_top_taper, &
     movmtn_plaunch, movmtn_psteer, movmtn_source, n_rdg_beta, pgwv, &
     pgwv_long, rdg_beta_cd_llb, tau_0_ubc, taubgnd, trpd_leewv_rdg_beta, &
     use_gw_rdg_beta, use_gw_rdg_gamma, use_gw_rdg_resid

  namelist /phys_ctl_nl/ cam_chempkg, cam_physpkg, &
     cam_snapshot_after_num, cam_snapshot_before_num, &
     cld_macmic_num_steps, convproc_do_aer, deep_scheme, do_clubb_sgs, &
     do_hb_above_clubb, eddy_scheme, history_aero_optics, history_aerosol, &
     history_amwg, history_budget, history_chemistry, history_chemspecies_srf, &
     history_clubb, history_dust, history_eddy, history_vdiag, history_waccm, &
     history_waccmx, macrop_scheme, microp_scheme, radiation_scheme, &
     shallow_scheme, srf_flux_avg, use_gw_convect_dp, use_gw_convect_sh, &
     use_gw_front, use_gw_front_igw, use_gw_movmtn_pbl, use_gw_oro, &
     use_hemco, use_hetfrz_classnuc, use_simple_phys, use_subcol_microp, &
     waccmx_opt
  

  unitn = getunit()
  open( unitn, file=trim(nlfile), status='old' )
  
  if ( trim(calculation_type) == 'camsnap') then
     call find_group_name(unitn, 'cam_initfiles_nl_camsnap', status=ierr)
     if (ierr == 0) then
        read(unitn, cam_initfiles_nl_camsnap, iostat=ierr)
        if (ierr /= 0) then
           call endrun(' ERROR reading namelist cam_initfiles')
        end if
     end if
  else if ( trim(calculation_type) == 'ERA5') then
     call find_group_name(unitn, 'cam_initfiles_nl_ERA5', status=ierr)
     if (ierr == 0) then
        read(unitn, cam_initfiles_nl_ERA5, iostat=ierr)
        if (ierr /= 0) then
           call endrun(' ERROR reading namelist cam_initfiles')
        end if
     end if
  else if ( trim(calculation_type) == 'xy') then
     call find_group_name(unitn, 'cam_initfiles_nl_xy', status=ierr)
     if (ierr == 0) then
        read(unitn, cam_initfiles_nl_xy, iostat=ierr)
        if (ierr /= 0) then
           call endrun(' ERROR reading namelist cam_initfiles')
        end if
     end if
  else if ( trim(calculation_type) == 'xympas') then
     call find_group_name(unitn, 'cam_initfiles_nl_xympas', status=ierr)
     if (ierr == 0) then
        read(unitn, cam_initfiles_nl_xympas, iostat=ierr)
        if (ierr /= 0) then
           call endrun(' ERROR reading namelist cam_initfiles')
        end if
     end if
  else
     STOP "NO vaid caculation type given"
  end if

     
  call find_group_name(unitn, 'gw_drag_nl', status=ierr)
  if (ierr == 0) then
     read(unitn, gw_drag_nl, iostat=ierr)
     if (ierr /= 0) then
        call endrun(' ERROR reading namelist gw_drag')
     end if
  end if
  call find_group_name(unitn, 'phys_ctl_nl', status=ierr)
  if (ierr == 0) then
     read(unitn, phys_ctl_nl, iostat=ierr)
     if (ierr /= 0) then
        call endrun(' ERROR reading namelist phys_ctl')
     end if
  end if
  close(unitn)
  call freeunit(unitn)
end subroutine atm_readnl

subroutine drv_readnl(nlfile ) !, bnd_topo, ncdata )

  ! File containing namelist input.
  character(len=*), intent(in) :: nlfile

  ! File names for code:
  ! character(len=*), intent(out) :: bnd_topo, ncdata

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: sub = 'gen_readnl'

  logical :: use_topo_file
  real(r8) ::  scale_dry_air_mass

  namelist /top_ctl_nl/ calculation_type,start_year, start_month, start_day, start_hour, nsteps, dt_step, &
       ncout_root,casename


  unitn = getunit()
  open( unitn, file=trim(nlfile), status='old' )

  call find_group_name(unitn, 'top_ctl_nl', status=ierr)
  if (ierr == 0) then
     read(unitn, top_ctl_nl, iostat=ierr)
     if (ierr /= 0) then
        call endrun(' ERROR reading namelist top level control')
     else
        write(*,*) "Read this calc type:", calculation_type
     end if
  end if

  close(unitn)
  call freeunit(unitn)
end subroutine drv_readnl
  
end module namelist_mod

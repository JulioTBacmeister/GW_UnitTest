  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   Huge messy import of all likely namelist params
  !   and then some, ... a lot more and then some
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use namelist_mod, only: genl_readnl
  ! --- cam_initfiles
  use namelist_mod, only: bnd_topo, ncdata, scale_dry_air_mass, use_topo_file
  ! --- gw_drag
  use namelist_mod, only: alpha_gw_movmtn, effgw_beres_dp, effgw_cm, effgw_movmtn_pbl, &
     effgw_rdg_beta, effgw_rdg_beta_max, effgw_rdg_resid, front_gaussian_width, &
     frontgfc, gw_apply_tndmax, gw_dc, gw_dc_long, gw_drag_file, &
     gw_drag_file_mm, gw_limit_tau_without_eff, gw_lndscl_sgh, gw_oro_south_fac, &
     gw_polar_taper, gw_prndl, gw_qbo_hdepth_scaling, gw_top_taper, &
     movmtn_plaunch, movmtn_psteer, movmtn_source, n_rdg_beta, pgwv, &
     pgwv_long, rdg_beta_cd_llb, tau_0_ubc, taubgnd, trpd_leewv_rdg_beta, &
     use_gw_rdg_beta, use_gw_rdg_gamma, use_gw_rdg_resid
   !-- phys_ctl
   use namelist_mod, only: cam_chempkg, cam_physpkg, &
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
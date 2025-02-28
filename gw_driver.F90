program gw_driver

  ! Uses
  
  use gw_utils, only: r8
  use gw_rdg
  use gw_common, only: gw_prof, gw_drag_prof, gw_common_init, GWband
  use ncread_mod
  use physconst, only: cpair, rair, gravit, zvir
  use coords_1d,  only: Coords1D
  use utils_mod, only: make_pressures, leovy_alpha !, make_geopht
  use physics_types,  only: physics_ptend,physics_ptend_init
  use ppgrid, only: pcnst,pcols, pver_in_ppgrid=>pver
  use geopotential, only: geopotential_t

  ! set_band_rdg required to get band construct into
  ! gw_rdg_calc_mod module.  
  use gw_rdg_calc_mod, only : gw_rdg_calc, set_band_rdg, set_vramp, report_from_within

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   Huge messy import of all likely namelist params
  !   hidden here
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  include 'namelist_imports.inc'

  implicit none
   
  ! Top level for gravity waves.
  integer, parameter :: ktop = 1

  ! File containing namelist input.
  character(len=80) :: nlfile
    
  ! Reporting of errors
  character(len=128) :: errstring
  
  ! Type of forcing data
  !character(len=128) :: ncdata_type

  integer ncid,status, dimid,varid  ! for netCDF data file

  integer :: ncol, nrdgs, pver, ntim, itime

  ! Ridge parameters
  real(r8)  , allocatable :: mxdis(:,:), angll(:,:), aniso(:,:), anixy(:,:), hwdth(:,:), clngt(:,:)
  real(r8)  , allocatable :: gbxar(:), isowgt(:), isovar(:), sgh(:)

  ! Meteorolgical data  and dims
  real(r8) , allocatable :: hyai(:), hybi(:), hyam(:), hybm(:)
  real(r8) , allocatable :: lon(:), lat(:)
  real(r8) , allocatable :: U(:,:,:), V(:,:,:), T(:,:,:), Q(:,:,:), PS(:,:)
  real(r8) , allocatable :: nm_(:,:,:), ni_(:,:,:), zm_(:,:,:), zi_(:,:,:), rhoi_(:,:,:)

  ! Additional fields that will be needed
  ! (Direct ... have dummy time dim)
  real(r8) , allocatable :: pint(:,:,:), pmid(:,:,:), piln(:,:,:), pmln(:,:,:)

  ! Additional fields
  ! (No dummy time dim)
  real(r8) , allocatable :: rhoi(:,:), nm(:,:), ni(:,:), zm(:,:), zi(:,:), dse(:,:)
  real(r8) , allocatable :: rair3D(:,:), zvir3D(:,:)

  ! Other Additional fields
  real(r8) , allocatable :: alpha(:), pref_edge(:)
  real(r8) , allocatable :: kvtt(:,:),flx_heat(:)

  real(r8) ::  prndl,dt
  
  ! Horrible, pointless construct 
  type(Coords1D) :: p               ! Pressure coordinates

  ! ptend derived typ ...
  type(physics_ptend) :: ptend

  ! A mid-scale "band" with only stationary waves (l = 0).
  type(GWBand) :: band_oro
  ! Horzontal wavelengths for bands [m].
  real(r8), parameter :: wavelength_mid = 1.e5_r8
  real(r8), parameter :: wavelength_long = 1.e6_r8
  ! fcrit2 for the mid-scale waves has been made a namelist variable to
  ! facilitate backwards compatibility with the CAM3 version of this
  ! parameterization.  In CAM3, fcrit2=0.5.
  real(r8) :: fcrit2 = 1.0   ! critical froude number squared

  logical   , allocatable :: lqptend(:)
  logical :: trpd_leewave
  integer :: lchnk, n_rdg

  nlfile='atm_in'
  !ncdata_type = 'ERA5_SE_IC'
  !ncdata_type = 'camsnap'

  write( *,* ) "Hello ... "
  call gw_rdg_readnl( nlfile )
  call genl_readnl( nlfile ) !, bnd_topo, ncdata )

  write(*,*) "ncdata :",ncdata
  write(*,*) "ncdata_type :",ncdata_type
  write(*,*) "bnd_topo :",bnd_topo
  write(*,*) " n_rdg_beta: ", n_rdg_beta
  write(*,*) " trpd_leewv_rdg_beta: ",  trpd_leewv_rdg_beta

  
  call ncread_topo( bnd_topo , mxdis, angll, aniso, anixy, hwdth, clngt, gbxar, isovar, isowgt, sgh )


  !! erafile = '/glade/campaign/cgd/amp/juliob/ERA5/ne30pg3/L93/2014/ERA5_x_ne30pg3_L93_rgC2_WO.2014-01-15-00000.nc'

  if ( trim(ncdata_type) == 'ERA5_SE_IC') then
     write(*,*) "ERA5 IC data"
     call ncread_era5_se_ic( ncdata , ncol, pver , ntim, &
          hyai , hybi , hyam , hybm , lon , lat , &
          PS, U , V , T  , Q )
  else if ( trim(ncdata_type) == 'camsnap') then
     write(*,*) "Wahhooo"
     call ncread_camsnap( ncdata , ncol, pver , ntim, &
          hyai , hybi , hyam , hybm , lon , lat , &
          PS, U , V , T  , Q , &
          zm_, zi_, nm_, ni_, rhoi_ )
  else
     write(*,*) "Wahhooo"
  end if

  pcols = ncol ! does this actaully go back to ppgrid? Yes.
  pver_in_ppgrid = pver

  itime=1

  write( *,* ) "Hello ... Back in driver"
  write(*,*) "Got ","PS  ",minval(PS), maxval(PS)  
  write(*,*) "Got ","U   ",minval(U), maxval(U)
  write(*,*) "Got ","V   ",minval(V), maxval(V)  
  write(*,*) "Got ","T   ",minval(T), maxval(T)
  write(*,*) "Got ","Q   ",minval(Q), maxval(Q)

  !allocate( pmid(ncol,pver) , pint(ncol, pver+1) )
  
  call make_pressures( ncol, pver, ntim, hyai, hybi, hyam, hybm, PS, pint, pmid )

  write(*,*) "Made ","pint   ",minval(pint), maxval(pint), shape(pint)
  write(*,*) "Made ","pmid   ",minval(pmid), maxval(pmid), shape(pmid)

  p = Coords1D(pint(:ncol,:,itime))
  write(*,*) "P%ifc "
  write(*,*) minval( p%ifc ),maxval( p%ifc ), shape(p%ifc)

  ! Will need these also
  allocate( piln(ncol,pver+1,ntim), pmln(ncol,pver,ntim) )
  piln = log( pint )
  pmln = log( pmid )

  allocate( pref_edge(pver+1)  )
  pref_edge(:) = 100000. * ( hyai(:) + hybi(:) )

  call leovy_alpha( pver, pref_edge, alpha )
  write(*,*) "Made ","alpha   ",minval(alpha), maxval(alpha), shape(alpha)

  call gw_common_init(pver,&
       tau_0_ubc, ktop, gravit, rair, alpha, & 
       gw_prndl, gw_qbo_hdepth_scaling, & 
       errstring)

  allocate( rhoi(ncol,pver+1) , ni(ncol, pver+1), nm(ncol, pver), dse(ncol,pver) )
  allocate( kvtt(ncol,pver+1) , flx_heat(ncol) )


  ! Profiles of background state variables
  call gw_prof(ncol, p, cpair, T(:,:,itime), rhoi, nm, ni)

  write(*,*) "Made ","nm     ",minval(nm), maxval(nm), shape(nm)
  write(*,*) "Made ","ni     ",minval(ni), maxval(ni), shape(ni)
  write(*,*) "Made ","rhoi   ",minval(rhoi), maxval(rhoi), shape(rhoi)

  if ( trim(ncdata_type) == 'camsnap') then
     write(*,*) "Override w/ camsnap "
     allocate( zm(ncol,pver) , zi(ncol,pver+1) )
     zm = zm_(:,:,itime)
     zi = zi_(:,:,itime)
     nm = nm_(:,:,itime)
     ni = ni_(:,:,itime)
     rhoi = rhoi_(:,:,itime)
     
  end if
  ! At this point if zm and zi have been allocated it means they have
  ! read in or calculated
  if ( (.not.allocated(zm)) .and. (.not.allocated(zi)) ) then
     write(*,*) "NEED GOEPHT"
     allocate( zm(ncol,pver) , zi(ncol,pver+1), rair3D(ncol,pver), zvir3D(ncol,pver) )
     rair3D(:,:) = rair
     zvir3D(:,:) = zvir
     call geopotential_t( ncol , pver ,                   &
       piln(:,:,itime)   , pmln(:,:,itime)   , &
       pint(:,:,itime)   , pmid(:,:,itime)   , P%del   , P%rdel  , &
       T(:,:,itime) , Q(:,:,itime)  , rair3D   , gravit , zvir3D   ,   &
       zi     , zm     )

     write(*,*)  "Made ","ZI   ",minval(zi), maxval(zi), shape(zi)
  end if

  dse = cpair*T(:,:,itime) + gravit*zm
  
  allocate( lqptend(pcnst))
  lqptend(:)=.TRUE.
  call physics_ptend_init(ptend, ncol , 'OGWD',.TRUE.,.TRUE.,.TRUE., lqptend )

  band_oro = GWBand(0, gw_dc, fcrit2, wavelength_mid)

  call set_vramp()
  call report_from_within()

  ! Write multiple records
  open(unit=20, file='GW.dat', form='unformatted', access='stream', status='replace', action='write')

  write(20) ncol,pver
  write(20) zm
  write(20) zi

  lchnk=1

#if 1
  ! Doing the calculation
  n_rdg=1
  dt = 86400._r8 /48._r8 
  call gw_rdg_calc( &
   'BETA ', ncol, lchnk, n_rdg, dt, &
   U(:,:,itime) , V(:,:,itime) , T(:,:,itime) , &
   P, piln(:,:,itime) , zm, zi, &
   nm, ni, rhoi, kvtt, Q(:,:,:), dse, &
   effgw_rdg_beta, effgw_rdg_beta_max, &
   effgw_rdg_resid, use_gw_rdg_resid, &
   hwdth, clngt, gbxar, &
   mxdis, angll, anixy, &
   isovar, isowgt, &
   rdg_beta_cd_llb, trpd_leewv_rdg_beta, &
   !++jtb  ptend defined in massively hacked physics_types
   ptend, flx_heat) 

#endif

  
end program gw_driver

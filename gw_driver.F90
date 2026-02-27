program gw_driver

  ! Uses
  
  use gw_utils, only: r8
  use gw_rdg
  use gw_common, only: gw_prof, gw_drag_prof, gw_common_init, GWband
  use ncread_mod
  use nc_flexout_mod
  use physconst, only: cpair, rair, gravit, zvir
  use coords_1d,  only: Coords1D
  use utils_mod, only: make_pressures, leovy_alpha, make_nc_filename
  use physics_types,  only: physics_ptend,physics_ptend_init
  use ppgrid, only: pcnst,pcols, pver_in_ppgrid=>pver
  use geopotential, only: geopotential_t
  use dt_bump, only: bump_ymdh
  use sphere_ops, only: sphere_fronto_wrapper

  ! set_band_rdg required to get band construct into
  use gw_rdg_calc_mod, only : gw_rdg_calc, set_band_rdg
  use gw_rdg_calc_mod, only : set_vramp_rdg => set_vramp
  use gw_rdg_calc_mod, only : report_from_within_rdg => report_from_within
  ! set_band_movmtn required to get band construct into
  use gw_movmtn_calc_mod, only : gw_movmtn_calc, set_band_movmtn
  use gw_movmtn_calc_mod, only : set_vramp_movmtn => set_vramp
  use gw_movmtn_calc_mod, only : report_from_within_movmtn => report_from_within
  ! set_band_movmtn required to get band construct into
  use gw_front_calc_mod, only : gw_front_calc, set_band_front
  use gw_front_calc_mod, only : set_vramp_front => set_vramp
  use gw_front_calc_mod, only : report_from_within_front => report_from_within

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
  character(len=256)  :: ncout
    
  ! Reporting of errors
  character(len=128) :: errstring
  
  ! Type of forcing data
  !character(len=128) :: ncdata_type

  integer ncid,status, dimid,varid  ! for netCDF data file

  integer :: ncol, nrdgs, pver, ntim, itime, nx, ny

  ! Ridge parameters
  real(r8)  , allocatable :: mxdis(:,:), angll(:,:), aniso(:,:), anixy(:,:), hwdth(:,:), clngt(:,:)
  real(r8)  , allocatable :: gbxar(:), isowgt(:), isovar(:), sgh(:)

  ! Meteorolgical data  and dims
  real(r8) , allocatable :: hyai(:), hybi(:), hyam(:), hybm(:)
  real(r8) , allocatable :: lon(:), lat(:), lon_R(:), lat_R(:)
  real(r8) , allocatable :: U(:,:,:), V(:,:,:), T(:,:,:), Q(:,:,:), ZETA(:,:,:), PS(:,:)
  real(r8) , allocatable :: nm_(:,:,:), ni_(:,:,:), zm_(:,:,:), zi_(:,:,:), rhoi_(:,:,:)
  real(r8) , allocatable :: pint_(:,:,:), piln_(:,:,:)

  ! Additional fields that will be needed
  ! (Direct ... have dummy time dim)
  real(r8) , allocatable :: pint(:,:,:), pmid(:,:,:), piln(:,:,:), pmln(:,:,:)
  real(r8) , allocatable :: TH(:,:,:), FRONTGF(:,:,:), FRONTGA(:,:,:), zgrid(:,:)

  ! Additional fields
  ! (No dummy time dim)
  real(r8) , allocatable :: rhoi(:,:), nm(:,:), ni(:,:), zm(:,:), zi(:,:), dse(:,:)
  real(r8) , allocatable :: rair3D(:,:), zvir3D(:,:)

  ! Other Additional fields
  real(r8) , allocatable :: alpha(:), pref_edge(:), ilev(:)
  real(r8) , allocatable :: kvtt(:,:),flx_heat(:)

  real(r8) ::  prndl,dt,time_val
  
  ! Horrible, pointless construct 
  type(Coords1D) :: PP               ! Pressure coordinates

  ! ptend derived typ ...
  type(physics_ptend) :: ptend

  ! Gravity wave "bands".
  type(GWBand) :: band_oro , band_movmtn , band_front
  ! Horzontal wavelengths for bands [m].
  real(r8), parameter :: wavelength_mid = 1.e5_r8
  real(r8), parameter :: wavelength_long = 1.e6_r8
  ! fcrit2 for the mid-scale waves has been made a namelist variable to
  ! facilitate backwards compatibility with the CAM3 version of this
  ! parameterization.  In CAM3, fcrit2=0.5.
  real(r8) :: fcrit2 = 1.0   ! critical froude number squared

  logical   , allocatable :: lqptend(:)
  logical :: trpd_leewave,llatlon,ldo_rdg
  integer :: lchnk, n_rdg, i,j,k,n
  integer :: yy,mm,dd,hh,ss,ii
  
  write( *,* ) "Hello ... "

  nlfile='drv_in'
  call drv_readnl( nlfile )

  nlfile='atm_in'
  call gw_rdg_readnl( nlfile )
  call atm_readnl( nlfile ) !, bnd_topo, ncdata )
  
  !ncout_root = &
  !      '/glade/derecho/scratch/juliob/GW_UnitTest/'//trim(calculation_type)// &
  !      '.h'
  ncout_root = &
        trim(ncout_root)//trim(casename)//'/'//trim(casename)//'.h'


  
  
  write(*,*) "scale_dry_air_mass : ",scale_dry_air_mass
  write(*,*) "use_topo_file : ",use_topo_file

  write(*,*) "calculation_type :",calculation_type
  write(*,*) "ncdata_root :",ncdata_root
  write(*,*) "ncdata_type :",ncdata_type
  write(*,*) "ncout_root :",ncout_root
  write(*,*) "bnd_topo :",bnd_topo
  write(*,*) " n_rdg_beta: ", n_rdg_beta
  write(*,*) " trpd_leewv_rdg_beta: ",  trpd_leewv_rdg_beta
  write(*,*) " .... "
  write(*,*) " Some non-orog parameters. "
  write(*,*) " movmtn_plaunch, movmtn_psteer, movmtn_source: ",movmtn_plaunch, movmtn_psteer, movmtn_source
  write(*,*) " alpha_gw_movmtn: ",alpha_gw_movmtn
  write(*,*) " use_gw_movmtn_pbl: ",use_gw_movmtn_pbl
  write(*,*) " use_gw_rdg_beta: ",use_gw_rdg_beta
  write(*,*) " use_gw_front: ",use_gw_front
  
  if ( ( trim(ncdata_type) == 'XY_DATA') .or. (trim(ncdata_type) == 'XYMPAS_DATA') ) then
     llatlon=.TRUE.
     call ncread_topo_latlon( bnd_topo , mxdis, angll, aniso, anixy, hwdth, clngt, gbxar, isovar, isowgt, sgh  )
  else
     llatlon=.FALSE.
     call ncread_topo( bnd_topo , mxdis, angll, aniso, anixy, hwdth, clngt, gbxar, isovar, isowgt, sgh )
  end if
  
  yy=start_year; mm=start_month; dd=start_day; hh=start_hour
  do ii=1,nsteps
     ss=hh*3600
     print '(I4,"-",I2.2,"-",I2.2," ",I2.2,":00")', yy,mm,dd,hh
     ncdata = make_nc_filename(ncdata_root, yy, mm, dd, ss)
     ncout  = make_nc_filename(ncout_root,  yy, mm, dd, ss)
     print '(A160)',trim(ncdata)
     print '(A160)',trim(ncout)
  
     !!call ncread_topo( bnd_topo , mxdis, angll, aniso, anixy, hwdth, clngt, gbxar, isovar, isowgt, sgh, llatlon=llatlon )
     
     
     itime=1
     time_val = yy * 10000._r8 + mm*100._r8 + dd*1._r8 + hh/24._r8

     if ( trim(ncdata_type) == 'ERA5_SE_IC') then
        write(*,*) "ERA5 IC data"
        call ncread_era5_se_ic( ncdata , ncol, pver , ntim, &
             hyai , hybi , hyam , hybm , lon , lat , &
             PS, U , V , T  , Q )

        pcols = ncol ! does this actaully go back to ppgrid? Yes.
        pver_in_ppgrid = pver

        call make_pressures( ncol, pver, ntim, hyai, hybi, hyam, hybm, PS, pint, pmid )

        write(*,*) "Made ","pint   ",minval(pint), maxval(pint), shape(pint)
        write(*,*) "Made ","pmid   ",minval(pmid), maxval(pmid), shape(pmid)

        ! Will need these also
        allocate( piln(ncol,pver+1,ntim), pmln(ncol,pver,ntim) )
        piln = log( pint )
        pmln = log( pmid )

        PP = Coords1D(pint(:ncol,:,itime))
        write(*,*) "PP%ifc "
        write(*,*) minval( PP%ifc ),maxval( PP%ifc ), shape(PP%ifc)


        allocate( pref_edge(pver+1)  )
        pref_edge(:) = 100000. * ( hyai(:) + hybi(:) )

        call leovy_alpha( pver, pref_edge, alpha )
        write(*,*) "Made ","alpha   ",minval(alpha), maxval(alpha), shape(alpha)

        call gw_common_init(pver,&
             tau_0_ubc, ktop, gravit, rair, alpha, & 
             gw_prndl, gw_qbo_hdepth_scaling, & 
             errstring)


        allocate( rhoi(ncol,pver+1) , ni(ncol, pver+1), nm(ncol, pver)  )
        ! Profiles of background state variables
        call gw_prof(ncol, PP, cpair, T(:,:,itime), rhoi, nm, ni)

        write(*,*) "NEED GOEPHT"
        allocate( zm(ncol,pver) , zi(ncol,pver+1), rair3D(ncol,pver), zvir3D(ncol,pver) )
        rair3D(:,:) = rair
        zvir3D(:,:) = zvir
        call geopotential_t( ncol , pver ,                   &
             piln(:,:,itime)   , pmln(:,:,itime)   , &
             pint(:,:,itime)   , pmid(:,:,itime)   , PP%del   , PP%rdel  , &
             T(:,:,itime) , Q(:,:,itime)  , rair3D   , gravit , zvir3D   ,   &
             zi     , zm     )

        write(*,*)  "Made ","ZI   ",minval(zi), maxval(zi), shape(zi)

        !-----------------------------
        ! Need to figure ZETA for ERA5
        ! Set dummy var=0 for now
        !------------------------------
        allocate( ZETA(ncol,pver,ntim)  )
        ZETA = 0._r8

     else if ( trim(ncdata_type) == 'XY_DATA') then
        write(*,*) "XY (regirdded) data"
        call ncread_xy_data( ncdata , ncol, ny, nx, pver , ntim, &
             hyai , hybi , hyam , hybm , lon , lat , &
             PS, U , V , T, Q, ZETA, lat_R, lon_R )

        pcols = ncol ! does this actaully go back to ppgrid? Yes.
        pver_in_ppgrid = pver

        call make_pressures( ncol, pver, ntim, hyai, hybi, hyam, hybm, PS, pint, pmid )

        write(*,*) "Made ","pint   ",minval(pint), maxval(pint), shape(pint)
        write(*,*) "Made ","pmid   ",minval(pmid), maxval(pmid), shape(pmid)

        ! Will need these also
        allocate( piln(ncol,pver+1,ntim), pmln(ncol,pver,ntim) )
        piln = log( pint )
        pmln = log( pmid )

        PP = Coords1D(pint(:ncol,:,itime))
        write(*,*) "PP%ifc "
        write(*,*) minval( PP%ifc ),maxval( PP%ifc ), shape(PP%ifc)


        allocate( pref_edge(pver+1)  )
        pref_edge(:) = 100000. * ( hyai(:) + hybi(:) )

        call leovy_alpha( pver, pref_edge, alpha )
        write(*,*) "Made ","alpha   ",minval(alpha), maxval(alpha), shape(alpha)

        call gw_common_init(pver,&
             tau_0_ubc, ktop, gravit, rair, alpha, & 
             gw_prndl, gw_qbo_hdepth_scaling, & 
             errstring)


        allocate( rhoi(ncol,pver+1) , ni(ncol, pver+1), nm(ncol, pver)  )
        ! Profiles of background state variables
        call gw_prof(ncol, PP, cpair, T(:,:,itime), rhoi, nm, ni)

        write(*,*) "NEED GOEPHT"
        allocate( zm(ncol,pver) , zi(ncol,pver+1), rair3D(ncol,pver), zvir3D(ncol,pver) )
        rair3D(:,:) = rair
        zvir3D(:,:) = zvir
        call geopotential_t( ncol , pver ,                   &
             piln(:,:,itime)   , pmln(:,:,itime)   , &
             pint(:,:,itime)   , pmid(:,:,itime)   , PP%del   , PP%rdel  , &
             T(:,:,itime) , Q(:,:,itime)  , rair3D   , gravit , zvir3D   ,   &
             zi     , zm     )

        allocate( TH(ncol,pver,ntim) ,  FRONTGF(ncol,pver,ntim) ,  FRONTGA(ncol,pver,ntim)  )

        write(*,*)  "Made ","ZI   ",minval(zi), maxval(zi), shape(zi)

        !TH(:,:,itime) = T(:,:,itime) * ( (100000._r8 / pmid)**(rair/cpair) )   ! + gravit * zm / cpair )
        TH  = T * ( (100000._r8 / pmid)**(rair/cpair) )   ! + gravit * zm / cpair )

        call sphere_fronto_wrapper( TH, U, V, lat_R, lon_R, FRONTGF, FRONTGA )

     else if ( trim(ncdata_type) == 'XYMPAS_DATA') then
        write(*,*) "XY (regirdded) data from MPAS"
        call ncread_xympas_data( ncdata , ncol, ny, nx, pver , ntim, &
             lon , lat , ilev, zgrid, &
             pint, U , V , T, Q, ZETA, lat_R, lon_R )

        write(*,*) " shape of ilev", shape(ilev)
        
        ! Will need these also, later
        allocate( piln(ncol,pver+1,ntim), pmln(ncol,pver,ntim) , pmid(ncol,pver,ntim) )
        ! Redundanf (w/ PINT) but saves a hassle
        allocate( PS(ncol,ntim) )

        pcols = ncol ! does this actaully go back to ppgrid? Yes.
        pver_in_ppgrid = pver

        
        do n=1,ntim
           PS(:,n) = pint(:,pver+1,n)
           do k=1,pver
              pmid(:,k,n) = 0.5*( pint(:,k+1,n) + pint(:,k,n) )
           end do
        end do

        write(*,*) "Made ","pint   ",minval(pint), maxval(pint), shape(pint)
        write(*,*) "Made ","pmid   ",minval(pmid), maxval(pmid), shape(pmid)


        piln = log( pint )
        pmln = log( pmid )

        PP = Coords1D(pint(:ncol,:,itime))
        write(*,*) "PP%ifc "
        write(*,*) minval( PP%ifc ),maxval( PP%ifc ), shape(PP%ifc)

        allocate( pref_edge(pver+1)  )
        pref_edge(:) = 100000. * exp( -ilev / 7000. )
        write(*,*) "Made ","pref_edge   ",minval(pref_edge), maxval(pref_edge)

        call leovy_alpha( pver, pref_edge, alpha )
        write(*,*) "Made ","alpha   ",minval(alpha), maxval(alpha), shape(alpha)

        call gw_common_init(pver,&
             tau_0_ubc, ktop, gravit, rair, alpha, & 
             gw_prndl, gw_qbo_hdepth_scaling, & 
             errstring)

        allocate( rhoi(ncol,pver+1) , ni(ncol, pver+1), nm(ncol, pver)  )
        ! Profiles of background state variables
        call gw_prof(ncol, PP, cpair, T(:,:,itime), rhoi, nm, ni)
        write(*,*)  "Made ","rhoi   ",minval(rhoi), maxval(rhoi), shape(rhoi)

        write(*,*) "NEED GOEPHT"
        allocate( zm(ncol,pver) , zi(ncol,pver+1), rair3D(ncol,pver), zvir3D(ncol,pver) )

        ! zi is relative to the local 'ground' elevation ...!!!!!
        do k=1,pver+1
           zi(:,k) = zgrid(:,k)-zgrid(:,pver+1)
        end do
        do k=1,pver
           zm(:,k) = 0.5*( zi(:,k+1) + zi(:,k) )
        end do

        write(*,*)  "Made ","ZI   ",minval(zi), maxval(zi), shape(zi)
        write(*,*)  "Made ","ZM   ",minval(zm), maxval(zm), shape(zm)

        allocate( TH(ncol,pver,ntim) ,  FRONTGF(ncol,pver,ntim) ,  FRONTGA(ncol,pver,ntim)  )


        !TH(:,:,itime) = T(:,:,itime) * ( (100000._r8 / pmid)**(rair/cpair) )   ! + gravit * zm / cpair )
        TH  = T * ( (100000._r8 / pmid)**(rair/cpair) )   ! + gravit * zm / cpair )

        call sphere_fronto_wrapper( TH, U, V, lat_R, lon_R, FRONTGF, FRONTGA )
        
     else if ( trim(ncdata_type) == 'camsnap') then
        write(*,*) "Wahhooo"
        call ncread_camsnap( ncdata , ncol, pver , ntim, &
             hyai , hybi , hyam , hybm , lon , lat , &
             PS, U , V , T  , Q , ZETA, &
             zm_, zi_, nm_, ni_, rhoi_, pint , piln )

        pcols = ncol ! does this actaully go back to ppgrid? Yes.
        pver_in_ppgrid = pver

        allocate( pmid(ncol,pver,ntim)  )
        allocate( zm(ncol,pver) , zi(ncol,pver+1) )
        allocate( rhoi(ncol,pver+1) , ni(ncol, pver+1), nm(ncol, pver)  )

        PP = Coords1D(pint(:ncol,:,itime))
        write(*,*) "PP%ifc "
        write(*,*) minval( PP%ifc ),maxval( PP%ifc ), shape(PP%ifc)

        zm = zm_(:,:,itime)
        zi = zi_(:,:,itime)
        nm = nm_(:,:,itime)
        ni = ni_(:,:,itime)
        rhoi = rhoi_(:,:,itime)

        do n=1,ntim
           do k=1,pver
              pmid(:,k,n) =  0.5 * ( pint(:,k+1,n) +  pint(:,k,n) )
           end do
        end do

        allocate( pref_edge(pver+1)  )
        pref_edge(:) = 100000. * ( hyai(:) + hybi(:) )

        call leovy_alpha( pver, pref_edge, alpha )
        write(*,*) "Made ","alpha   ",minval(alpha), maxval(alpha), shape(alpha)

        call gw_common_init(pver,&
             tau_0_ubc, ktop, gravit, rair, alpha, & 
             gw_prndl, gw_qbo_hdepth_scaling, & 
             errstring)

     else
        write(*,*) "  You didn't give any valid calculation type/ncdata_type ..."
        STOP
     end if

     allocate( dse(ncol,pver) )
     dse = cpair*T(:,:,itime) + gravit*zm
  
     allocate( lqptend(pcnst))
     lqptend(:)=.TRUE.
     call physics_ptend_init(ptend, ncol , 'OGWD',.TRUE.,.TRUE.,.TRUE., lqptend )

     band_oro = GWBand(0, gw_dc, fcrit2, wavelength_mid)
     band_movmtn = GWBand(0, gw_dc, fcrit2, wavelength_mid)


     write(*,*) " PGWV for fronts ",pgwv
     band_front = GWBand(pgwv, gw_dc, fcrit2, wavelength_mid)

     call set_band_rdg(band_oro)
     call set_vramp_rdg()
     call report_from_within_rdg()

     call set_band_movmtn(band_movmtn)
     call set_vramp_movmtn()
     call report_from_within_movmtn()

     call set_band_front(band_front)
     call set_vramp_front()
     call report_from_within_front()

     if ( ( trim(ncdata_type) == 'XY_DATA') .or. &
          ( trim(ncdata_type) == 'XYMPAS_DATA')) then
        write(*,*) " Adding nx ny to dims ",ny,nx
        call ncfile_init_col( ncout , ncol, pver, time_val, ny, nx, lat_R, lon_R)
     else
        call ncfile_init_col( ncout , ncol, pver, time_val)
     end if

     call ncfile_set_globals(ncdata=ncdata, calculation_type=calculation_type )

     if ( trim(ncdata_type) .ne. 'XYMPAS_DATA' ) then
        call ncfile_put_col1d_notime('hyam', hyam, '1', 'hybrid midl a-coeff' )
        call ncfile_put_col1d_notime('hybm', hybm, '1', 'hybrid midl b-coeff' )
        call ncfile_put_col1d_notime('hyai', hyai, '1', 'hybrid intf a-coeff' )
        call ncfile_put_col1d_notime('hybi', hybi, '1', 'hybrid intf b-coeff' )
     else
        call ncfile_put_col1d_notime('ilev', ilev, '1', 'hybrid midl a-coeff' )
        call ncfile_put_col3d('zgrid' , zgrid, itime, 'm', 'real height at intfcs' )
     end if
     
     call ncfile_put_col2d_notime('lat', lat, 'deg', 'latitude' )
     call ncfile_put_col2d_notime('lon', lon, 'deg', 'longitude' )
     call ncfile_put_col2d_notime('SGH', sgh, 'm', 'topography std' )
     call ncfile_put_col2d('PS' , PS(:,itime), itime, 'Pa', 'surface pressure' )
     call ncfile_put_col3d('ZM' , zm, itime, 'm', 'height above ground at midlayer' )
     call ncfile_put_col3d('ZI' , zi, itime, 'm', 'height above ground at intfcs' )
     call ncfile_put_col3d('PMID' , pmid(:,:,itime), itime, 'Pa', 'pressure at midlayer' )
     call ncfile_put_col3d('PINT' , pint(:,:,itime), itime, 'Pa', 'pressure at intfcs' )
     call ncfile_put_col3d('U' , U(:,:,itime), itime, 'ms-1', 'zonal wind' )
     call ncfile_put_col3d('V' , V(:,:,itime), itime, 'ms-1', 'meridional wind' )
     call ncfile_put_col3d('T' , T(:,:,itime), itime, 'K', 'temperature' )

     if ( allocated( TH ) ) then
        call ncfile_put_col3d('TH' ,TH(:,:,itime), itime, 'K', 'potential temp' )
     end if
     if ( allocated( FRONTGF ) ) then
        call ncfile_put_col3d('FRONTGF' , FRONTGF(:,:,itime), itime, 'quaqu', 'frontogenesis' )
     end if
     if ( allocated( FRONTGA ) ) then
        call ncfile_put_col3d('FRONTGA' , FRONTGA(:,:,itime), itime, 'quaqu', 'frontogenesis' )
     end if

     lchnk=1


     write(*,*) "At input ","U      ",minval(U),    maxval(U),    shape(U)
     write(*,*) "At input ","V      ",minval(V),    maxval(V),    shape(V)
     write(*,*) "At input ","T      ",minval(T),    maxval(T),    shape(T)
     write(*,*) "At input ","Q      ",minval(Q),    maxval(Q),    shape(Q)
     write(*,*) "At input ","zm     ",minval(zm),   maxval(zm),   shape(zm)
     write(*,*) "At input ","zi     ",minval(zi),   maxval(zi),   shape(zi)
     write(*,*) "At input ","nm     ",minval(nm),   maxval(nm),   shape(nm)
     write(*,*) "At input ","ni     ",minval(ni),   maxval(ni),   shape(ni)
     write(*,*) "At input ","rhoi   ",minval(rhoi), maxval(rhoi), shape(rhoi)
     write(*,*) "At input ","pint   ",minval(pint), maxval(pint), shape(pint)
     write(*,*) "At input ","piln   ",minval(piln), maxval(piln), shape(piln)
     write(*,*) "At input ","PP%ifc ",minval(PP%ifc), maxval(PP%ifc), shape(PP%ifc)

  
     allocate( kvtt(ncol,pver+1) , flx_heat(ncol) )
     dt = 86400._r8 /48._r8 
     ! "molecular diffusivity" ... init(never changes) =>0.
     kvtt = 0._r8

     if (use_gw_rdg_beta) then
        n_rdg=1
        call gw_rdg_calc( &
             'BETA ', ncol, lchnk, n_rdg, dt, &
             U(:,:,itime) , V(:,:,itime) , T(:,:,itime) , &
             PP , piln(:,:,itime) , zm, zi, &
             nm, ni, rhoi, kvtt, Q(:,:,:), dse, &
             effgw_rdg_beta, effgw_rdg_beta_max, &
             effgw_rdg_resid, use_gw_rdg_resid, &
             hwdth, clngt, gbxar, &
             mxdis, angll, anixy, &
             isovar, isowgt, &
             rdg_beta_cd_llb, trpd_leewv_rdg_beta, &
                                !++jtb  ptend defined in massively hacked physics_types
             ptend, flx_heat) 
     end if


  
     if (use_gw_movmtn_pbl) then
        call gw_movmtn_calc( &
             ncol, lchnk, dt, pref_edge, &
             U(:,:,itime) , V(:,:,itime) , T(:,:,itime) , &
             ZETA(:,:,itime) , &
             PP , piln(:,:,itime) , zm, zi, &
             nm, ni, rhoi, kvtt, Q(:,:,:), dse, &
             effgw_movmtn_pbl, &
                                !++jtb  ptend defined in massively hacked physics_types
             ptend, flx_heat) 
     end if

     if (use_gw_front) then
        call gw_front_calc( &
             ncol, lchnk, dt, pref_edge, &
             U(:,:,itime) , V(:,:,itime) , T(:,:,itime) , &
             FRONTGF(:,:,itime) , &
             PP , piln(:,:,itime) , zm, zi, &
             nm, ni, rhoi, kvtt, Q(:,:,:), dse, &
             effgw_cm, &
             frontgfc, taubgnd, front_gaussian_width, &
                                !++jtb  ptend defined in massively hacked physics_types
             ptend, flx_heat) 
     end if

     close( unit=20 )
     call ncfile_close()
     
     call bump_ymdh(yy,mm,dd,hh, dt_step,  yy,mm,dd,hh , use_leap=.false. )                      ! -> 2004-02-29 02 (leap)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !  Need to deallocate !!!!                   !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call free_fields ()
  
  end do ! end of "time loop" ...


contains

  !===========================================================
subroutine free_fields( )

  use iso_fortran_env, only : real64
  implicit none

  !---- Local variables ----
  integer :: ierr

  !---- Deallocate all safely ----
  if (allocated(hyai))       deallocate(hyai)
  if (allocated(hybi))       deallocate(hybi)
  if (allocated(hyam))       deallocate(hyam)
  if (allocated(hybm))       deallocate(hybm)

  if (allocated(lon))        deallocate(lon)
  if (allocated(lat))        deallocate(lat)
  if (allocated(lon_R))      deallocate(lon_R)
  if (allocated(lat_R))      deallocate(lat_R)

  if (allocated(U))          deallocate(U)
  if (allocated(V))          deallocate(V)
  if (allocated(T))          deallocate(T)
  if (allocated(TH))         deallocate(TH)
  if (allocated(Q))          deallocate(Q)
  if (allocated(ZETA))       deallocate(ZETA)
  if (allocated(PS))         deallocate(PS)

  if (allocated(FRONTGF))    deallocate(FRONTGF)
  if (allocated(FRONTGA))    deallocate(FRONTGA)

  if (allocated(zgrid))      deallocate(zgrid)
  if (allocated(ilev))       deallocate(ilev)

  
  if (allocated(nm_))        deallocate(nm_)
  if (allocated(ni_))        deallocate(ni_)
  if (allocated(zm_))        deallocate(zm_)
  if (allocated(zi_))        deallocate(zi_)
  if (allocated(rhoi_))      deallocate(rhoi_)
  if (allocated(pint_))      deallocate(pint_)
  if (allocated(piln_))      deallocate(piln_)

  if (allocated(pint))       deallocate(pint)
  if (allocated(pmid))       deallocate(pmid)
  if (allocated(piln))       deallocate(piln)
  if (allocated(pmln))       deallocate(pmln)

  if (allocated(rhoi))       deallocate(rhoi)
  if (allocated(nm))         deallocate(nm)
  if (allocated(ni))         deallocate(ni)
  if (allocated(zm))         deallocate(zm)
  if (allocated(zi))         deallocate(zi)
  if (allocated(dse))        deallocate(dse)

  if (allocated(rair3D))     deallocate(rair3D)
  if (allocated(zvir3D))     deallocate(zvir3D)

  if (allocated(alpha))      deallocate(alpha)
  if (allocated(pref_edge))  deallocate(pref_edge)

  if (allocated(kvtt))       deallocate(kvtt)
  if (allocated(flx_heat))   deallocate(flx_heat)
  if (allocated(lqptend))    deallocate(lqptend)

end subroutine free_fields
!===========================================================

  
end program gw_driver

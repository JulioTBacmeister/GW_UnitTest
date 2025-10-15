##########
# get validation data
##########
import sys
import os
import glob

workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} ")

# Own local packages
from myPythonTools.Utils import AveragingUtils as Av
from myPythonTools.Utils import numerical_utils as nuti
from PyRegridding.Utils import MakePressures as MkP
from PyRegridding.Utils import GridUtils as GrU
from PyRegridding.Regridder import var_A_x_B as vAB

import numpy as np
import xarray as xr

import importlib

importlib.reload(vAB)

def data(fld,season=None ,months=-999,**kwargs):

    ##################################
    # Parse out logic etc of kwargs
    ##################################
    if 'mgrid' in kwargs: 
        if kwargs['mgrid'] in [False,None]:
            regrid_x_mgrid = False
        else:
            regrid_x_mgrid = True
    else:
        regrid_x_mgrid = False

    ##################################################################
    # A block for 'hourly' ERA5. Not even sure it still works.
    ##################################################################
    if 'ERA5hourly' in kwargs: 
        era5dir = "/glade/campaign/collections/rda/data/ds633.6/e5.oper.an.ml/"
        year,month,day,hour = kwargs['year'],kwargs['month'],kwargs['day'],kwargs['hour']
        hour0=(hour//6)*6
        hour1=hour0+5
        ymdh0=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+str(hour0).zfill(2)
        ymdh1=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+str(hour1).zfill(2)
        ymdh=ymdh0+'_'+ymdh1
        monStr=str( year ).zfill(4)+str(month).zfill(2)
        wrkdir=era5dir+monStr+"/"

        spfile= wrkdir + 'e5.oper.an.ml.128_134_sp.regn320sc.'+ymdh+'.nc'
        tfile = wrkdir + 'e5.oper.an.ml.0_5_0_0_0_t.regn320sc.'+ymdh+'.nc'
        qfile = wrkdir + 'e5.oper.an.ml.0_5_0_1_0_q.regn320sc.'+ymdh+'.nc'
        ufile = wrkdir + 'e5.oper.an.ml.0_5_0_2_2_u.regn320uv.'+ymdh+'.nc'
        vfile = wrkdir + 'e5.oper.an.ml.0_5_0_2_3_v.regn320uv.'+ymdh+'.nc'
        wfile = wrkdir + 'e5.oper.an.ml.0_5_0_2_8_w.regn320sc.'+ymdh+'.nc'

        ihour = hour-hour0
        if (fld == 'U'):
            fileN = wrkdir + 'e5.oper.an.ml.0_5_0_2_2_u.regn320uv.'+ymdh+'.nc'
            print(f'Reading {fld} in {fileN} at time index={ihour}')
            Dc = xr.open_dataset( fileN )
            aa = Dc.U[ihour,:,:,:].values
        if (fld == 'V'):
            fileN = wrkdir + 'e5.oper.an.ml.0_5_0_2_3_v.regn320uv.'+ymdh+'.nc'
            print(f'Reading {fld} in {fileN} at time index={ihour}')
            Dc = xr.open_dataset( fileN )
            aa = Dc.V[ihour,:,:,:].values
        if (fld == 'T'):
            fileN = wrkdir + 'e5.oper.an.ml.0_5_0_0_0_t.regn320sc.'+ymdh+'.nc'
            print(f'Reading {fld} in {fileN} at time index={ihour}')
            Dc = xr.open_dataset( fileN )
            aa = Dc.T[ihour,:,:,:].values
        if (fld == 'Q'):
            fileN = wrkdir + 'e5.oper.an.ml.0_5_0_1_0_q.regn320sc.'+ymdh+'.nc'
            print(f'Reading {fld} in {fileN} at time index={ihour}')
            Dc = xr.open_dataset( fileN )
            aa = Dc.Q[ihour,:,:,:].values

        plev = 100_000. * Dc.b_model.values +  Dc.a_model.values
        zlev = -7.0*np.log( plev /100_000. )
        #----- Pack output into a 'dict'
        dic={'aa':aa,
             'lev':zlev,
             'lat':Dc.latitude.values,
             'lon':Dc.longitude.values,
             'hyam':Dc.a_model.values,
             'hyai':Dc.a_half.values,
             'hybm':Dc.b_model.values,
             'hybi':Dc.b_half.values,
             'date':[year,month,day,hour],
             'years':ymdh, 'data_path':fileN,'data_source':'ERA5','data_exists':True, 'rcode':0}
        
        return dic
    ##################################################################
    # This ENDS the 'hourly' ERA block
    ##################################################################

    ADFobsdir = '/glade/campaign/cgd/amp/juliob/ADF_obs' #'/glade/work/nusbaume/SE_projects/model_diagnostics/ADF_obs/'
    AMWGobsdir = '/glade/campaign/cgd/amp/juliob/amwg_dev/obs_data/'
    
    if (fld in ('U','V','T','Q','PS',) ): #'OMEGA') ):

        if 'ERA5native' in kwargs:
            ERA5native=kwargs['ERA5native']
        else:
            ERA5native=False
            
        if ('Years' in kwargs) and ('Climo' not in kwargs):
            yearsA = kwargs['Years']
            if ( yearsA == '*' ):
                if ERA5native==True:
                    ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/native_monthly_climo/'
                else:
                    ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/fv1x1_monthly_climo/ZYX/'
                yearsA = '1979_2022'
            elif ( yearsA == '1979_2022' ):
                if ERA5native==True:
                    ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/native_monthly_climo/'
                else:
                    ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/fv1x1_monthly_climo/ZYX/'
            else: 
                if ERA5native==True:
                    ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/native_monthly/'
                else:
                    ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/fv1x1_monthly/ZYX/'
        elif 'Climo' in kwargs:
            ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/fv1x1_monthly_climo/ZYX/'
            if 'Years' in kwargs:
                yearsA = kwargs['Years']
            else:
                yearsA = '1979_2022'
        elif 'Timecube' in kwargs:
            yearsA = '*'
            if ERA5native==True:
                ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/native_monthly/'
            else:
                ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/fv1x1_monthly/ZYX/'
        else:
            yearsA = '1979_2022'
            if ERA5native==True:
                ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/native_monthly_climo/'
            else:
                ERA5dir = '/glade/campaign/cgd/amp/juliob/ERA5/fv1x1_monthly_climo/ZYX/'
        # get ERA5 
        if ERA5native==True:
            path_C = ERA5dir + 'ERA5.native.time.' + yearsA +'-*.nc'
        else:
            path_C = ERA5dir + 'ERA5.fv1x1_zERA.time.' + yearsA +'-*.nc'

        print( f" Trying validation data from {path_C}")

        Dc = xr.open_mfdataset( path_C ,data_vars='different', coords='different' )
        print( f" Validation data {fld} has dims {Dc[fld].dims}")
        print( f" Validation data read from {path_C}")
        
        if (Dc.hybm.dims[0] == 'time'):
            hybm = Dc.hybm.values[0,:]
        else:
            hybm = Dc.hybm.values
        if (Dc.hyam.dims[0] == 'time'):
            hyam = Dc.hyam.values[0,:]
        else:
            hyam = Dc.hyam.values
        if (Dc.hybi.dims[0] == 'time'):
            hybi = Dc.hybi.values[0,:]
        else:
            hybi = Dc.hybi.values
        if (Dc.hyai.dims[0] == 'time'):
            hyai = Dc.hyai.values[0,:]
        else:
            hyai = Dc.hyai.values

        lev = 1_000.*hyam[:] +1_000.*hybm[:] 
        if ('longitude' in Dc):
            lon =Dc.longitude.values
        elif ('lon' in Dc):
            lon =Dc.lon.values
        if ('latitude' in Dc):
            lat =Dc.latitude.values
        elif ('lat' in Dc):
            lat =Dc.lat.values

        if 'Timecube' in kwargs:
            aa=Dc[fld].values
            if 'zlev' in kwargs:
                lev = -7. * np.log( lev /1_000. )
            #----- Pack output into a 'dict'
            dic={'aa':aa,'lev':lev,'lat':lat,'lon':lon,
                 'hyai':hyai,'hyam':hyam,'hybi':hybi,'hybm':hybm,
                 'years':yearsA, 'data_path':path_C,'data_source':'ERA5',
                 'data_Hgrid':'ERA5','data_Vgrid':'ERA5','rcode':0}
            return dic

        
        aa = Av.Seasonal( ds=Dc, season=season , fld=fld)

        if 'zlev' in kwargs:
            lev = -7. * np.log( lev /1_000. )

        #----- Pack output into a 'dict'
        dic={'aa':aa,'lev':lev,'lat':lat,'lon':lon,
             'hyai':hyai,'hyam':hyam,'hybi':hybi,'hybm':hybm,
             'years':yearsA, 'data_exists':True, 'data_path':path_C,'data_source':'ERA5',
             'data_Hgrid':'ERA5','data_Vgrid':'ERA5','rcode':0}
        
    # /glade/work/nusbaume/SE_projects/model_diagnostics/ADF_obs/CERES_EBAF_Ed4.1_2001-2020.nc
    # need to map names from cesm to obs
       
    elif (fld in ('PRECT',) ):
        if ('data_source' in kwargs):
            data_source = kwargs['data_source']
        else:
            data_source = 'GPCP'
            
        if (data_source=='GPCP'):
            path_C = f'/glade/campaign/cgd/amp/juliob/cgd_amp/GPCP/gpcp.mon.mean.197901-201607.nc' #   /glade/work/nusbaume/SE_projects/model_diagnostics/ADF_obs/CERES_EBAF_Ed4.1_2001-2020.nc'
            Dc = xr.open_mfdataset( path_C ,data_vars='different', coords='different' )
            Dc = Dc.where(Dc.time.dt.year <= 2015, drop=True)
            print(f'Last time in GPCP {Dc.time[-1].values}')
            lon =Dc.lon.values
            lat =Dc.lat.values
            lev =np.asarray( [9999.] )
            aa = Av.Seasonal( ds=Dc, season=season , fld='precip' )
        else:
            dic = {'aa':-999e10, 'data_exists':False, 'rcode':-99}
            return dic

        dic={'aa':aa,'lev':lev,'lat':lat,'lon':lon,
             'years':'1979-2015', 'data_exists':True, 'data_path':path_C,'data_source':data_source,
             'data_Hgrid':None,'data_Vgrid':None,'rcode':0}


    elif (fld in ('SWCF',) ):
        path_C = f'{ADFobsdir}/CERES_EBAF_Ed4.1_2001-2020.nc' #   /glade/work/nusbaume/SE_projects/model_diagnostics/ADF_obs/CERES_EBAF_Ed4.1_2001-2020.nc'
        Dc = xr.open_mfdataset( path_C ,data_vars='different', coords='different' )
        lon =Dc.lon.values
        lat =Dc.lat.values
        lev =np.asarray( [9999.] )
        aa = Av.Seasonal( ds=Dc, season=season , fld='toa_cre_sw_mon' )
        dic={'aa':aa,'lev':lev,'lat':lat,'lon':lon,
             'years':'2001-2020', 'data_exists':True, 'data_path':path_C,'data_source':'CERES-EBAF',
             'data_Hgrid':None,'data_Vgrid':None,'rcode':0}


    elif (fld in ('SURFACE_STRESS',) ):
        path_C = f'{ADFobsdir}/LARYEA_climo.nc' #   /glade/work/nusbaume/SE_projects/model_diagnostics/ADF_obs/CERES_EBAF_Ed4.1_2001-2020.nc'
        Dc = xr.open_mfdataset( path_C ,data_vars='different', coords='different' )
        lon =Dc.lon.values
        lat =Dc.lat.values
        lev =np.asarray( [9999.] )
        aa = Av.Seasonal( ds=Dc, season=season , fld='STRESS_MAG' )
        dic={'aa':aa,'lev':lev,'lat':lat,'lon':lon,
             'years':'1979-2000', 'data_exists':True, 'data_path':path_C,'data_source':'Large-Yeager',
             'data_Hgrid':None,'data_Vgrid':None,'rcode':0}

    elif (fld in ('WIND_STRESS_CURL',) ):
        path_C = f'{ADFobsdir}/LARYEA_climo.nc' #   /glade/work/nusbaume/SE_projects/model_diagnostics/ADF_obs/CERES_EBAF_Ed4.1_2001-2020.nc'
        Dc = xr.open_mfdataset( path_C ,data_vars='different', coords='different' )
        lon =Dc.lon.values
        lat =Dc.lat.values
        lev =np.asarray( [1000.] )
        taux = Av.Seasonal( ds=Dc, season=season , fld='TAUX' )
        tauy = Av.Seasonal( ds=Dc, season=season , fld='TAUY' )
        aa=nuti.Sphere_Curl2( f_x=taux, f_y=tauy , lat=lat , lon=lon , wrap=True )
        dic={'aa':aa,'lev':lev,'lat':lat,'lon':lon,
             'years':'1979-2000', 'data_exists':True, 'data_path':path_C,'data_source':'Large-Yeager',
             'data_Hgrid':None,'data_Vgrid':None,'rcode':0}

    elif (fld in ('STRESS_MAG','TAUX','TAUY',) ):
        path_C = f'{ADFobsdir}/LARYEA_climo.nc' #   /glade/work/nusbaume/SE_projects/model_diagnostics/ADF_obs/CERES_EBAF_Ed4.1_2001-2020.nc'
        Dc = xr.open_mfdataset( path_C ,data_vars='different', coords='different' )
        lon =Dc.lon.values
        lat =Dc.lat.values
        lev =np.asarray( [1000.] )
        aa = Av.Seasonal( ds=Dc, season=season , fld=fld )
        if( fld in ('TAUX','TAUY',) ):
            aa = -1.0*aa
        dic={'aa':aa,'lev':lev,'lat':lat,'lon':lon,
             'years':'1979-2000', 'data_exists':True, 'data_path':path_C,'data_source':'Large-Yeager',
             'data_Hgrid':None,'data_Vgrid':None,'rcode':0}

    elif (fld in ('VQ') ):
        SuperDir = "/glade/campaign/cgd/amp/juliob/ERA5"
        Dst='fv1x1'
        years=[1991,1992,1993,1994]
        path=[]
        for year in years:
            regrd_monthly_dir = f"{SuperDir}/Fluxes/VQ/monthly/{Dst}-{str(year).zfill(4)}/"
            regrd_monthly_fil = f"{regrd_monthly_dir}/e5.oper.an.ml.VQ.{Dst}.{str(year).zfill(4)}-*.nc"
            newpath=sorted( glob.glob( regrd_monthly_fil ) )
            path = path + newpath

        Dc = xr.open_mfdataset( path ,data_vars='different', coords='different' )
        aa = Av.Seasonal( ds=Dc, season=season , fld=fld)
        lev = 1.*Dc.a_model.values +100_000.*Dc.b_model.values 
        lat=Dc.latitude
        lon=Dc.longitude
        print( Dc.VQ.shape )
        print(list(Dc.data_vars))
        yearsA='1991-1994'
        if 'zlev' in kwargs:
            lev = -7. * np.log( lev /100_000. )
        #----- Pack output into a 'dict'
        dic={'aa':aa,'lev':lev,'lat':lat,'lon':lon,
             'years':yearsA, 'data_exists':True, 'data_path':path,'data_source':'ERA5',
             'data_Hgrid':'fv1x1','data_Vgrid':'ERA5','rcode':0}
 

    else:
        dic = {'aa':-999e10, 'data_exists':False, 'rcode':-99}
        
            
    return dic  # aa,lev,lat,lon

def to_mgrid(aa=None,mgrid=None,obsgrid=None):

    print(f" In to_mgrid {mgrid.keys()}")  # Output will be dict_keys(['a', 'b', 'c'])
    
    #-----------------------------------------
    # Make 3D prssure field on Obs grid
    #-----------------------------------------
    pmid_o, pint_o, delp_o \
        = MkP.Pressure (am=obsgrid['hyam'],
                        bm=obsgrid['hybm'] ,
                        ai=obsgrid['hyai'] ,
                        bi=obsgrid['hybi'] ,
                        ps=obsgrid['ps'] ,
                        p_00=100_000. ,
                        Gridkey = 'zyx' )
    print(f" In to_mgrid shap pmid_o {np.shape(pmid_o)}")

    """
    if ( Dst == 'ne30pg3'):
    griddir = "/glade/work/juliob/ERA5-proc/ERA5interp/grids/"
    wgts_file_Con = griddir + "ERA5_ne30pg3_Conserv_wgts.nc"
    write_weights = False 
    read_weights = True 
    """

    pmid_o_Hm , lat_o_Hm , lon_o_Hm \
        = vAB.Hregrid(avar=pmid_o,
                        agrid=obsgrid['hgrid'],
                        akey='zyx',
                        bgrid=mgrid['hgrid'] )

    print(f" In to_mgrid shap pmid_o_Hm {np.shape(pmid_o_Hm)}")
 
    aa_Hm , lat_o_Hm , lon_o_Hm \
        = vAB.Hregrid(avar=aa,
                        agrid=obsgrid['hgrid'],
                        akey='zyx',
                        bgrid=mgrid['hgrid'] )

    print(f" In to_mgrid shap aa_Hm {np.shape(aa_Hm)}")
    print(f" In to_mgrid shap mgrid ps {np.shape(mgrid['ps'])}")

    
    pmid_m, pint_m, delp_m \
        = MkP.Pressure (am=mgrid['hyam'],
                        bm=mgrid['hybm'] ,
                        ai=mgrid['hyai'] ,
                        bi=mgrid['hybi'] ,
                        ps=mgrid['ps'] ,
                        p_00=100_000. ,
                        Gridkey = 'zyx' )

    p_00 = 100_000. # Here we just use the sensible value of p_00
    lnpmid_o_Hm   = -7_000. * np.log( pmid_o_Hm / p_00 )
    lnpmid_m = -7_000. * np.log( pmid_m / p_00 )
    
    aa_m = vAB.VertRG( a_x = aa_Hm ,
                zSrc = lnpmid_o_Hm ,
                zDst = lnpmid_m ,
                Gridkey ='zyx' ,
                kind = 'quadratic' ) #linea


    GridInfo = GrU.gridInfo( grid=mgrid['hgrid'] )
    lat,lon = GrU.latlon( scrip=GridInfo['scrip'],Hkey=GridInfo['Hkey'] )

    lev = 1_000.*mgrid['hyam'][:] +1_000.*mgrid['hybm'][:] 

    oDict = {'aa':aa_m, 'lev':lev, 'lat':lat,'lon':lon }
    
    return oDict





   

from datetime import date
import numpy as np

# Some useful packages 
import importlib
import copy
import time as ttime
import cftime

def MonthsSeason( season ):
    if (season.lower()=='ann'):
        monthsx=[1,2,3,4,5,6,7,8,9,10,11,12]
    if (season.lower()=='djf'):
        monthsx=[12,1,2]
    if (season.lower()=='jfm'):
        monthsx=[1,2,3]
    if (season.lower()=='jja'):
        monthsx=[6,7,8]
    if (season.lower()=='jas'):
        monthsx=[7,8,9]
    if (season.lower()=='mam'):
        monthsx=[3,4,5]
    if (season.lower()=='son'):
        monthsx=[9,10,11]
    if (season.lower()=='jan'):
        monthsx=[1]
    if (season.lower()=='feb'):
        monthsx=[2]
    if (season.lower()=='mar'):
        monthsx=[3]
    if (season.lower()=='apr'):
        monthsx=[4]
    if (season.lower()=='may'):
        monthsx=[5]
    if (season.lower()=='jun'):
        monthsx=[6]
    if (season.lower()=='jul'):
        monthsx=[7]
    if (season.lower()=='aug'):
        monthsx=[8]
    if (season.lower()=='sep'):
        monthsx=[9]
    if (season.lower()=='oct'):
        monthsx=[10]
    if (season.lower()=='nov'):
        monthsx=[11]
    if (season.lower()=='dec'):
        monthsx=[12]
    return monthsx

def ListMatch( list1 , list2 ):
    
    theMonth=[]  # This can be any 'element' 
    IndexIn1=[]
    IndexIn2=[]
    # Check each element in list1
    icommon=0
    for idx1, elem1 in enumerate(list1):
        for idx2, elem2 in enumerate(list2):
            if elem1 == elem2:
                theMonth.append(elem1)
                IndexIn1.append(idx1)
                IndexIn2.append(idx2)
                icommon += 1
    
    return IndexIn1, IndexIn2, theMonth 

def MonthsInDataset( ds ):

    print( f" In MonthsInDataset function ")
    if ('time_bnds' in ds):        
        
        if (ds['time_bnds'].dtype == 'datetime64[ns]' ):
            # normally the first element of time_bnds
            # gives the month you want.
            time_array = ds['time_bnds'].values[:,0]
            months=[]
            years=[]
            for ixtime in time_array:
                years.append( ixtime.astype('datetime64[Y]' ).astype(int) + 1970 )                
                months.append( ixtime.astype('datetime64[M]' ).astype(int) %12 +1 )                
        else:
            print(f" time_bnds present in DataSet " )
            time_bnds = ds['time_bnds']
    
            #print(np.shape(time))
            #print(time_bnds[0].values[0])
            time0=[]
            months=[]
            years=[]
    
            for ixtime in time_bnds:
                time0.append( ixtime.values[0] )
                years.append( ixtime.values[0].year )
                months.append( ixtime.values[0].month )
            
    elif ('time_bounds' in ds):        
        
        if (ds['time_bounds'].dtype == 'datetime64[ns]' ):
            # normally the first element of time_bnds
            # gives the month you want.
            time_array = ds['time_bounds'].values[:,0]
            months=[]
            years=[]
            for ixtime in time_array:
                years.append( ixtime.astype('datetime64[Y]' ).astype(int) + 1970 )                
                months.append( ixtime.astype('datetime64[M]' ).astype(int) %12 +1 )                
        else:
            print(f" time_bounds (post cam6_3_153 ... ) present in DataSet " )
            time_bnds = ds['time_bounds']
    
            #print(np.shape(time))
            #print(time_bnds[0].values[0])
            time0=[]
            months=[]
            years=[]
    
            for ixtime in time_bnds:
                time0.append( ixtime.values[0] )
                years.append( ixtime.values[0].year )
                months.append( ixtime.values[0].month )
            
    elif ('time' in ds) :
        if (ds['time'].dtype == 'datetime64[ns]' ):
            print(f" time of type datetime64[ns] present in DataSet " )
            years = ds['time'].values.astype('datetime64[Y]' ).astype(int)
            months = ds['time'].values.astype('datetime64[M]' ).astype(int) %12 +1
        elif (ds['time'].dtype == 'int' ):
            print(f" time of type 'int' present in DataSet " )
            #if (np.shape( ds['time'].shape )==(1,)): 
            if ('title' in ds.attrs):
                # if (('CERES EBAF' in ds.title) and ('Climatology' in ds.title)):
                if (('Climatology' in ds.title)):
                    print( f'{ds.title}' )
                    months = ds['time']
                    years  = np.zeros( 12 )+9999
        else:
            # Raise an error if no valid dtype is found
            raise ValueError(f"Unsupported time array dtype: {ds['time'].dtype}")


    return months,years

def parse_dims( dims ):

    ndims = len( dims )
    n=0
    if (dims[n]=='time'):
        gridKey = 't'
    elif (dims[n] in ('lev','level')):
        gridKey = 'z'
    elif (dims[n] in ('lat','latitude')):
        gridKey = 'y'
    elif (dims[n]=='ncol'):
        gridKey = 'c'
        
    n=1
    if (dims[n]=='time'):
        gridKey = gridKey + 't'
    elif (dims[n] in ('lev','level')):
        gridKey = gridKey + 'z'
    elif (dims[n] in ('lat','latitude')):
        gridKey = gridKey + 'y'
    elif (dims[n]=='ncol'):
        gridKey = gridKey + 'c'

    if (ndims > 2):
        n=2
        if (dims[n]=='time'):
            gridKey = gridKey + 't'
        elif (dims[n] in ('lev','level')):
            gridKey = gridKey + 'z'
        elif (dims[n] in ('lat','latitude')):
            gridKey = gridKey + 'y'
        elif (dims[n] in ('lon','longitude')):
            gridKey = gridKey + 'x'
        elif (dims[n]=='ncol'):
            gridKey = gridKey + 'c'

    if (ndims > 3):
        n=3
        if (dims[n]=='time'):
            gridKey = gridKey + 't'
        elif (dims[n] in ('lev','level')):
            gridKey = gridKey + 'z'
        elif (dims[n] in ('lat','latitude')):
            gridKey = gridKey + 'y'
        elif (dims[n] in ('lon','longitude')):
            gridKey = gridKey + 'x'
        elif (dims[n]=='ncol'):
            gridKey = gridKey + 'c'
    
    return gridKey

def Seasonal( ds, season, **kwargs ):
    
    #----------------------------------------------------
    # Fist thing to do is trim the Dataset to ensure/test
    # 'seasonal contiguity' if desired. Note, that
    # unless something is done in except block below,
    # code will just print an error message and go on 
    # to clauclate seasonal mean anyway
    #----------------------------------------------------
    if 'contiguous' in kwargs:
        contiguous = kwargs['contiguous']
    else:
        contiguous=True
    if 'verbose' in kwargs:
        verbose = kwargs['verbose']
    else:
        verbose=False

    if ( (season=='djf')and(contiguous==True) ):
        try:
            ds=contiguous_DJF( ds )
            print("Ensured contiguous DJF")
        except ValueError as e:
            print(e)

    months,years = MonthsInDataset( ds )
    
    imos = MonthsSeason( season=season )
    imonths=np.asarray(months)
    Iseason = ListMatch( list1=imonths, list2=imos )
    
    print( "Indices of months", Iseason[0] )
    nmos=len(Iseason[0]) 
    print( f" length of Iseason {nmos} ")

    months_in_av = np.zeros( nmos, dtype=int )
    years_in_av  = np.zeros( nmos, dtype=int )
    for n in np.arange( nmos ):
        t = Iseason[0][n]
        months_in_av[n] = months[t]
        years_in_av[n]  = years[t]

    #--------------------------------------------
    # Now get varibale from Dataset and try to 
    # figure out its shape
    #-------------------------------------------
    if 'fld' in kwargs:
        fld = kwargs['fld']
        try:
            A = ds[fld]
        except:
            print(f'{fld} not dataset')
            A_mmm=-999.e9
            if kwargs['return_time']==True:
                return A_mmm,years_in_av,months_in_av
            else:
                return A_mmm
        
    elif 'data' in kwargs:
        A = kwargs['data']
    else:
        assert A is not None, "A needs to given "
    
    if 'dims' in kwargs:
        VarDims = kwargs['dims']
    else:
        VarDims = parse_dims(ds[fld].dims) #'tzyx'
        print(f" dims key parsed to {VarDims} ")

    if (VarDims=='tzyx'):
        nt,nz,ny,nx = np.shape( A )
        A_mmm = np.zeros( (nz,ny,nx) )
        for n in np.arange( nmos ):
            t = Iseason[0][n]
            print( t ,end=',')
            A_mmm = A_mmm + A[ t,:,:,:].values/nmos
            
        print(" time _X_ 3D lev-lat-lon variable " )
        
    if (VarDims=='tyx'):
        A_mmm = np.average( A[ Iseason[0] ,:,:], axis=0 )
        print(" time _X_ 2D lat-lon variable " )

    if (VarDims=='tzc'):
        nt,nz,ncol = np.shape( A )
        A_mmm = np.zeros( (nz,ncol ) )
        for n in np.arange( nmos ):
            t = Iseason[0][n]
            print( t ,end=',')
            A_mmm = A_mmm + A[ t,:,:].values/nmos
            
        print(" time _X_ 3D lev-ncol variable " )
        
    if (VarDims=='tc'):
        A_mmm = np.average( A[ Iseason[0] ,: ], axis=0 )
        print(" time _X_ 2D ncol variable " )
    
    if ('return_time' in kwargs):
        if kwargs['return_time']==True:
            return A_mmm,years_in_av,months_in_av
    else:
        return A_mmm

def contiguous_DJF(ds):
    """
    Finds a contiguous sequence of months from December to February (DJF) in the input xarray dataset.
    
    Parameters:
        months (list or array): A list or array of month numbers (1-12).
    
    Returns:
        array: The contiguous DJF months if found.
    
    Raises:
        ValueError: If December or February is not found, or if a contiguous DJF period is not possible.
    """

    months,years = MonthsInDataset( ds )
    monthsN = np.array(months)
    
    # Find the first occurrence of December (12)
    indices = np.where(monthsN == 12)[0]
    if indices.size > 0:
        firstDec = indices[0]
    else:
        raise ValueError("December not found in the dataset.")
    
    # Find the last occurrence of February (2)
    indices = np.where(monthsN == 2)[0]
    if indices.size > 0:
        lastFeb = indices[-1]
    else:
        raise ValueError("February not found in the dataset.")
    
    # Ensure that December occurs before February for a valid DJF period
    if firstDec < lastFeb:
        print( f"Succesfully trimming dataset from {firstDec} to {lastFeb}" )
        print( f"{ds.time[firstDec].values}" )
        print( f"{ds.time[lastFeb].values}" )
        ds_djf = ds.isel(time=slice(firstDec, lastFeb+1))  #[firstDec:lastFeb+1]
        return ds_djf
    else:
        raise ValueError("Contiguous DJF not possible")










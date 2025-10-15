##########################################
#
##########################################
import sys

workdir_ = '/glade/work/juliob/'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} ")


# The usual
from datetime import date
import numpy as np
import xarray as xr

# Some other useful packages 
import importlib
import copy
import time
import cftime
import glob
import os

def extract_save( X1=None,X2=None,X3=None,exp_X=None, read=None ):
    ############################################################
    # Bespoke for analysis of moving mountain development runs.
    # No point in generalizing here
    ############################################################
    # X1 contains state vars
    # X2 contains tendencies
    # X3 contains CLUBB vars or potential movmtn forcers

    if ((read==None) or (read==False)):
        dic_X={}
        
        scaletend=86_400.
        scaleprec=86_400. * 1000.
        
        fld='lon'
        dic_X[fld] = X1[fld].values
        print(f" Added {fld} ")
        fld='lat'
        dic_X[fld] = X1[fld].values
        print(f" Added {fld} ")
        fld='lev'
        dic_X[fld] = X1[fld].values
        print(f" Added {fld} ")
    
        vars_in_X1 = ['U','V','T',]
        vars_in_X2 = ['UTEND_GWDTOT', 'UTEND_PHYSTOT', 'UTEND_CORE', 'Nudge_U' ,  ]
        scale_in_X2 = [ scaletend,     scaletend,       scaletend,    scaletend,  ]
        vars_in_X3  = ['UPWP_CLUBB', 'VPWP_CLUBB', 'WP2_CLUBB', 'THLP2_CLUBB', 'STEND_CLUBB', 'PRECC', 'PRECL',  ]
        scale_in_X3 = [   1.0      ,    1.0   ,      1.0,          1.0,           1.0,        scaleprec, scaleprec, ] 
        
        fld='U'
        if (fld in X1):
            dic_X[fld] = X1[fld].values
            print(f" Added {fld} ")
        fld='V'
        if (fld in X1):
            dic_X[fld] = X1[fld].values
            print(f" Added {fld} ")
        fld='T'
        if (fld in X1):
            dic_X[fld] = X1[fld].values
            print(f" Added {fld} ")
    
    
        
        fld='UTEND_GWDTOT'
        if (fld in X2):
            dic_X[fld] = X2[fld].values  * scaletend
            print(f" Added {fld} ")
        fld='UTEND_PHYSTOT'
        if (fld in X2):
            dic_X[fld] = X2[fld].values  * scaletend
            print(f" Added {fld} ")
        fld='UTEND_CORE'
        if (fld in X2):
            dic_X[fld] = X2[fld].values  * scaletend
            print(f" Added {fld} ")
        fld='Nudge_U'
        if (fld in X2):
            dic_X[fld] = X2[fld].values  * scaletend
            print(f" Added {fld} ")
    
        
        fld='UPWP_CLUBB'
        if (fld in X3):
            dic_X[fld] = X3[fld].values
            print(f" Added {fld} ")
        fld='VPWP_CLUBB'
        if (fld in X3):
            dic_X[fld] = X3[fld].values
            print(f" Added {fld} ")
        fld='WP2_CLUBB'
        if (fld in X3):
            dic_X[fld] = X3[fld].values
            print(f" Added {fld} ")
        fld='THLP2_CLUBB'
        if (fld in X3):
            dic_X[fld] = X3[fld].values
            print(f" Added {fld} ")
        fld='STEND_CLUBB'
        if (fld in X3):
            dic_X[fld] = X3[fld].values
            print(f" Added {fld} ")
        fld='PRECC'
        if (fld in X3):
            dic_X[fld] = X3[fld].values  * scaleprec
            print(f" Added {fld} ")
        fld='PRECL'
        if (fld in X3):
            dic_X[fld] = X3[fld].values  * scaleprec
            print(f" Added {fld} ")
    
        savedir = '/glade/derecho/scratch/juliob/numpy_Saves'
        os.makedirs( savedir , exist_ok=True )
        savefil = f"{savedir}/{exp_X}_hiFrequency_extract.npz"
        print( savefil )
        np.savez( savefil , **dic_X)
    else:
        savedir = '/glade/derecho/scratch/juliob/numpy_Saves'
        savefil = f"{savedir}/{exp_X}_hiFrequency_extract.npz"
        print( savefil )
        dic_X = np.load( savefil )
    
    return dic_X


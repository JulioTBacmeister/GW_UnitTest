##############################################
# Run like
#    python setup_gw_build.py build_ext --inplace
###############################################
from numpy.distutils.core import setup, Extension

extension_mod = Extension(
    name='fortran_gwave_module',
    sources=[
        'share/shr_kind_mod.F90', 
        'gw_utils.F90',
        'utils/cam_abortutils.F90',
        'share/shr_const_mod.F90',
        'utils/coords_1d.F90',
        'utils/linear_1d_operators.F90',        
        'main_routine.F90',
    ],
    extra_f90_compile_args=['-O0']  # No optimization, easy debugging
)

setup(
    name='fortran_gwave_module',
    ext_modules=[extension_mod]
)

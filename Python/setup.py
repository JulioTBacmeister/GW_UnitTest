from numpy.distutils.core import setup, Extension

extension_mod = Extension(
    name='fortran_array_module',
    sources=['mod_array.F90', 'main_routine.F90'],
    extra_f90_compile_args=['-O0']  # No optimization, easy debugging
)

setup(
    name='fortran_array_module',
    ext_modules=[extension_mod]
)

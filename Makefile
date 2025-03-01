# Compiler and flags
FC = gfortran
FFLAGS = -O0 -fPIC -Wall          # Debug-friendly flags
LDFLAGS = -shared
NETCDF_INC = $(shell nf-config --fflags)
#NETCDF_LIBS = $(shell nf-config --flibs)
NETCDF_LIBS = -L/glade/u/apps/opt/conda/envs/npl-2024b/lib -lnetcdff -lnetcdf
TARGET_LIB = libfortran_gwave_module.so
PYTHON_MODULE = gwave_wrapper.cpython-*.so
EXEC = gw_driver.x

# Source files for the shared library
SRCS = \
    share/shr_kind_mod.F90 \
    share/shr_const_mod.F90 \
    share/shr_string_mod.F90 \
    share/shr_file_mod.F90 \
    share/shr_nl_mod.F90 \
    utils/physconst.F90 \
    utils/ncread_mod.F90 \
    utils/units.F90 \
    utils/namelist_utils.F90 \
    utils/cam_abortutils.F90 \
    utils/namelist_mod.F90 \
    utils/coords_1d.F90 \
    utils/linear_1d_operators.F90 \
    utils/vdiff_lu_solver.F90 \
    utils/interpolate_data.F90 \
    utils/utils_mod.F90 \
    src/ppgrid.F90 \
    src/physics_types.F90 \
    src/geopotential.F90 \
    src/gw_utils.F90 \
    src/gw_diffusion.F90 \
    src/gw_common.F90 \
    src/gw_rdg.F90 \
    src/gw_rdg_calc_mod.F90

OBJS = $(SRCS:.F90=.o)

# === Default target ===
all: $(TARGET_LIB)

# === Shared library ===
$(TARGET_LIB): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(NETCDF_LIBS)

%.o: %.F90
	$(FC) $(FFLAGS) $(NETCDF_INC) -c $< -o $@

# === Python wrapper using f2py ===
python_module: $(TARGET_LIB) gw_interface.F90
	f2py -c gw_interface.F90 -m gwave_wrapper -L. -lfortran_gwave_module --fcompiler=gnu95

# === Fortran executable (standalone) ===
$(EXEC): $(TARGET_LIB) gw_driver.F90
	$(FC) -o $@ gw_driver.F90 -L. -lfortran_gwave_module  -Wl,-rpath=. $(NETCDF_LIBS)

# === Run the Fortran executable ===
run: $(EXEC)
	./$(EXEC)

# === Clean targets ===
clean:
	rm -f $(OBJS) $(TARGET_LIB) $(EXEC) $(PYTHON_MODULE) *.mod GW.dat

clean_python:
	rm -f $(PYTHON_MODULE)

clean_all: clean clean_python

# Phony targets (always run)
.PHONY: all clean clean_python clean_all python_module run

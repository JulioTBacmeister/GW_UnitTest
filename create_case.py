#!/usr/bin/env python
import os
import sys
import shutil
import subprocess
import glob
import argparse
import yaml
import re

#--------------------------------------------------------
# Minimal call(s)
# if YAML control file exists
#    ./create_case.py
# else if don't want to specify run type
#    ./create_case.py --casename YourCase
# else
#    ./create_case.py --casename YourCase --type RunType
#---------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate topo cases.")
    parser.add_argument("--output-dir", type=str, default="./cases", dest="output_dir", help="Path to mother of casedirs ")
    parser.add_argument("--casename", type=str, default="exp00", help="casename ")
    parser.add_argument("--type", type=str, default="xy", help="type of calc - xy, camsnap, or ERA5")
    parser.add_argument("--config", type=str, default="create_case.yaml", help="Path to YAML configuration file (default: create_case.yaml)")
    parser.add_argument("--clean-case", action="store_true", dest="clean_case", help="Delete existing case directory if it exists")
    parser.add_argument("--run", action="store_true", dest="run_case", help="If set runs case, if not just builds")
    parser.add_argument("--build", action="store_true", dest="build_case", help="If set builds case")

    parser.add_argument("--ncdata_root", type=str, dest="ncdata_root", help="Path to ncdata ")
    parser.add_argument("--bnd_topo", type=str, dest="bnd_topo", help="topo file ")


    return parser.parse_args()

def load_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def clean_case_directory(case_dir):
    if os.path.exists(case_dir):
        print(f"Cleaning existing case directory: {case_dir}")
        shutil.rmtree(case_dir)

def update_or_add_namelist_value(filename, key, new_value):
    """
    Update a Fortran namelist assignment if it exists, otherwise add it
    before the namelist terminator ('/').

    Preserves spacing and comments.
    """
    pattern = re.compile(rf"^(\s*{key}\s*=\s*)'([^']*)'(.*)$")
    out_lines = []
    found = False

    with open(filename, "r") as f:
        for line in f:
            m = pattern.match(line)
            if m:
                prefix, old_val, suffix = m.groups()
                out_lines.append(f"{prefix}'{new_value}'{suffix}\n")
                found = True
            else:
                out_lines.append(line)

    # If key was not found, insert it before '/'
    if not found:
        for i, line in reversed(list(enumerate(out_lines))):
            if line.strip().startswith("/"):
                indent = "  "  # or infer from file
                out_lines.insert(i, f"{indent}{key} = '{new_value}'\n")
                break
        else:
            raise RuntimeError("No namelist terminator '/' found")

    with open(filename, "w") as f:
        f.writelines(out_lines)
        
def update_namelist_value(filename, key, new_value):
    """
    Updates a Fortran namelist assignment like:
        key = 'old'
    and replaces it with:
        key = 'new_value'
    
    Preserves spacing and comments.
    """
    print( filename )
    
    pattern = re.compile(rf"^(\s*{key}\s*=\s*)'([^']*)'(.*)$")
    out_lines = []

    with open(filename, 'r') as f:
        for line in f:
            m = pattern.match(line)
            if m:
                # reconstruct the line, preserving whitespace & trailing text
                prefix, old_val, suffix = m.groups()
                new_line = f"{prefix}'{new_value}'{suffix}\n"
                out_lines.append(new_line)
            else:
                out_lines.append(line)

    with open(filename, 'w') as f:
        f.writelines(out_lines)

def read_namelist_value(filename, key):
    """
    Updates a Fortran namelist assignment like:
        key = 'old'
    and replaces it with:
        key = 'new_value'
    
    Preserves spacing and comments.
    """
    print( filename )
    
    pattern = re.compile(rf"^(\s*{key}\s*=\s*)'([^']*)'(.*)$")
    out_lines = []

    with open(filename, 'r') as f:
        for line in f:
            m = pattern.match(line)
            if m:
                # reconstruct the line, preserving whitespace & trailing text
                prefix, current_val, suffix = m.groups()

    return current_val

def main():

    # Parse command-line arguments
    args = parse_arguments()

    # Load settings from YAML if the file exists
    config = {}
    try:
        config = load_config(args.config)
    except FileNotFoundError:
        print(f"Configuration file '{args.config}' not found. Using defaults where applicable.")

    # Override YAML settings with explicit command-line arguments
    odir = args.output_dir or config.get("output_dir")
    casename = args.casename or config.get("casename")
    runtype = args.type or config.get("type")
    clean_case = args.clean_case if args.clean_case else config.get("clean_case", False)
    run_case = args.run_case if args.run_case else config.get("run_case", False)
    build_case = args.build_case if args.build_case else config.get("build_case", False)

    if ( runtype == 'xy'):
        bnd_topo_def = '/glade/work/juliob/bndtopo/fv1x1_gmted2010_modis_bedmachine_nc3000_Laplace0100_noleak_greenlndantarcsgh30fac2.50_20251009.nc'
        ncdata_root_def  = '/glade/derecho/scratch/juliob/archive//c153_topfix_ne240pg3_FMTHIST_xic_x02/atm/fv1x1/c153_topfix_ne240pg3_FMTHIST_xic_x02.cam.h1i'
        ncdata_type = 'XY_DATA'
    elif ( runtype == 'xympas'):
        bnd_topo_def = '/glade/work/juliob/bndtopo/fv1x1_gmted2010_modis_bedmachine_nc3000_Laplace0100_noleak_greenlndantarcsgh30fac2.50_20251009.nc'
        ncdata_root_def  = '/glade/derecho/scratch/juliob/archive/c124_dyamond1/atm/hist/DynVars_dyamond_fv1x1'
        ncdata_type = 'XYMPAS_DATA'
    elif ( runtype == 'camsnap'):
        bnd_topo_def  = '/glade/work/juliob/bndtopo/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_noleak_20240720.nc'
        ncdata_root_def	 = '/glade/derecho/scratch/juliob/archive/ndg_x02_ne30pg3_fmt_c64109/atm/hist/ndg_x02_ne30pg3_fmt_c64109.cam.h2i' 
        ncdata_type = 'camsnap'
    elif ( runtype == 'ERA5'):
        bnd_topo_def  = '/glade/work/juliob/bndtopo/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_noleak_20240720.nc'
        ncdata_root_def	 = '/glade/campaign/cgd/amp/juliob/ERA5/ne30pg3/L93/2014/ERA5_x_ne30pg3_L93_rgC2_WO'  
        ncdata_type = 'ERA5_SE_IC'
    else:
        bnd_topo_def  = 'None'
        ncdata_root_def	 = 'None'  
        
    if args.ncdata_root is not None:
        ncdata_root = args.ncdata_root
    elif config.get("ncdata_root") is not None:
        ncdata_root = config["ncdata_root"]
    else:
        ncdata_root = ncdata_root_def
    print( f" ncdata_root , {ncdata_root}" )
        
    if args.bnd_topo is not None:
        bnd_topo = args.bnd_topo
    elif config.get("bnd_topo") is not None:
        bnd_topo = config["ncdata_root"]
    else:
        bnd_topo = bnd_topo_def
    
    casename = f'{runtype}-{casename}'
    case_dir = os.path.join(odir, casename )

    # Clean the case directory if the flag is set
    if clean_case==True :
        clean_case_directory(case_dir)

    # Create directories
    os.makedirs(os.path.join(case_dir, "src"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "share"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "utils"), exist_ok=True)

    # Copy files
    for file in ["Makefile", "gw_driver.F90","namelist_imports.inc", "atm_in", "drv_in"]:
        shutil.copy(file, case_dir)
    shutil.copy(__file__, os.path.join(case_dir, "create_case.py"))
    
    # Copy all .F90 files to the target directories
    for subdir in ["src","share","utils"]:
        f90_files = glob.glob(f"{subdir}/*.F90")  # This will match all .F90 files in the current directory
        for f90_file in f90_files:
            shutil.copy(f90_file, f"{case_dir}/{subdir}/")

    # Move to case directory
    os.chdir(case_dir)

    
    #----------------------------------------------------------
    # Load gcc and compile.  Note these command are glommed 
    # together because otherwise gcc's settting of FC=gfortran
    # is not passed to the make commands ...
    #----------------------------------------------------------
    if  ( (run_case==True) or (build_case==True) ):
        command=f"module load gcc; make clean; make gw_driver.x"
        result = subprocess.run(command, shell=True, executable="/bin/tcsh", capture_output=True, text=True)
        print(result.stdout)
    else:
        command=f"module load gcc; make clean"
        result = subprocess.run(command, shell=True, executable="/bin/tcsh", capture_output=True, text=True)
        print(result.stdout)

    #print(result.stderr)
    print("Return code:", result.returncode)
    
    #Modify namelist drv_in
    update_namelist_value( f"drv_in", "calculation_type", runtype)
    update_namelist_value( f"drv_in", "casename", casename)
    
    update_or_add_namelist_value( f"drv_in", "ncdata_root", ncdata_root)
    update_or_add_namelist_value( f"drv_in", "bnd_topo", bnd_topo)
    update_or_add_namelist_value( f"drv_in", "ncdata_type", ncdata_type)

    ncout_root = read_namelist_value('drv_in', 'ncout_root' ) 
    os.makedirs(os.path.join(ncout_root, casename ), exist_ok=True)
    
    if ( (result.returncode == 0) and (run_case==True) ):
        command=f"pwd; ./gw_driver.x"
        result = subprocess.run(command, shell=True, executable="/bin/tcsh", capture_output=True, text=True)
        print(result.stdout)
    else:
        print(result.stderr)
    
if __name__ == "__main__":
    main()

    

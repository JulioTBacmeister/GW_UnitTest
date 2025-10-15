#!/usr/bin/env python
import os
import sys
import shutil
import subprocess
import glob
import argparse
import yaml

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate topo cases.")
    #parser.add_argument("--smoothing_scale", type=int, help="Smoothing scale")
    parser.add_argument("--output_dir", type=str, default="./cases", help="Path to mother of casedirs ")
    parser.add_argument("--casename", type=str, default="exp00", help="casename ")
    parser.add_argument("--config", type=str, default="create_topo.yaml", help="Path to YAML configuration file (default: create_topo.yaml)")
    parser.add_argument("--clean_case", action="store_true", help="Delete existing case directory if it exists")
    return parser.parse_args()

def load_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def clean_case_directory(case_dir):
    if os.path.exists(case_dir):
        print(f"Cleaning existing case directory: {case_dir}")
        shutil.rmtree(case_dir)


def create_command(odir=None, casename=None ):
    """
    Creates a command for subprocess.run, in this list form:
 
    command = [
        "./gw_driver.x",
        f"--output_dir={odir}",
        f"--casename={casename}",
         ]
    """
    #----------------------------------------------------
    # Some of these 'optional' arguments are not really
    # optional - ogrid, cstopo, scrip and smoothing_scale.
    # Sorry.
    #-----------------------------------------------------
    command = [  "./cube_to_target" ]
    #------------------------------------------------------------
    # We are to rely on Python's 'truthy' vs 'falsy' idea below,
    # to decide whether an input is 'provided'. If this doesn't 
    # work can always go to 'if var is not None:'
    #------------------------------------------------------------
    if odir: 
        command.append(f"--output_dir={odir}")
    else:
        raise ValueError("output directory is required but was not provided.")

    print("created command line", flush=True)
    print( command, flush=True )
    
    return command

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
    clean_case = args.clean_case if args.clean_case else config.get("clean_case", False)

    case_dir = os.path.join(odir, casename )

    # Clean the case directory if the flag is set
    if clean_case==True :
        clean_case_directory(case_dir)

    # Create directories
    os.makedirs(os.path.join(case_dir, "src"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "share"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "utils"), exist_ok=True)

    # Copy files
    for file in ["Makefile", "gw_driver.F90","namelist_imports.inc", "atm_in"]:
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
    command=f"module load gcc; make clean; make gw_driver.x"
    result = subprocess.run(command, shell=True, executable="/bin/tcsh", capture_output=True, text=True)
    print(result.stdout)

    #print(result.stderr)
    print("Return code:", result.returncode)

    if (result.returncode == 0):
        command=f"pwd; ./gw_driver.x"
        result = subprocess.run(command, shell=True, executable="/bin/tcsh", capture_output=True, text=True)
        print(result.stdout)
    else:
        print(result.stderr)
    
    """
    #------------------------------------------------------------------------------------
    # Get information about destination grid 'ogrid' 
    #------------------------------------------------------------------------------------
    grid_info = gridInfo( ogrid )
    scrip = grid_info['scrip']
    scrip_gll = grid_info['scrip_gll']
    yfac = grid_info['yfac']
    
    #------------------------------------------------------------------------------------
    # Specify cubed-sphere 'intermediate' topography. Currently, code is configured 
    # to work with a 3000x3000x6 grid, which corresponds to a grid size of around 3km.
    # Clearly, this is a problem at k-scale resolutions.
    #------------------------------------------------------------------------------------
    cstopo = "/glade/work/juliob/Topo/CubeData/gmted2010_modis_bedmachine-ncube3000-220518.nc"

    cog = f"Co{int(smoothing_scale):03d}"
    smtopo = f"topo_smooth_gmted2010_bedmachine_nc3000_{cog}.nc"

    print("Here you are")
    print(f"Smoothing scale: {cog}")
    print(f"Smooth topo file: {smtopo}")

    # Create command line that runs cube_to_target
    command = create_command(ogrid=ogrid, 
                             cstopo=cstopo, 
                             smoothing_scale=smoothing_scale, 
                             scrip=scrip, 
                             scrip_gll=scrip_gll, 
                             yfac=yfac, 
                             development_diags=True )
        
    subprocess.run(command, check=True)
    """
if __name__ == "__main__":
    main()

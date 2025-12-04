# GW_UnitTest
Test harness for GW codes

1) Fortran infrastructure from cam (more or less tag=cam6_4_055) to allow gravity wave codes to be run in
   Unit test mode. Code currently reproduces results from full model to within round-off
   when forced with 'snapshot' of physics state right before calls to GW code. These fields are found
   ncdata file specified in drv_in namelist file.

   A Makefile is present. To build the test, run 'module load gcc' (once per login), and then simply type

     	    make gw_driver.x

   This builds the executable. Executable produces netcdf output files that can read by notebooks in
   'Python' subdirectory - see ./Python/AnaGenl.ipynb .

   Infrastructure in Python has also been added to create cases - see create_case.py in root dir .

2) Names of cam files have not been changed, although many have been dramatically hacked to simplify build.
   Original gravity wave source code which forms the basis for this Unit test are in subdir:

   	    ./original_gw_codes
	    
   Modified or adapted GW Codes unsed in the Unit test are in subdirectory

   	    ./src

   These codes were ported with minimal possible changes

3) Some new files have been created that are not in cam code base.

4) At some point a guide to the provenance of each file will be provided.

5) Added ChangeLog 9/2/2025. Still bit-for-bit* with cam6_4_109.

6) Major add-ons 10/5/2025. See ChangeLog.


*NOTE (11/30/25) - I've realized these results are NOT bit-for-bit with model, just equal at round-off level
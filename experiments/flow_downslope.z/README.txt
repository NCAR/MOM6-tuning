# MOM6 flow_downslope expirement for precision tuning.
  Target module : MOM_CoriolisAdv.F90

## Files:

### MOM6 input files:
 - MOM_input : Default MOM6 parameters.
 - MOM_override: Parameter definitions in addition to (or overriding) MOM_input parameters.
 - diag_table : Configures diagnostics output recorded in netcdf format.
 - input.nml : Some general MOM6 and FMS namelist parameters.

### Cheyenne-specific files used for running Precimonious:

 * Note: Cheyenne is an NCAR supercomputer. To run this experiment on a different machine, make
    sure to copy cheyenne directory and modify the copied files accordingly.

 - cheyenne/build_MOM6_GNU.sh : Bash script to build MOM6 using the GNU compiler
 - cheyenne/build_MOM6_ROSE_depends.sh : Bash script to build ROSE dependencies, e.g., rmod files,
    for a given target module.
 - cheyenne/GNU_env.sh : Bash script to switch to GNU environment. Running this script should
    activate GNU as the default fortran compiler. The script should also update environment
    variables, e.g., PATH, LD_LIBRARY_PATH, etc., accordingly.
 - cheyenne/ROSE_env.sh : Bash script to switch to ROSE environment. Running this script should
    activate the ROSE compiler. The script should also update environment variables, e.g., PATH,
    LD_LIBRARY_PATH, etc., accordingly.
 - cheyenne/setup.ini : Precimonious setup.ini file to be passed to both prec_preprocess.py and
    prec_search.py scripts. This file configures different stages of the Precimonious execution.
    Users running on a different machine should ideally need to change only the [machine] section
    in this file.

### User config file:

 - user_cfg_CorAdCalc.txt : a predefined user config file containing only the local variables of
    CorAdCalc subroutine in MOM_CoriolisAdv module

## Instructions:

### Step 0: If you are not running on the cheyenne supercomputer, copy the cheyenne directory and
 modify all the files under the new copy of the directory according to your machine specifications,
 paths, modules, etc.

### Step 1: Make sure that your copies of the bash scripts work successfully. If you are on
 cheyenne, for example, run the following scripts to confirm that the scripts can successfully
 build ROSE dependenciesi and build the MOM6 executable:

    source ./cheyenne/ROSE_env.sh
    ./cheyenne/build_MOM6_ROSE_depends.sh MOM_CoriolisAdv.o

    source ./cheyenne/GNU_env.sh
    ./cheyenne/build_MOM6_GNU.sh MOM6

 If the above commands ran successfully, run them again to make sure that scripts can indeed make
 the switch from ROSE environment to GNU and visa versa.

### Step 2: Run the Precimonious preprocessor:

    /glade/work/altuntas/ROSE/precimonious-w-rose/scripts/prec_preprocess.py -s [YOUR-MACH-DIR]/setup.ini

### Step 3: Run the Precimonious search script:

 Note: before running this step, make sure that MPI is active and multiple cores are available.

    /glade/work/altuntas/ROSE/precimonious-w-rose/scripts/prec_search.py -s [YOUR-MACH-DIR]/setup.ini -c user_cfg_CorAdCalc.txt



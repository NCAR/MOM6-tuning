#!/bin/bash

REPO_ROOT=/root/MOM6-tuning

FMS_BLD_DIR=${REPO_ROOT}/build/rose/shared/repro
MOM6_BLD_DIR=${REPO_ROOT}/build/rose/ocean_only/repro
ROSE_MKMF_TEMPLATE=${REPO_ROOT}/src/mkmf/templates/linux-ubuntu-trusty-rose.mk

TARGET_MOM6_OBJ=$1

set -v 
set -e # abort if any command fails

# can we put in the if/else check from the GNU.sh script to see if FMS is built already?

# build FMS:
mkdir -p $FMS_BLD_DIR
cd $FMS_BLD_DIR
rm -f path_names *.?90 *.mod *.rmod *.o
../../../../src/mkmf/bin/list_paths ../../../../src/FMS
../../../../src/mkmf/bin/mkmf -t $ROSE_MKMF_TEMPLATE -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD" path_names
make clean
set +e # do not abort if FMS build fails
make NETCDF=3 REPRO=1 libfms.a

set -e # abort if any command fails

# build MOM6:
mkdir -p $MOM6_BLD_DIR
cd $MOM6_BLD_DIR
rm -f path_names *.?90 *.mod *.rmod *.o
../../../../src/mkmf/bin/list_paths ./ ../../../../src/MOM6/{config_src/dynamic,config_src/solo_driver,src/{*,*/*}}/
../../../../src/mkmf/bin/mkmf -t $ROSE_MKMF_TEMPLATE -o '-I../../shared/repro' -p MOM6 -l '-L../../shared/repro -lfms' -c '-Duse_libMPI -Duse_netCDF -DSPMD' path_names
make NETCDF=3 REPRO=1 $TARGET_MOM6_OBJ

echo "Successfully generated the ROSE-specific dependencies"

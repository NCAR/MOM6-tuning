#!/bin/bash

REPO_ROOT=/root/MOM6-tuning

FMS_BLD_DIR=${REPO_ROOT}/build/intel/shared/repro
MOM6_BLD_DIR=${REPO_ROOT}/build/intel/ocean_only/repro
INTEL_MKMF_TEMPLATE=${REPO_ROOT}/src/mkmf/templates/linux-ubuntu-trusty-intel.mk

TARGET_MOM6_OBJ=$1

set -v # print shell input lines as they are read
set -e # abort if any command fails

# build FMS:
if test -f "$FMS_BLD_DIR/libfms.a"; then
    echo "FMS already built."
else
  mkdir -p $FMS_BLD_DIR
  cd $FMS_BLD_DIR
  rm -f path_names *.?90 *.mod *.rmod *.o
  ../../../../src/mkmf/bin/list_paths ../../../../src/FMS
  ../../../../src/mkmf/bin/mkmf -t $INTEL_MKMF_TEMPLATE -p libfms.a -c "-Duse_libMPI -Duse_netCDF" path_names
  make clean
  make NETCDF=3 REPRO=1 INCLUDES=-I/usr/include/x86_64-linux-gnu/ libfms.a -j
fi

# build MOM6:
mkdir -p $MOM6_BLD_DIR
cd $MOM6_BLD_DIR
rm -f path_names *.?90 *.mod *.rmod *.o
../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/MOM6/{config_src/dynamic,config_src/solo_driver,src/{*,*/*}}/
../../../../src/mkmf/bin/mkmf -t $INTEL_MKMF_TEMPLATE -o '-I../../shared/repro' -p 'MOM6 -L../../shared/repro -lfms' -c '-Duse_libMPI -Duse_netCDF -DSPMD' path_names
make NETCDF=3 REPRO=1 $TARGET_MOM6_OBJ -j

echo "Successfully built MOM6 with Intel Compiler."

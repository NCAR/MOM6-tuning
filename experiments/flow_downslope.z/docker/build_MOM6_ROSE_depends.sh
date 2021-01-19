#!/bin/bash

REPO_ROOT=/root/MOM6-tuning

FMS_BLD_DIR=${REPO_ROOT}/build/rose/shared/repro
MOM6_BLD_DIR=${REPO_ROOT}/build/rose/ocean_only/repro
ROSE_MKMF_TEMPLATE=${REPO_ROOT}/src/mkmf/templates/linux-ubuntu-trusty-rose.mk

TARGET_MOM6_OBJ=$1

set -v 

# build FMS:
if test -f "$FMS_BLD_DIR/generated_dependencies_flag"; then
    echo "Already built rose dependencies."
else
  mkdir -p $FMS_BLD_DIR
  cd $FMS_BLD_DIR
  rm -f path_names *.?90 *.mod *.rmod *.o
  ../../../../src/mkmf/bin/list_paths ../../../../src/FMS
  ../../../../src/mkmf/bin/mkmf -t $ROSE_MKMF_TEMPLATE -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD" path_names
  make clean

  # make command expected to fail; we don't need the final libfms.a, we just need the rose
  # compiler to generate its versions of all of the prerequisite modules
  # set +e to prevent script from aborting when make fails
  set +e
  make NETCDF=3 REPRO=1 libfms.a

  # mark dependencies as having been generated
  touch generated_dependencies_flag
fi

# build MOM6:
mkdir -p $MOM6_BLD_DIR
cd $MOM6_BLD_DIR
rm -f path_names *.?90 *.mod *.rmod *.o
../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/MOM6/{config_src/dynamic,config_src/solo_driver,src/{*,*/*}}/
#../../../../src/mkmf/bin/mkmf -t $ROSE_MKMF_TEMPLATE -o '-I../../shared/repro' -p MOM6 -l '-L../../shared/repro -lfms' -c '-Duse_libMPI -Duse_netCDF -DSPMD' path_names
../../../../src/mkmf/bin/mkmf -t $ROSE_MKMF_TEMPLATE -o '-I../../shared/repro' -p 'MOM6 -L../../shared/repro -lfms' -c '-Duse_libMPI -Duse_netCDF -DSPMD' path_names
# command taken from wiki instructions for ubuntu; removed -l flag, moved MOM6 inside quotes
make NETCDF=3 REPRO=1 $TARGET_MOM6_OBJ -j #(this flag is from the wiki instructions for ubuntu...)

echo "Successfully generated the ROSE-specific dependencies"

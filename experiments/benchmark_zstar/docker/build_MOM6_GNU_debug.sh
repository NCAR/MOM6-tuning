#!/bin/bash

REPO_ROOT=/root/MOM6-tuning

FMS_BLD_DIR=${REPO_ROOT}/build/gnu/shared/debug
MOM6_BLD_DIR=${REPO_ROOT}/build/gnu/ocean_only/debug
GNU_MKMF_TEMPLATE=${REPO_ROOT}/src/mkmf/templates/linux-ubuntu-trusty-gnu.mk

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
  ../../../../src/mkmf/bin/mkmf -t $GNU_MKMF_TEMPLATE -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD" path_names
  make clean
  make NETCDF=3 DEBUG=1 libfms.a -j #(this flag is from the wiki instructions for ubuntu...)
fi

# build MOM6:
mkdir -p $MOM6_BLD_DIR
cd $MOM6_BLD_DIR
rm -f path_names *.?90 *.mod *.rmod *.o
../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/MOM6/{config_src/dynamic,config_src/solo_driver,src/{*,*/*}}/
#../../../../src/mkmf/bin/mkmf -t $GNU_MKMF_TEMPLATE -o '-I../../shared/repro' -p MOM6 -l '-L../../shared/repro -lfms' -c '-Duse_libMPI -Duse_netCDF -DSPMD' path_names
../../../../src/mkmf/bin/mkmf -t $GNU_MKMF_TEMPLATE -o '-I../../shared/debug' -p 'MOM6 -L../../shared/debug -lfms' -c '-Duse_libMPI -Duse_netCDF -DSPMD' path_names
# command taken from wiki instructions for ubuntu; removed -l flag, moved MOM6 inside quotes
make NETCDF=3 DEBUG=1 $TARGET_MOM6_OBJ -j #(this flag is from the wiki instructions for ubuntu...)

echo "Successfully built MOM6 with GNU."

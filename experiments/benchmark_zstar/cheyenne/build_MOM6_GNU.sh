#!/bin/bash

MOM6_TUNING_ROOT=/glade/work/jdvanover/MOM6-tuning

FMS_BLD_DIR="${MOM6_TUNING_ROOT}"/build/gnu/shared/repro
MOM6_BLD_DIR="${MOM6_TUNING_ROOT}"/build/gnu/ocean_only/repro
GNU_MKMF_TEMPLATE="${MOM6_TUNING_ROOT}"/src/mkmf/templates/cheyenne-gnu.mk

TARGET_MOM6_OBJ=$1

set -v 
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
  make NETCDF=3 REPRO=1 libfms.a
fi

# build MOM6:
mkdir -p $MOM6_BLD_DIR
cd $MOM6_BLD_DIR
rm -f path_names *.?90 *.mod *.rmod *.o
../../../../src/mkmf/bin/list_paths ./ ../../../../src/MOM6/{config_src/dynamic,config_src/solo_driver,src/{*,*/*}}/
../../../../src/mkmf/bin/mkmf -t $GNU_MKMF_TEMPLATE -o '-I../../shared/repro' -p MOM6 -l '-L../../shared/repro -lfms' -c '-Duse_libMPI -Duse_netCDF -DSPMD' path_names
make NETCDF=3 REPRO=1 $TARGET_MOM6_OBJ

echo "Successfully built MOM6 with GNU."

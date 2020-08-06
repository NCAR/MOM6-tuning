#!/bin/sh

FMS_BLD_DIR=/glade/work/altuntas/mom6.standalone/MOM6-tuning/build/intel/shared/repro
MOM6_BLD_DIR=/glade/work/altuntas/mom6.standalone/MOM6-tuning/build/intel/ocean_only/repro
INTEL_MKMF_TEMPLATE=/glade/work/altuntas/mom6.standalone/MOM6-tuning/src/mkmf/templates/cheyenne-intel.mk

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
  ../../../../src/mkmf/bin/mkmf -t $INTEL_MKMF_TEMPLATE -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD" path_names
  make clean
  make NETCDF=3 REPRO=1 libfms.a
fi

# build MOM6:
mkdir -p $MOM6_BLD_DIR
cd $MOM6_BLD_DIR
rm -f path_names *.?90 *.mod *.rmod *.o
../../../../src/mkmf/bin/list_paths ./ ../../../../src/MOM6/{config_src/dynamic,config_src/solo_driver,src/{*,*/*}}/
../../../../src/mkmf/bin/mkmf -t $INTEL_MKMF_TEMPLATE -o '-I../../shared/repro' -p MOM6 -l '-L../../shared/repro -lfms' -c '-Duse_libMPI -Duse_netCDF -DSPMD' path_names
make NETCDF=3 REPRO=1 $TARGET_MOM6_OBJ

echo "Successfully built MOM6 with INTEL."

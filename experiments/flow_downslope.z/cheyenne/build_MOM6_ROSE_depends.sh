#!/bin/sh

FMS_BLD_DIR=/glade/work/altuntas/mom6.standalone/MOM6-tuning/build/rose/shared/repro
MOM6_BLD_DIR=/glade/work/altuntas/mom6.standalone/MOM6-tuning/build/rose/ocean_only/repro
ROSE_MKMF_TEMPLATE=/glade/work/altuntas/mom6.standalone/MOM6-tuning/src/mkmf/templates/cheyenne-rose.mk

TARGET_MOM6_OBJ=$1

set -v 
set -e # abort if any command fails

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

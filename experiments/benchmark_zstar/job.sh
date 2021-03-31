#!/bin/bash
#PBS -N mom6_job
#PBS -A UCDV0023
#PBS -l walltime=12:00:00
#PBS -q regular
#PBS -j oe
#PBS -k eod
#PBS -l select=3:ncpus=36:mpiprocs=36
#PBS -l inception=login

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

source /glade/u/home/jdvanover/precimonious-w-rose/scripts/cheyenne/activate_ROSE.sh

### Run the executable
python3 /glade/u/home/jdvanover/precimonious-w-rose/scripts/prose_full.py -s cheyenne/setup.ini

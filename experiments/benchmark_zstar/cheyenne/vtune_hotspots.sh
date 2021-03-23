#!/bin/bash
#PBS -N mom6_job
#PBS -A UCDV0023
#PBS -l walltime=00:20:00
#PBS -q regular
#PBS -j oe
#PBS -k eod
#PBS -l select=1:ncpus=36:mpiprocs=36

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

### Run the executable
rm -rf result result.txt
amplxe-cl --collect=hotspots --result-dir=result --return-app-exitcode -- mpiexec "${1}"
err="${?}"
if [ "${err}" -eq 0 ]; then
    amplxe-cl --report=gprof-cc --result-dir=result --format=text --report-output=result.txt
    pushd prose_logs > /dev/null 2>&1
    dest="$(printf "./prose_logs/%03d/report.txt" $(($(find -type d -regex "./[0-9]+$" | wc -l) - 2)))"
    popd > /dev/null 2>&1
    cp result.txt "${dest}"
    exit 0
else
    exit 1
fi


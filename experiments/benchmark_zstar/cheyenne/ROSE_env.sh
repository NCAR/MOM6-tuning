#!/usr/bin/bash

# variables for this script
DEPENDENCIES_ROOT=/glade/work/altuntas/ROSE
ROSE_INSTALL=/glade/work/jdvanover/ROSE/install
PROSE_REPO=/glade/u/home/jdvanover/precimonious-w-rose

source /etc/profile.d/modules.sh

module purge
module reset
module load python
ncar_pylib

# set environment variables
export LD_LIBRARY_PATH="${ROSE_INSTALL}"/lib:"${DEPENDENCIES_ROOT}"/jdk1.8.0_241/jre/lib/amd64/server:"${DEPENDENCIES_ROOT}"/flex-2.6.4/install/lib:"${DEPENDENCIES_ROOT}"/boost/1_67_0/install/lib:"${DEPENDENCIES_ROOT}"/gcc/7.4.0/install/lib64:$LD_LIBRARY_PATH
export PATH="${PROSE_REPO}"/scripts:"${ROSE_INSTALL}"/bin:"${DEPENDENCIES_ROOT}"/jdk1.8.0_241/bin/:"${DEPENDENCIES_ROOT}"/jdk1.8.0_241/jre/bin/:"${DEPENDENCIES_ROOT}"/automake-1.16.2/install/bin/:"${DEPENDENCIES_ROOT}"/flex-2.6.4/install/bin:"${DEPENDENCIES_ROOT}"/gettext-0.19.7/install/bin:"${DEPENDENCIES_ROOT}"/gcc/7.4.0/install/bin:$PATH
install/bin:/glade/work/altuntas/ROSE/gcc/7.4.0/install/bin:$PATH

export JRE_HOME="${DEPENDENCIES_ROOT}"/jdk1.8.0_241/jre
export JAVA_BINDIR="${DEPENDENCIES_ROOT}"/jdk1.8.0_241/bin
export JAVA_HOME="${DEPENDENCIES_ROOT}"/jdk1.8.0_241
export SDK_HOME="${DEPENDENCIES_ROOT}"/jdk1.8.0_241
export JDK_HOME="${DEPENDENCIES_ROOT}"/jdk1.8.0_241
export JAVA_ROOT="${DEPENDENCIES_ROOT}"/jdk1.8.0_241

export ROSE_EXE_PATH="${ROSE_INSTALL}"/bin
export ROSE_PLUGIN_PATH="${PROSE_REPO}"/plugins
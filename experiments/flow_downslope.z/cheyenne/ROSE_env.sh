#!/usr/bin/sh

source /etc/profile.d/modules.sh

module purge
module load python

export LD_LIBRARY_PATH=/glade/work/altuntas/ROSE/rose/install/lib:/glade/work/altuntas/ROSE/jdk1.8.0_241/jre/lib/amd64/server:/glade/work/altuntas/ROSE/flex-2.6.4/install/lib:/glade/work/altuntas/ROSE/boost/1_67_0/install/lib:/glade/work/altuntas/ROSE/gcc/7.4.0/install/lib64:$LD_LIBRARY_PATH

export PATH=/glade/work/altuntas/ROSE/rose/install/bin:/glade/work/altuntas/ROSE/jdk1.8.0_241/bin/:/glade/work/altuntas/ROSE/jdk1.8.0_241/jre/bin/:/glade/work/altuntas/ROSE/automake-1.16.2/install/bin/:/glade/work/altuntas/ROSE/flex-2.6.4/install/bin:/glade/work/altuntas/ROSE/gettext-0.19.7/install/bin:/glade/work/altuntas/ROSE/gcc/7.4.0/install/bin:$PATH

export JRE_HOME=/glade/work/altuntas/ROSE/jdk1.8.0_241/jre
export JAVA_BINDIR=/glade/work/altuntas/ROSE/jdk1.8.0_241/bin
export JAVA_HOME=/glade/work/altuntas/ROSE/jdk1.8.0_241
export SDK_HOME=/glade/work/altuntas/ROSE/jdk1.8.0_241
export JDK_HOME=/glade/work/altuntas/ROSE/jdk1.8.0_241
export JAVA_ROOT=/glade/work/altuntas/ROSE/jdk1.8.0_241

export ROSE_EXE_PATH=/glade/work/altuntas/ROSE/rose/install/bin
export ROSE_PLUGIN_PATH=/glade/work/altuntas/ROSE/precimonious-w-rose/plugins

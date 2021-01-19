# Template for the GNU Compiler Collection on Trusty version of Ubuntu Linux systems (used by Travis-CI)
#
# Typical use with mkmf
# mkmf -t linux-ubuntu-trusty-gnu.mk -c"-Duse_libMPI -Duse_netCDF" path_names /usr/local/include

############
# Commands Macors
FC = identityTranslator
CC = mpicc
LD = mpif90 $(MAIN_PROGRAM)

#######################
# Build target macros
#
# Macros that modify compiler flags used in the build.  Target
# macrose are usually set on the call to make:
#
#    make REPRO=on NETCDF=3
#
# Most target macros are activated when their value is non-blank.
# Some have a single value that is checked.  Others will use the
# value of the macro in the compile command.

DEBUG =              # If non-blank, perform a debug build (Cannot be
                     # mixed with REPRO or TEST)

REPRO =              # If non-blank, erform a build that guarentees
                     # reprodicuibilty from run to run.  Cannot be used
                     # with DEBUG or TEST

TEST  =              # If non-blank, use the compiler options defined in
                     # the FFLAGS_TEST and CFLAGS_TEST macros.  Cannot be
                     # use with REPRO or DEBUG

VERBOSE =            # If non-blank, add additional verbosity compiler
                     # options

OPENMP =             # If non-blank, compile with openmp enabled

NO_OVERRIDE_LIMITS = # If non-blank, do not use the -qoverride-limits
                     # compiler option.  Default behavior is to compile
                     # with -qoverride-limits.

NETCDF =             # If value is '3' and CPPDEFS contains
                     # '-Duse_netCDF', then the additional cpp macro
                     # '-Duse_LARGEFILE' is added to the CPPDEFS macro.

INCLUDES =           # A list of -I Include directories to be added to the
                     # the compile command.

SSE =                # The SSE options to be used to compile.  If blank,
                     # than use the default SSE settings for the host.
                     # Current default is to use SSE2.

COVERAGE =           # Add the code coverage compile options.

# Need to use at least GNU Make version 3.81
need := 3.81
ok := $(filter $(need),$(firstword $(sort $(MAKE_VERSION) $(need))))
ifneq ($(need),$(ok))
$(error Need at least make version $(need).  Load module gmake/3.81)
endif

# REPRO, DEBUG and TEST need to be mutually exclusive of each other.
# Make sure the user hasn't supplied two at the same time
ifdef REPRO
ifneq ($(DEBUG),)
$(error Options REPRO and DEBUG cannot be used together)
else ifneq ($(TEST),)
$(error Options REPRO and TEST cannot be used together)
endif
else ifdef DEBUG
ifneq ($(TEST),)
$(error Options DEBUG and TEST cannot be used together)
endif
endif

# modified by Jackson to be a single job so that files can be written to
# MAKEFLAGS += --jobs=$(shell grep '^processor' /proc/cpuinfo | wc -l)
MAKEFLAGS += --jobs=1

# Required Preprocessor Macros:
CPPDEFS += -Duse_netCDF

# Additional Preprocessor Macros needed due to  Autotools and CMake
# commented out by jackson b/c this flag wasn't present in the calls to the compiler made by cheyenne-rose.mk
#CPPDEFS += -DHAVE_SCHED_GETAFFINITY

# Macro for Fortran preprocessor
FPPFLAGS := $(INCLUDES)
# Fortran Compiler flags for the NetCDF library
FPPFLAGS += $(shell nf-config --fflags)

# Base set of Fortran compiler flags
FFLAGS := -fcray-pointer -fdefault-double-8 -fdefault-real-8 -Waliasing -ffree-line-length-none -fno-range-check

# Flags based on perforance target (production (OPT), reproduction (REPRO), or debug (DEBUG)
FFLAGS_OPT = -O3
FFLAGS_REPRO = -O2 -fbounds-check
FFLAGS_DEBUG = -O0 -g -W -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow

# Flags to add additional build options
FFLAGS_OPENMP = -fopenmp
FFLAGS_VERBOSE =
FFLAGS_COVERAGE =

# Macro for C preprocessor
CPPFLAGS = $(INCLUDES)
# C Compiler flags for the NetCDF library
CPPFLAGS += $(shell nf-config --cflags)

# Base set of C compiler flags
CFLAGS := -D__IFC

# Flags based on perforance target (production (OPT), reproduction (REPRO), or debug (DEBUG)
CFLAGS_OPT = -O2
CFLAGS_REPRO = -O2
CFLAGS_DEBUG = -O0 -g

# Flags to add additional build options
CFLAGS_OPENMP = -fopenmp
CFLAGS_VERBOSE =
CFLAGS_COVERAGE =

# Optional Testing compile flags.  Mutually exclusive from DEBUG, REPRO, and OPT
# *_TEST will match the production if no new option(s) is(are) to be tested.
FFLAGS_TEST = $(FFLAGS_OPT)
CFLAGS_TEST = $(CFLAGS_OPT)

# Linking flags
LDFLAGS :=
LDFLAGS_OPENMP := -fopenmp
LDFLAGS_VERBOSE :=
LDFLAGS_COVERAGE :=

# Start with a blank LIBS
LIBS =
# NetCDF library flags
LIBS += $(shell nf-config --flibs)

# Get compile flags based on target macros.
ifdef REPRO
CFLAGS += $(CFLAGS_REPRO)
FFLAGS += $(FFLAGS_REPRO)
else ifdef DEBUG
CFLAGS += $(CFLAGS_DEBUG)
FFLAGS += $(FFLAGS_DEBUG)
else ifdef TEST
CFLAGS += $(CFLAGS_TEST)
FFLAGS += $(FFLAGS_TEST)
else
CFLAGS += $(CFLAGS_OPT)
FFLAGS += $(FFLAGS_OPT)
endif

ifdef OPENMP
CFLAGS += $(CFLAGS_OPENMP)
FFLAGS += $(FFLAGS_OPENMP)
LDFLAGS += $(LDFLAGS_OPENMP)
endif

ifdef SSE
CFLAGS += $(SSE)
FFLAGS += $(SSE)
endif

ifdef NO_OVERRIDE_LIMITS
FFLAGS += $(FFLAGS_OVERRIDE_LIMITS)
endif

ifdef VERBOSE
CFLAGS += $(CFLAGS_VERBOSE)
FFLAGS += $(FFLAGS_VERBOSE)
LDFLAGS += $(LDFLAGS_VERBOSE)
endif

ifeq ($(NETCDF),3)
  # add the use_LARGEFILE cppdef
  ifneq ($(findstring -Duse_netCDF,$(CPPDEFS)),)
    CPPDEFS += -Duse_LARGEFILE
  endif
endif

ifdef COVERAGE
ifdef BUILDROOT
PROF_DIR=-prof-dir=$(BUILDROOT)
endif
CFLAGS += $(CFLAGS_COVERAGE) $(PROF_DIR)
FFLAGS += $(FFLAGS_COVERAGE) $(PROF_DIR)
LDFLAGS += $(LDFLAGS_COVERAGE) $(PROF_DIR)
endif

LDFLAGS += $(LIBS)

################### Added by jackson ########################

# these flags were present in cheyenne-rose.mk
FFLAGS += -DOVERLOAD_R8 -DOVERLOAD_R4

# ROSE flags:
FFLAGS += -rose:plugin_lib /root/precimonious-w-rose/plugins/ProsePlugin.so -rose:plugin_action prose-generate-graph -rose:plugin_arg_prose-generate-graph /root/MOM6-tuning/experiments/flow_downslope.z/
FFLAGS += -rose:skip_syntax_check
FFLAGS += -rose:skipfinalCompileStep
FFLAGS += -rose:cray_pointer_support
FFLAGS += -DROSEPREP

# ROSE+MPI flags:
#MPIPATH = /glade/work/altuntas/ROSE/openmpi-4.0.3/install/
MPIPATH = /usr/lib/x86_64-linux-gnu/openmpi
FFLAGS += -I${MPIPATH}/include -pthread -I${MPIPATH}/lib -Wl,-rpath -Wl,${MPIPATH}/lib -Wl,--enable-new-dtags -L${MPIPATH}/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi 
CFLAGS +=  -I${MPIPATH}/include -pthread -Wl,-rpath -Wl,${MPIPATH}/lib -Wl,--enable-new-dtags -L${MPIPATH}/lib -lmpi

# ROSE+NETCDF flags:
#NETCDFPATH = /glade/u/apps/ch/opt/netcdf/4.7.3/gnu/7.4.0/
#FFLAGS += -lnetcdff -lnetcdf -L${NETCDFPATH}/lib -I${NETCDFPATH}/include -Wl,-rpath,${NETCDFPATH}/lib
#CFLAGS += -lnetcdff -lnetcdf -L${NETCDFPATH}/lib -I${NETCDFPATH}/include -Wl,-rpath,${NETCDFPATH}/lib
NETCDF_LIB = /usr/lib/x86_64-linux-gnu/
NETCDF_INC = /usr/include/
FFLAGS += -lnetcdff -lnetcdf -L${NETCDF_LIB} -I${NETCDF_INC} -Wl,-rpath,${NETCDF_LIB}
CFLAGS += -lnetcdff -lnetcdf -L${NETCDF_LIB} -I${NETCDF_INC} -Wl,-rpath,${NETCDF_LIB}

#############################################################

#---------------------------------------------------------------------------
# you should never need to change any lines below.

# see the MIPSPro F90 manual for more details on some of the file extensions
# discussed here.
# this makefile template recognizes fortran sourcefiles with extensions
# .f, .f90, .F, .F90. Given a sourcefile <file>.<ext>, where <ext> is one of
# the above, this provides a number of default actions:

# make <file>.opt	create an optimization report
# make <file>.o		create an object file
# make <file>.s		create an assembly listing
# make <file>.x		create an executable file, assuming standalone
#			source
# make <file>.i		create a preprocessed file (for .F)
# make <file>.i90	create a preprocessed file (for .F90)

# The macro TMPFILES is provided to slate files like the above for removal.

RM = rm -f
TMPFILES = .*.m *.B *.L *.i *.i90 *.l *.s *.mod *.opt

################### Added by jackson ########################
#ROSE files:
TMPFILES += *.rmod *_postprocessed.?90 rose_*.?90
#############################################################


.SUFFIXES: .F .F90 .H .L .T .f .f90 .h .i .i90 .l .o .s .opt .x

.f.L:
	$(FC) $(FFLAGS) -c -listing $*.f
.f.opt:
	$(FC) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.f
.f.l:
	$(FC) $(FFLAGS) -c $(LIST) $*.f
.f.T:
	$(FC) $(FFLAGS) -c -cif $*.f
.f.o:
	$(FC) $(FFLAGS) -c $*.f
.f.s:
	$(FC) $(FFLAGS) -S $*.f
.f.x:
	$(FC) $(FFLAGS) -o $*.x $*.f *.o $(LDFLAGS)
.f90.L:
	$(FC) $(FFLAGS) -c -listing $*.f90
.f90.opt:
	$(FC) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.f90
.f90.l:
	$(FC) $(FFLAGS) -c $(LIST) $*.f90
.f90.T:
	$(FC) $(FFLAGS) -c -cif $*.f90
.f90.o:
	$(FC) $(FFLAGS) -c $*.f90
.f90.s:
	$(FC) $(FFLAGS) -c -S $*.f90
.f90.x:
	$(FC) $(FFLAGS) -o $*.x $*.f90 *.o $(LDFLAGS)
.F.L:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -listing $*.F
.F.opt:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.F
.F.l:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $(LIST) $*.F
.F.T:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -cif $*.F
.F.f:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -EP $*.F > $*.f
.F.i:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -P $*.F
.F.o:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $*.F
.F.s:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -S $*.F
.F.x:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -o $*.x $*.F *.o $(LDFLAGS)
.F90.L:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -listing $*.F90
.F90.opt:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -opt_report_level max -opt_report_phase all -opt_report_file $*.opt $*.F90
.F90.l:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $(LIST) $*.F90
.F90.T:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -cif $*.F90
.F90.f90:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -EP $*.F90 > $*.f90
.F90.i90:
	$(FC) $(CPPDEFS) $(FPPFLAGS) -P $*.F90
.F90.o:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c $*.F90
.F90.s:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -c -S $*.F90
.F90.x:
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) -o $*.x $*.F90 *.o $(LDFLAGS)

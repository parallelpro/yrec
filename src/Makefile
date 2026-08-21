# YREC Makefile
#
# Usage:
# Several build targets are supported.
#   release : The default target, maximum performance, warning and optimization flags only
#   devel   : standard development flags, display all warnings, trap exceptions
#   debug   : Default flags + debug flags, address sanitizer enabled
#   prof    : Allows use of gprof profiling tool + devel flags + optimization
#
# Several compilers are supported.
#   gfortran : GNU Fortran compiler, the default
#   ifort    : legacy Intel Fortran compiler
#   ifx      : newer Intel compiler may also be used by defining FC=<compiler> on
# the command line when invoking make.
#
#  make [FC=<compiler>] <target>
#
# Examples:
# gfortran, optimized release build:
#    make
# gfortran, debug build
#    make debug
# gfortran, profiling build
#    make prof
# ifort, development build:
#    make FC=ifort devel
# ifx, release build:
#    make FC=ifx
#
# To remove executable and all intermediate object files:
#    make clean

# The name for the executable
NAME = yrec

# GNU Make picks the default goal (what bare `make`/`make -j4` builds) as
# the target of the *first* rule with a real target list in this file --
# not necessarily the `default:` target below. The MODULE_OBJS order-only
# rule further down expands to a long list of plain object-file targets,
# and would otherwise become that first rule, silently hijacking the
# default goal to a single arbitrary .o file instead of linking $(NAME).
# Pin it explicitly so bare `make` always does what the `default:` target
# below says.
.DEFAULT_GOAL := default

# The default path into which the 'install' target will try to copy the executable.
PREFIX ?= /usr/local/bin

# Compose list of all object files to produce. Sources are organized into
# subdirectories by physics domain (eos/, kap/, rotation/, etc.) -- find
# every .f90 file recursively and keep each .o alongside its source.
SRCS := $(shell find . -name '*.f90' | sed 's|^\./||')
OBJS := $(patsubst %.f90,%.o,$(SRCS))

# Files that define a Fortran MODULE must be compiled before anything that
# `use`s them, so the compiler can find the resulting .mod file. GNU Make
# has no built-in Fortran module dependency scanning, so this list is
# maintained by hand -- add a source here whenever it gains a `module`
# statement. The order-only prerequisite below (declared once all module
# sources are known) ensures every module object is built first, safely
# under parallel (-j) builds, without forcing spurious rebuilds.
MODULE_SRCS := const/const_lib.f90 const/luout_lib.f90 const/intpar_lib.f90 state/oldmod_lib.f90 state/scrtch_lib.f90 state/turnover_lib.f90 state/light_burn_lib.f90 state/engeb_diag_lib.f90 state/fluxes_lib.f90 state/mdphy_lib.f90 state/envelope_comp_lib.f90 state/temp_lib.f90 state/envstruct_lib.f90 state/pulse_diag_lib.f90 state/atmstruct_lib.f90 state/temp2_lib.f90 state/run_diag_lib.f90 state/rotdiff_lib.f90 state/atm_table_lib.f90 state/mhd_eos_lib.f90 state/opacity_table_lib.f90 nuclear/nuclear_lib.f90 numerics/numerics_lib.f90
MODULE_OBJS := $(patsubst %.f90,%.o,$(MODULE_SRCS))

$(filter-out $(MODULE_OBJS),$(OBJS)): | $(MODULE_OBJS)

# nuclear_lib.f90's rates/eqburn routines `use` const_lib/oldmod_lib --
# the block rule above only orders non-module objects after the full
# set of module objects, not module objects relative to each other.
nuclear/nuclear_lib.o: | const/const_lib.o state/oldmod_lib.o

# Flags for GNU Fortran (gfortran)
GNU_BASE_FFLAGS = -ffixed-line-length-132 -Wall -cpp -DGIT_HASH="'$(REPO_HASH)'"
GNU_DEV_FFLAGS = -fcheck=all -ffpe-trap=invalid,zero,overflow
GNU_DBG_FFLAGS = -g -fsanitize=undefined -fsanitize=address -fstack-protector-strong -O0
GNU_PROF_FFLAGS = -O3 -g -pg
GNU_REL_FFLAGS = -O3

# Flags for Intel Fortran (ifort)
IFORT_BASE_FFLAGS = -132 -traceback -warn all -nogen-interfaces -fpp -DGIT_HASH="'$(REPO_HASH)'"
IFORT_DEV_FFLAGS = -check all -fpe:0
IFORT_DBG_FFLAGS = -g
IFORT_REL_FFLAGS = -O3

# Flags for newer Intel Fortran compiler (ifx)
IFX_BASE_FFLAGS = -132 -traceback -warn all -nogen-interfaces -fpp -DGIT_HASH="'$(REPO_HASH)'"
IFX_DEV_FFLAGS = -check all -fpe:0
IFX_DBG_FFLAGS = -g
IFX_REL_FFLAGS = -O3

# Assign flags for the requested compiler
ifeq ("$(FC)", "ifort")
	FFLAGS = $(IFORT_BASE_FFLAGS)
	DEV_FFLAGS = $(IFORT_DEV_FFLAGS)
	DBG_FFLAGS = $(IFORT_DBG_FFLAGS)
	REL_FFLAGS = $(IFORT_REL_FFLAGS)
else ifeq ("$(FC)", "ifx")
	FFLAGS = $(IFX_BASE_FFLAGS)
	DEV_FFLAGS = $(IFX_DEV_FFLAGS)
	DBG_FFLAGS = $(IFX_DBG_FFLAGS)
	REL_FFLAGS = $(IFX_REL_FFLAGS)
else
	FC = gfortran
	FFLAGS = $(GNU_BASE_FFLAGS)
	DEV_FFLAGS = $(GNU_DEV_FFLAGS)

	DBG_FFLAGS = $(GNU_DBG_FFLAGS)
	DBG_LDFLAGS = $(GNU_DBG_FFLAGS)

	PROF_FFLAGS = $(GNU_PROF_FFLAGS)
	PROF_LDFLAGS = $(GNU_PROF_FFLAGS)

	REL_FFLAGS = $(GNU_REL_FFLAGS)
endif

# Test for presence of git command for version hash collection.
GIT_AVAIL := $(shell command -v git)
ifdef GIT_AVAIL
REPO_DIRTY := $(shell command git diff --name-only --quiet .; echo $$?)
SHORT_HASH := $(shell command git rev-parse --short HEAD)
ifeq ($(REPO_DIRTY),1)
  REPO_HASH=$(SHORT_HASH)-dirty
else
  REPO_HASH=$(SHORT_HASH)
endif
else
  REPO_HASH=unknown
endif

# Explicit rule for free-form modern Fortran sources (.f90); GNU Make's
# built-in implicit rules don't reliably cover this suffix.
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Targets
default: FFLAGS += $(REL_FFLAGS)
default: $(NAME)

devel: FFLAGS += $(DEV_FFLAGS)
devel: $(NAME)

debug: FFLAGS += $(DEV_FFLAGS) $(DBG_FFLAGS)
debug: LDFLAGS += $(DBG_LDFLAGS)
debug: $(NAME)

release: FFLAGS += $(REL_FFLAGS)
release: $(NAME)

prof: FFLAGS += $(DEV_FFLAGS) $(PROF_FFLAGS)
prof: LDFLAGS += $(PROF_LDFLAGS)
prof: $(NAME)

$(NAME): $(OBJS)
	$(FC) $(LDFLAGS) -o $(NAME) $(OBJS)

# Remove all object files and executable
clean:
	@rm -f $(NAME)
	@rm -f $(OBJS)

install:
	install -d $(PREFIX)
	install -m 755 $(NAME) $(PREFIX)

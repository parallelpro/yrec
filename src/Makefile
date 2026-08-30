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
# Per-domain test/ directories hold standalone test PROGRAMs (phase
# three, ROADMAP.md stage 2) and are excluded here so the yrec link
# never sees a second `program` unit; they build via their own targets
# below.
SRCS := $(shell find . -name '*.f90' -not -path '*/test/*' | sed 's|^\./||')
OBJS := $(patsubst %.f90,%.o,$(SRCS))

# Fortran module and include dependencies are GENERATED -- deps.mk
# defines MODULE_SRCS (every module-defining source), an order-only
# rule for every object on the module objects whose .mod files it
# reads, and real prerequisites on include files. Regenerate with
#   python3 tools/gen_makefile_deps.py
# whenever a module, use statement, or include is added or removed
# (CI verifies with --check). The blanket order-only rule below stays
# as coarse belt-and-braces on top of deps.mk's precise rules.
include deps.mk

# the default output column selections are GENERATED from the
# authoritative defaults/{history,profile}_columns.list (uncommented
# names, file order) -- edit a .list and make reruns the generator,
# and deps.mk makes io/yrec_output.o depend on the include
io/default_columns.inc: defaults/history_columns.list \
                        defaults/profile_columns.list \
                        tools/gen_default_columns.py
	python3 tools/gen_default_columns.py
MODULE_OBJS := $(patsubst %.f90,%.o,$(MODULE_SRCS))

$(filter-out $(MODULE_OBJS),$(OBJS)): | $(MODULE_OBJS)

# Flags for GNU Fortran (gfortran)
GNU_BASE_FFLAGS = $(HDF5_FFLAGS) -ffixed-line-length-132 -ffree-line-length-none -finit-local-zero -Wall -cpp -DGIT_HASH="'$(REPO_HASH)'" -fPIC

# Optional HDF5 support (GSM / GYRE-HDF5 pulse output): build with
# `make USE_HDF5=1`. HDF5_DIR defaults to the MESA SDK, which ships
# static HDF5 libraries and gfortran modules; point it elsewhere for
# a different HDF5 installation (needs hdf5.mod/h5lt.mod +
# libhdf5*_fortran). Without USE_HDF5, write_gsm_pulse compiles to a
# stub that reports the missing capability.
ifdef USE_HDF5
HDF5_DIR ?= /Applications/mesasdk
HDF5_FFLAGS = -DYREC_USE_HDF5 -I$(HDF5_DIR)/include
HDF5_LIBS = $(HDF5_DIR)/lib/libhdf5hl_fortran.a $(HDF5_DIR)/lib/libhdf5_hl.a \
	$(HDF5_DIR)/lib/libhdf5_fortran.a $(HDF5_DIR)/lib/libhdf5.a -lz -ldl -lm
else
HDF5_FFLAGS =
HDF5_LIBS =
endif
# Rebuild the GSM writer whenever the HDF5 mode flips (make has no
# flag tracking; without this a stale HDF5-flavored object links into
# a non-HDF5 build and vice versa).
HDF5_MODE_STAMP = .hdf5_mode_$(if $(USE_HDF5),on,off)
$(HDF5_MODE_STAMP):
	@rm -f .hdf5_mode_on .hdf5_mode_off
	@touch $@
io/write_gsm_pulse.o: $(HDF5_MODE_STAMP)
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

# Per-domain static libraries (phase three, ROADMAP.md stage 2),
# following MESA's one-library-per-module convention. `make libs`
# builds all of them; the yrec executable itself still links the flat
# object list (the transition to library-based linking can go domain
# by domain). Standalone per-domain test programs link against these.
DOMAIN_DIRS := eos kap atm net wind mixing rotation numerics \
	state setup io core util

define LIB_template
lib$(1).a: $$(filter $(1)/%,$$(OBJS))
	ar rcs $$@ $$^
endef
$(foreach d,$(DOMAIN_DIRS),$(eval $(call LIB_template,$(d))))

LIBS := $(foreach d,$(DOMAIN_DIRS),lib$(d).a)

.PHONY: libs
libs: $(LIBS)

# Library-linked executable (2026, ROADMAP "Library-based yrec link"):
# links core/main.o against the per-domain archives instead of the
# flat object list. Proves every archive is complete and the domain
# layout self-consistent at link time (an object missing from its
# domain library surfaces as an undefined symbol here, which the flat
# link would silently paper over). The archive list is repeated so
# single-pass archive resolution converges regardless of dependency
# order. The default `yrec` binary stays flat-linked (byte-stable
# object order); CI builds this target as the layout check.
YREC_LIBS_ORDER := libstate.a libio.a libnumerics.a libeos.a libkap.a \
	libatm.a libnet.a libmixing.a librotation.a libwind.a \
	libsetup.a libutil.a libcore.a
yrec_libs: core/main.o $(LIBS)
	$(FC) $(LDFLAGS) -o $@ core/main.o \
	$(YREC_LIBS_ORDER) $(YREC_LIBS_ORDER) $(YREC_LIBS_ORDER) \
	$(HDF5_LIBS)

# Standalone per-domain test programs (phase three, ROADMAP.md stage
# 2). Each links the full object set minus core/main.o (the flat-link
# transitional form; library-based linking can follow the yrec
# executable's own transition). Compiled with the release flags so
# pinned expected outputs are flag-stable.
test_eos: FFLAGS += $(REL_FFLAGS)
test_eos: eos/test/test_eos.o $(filter-out core/main.o,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(HDF5_LIBS)

test_kap: FFLAGS += $(REL_FFLAGS)
test_kap: kap/test/test_kap.o $(filter-out core/main.o,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(HDF5_LIBS)

test_atm: FFLAGS += $(REL_FFLAGS)
test_atm: atm/test/test_atm.o $(filter-out core/main.o,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(HDF5_LIBS)

test_reentry: FFLAGS += $(REL_FFLAGS)
test_reentry: core/test/test_reentry.o $(filter-out core/main.o,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(HDF5_LIBS)

test_net: FFLAGS += $(REL_FFLAGS)
test_net: net/test/test_net.o $(filter-out core/main.o,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^ $(HDF5_LIBS)

.PHONY: tests
tests: test_eos test_kap test_atm test_net test_reentry

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
	$(FC) $(LDFLAGS) -o $(NAME) $(OBJS) $(HDF5_LIBS)

# libyrec (2026): the embeddable engine as a shared library -- every
# object except the CLI wrapper. Entry point: core/yrec_capi.f90's
# C-interoperable yrec_run, consumed by the pyyrec ctypes package.
# (-fPIC is in GNU_BASE_FFLAGS: a no-op on Darwin, required on Linux.)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
SHLIB := libyrec.dylib
else
SHLIB := libyrec.so
endif
$(SHLIB): $(OBJS)
	$(FC) $(LDFLAGS) -shared -o $(SHLIB) $(filter-out core/main.o,$(OBJS)) $(HDF5_LIBS)

lib: $(SHLIB)
.PHONY: lib

# Remove all object files, libraries, and executable.
# *.mod removal is load-bearing: a module .mod that survives clean can
# be picked up by a parallel (-j) module-phase compile BEFORE its
# regenerating module object builds, silently baking stale derived-type
# offsets into the consumer (2026: run_log.o read a shifted
# star%ctrl%report_solver_diagnostics this way).
clean:
	@rm -f $(NAME)
	@rm -f *.mod
	@rm -f $(OBJS)
	@rm -f $(LIBS)
	@rm -f libyrec.dylib libyrec.so
	@rm -f */test/*.o test_eos test_kap test_atm test_net test_reentry

install:
	install -d $(PREFIX)
	install -m 755 $(NAME) $(PREFIX)

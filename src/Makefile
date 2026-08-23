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

# Files that define a Fortran MODULE must be compiled before anything that
# `use`s them, so the compiler can find the resulting .mod file. GNU Make
# has no built-in Fortran module dependency scanning, so this list is
# maintained by hand -- add a source here whenever it gains a `module`
# statement. The order-only prerequisite below (declared once all module
# sources are known) ensures every module object is built first, safely
# under parallel (-j) builds, without forcing spurious rebuilds.
MODULE_SRCS := const/phys_const_lib.f90 const/controls_lib.f90 const/controls_reset_lib.f90 const/const_lib.f90 const/luout_lib.f90 const/intpar_lib.f90 atm/envstruct_lib.f90 atm/atmstruct_lib.f90 eos/scv_eos_lib.f90 eos/eos_mixture_lib.f90 kap/conductive_table_lib.f90 state/star_info_lib.f90 atm/atm_table_lib.f90 eos/mhd_eos_lib.f90 kap/opacity_table_lib.f90 eos/opal_eos_lib.f90 eos/yale_eos_lib.f90 eos/eos_lib.f90 kap/kap_lib.f90 numerics/numerics_lib.f90 net/net_lib.f90 core/burn_lib.f90 atm/atm_lib.f90 core/envint_lib.f90 core/yrec_reset.f90 core/yrec_capi.f90 io/yrec_output.f90
MODULE_OBJS := $(patsubst %.f90,%.o,$(MODULE_SRCS))

$(filter-out $(MODULE_OBJS),$(OBJS)): | $(MODULE_OBJS)

# Module objects that `use` other module objects need an explicit
# order-only prerequisite too -- the block rule above only orders
# non-module objects after the full set of module objects, not module
# objects relative to each other, and MODULE_SRCS list position alone
# does not guarantee build order under parallel (-j) builds.
const/const_lib.o: | const/phys_const_lib.o const/controls_lib.o
# parmin embeds the generated new-style inlist reader
core/parmin.o: core/inlist_new_decl.inc core/inlist_new_read.inc
core/yrec_capi.o: | state/star_info_lib.o
io/yrec_output.o: | const/const_lib.o const/luout_lib.o state/star_info_lib.o
numerics/numerics_lib.o: | const/const_lib.o const/intpar_lib.o \
	const/luout_lib.o
state/star_info_lib.o: | eos/scv_eos_lib.o eos/eos_mixture_lib.o \
	kap/conductive_table_lib.o
net/net_lib.o: | const/const_lib.o state/star_info_lib.o \
	const/luout_lib.o numerics/numerics_lib.o
eos/eos_lib.o: | const/const_lib.o eos/opal_eos_lib.o \
	eos/scv_eos_lib.o eos/eos_mixture_lib.o \
	eos/yale_eos_lib.o const/luout_lib.o
kap/kap_lib.o: | const/const_lib.o state/star_info_lib.o \
	kap/opacity_table_lib.o \
	const/luout_lib.o
const/controls_reset_lib.o: | const/controls_lib.o const/phys_const_lib.o
core/burn_lib.o: | net/net_lib.o state/star_info_lib.o \
	const/const_lib.o const/luout_lib.o
core/yrec_reset.o: | state/star_info_lib.o const/controls_reset_lib.o
atm/atm_lib.o: | atm/atm_table_lib.o atm/atmstruct_lib.o \
	state/star_info_lib.o \
	const/const_lib.o atm/envstruct_lib.o \
	eos/eos_lib.o const/intpar_lib.o kap/kap_lib.o \
	const/luout_lib.o numerics/numerics_lib.o
core/envint_lib.o: | atm/atm_lib.o state/star_info_lib.o \
	eos/eos_lib.o kap/kap_lib.o numerics/numerics_lib.o \
	atm/atm_table_lib.o atm/atmstruct_lib.o atm/envstruct_lib.o \
	const/const_lib.o const/luout_lib.o const/intpar_lib.o

# Flags for GNU Fortran (gfortran)
GNU_BASE_FFLAGS = -ffixed-line-length-132 -ffree-line-length-none -finit-local-zero -Wall -cpp -DGIT_HASH="'$(REPO_HASH)'" -fPIC
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
	const state setup io core util

define LIB_template
lib$(1).a: $$(filter $(1)/%,$$(OBJS))
	ar rcs $$@ $$^
endef
$(foreach d,$(DOMAIN_DIRS),$(eval $(call LIB_template,$(d))))

LIBS := $(foreach d,$(DOMAIN_DIRS),lib$(d).a)

.PHONY: libs
libs: $(LIBS)

# Standalone per-domain test programs (phase three, ROADMAP.md stage
# 2). Each links the full object set minus core/main.o (the flat-link
# transitional form; library-based linking can follow the yrec
# executable's own transition). Compiled with the release flags so
# pinned expected outputs are flag-stable.
test_eos: FFLAGS += $(REL_FFLAGS)
test_eos: eos/test/test_eos.o $(filter-out core/main.o,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^

test_kap: FFLAGS += $(REL_FFLAGS)
test_kap: kap/test/test_kap.o $(filter-out core/main.o,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^

test_atm: FFLAGS += $(REL_FFLAGS)
test_atm: atm/test/test_atm.o $(filter-out core/main.o,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^

test_reentry: FFLAGS += $(REL_FFLAGS)
test_reentry: core/test/test_reentry.o $(filter-out core/main.o,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^

test_net: FFLAGS += $(REL_FFLAGS)
test_net: net/test/test_net.o $(filter-out core/main.o,$(OBJS))
	$(FC) $(LDFLAGS) -o $@ $^

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
	$(FC) $(LDFLAGS) -o $(NAME) $(OBJS)

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
	$(FC) $(LDFLAGS) -shared -o $(SHLIB) $(filter-out core/main.o,$(OBJS))

lib: $(SHLIB)
.PHONY: lib

# Remove all object files, libraries, and executable
clean:
	@rm -f $(NAME)
	@rm -f $(OBJS)
	@rm -f $(LIBS)
	@rm -f libyrec.dylib libyrec.so

install:
	install -d $(PREFIX)
	install -m 755 $(NAME) $(PREFIX)

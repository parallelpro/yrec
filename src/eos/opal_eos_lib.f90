!----------------------------------------------------------------------
! opal_eos_lib
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Bundles the OPAL 1995/2001/2006 equation-of-
! state interpolation tables -- former common/a/, aa/, b/, bb/, e/,
! ee/, eee/, rmpopeos/ (1995), aeos/, aaeos/, beos/, bbeos/, eeos/,
! eeeos/, eeeeos/, rmpopeos01/ (2001), aeos06/, aaeos06/, beos06/,
! bbeos06/, eeeos06/, eeeeos06/ (2006), and the single common/lreadco/
! flag shared by name across all three years' readers -- into
! per-vintage derived-type instances.
!
! Readability W3 (2026): each vintage's mutable state is one instance
! (opal95 / opal01 / opal06) whose members carry no year suffix; the
! interpolation routines receive the instance as their first argument
! ("v"). The read-only tables that used to sit inside the state as
! declaration-time initializers (x_grid, eos_var_order, t6_index_lo,
! density_index_edge*) are named constants below, one per vintage,
! transcribed verbatim from eos/blkdta000.f90 (1995 -- since deleted),
! eos/opal/esac01.f90 + eos/opal/readcoeos01.f90 (2001),
! eos/opal/esac06.f90 (2006), and eos/opal/readco.f90 (1995's rmpopeos
! ramp table); every repeat-count list's element count was
! independently verified against its table's nr/nt dimension. The
! 2001/2006 x_grid constants are intentionally written WITHOUT a d0
! suffix (0.0, 0.2, 0.4, 0.6, 0.8), exactly matching the original
! esac01.f/esac06.f DATA statements: those literals are parsed as
! single precision and then widened to double, which is NOT
! bit-identical to the correctly-rounded double value -- preserved
! verbatim rather than "fixed", per this project's policy on
! pre-existing behavior.
!
! Two state types, not one: the 2001 and 2006 interpolators
! (esac01/esac06) clamp every stencil to the table (t6_order_hi /
! density_order_hi), so they share opal_eos_vintage, dimensioned by the
! larger 2006 table. The 1995 interpolator (esac.f90) does not: its
! boundary-sum stencils run to t6_index_1+3 = nt+1 and to x_index_lo+3
! = 6 (> n_eos_mx) whenever T6 lies in the lowest table interval
! (5000 K <= T < 5500 K, reachable in any 1995-EOS run), and the
! interpolation that follows reads t6_grid(nt+1), i.e. the member that
! happens to follow t6_grid. Re-dimensioning the 1995 arrays would
! therefore change 1995 results, so opal_eos_vintage95 keeps the 1995
! dimensions and member order exactly as they were laid out before.
module opal_eos_lib
      implicit none
      integer, parameter :: n_eos_mx = 5, n_eos_mv = 10
      integer, parameter :: n_eos95_nr = 77, n_eos95_nt = 56
      integer, parameter :: n_eos01_nr = 169, n_eos01_nt = 191
      integer, parameter :: n_eos06_nr = 169, n_eos06_nt = 197
! Largest 2001/2006 table dimensions: the physical size of every
! opal_eos_vintage table array (the 2001 loops run to n_eos01_*).
      integer, parameter :: n_eos_nr_max = n_eos06_nr, n_eos_nt_max = n_eos06_nt
! Vintage ids: indices of opal_eos%readco_init_flag / table_metal_fraction.
      integer, parameter :: n_opal_vintages = 3
      integer, parameter :: iv_opal95 = 1, iv_opal01 = 2, iv_opal06 = 3
! Variable numbers of the OPAL EOS results (the slots of
! opal95%eos_output / opal06%eos_output and the arguments of their
! eos_index_inverse): the 1995 and 2006 tables carry the same ten
! variables in the same order -- this block is the layout descriptor
! of both vintages' eos_output. Units as documented in esac.f90
! (pressure in 1e12 dyn/cm**2, energy in 1e12 erg/g, ...).
      integer, parameter :: i_opal_p = 1            ! pressure
      integer, parameter :: i_opal_e = 2            ! internal energy
      integer, parameter :: i_opal_s = 3            ! entropy
      integer, parameter :: i_opal_dedrho = 4       ! dE/dRho at constant T6
      integer, parameter :: i_opal_cv = 5           ! dE/dT6 at constant volume
      integer, parameter :: i_opal_chi_rho = 6      ! dlogP/dlogRho at constant T6
      integer, parameter :: i_opal_chi_t = 7        ! dlogP/dlogT6 at constant Rho
      integer, parameter :: i_opal_gamma1 = 8       ! gamma1
      integer, parameter :: i_opal_gamma2_ratio = 9 ! gamma2/(gamma2-1)
      integer, parameter :: i_opal_gamma3m1 = 10    ! gamma3-1
! The 2001 tables (opal01%eos_output / eos_index_inverse) carry no
! dE/dRho variable, so every slot after entropy is one lower -- the
! layout descriptor of the 2001 eos_output.
      integer, parameter :: i_opal01_p = 1
      integer, parameter :: i_opal01_e = 2
      integer, parameter :: i_opal01_s = 3
      integer, parameter :: i_opal01_cv = 4
      integer, parameter :: i_opal01_chi_rho = 5
      integer, parameter :: i_opal01_chi_t = 6
      integer, parameter :: i_opal01_gamma1 = 7
      integer, parameter :: i_opal01_gamma2_ratio = 8
      integer, parameter :: i_opal01_gamma3m1 = 9
! Density returned by rhoofp/rhoofp01/rhoofp06 when the (P,T6) point
! is not in the table (the oeqos* callers test for .le. -998).
      double precision, parameter :: opal_rho_not_found = -999.0d0
! Value of the table_loaded_flag / readco_init_flag members once the
! corresponding tables have been read.
      integer, parameter :: opal_flag_set = 12345678
! Loop variable for the repeat-count array constructors used below to
! transcribe the original DATA statements' "N*value" runs -- must be a
! separately-declared integer in scope for the implied-do to resolve
! under `implicit none`, even though it is otherwise unused. Private
! so it can't collide with any `use`ing file's own local `i`.
      integer, private :: i

! ==================== read-only per-vintage tables ====================
! Hydrogen-fraction grid of each vintage's table set (former common/a/
! x_grid, common/aeos/ x_grid_01, common/aeos06/ x_grid_06).
      double precision, parameter :: opal95_x_grid(n_eos_mx) = &
           [0.0d0,0.2d0,0.4d0,0.6d0,0.8d0]
! NOTE: no D-suffix, as in the original (esac01.f/esac06.f) -- preserved
! verbatim: parsed as single-precision constants and widened to double
! precision (0.2/0.4/0.6/0.8 are not exactly representable in binary
! floating point, so these are NOT bit-identical to the correctly-
! rounded double values).
      double precision, parameter :: opal01_x_grid(n_eos_mx) = [0.0, 0.2, 0.4, 0.6, 0.8]
      double precision, parameter :: opal06_x_grid(n_eos_mx) = [0.0, 0.2, 0.4, 0.6, 0.8]
! Table-column -> eos_output-slot permutation (the identity for all
! three vintages; former common/b/ eos_var_order and its _01/_06 twins).
      integer, parameter :: opal_eos_var_order(n_eos_mv) = [1,2,3,4,5,6,7,8,9,10]
! Row of the lowest tabulated T6 at each density column (former
! common/b/ t6_index_lo and its _01/_06 twins).
      integer, parameter :: opal95_t6_index_lo(n_eos95_nr) = &
           [(56,i=1,37),(54,i=1,7), &
            39,37,36,34,33,31,31,30,29,28,27,26,26,25,24,23,22,21,21,20, &
            19,18,17,16,15,15,13,13,11,11,9,9,7]
      integer, parameter :: opal01_t6_index_lo(n_eos01_nr) = &
           [(191,i=1,92),190,189,188,187,186,185,184,174,(134,i=1,4), &
            (133,i=1,3),(132,i=1,2),98,92,(85,i=1,2),(77,i=1,2),71, &
            (63,i=1,3),(59,i=1,2),53,51,(46,i=1,2),(44,i=1,9),(38,i=1,3), &
            (33,i=1,6),(29,i=1,16),27,26,25,23,22,20,19,18,17,16]
      integer, parameter :: opal06_t6_index_lo(n_eos06_nr) = &
           [(197,i=1,87),(191,i=1,7),190,(189,i=1,2),185,179,170, &
            (149,i=1,2),133,125,123,122,120,115,113,107,102,(80,i=1,2), &
            72,68,66,64,62,56,54,52,51,(50,i=1,2),49,47,(45,i=1,2),43,42, &
            (40,i=1,28),39,37,36,35,34,32,31,30,29,27,26]
! Column of the highest tabulated density at each T6 row (former
! common/rmpopeos/ density_index_edge_at_t, common/rmpopeos01/
! density_index_edge_at_t_01, common/beos06/ density_index_edge_06).
      integer, parameter :: opal95_density_index_edge(n_eos95_nt) = &
           [(77,i=1,7),(76,i=1,2),(74,i=1,2),(72,i=1,2),(70,i=1,2), &
            68,67,66,65,64,63,61,60,59,58,57,55,54,53,52,51, &
            (49,i=1,2),48,(47,i=1,2),46,(45,i=1,2),(44,i=1,15),(37,i=1,2)]
      integer, parameter :: opal01_density_index_edge(n_eos01_nt) = &
           [(169,i=1,16),168,167,166,165,(164,i=1,2),163,(162,i=1,2), &
            161,160,(159,i=1,2),(143,i=1,4),(137,i=1,5),(134,i=1,6), &
            (125,i=1,2),(123,i=1,5),(122,i=1,2),(121,i=1,6),(119,i=1,4), &
            (116,i=1,8),(115,i=1,9),(113,i=1,5),(111,i=1,7),(110,i=1,6), &
            (109,i=1,34),107,104,(100,i=1,40),(99,i=1,10),98,97,96,95,94,93,92]
      integer, parameter :: opal06_density_index_edge(n_eos06_nt) = &
           [(169,i=1,26),168,(167,i=1,2),166,165,164,(163,i=1,2),162, &
            161,160,(159,i=1,2),158,(130,i=1,2),129,(128,i=1,2), &
            (126,i=1,2),(125,i=1,2),124,122,121,(120,i=1,2),(119,i=1,2), &
            (118,i=1,6),(117,i=1,2),(116,i=1,2),(115,i=1,2),(114,i=1,4), &
            (113,i=1,8),(111,i=1,22),(110,i=1,5),(109,i=1,6),(108,i=1,2), &
            (107,i=1,5),(106,i=1,2),105,(104,i=1,2),(103,i=1,8), &
            (102,i=1,16),(100,i=1,21),(99,i=1,9),(98,i=1,6),(97,i=1,4), &
            95,94,(87,i=1,6)]

! ==================== mutable state ====================
! quad's per-slot Lagrange-coefficient cache (recomputed when
! recompute_flag==0, read on every other call) -- promoted from SAVEd
! function locals in the save-migration campaign; one per vintage.
      type, public :: opal_quad_cache
            double precision :: x12_inv(30), x13_inv(30), x23_inv(30), &
                 x1_squared(30), x1_plus_x2(30)
      end type opal_quad_cache

! State of one OPAL 2001 or 2006 table set (see the header note on why
! 1995 has its own type). Table arrays are dimensioned by the 2006
! table; the 2001 code addresses only its first n_eos01_nr/n_eos01_nt
! rows/columns.
      type, public :: opal_eos_vintage
            type(opal_quad_cache) :: quad
! former common/aeos/ and aeos06/
           double precision :: eos_table(n_eos_mx,n_eos_mv,n_eos_nt_max,n_eos_nr_max), &
                t6_list(n_eos_nr_max,n_eos_nt_max), density_grid(n_eos_nr_max), &
                t6_grid(n_eos_nt_max), x_interp_result(n_eos_nt_max,n_eos_nr_max), &
                x_interp_result_alt(n_eos_nt_max,n_eos_nr_max), &
                x_grid_spacing_inv(n_eos_mx), t6_grid_spacing_inv(n_eos_nt_max), &
                density_grid_spacing_inv(n_eos_nr_max)
! x_loop_index is a plain local (not a module field): it is used as a
! DO-loop control variable in esac*.f90/readco*.f90, which the Fortran
! standard requires to be a named scalar variable, not a derived-type
! component.
           integer :: x_index_lo
! former common/aaeos/ and aaeos06/ (xxh is not referenced anywhere)
           double precision :: rho_interp_hi(4), rho_interp_lo(4), xxh
! former common/beos/ and beos06/
           double precision :: z_table(n_eos_mx)
           integer :: eos_index_inverse(10)
! former common/bbeos/ and bbeos06/
           integer :: density_index_1, density_index_2, density_index_3, &
                density_index_4, t6_index_1, t6_index_2, t6_index_3, &
                t6_index_4, t6_interp_order, density_interp_order
! former common/eeos/ and eeos06/ -- the latter relocated here from
! atm_table_lib in 2026 (ROADMAP.md stage 1); its former non-eos reader
! (a direct eos_output(8)/(9) read in the turnover-timescale code,
! since removed) was replaced by eos_lib's eos_get_gamma1.
           double precision :: esact, eos_output(n_eos_mv)
! former common/eeeos/ and eeeos06/
           double precision :: x_interp_workspace(n_eos_mx,n_eos_nt_max,n_eos_nr_max), &
                x_grid_copy(n_eos_mx)
! former common/eeeeos/ and eeeeos06/ (amu_grid and log10_ne_grid are
! read from the 2006 tables only)
           double precision :: moles_per_gram_table(n_eos_mx), &
                hydrogen_fraction_header(n_eos_mx), &
                mean_molecular_weight_header(n_eos_mx), &
                amu_grid(n_eos_nr_max,n_eos_nt_max), &
                log10_ne_grid(n_eos_nr_max,n_eos_nt_max), &
                density_grid_table(n_eos_mx,n_eos_nr_max), &
                species_fraction_header(n_eos_mx,6), &
                log10_r_value(n_eos_nr_max,n_eos_nt_max)
           integer :: temperature_count_used(n_eos_mx,n_eos_nr_max)
! former common/rmpopeos01/ (2001 only: eqbound06 works from the
! density_index_edge table directly)
           double precision :: density_edge_at_t(n_eos_nt_max)
           integer :: t_row_index
      end type opal_eos_vintage

! State of the OPAL 1995 table set: the same members as
! opal_eos_vintage (minus the 2006-only amu/log10_ne grids) at the
! 1995 dimensions, in the member order the 1995 stencil overruns rely
! on (header note).
      type, public :: opal_eos_vintage95
            type(opal_quad_cache) :: quad
! former common/a/
           double precision :: eos_table(n_eos_mx,n_eos_mv,n_eos95_nt,n_eos95_nr), &
                t6_list(n_eos95_nr,n_eos95_nt), density_grid(n_eos95_nr), &
                t6_grid(n_eos95_nt), x_interp_result(n_eos95_nt,n_eos95_nr), &
                x_interp_result_alt(n_eos95_nt,n_eos95_nr), &
                x_grid_spacing_inv(n_eos_mx), t6_grid_spacing_inv(n_eos95_nt), &
                density_grid_spacing_inv(n_eos95_nr)
           integer :: x_index_lo
! former common/aa/ (xxh is not referenced anywhere)
           double precision :: rho_interp_hi(4), rho_interp_lo(4), xxh
! former common/b/
           double precision :: z_table(n_eos_mx)
           integer :: eos_index_inverse(10)
! former common/bb/
           integer :: density_index_1, density_index_2, density_index_3, &
                density_index_4, t6_index_1, t6_index_2, t6_index_3, &
                t6_index_4, t6_interp_order, density_interp_order
! former common/e/
           double precision :: esact, eos_output(n_eos_mv)
! former common/ee/
           double precision :: x_interp_workspace(n_eos_mx,n_eos95_nt,n_eos95_nr), &
                x_grid_copy(n_eos_mx)
! former common/eee/
           double precision :: moles_per_gram_table(n_eos_mx), &
                hydrogen_fraction_header(n_eos_mx), &
                mean_molecular_weight_header(n_eos_mx), &
                density_grid_table(n_eos_mx,n_eos95_nr), &
                species_fraction_header(n_eos_mx,6), &
                log10_r_value(n_eos95_nr,n_eos95_nt)
           integer :: temperature_count_used(n_eos_mx,n_eos95_nr)
! former common/rmpopeos/
           double precision :: density_edge_at_t(n_eos95_nt)
           integer :: t_row_index
      end type opal_eos_vintage95

! Flags shared by, or indexed by, vintage. Kept apart from the vintage
! instances so that those (48 MB apiece) stay initializer-free and
! zero-filled at load rather than becoming initialized data.
      type, public :: opal_eos_state
! former common/lreadco/: shared by COMMON block name across all three
! years' readers/interpolators in the original.
           integer :: table_loaded_flag
! save-migration campaign (2026): lazy-load guards and first-load
! table-Z memory, promoted from SAVEd locals; indexed by iv_opal95/01/06.
            integer :: readco_init_flag(n_opal_vintages) = 0
            double precision :: table_metal_fraction(n_opal_vintages) = 0.0d0
      end type opal_eos_state

      type(opal_eos_state), public :: opal_eos
      type(opal_eos_vintage95), public :: opal95
      type(opal_eos_vintage), public :: opal01, opal06

end module opal_eos_lib

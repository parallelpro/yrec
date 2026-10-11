!----------------------------------------------------------------------
! opal_eos_tables
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
! esac01 + readcoeos01 (2001),
! esac06 (2006), and readco (1995's rmpopeos
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
! larger 2006 table. The 1995 interpolator (esac) does not: its
! boundary-sum stencils run to t6_index_1+3 = nt+1 and to x_index_lo+3
! = 6 (> n_eos_mx) whenever T6 lies in the lowest table interval
! (5000 K <= T < 5500 K, reachable in any 1995-EOS run), and the
! interpolation that follows reads t6_grid(nt+1), i.e. the member that
! happens to follow t6_grid. Re-dimensioning the 1995 arrays would
! therefore change 1995 results, so opal_eos_vintage95 keeps the 1995
! dimensions and member order exactly as they were laid out before.
module opal_eos_tables
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use math_lib
      implicit none
      private
      public :: find_opal_ramp, gmass01, i_opal01_chi_rho, i_opal01_chi_t, &
           i_opal01_cv, i_opal01_e, i_opal01_gamma1, i_opal01_gamma2_ratio, &
           i_opal01_gamma3m1, i_opal01_p, i_opal01_s, i_opal_chi_rho, &
           i_opal_chi_t, i_opal_cv, i_opal_dedrho, i_opal_e, i_opal_gamma1, &
           i_opal_gamma2_ratio, i_opal_gamma3m1, i_opal_p, i_opal_s, &
           interpolate_opal_t6_rho, iv_opal01, iv_opal06, iv_opal95, &
           n_eos01_nr, n_eos01_nt, n_eos06_nr, n_eos06_nt, n_eos95_nr, &
           n_eos95_nt, n_eos_mv, n_eos_mx, opal01, opal01_density_index_edge, &
           opal01_t6_index_lo, opal01_x_grid, opal06, &
           opal06_density_index_edge, opal06_t6_index_lo, opal06_x_grid, &
           opal95, opal95_density_index_edge, opal95_t6_index_lo, &
           opal95_x_grid, opal_eos, opal_eos_var_order, opal_eos_vintage, &
           opal_eos_vintage95, opal_flag_set, opal_rho_not_found, quad
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
! of both vintages' eos_output. Units as documented in esac
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
! is not in the table (the compute_opal95_eos* callers test for .le. -998).
      real(dp), parameter :: opal_rho_not_found = -999.0d0
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
      real(dp), parameter :: opal95_x_grid(n_eos_mx) = &
           [0.0d0,0.2d0,0.4d0,0.6d0,0.8d0]
! NOTE: no D-suffix, as in the original (esac01.f/esac06.f) -- preserved
! verbatim: parsed as single-precision constants and widened to double
! precision (0.2/0.4/0.6/0.8 are not exactly representable in binary
! floating point, so these are NOT bit-identical to the correctly-
! rounded double values).
      real(dp), parameter :: opal01_x_grid(n_eos_mx) = [0.0, 0.2, 0.4, 0.6, 0.8]
      real(dp), parameter :: opal06_x_grid(n_eos_mx) = [0.0, 0.2, 0.4, 0.6, 0.8]
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
      type :: opal_quad_cache
            real(dp) :: x12_inv(30), x13_inv(30), x23_inv(30), &
                 x1_squared(30), x1_plus_x2(30)
      end type opal_quad_cache

! State of one OPAL 2001 or 2006 table set (see the header note on why
! 1995 has its own type). Table arrays are dimensioned by the 2006
! table; the 2001 code addresses only its first n_eos01_nr/n_eos01_nt
! rows/columns.
      type :: opal_eos_vintage
            type(opal_quad_cache) :: quad
! former common/aeos/ and aeos06/
           real(dp) :: eos_table(n_eos_mx,n_eos_mv,n_eos_nt_max,n_eos_nr_max), &
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
           real(dp) :: rho_interp_hi(4), rho_interp_lo(4), xxh
! former common/beos/ and beos06/
           real(dp) :: z_table(n_eos_mx)
           integer :: eos_index_inverse(10)
! former common/bbeos/ and bbeos06/
           integer :: density_index_1, density_index_2, density_index_3, &
                density_index_4, t6_index_1, t6_index_2, t6_index_3, &
                t6_index_4, t6_interp_order, density_interp_order
! former common/eeos/ and eeos06/ -- the latter relocated here from
! atm_tables in 2026 (docs/ROADMAP.md stage 1); its former non-eos reader
! (a direct eos_output(8)/(9) read in the turnover-timescale code,
! since removed) was replaced by eos_lib's eos_get_gamma1.
           real(dp) :: esact, eos_output(n_eos_mv)
! former common/eeeos/ and eeeos06/
           real(dp) :: x_interp_workspace(n_eos_mx,n_eos_nt_max,n_eos_nr_max), &
                x_grid_copy(n_eos_mx)
! former common/eeeeos/ and eeeeos06/ (amu_grid and log10_ne_grid are
! read from the 2006 tables only)
           real(dp) :: moles_per_gram_table(n_eos_mx), &
                hydrogen_fraction_header(n_eos_mx), &
                mean_molecular_weight_header(n_eos_mx), &
                amu_grid(n_eos_nr_max,n_eos_nt_max), &
                log10_ne_grid(n_eos_nr_max,n_eos_nt_max), &
                density_grid_table(n_eos_mx,n_eos_nr_max), &
                species_fraction_header(n_eos_mx,6), &
                log10_r_value(n_eos_nr_max,n_eos_nt_max)
           integer :: temperature_count_used(n_eos_mx,n_eos_nr_max)
! former common/rmpopeos01/ (2001 only: find_opal2006_ramp works from the
! density_index_edge table directly)
           real(dp) :: density_edge_at_t(n_eos_nt_max)
           integer :: t_row_index
      end type opal_eos_vintage

! State of the OPAL 1995 table set: the same members as
! opal_eos_vintage (minus the 2006-only amu/log10_ne grids) at the
! 1995 dimensions, in the member order the 1995 stencil overruns rely
! on (header note).
      type :: opal_eos_vintage95
            type(opal_quad_cache) :: quad
! former common/a/
           real(dp) :: eos_table(n_eos_mx,n_eos_mv,n_eos95_nt,n_eos95_nr), &
                t6_list(n_eos95_nr,n_eos95_nt), density_grid(n_eos95_nr), &
                t6_grid(n_eos95_nt), x_interp_result(n_eos95_nt,n_eos95_nr), &
                x_interp_result_alt(n_eos95_nt,n_eos95_nr), &
                x_grid_spacing_inv(n_eos_mx), t6_grid_spacing_inv(n_eos95_nt), &
                density_grid_spacing_inv(n_eos95_nr)
           integer :: x_index_lo
! former common/aa/ (xxh is not referenced anywhere)
           real(dp) :: rho_interp_hi(4), rho_interp_lo(4), xxh
! former common/b/
           real(dp) :: z_table(n_eos_mx)
           integer :: eos_index_inverse(10)
! former common/bb/
           integer :: density_index_1, density_index_2, density_index_3, &
                density_index_4, t6_index_1, t6_index_2, t6_index_3, &
                t6_index_4, t6_interp_order, density_interp_order
! former common/e/
           real(dp) :: esact, eos_output(n_eos_mv)
! former common/ee/
           real(dp) :: x_interp_workspace(n_eos_mx,n_eos95_nt,n_eos95_nr), &
                x_grid_copy(n_eos_mx)
! former common/eee/
           real(dp) :: moles_per_gram_table(n_eos_mx), &
                hydrogen_fraction_header(n_eos_mx), &
                mean_molecular_weight_header(n_eos_mx), &
                density_grid_table(n_eos_mx,n_eos95_nr), &
                species_fraction_header(n_eos_mx,6), &
                log10_r_value(n_eos95_nr,n_eos95_nt)
           integer :: temperature_count_used(n_eos_mx,n_eos95_nr)
! former common/rmpopeos/
           real(dp) :: density_edge_at_t(n_eos95_nt)
           integer :: t_row_index
      end type opal_eos_vintage95

! Flags shared by, or indexed by, vintage. Kept apart from the vintage
! instances so that those (48 MB apiece) stay initializer-free and
! zero-filled at load rather than becoming initialized data.
      type :: opal_eos_state
! whether each vintage's table is loaded, indexed by iv_opal95/01/06
! (2026-10: one flag shared by all three, from the original common
! /lreadco/, so loading one vintage kept the next from reading its own
! table in the same process)
           integer :: table_loaded_flag(n_opal_vintages) = 0
! the unit eos_init opened the run's OPAL table on (2026-10:
! read_controls opened it on ctrl%iopale)
           integer :: table_unit = 0
! save-migration campaign (2026): lazy-load guards and first-load
! table-Z memory, promoted from SAVEd locals; indexed by iv_opal95/01/06.
            integer :: readco_init_flag(n_opal_vintages) = 0
            real(dp) :: table_metal_fraction(n_opal_vintages) = 0.0d0
      end type opal_eos_state

      type(opal_eos_state) :: opal_eos
      type(opal_eos_vintage95) :: opal95
      type(opal_eos_vintage) :: opal01, opal06


contains

!----------------------------------------------------------------------
! interpolate_opal_t6_rho (formerly t6rint_core)
!----------------------------------------------------------------------
! Readability W3 (2026): the T6/density interpolation shared by the
! OPAL 2001 and 2006 EOS (t6rinteos01 / t6rinteos06). Those two
! routines were identical statement for statement apart from the state
! instance they addressed, the RECURSIVE attribute of the 2006 copy,
! and the text of the out-of-range message; the arithmetic lives here
! once, and each vintage's wrapper keeps its own message. Not used by
! the 1995 EOS (t6rinterp orders its esactq2 block differently).
!
! Interpolates the already X-interpolated table slice
! (v%x_interp_result) in density (slr) and T6 (slt) with the stencil
! indices esac01/esac06 stored in v (density_index_*, t6_index_*,
! *_interp_order), leaving the result in v%esact. The caller tests
! v%esact against the out-of-range sentinel.
subroutine interpolate_opal_t6_rho(v, slr, slt)


      type(opal_eos_vintage), intent(inout) :: v

      real(dp), intent(in) :: slr, slt

! --- locals ---
      integer :: hi_loop_count, recompute_flag, cache_slot, t6_grid_idx
      real(dp) :: esactq, esact2, esactq2, dix, dix2

      hi_loop_count = 0
      recompute_flag = 0

      do t6_grid_idx = v%t6_index_1, v%t6_index_1 + v%t6_interp_order
         cache_slot = 1
         hi_loop_count = hi_loop_count + 1
         v%rho_interp_lo(hi_loop_count) = quad(v%quad, recompute_flag, cache_slot, &
              slr, v%x_interp_result(t6_grid_idx,v%density_index_1), &
              v%x_interp_result(t6_grid_idx,v%density_index_2), &
              v%x_interp_result(t6_grid_idx,v%density_index_3), &
              v%density_grid(v%density_index_1), v%density_grid(v%density_index_2), &
              v%density_grid(v%density_index_3))
         if (v%density_interp_order == 3) then
            cache_slot = 2
            v%rho_interp_hi(hi_loop_count) = quad(v%quad, recompute_flag, cache_slot, &
                 slr, v%x_interp_result(t6_grid_idx,v%density_index_2), &
                 v%x_interp_result(t6_grid_idx,v%density_index_3), &
                 v%x_interp_result(t6_grid_idx,v%density_index_4), &
                 v%density_grid(v%density_index_2), v%density_grid(v%density_index_3), &
                 v%density_grid(v%density_index_4))
         end if
         recompute_flag = 1
      end do

      recompute_flag = 0
      cache_slot = 1
! ..... eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      v%esact = quad(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_lo(1), &
           v%rho_interp_lo(2), v%rho_interp_lo(3), v%t6_grid(v%t6_index_1), &
           v%t6_grid(v%t6_index_2), v%t6_grid(v%t6_index_3))
      if (v%density_interp_order == 3) then
! .....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
         esactq = quad(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_hi(1), &
              v%rho_interp_hi(2), v%rho_interp_hi(3), v%t6_grid(v%t6_index_1), &
              v%t6_grid(v%t6_index_2), v%t6_grid(v%t6_index_3))
      end if
      if (v%t6_interp_order == 3) then
! .....    eos(i) in lower-left 3x3.
         esact2 = quad(v%quad, recompute_flag, cache_slot, slt, v%rho_interp_lo(2), &
              v%rho_interp_lo(3), v%rho_interp_lo(4), v%t6_grid(v%t6_index_2), &
              v%t6_grid(v%t6_index_3), v%t6_grid(v%t6_index_4))
! .....    eos(i) smoothed in left 3x4
         dix = (v%t6_grid(v%t6_index_3) - slt)*v%t6_grid_spacing_inv(v%t6_index_3)
         v%esact = v%esact*dix + esact2*(1.0d0 - dix)
! endif   ! moved to loc a
         if (v%density_interp_order == 3) then

! .....     eos(i) in upper-right 3x3.
            esactq2 = quad(v%quad, recompute_flag, cache_slot, slt, &
                 v%rho_interp_hi(2), v%rho_interp_hi(3), v%rho_interp_hi(4), &
                 v%t6_grid(v%t6_index_2), v%t6_grid(v%t6_index_3), &
                 v%t6_grid(v%t6_index_4))
            esactq = esactq*dix + esactq2*(1.0d0 - dix)
         end if
      end if  ! loc a
!
      if (v%density_interp_order == 3) then
         dix2 = (v%density_grid(v%density_index_3) - slr)* &
              v%density_grid_spacing_inv(v%density_index_3)
         if (v%t6_interp_order == 3) then
! .....        eos(i) smoothed in both log(T6) and log(R)
            v%esact = v%esact*dix2 + esactq*(1.0d0 - dix2)
         end if
      end if

end subroutine interpolate_opal_t6_rho

!----------------------------------------------------------------------
! find_opal_ramp (formerly eqbound_core)
!----------------------------------------------------------------------
! Readability W3 (2026): the table-edge search and ramp-factor
! calculation shared by find_opal95_ramp (OPAL 1995) and find_opal2001_ramp (OPAL
! 2001). Those two routines were identical statement for statement
! from the T6 row search onward, differing only in the table they
! address and its row count nt; each wrapper keeps its own density
! pre-check (5.0 for 1995, 7.0 for 2001) and then calls this. Not used
! by find_opal2006_ramp, whose edge search is a different algorithm.
!
! Arguments: the vintage's row/column counts nt/nr, t6_grid (nt
! descending T6 rows), the caller's cached row index t_row_index
! (updated here), density_edge_at_t (table-edge density per row),
! density_index_edge (column of that edge per row) and density_grid
! (the nr density columns), then the point (t6, log10_density) and the
! find_opal95_ramp outputs.
subroutine find_opal_ramp(nt, nr, t6_grid, t_row_index, density_edge_at_t, &
     density_index_edge, density_grid, t6, log10_density, ramp_factor, &
     in_opal_table, needs_ramp, ierr)


      integer, intent(in) :: nt, nr
      real(dp), intent(in) :: t6_grid(:), density_edge_at_t(:)
      integer, intent(in) :: density_index_edge(:)
      real(dp), intent(in) :: density_grid(:)
      integer, intent(inout) :: t_row_index

      real(dp), intent(in) :: t6, log10_density
      real(dp), intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp
! --- locals ---
      real(dp) :: table_edge_density, ramp_start_density
      real(dp) :: t6_top_of_table
      integer :: t6_scan_idx

      integer, intent(out) :: ierr

      ierr = 0

!     find nearest table element in t.
      if (t6 < t6_grid(t_row_index)) then
         do t6_scan_idx = t_row_index+1, nt
            if (t6 >= t6_grid(t6_scan_idx)) then
               t_row_index = t6_scan_idx - 1
               exit
            end if
         end do
         if (t6_scan_idx > nt) then
         t6_top_of_table = t6_grid(nt)
!        caller should have stopped outside table bounds; error exit
         write(*,5) t6, t6_top_of_table, t_row_index
    5    format(' ERROR IN OPAL EOS: OUTSIDE TABLE IN T6',2F10.6,I5)
         ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
         end if
      else
         do t6_scan_idx = t_row_index, 1, -1
            if (t6 <= t6_grid(t6_scan_idx)) then
               t_row_index = t6_scan_idx
               exit
            end if
         end do
         if (t6_scan_idx < (1)) then
         t6_top_of_table = t6_grid(1)
!        caller should have stopped outside table bounds; error exit
         write(*,5) t6, t6_top_of_table, t_row_index
         ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
         end if
      end if
!     define table edge in rho by linear interpolation.
      table_edge_density = density_edge_at_t(t_row_index+1)
      table_edge_density = log10(table_edge_density)
!     define beginning of ramp in the same fashion -
!     ramp is defined as one table element wide.
      ramp_start_density = density_grid(density_index_edge(t_row_index+1)-1)
      ramp_start_density = log10(ramp_start_density)
!     check if within table bounds in rho
      if (log10_density > table_edge_density) then
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         return
      end if

!     If we get here, the point is in the table.
      in_opal_table = .true.

!     Now we check if ramping is needed.
!     First we check if ramping in temperature is needed.
      if (t6 <= t6_grid(nt-1)) then
         needs_ramp = .true.
         ramp_factor = (t6-t6_grid(nt))/(t6_grid(nt-1)-t6_grid(nt))
      else if (t6 >= t6_grid(2)) then
         needs_ramp = .true.
         ramp_factor = (t6_grid(1)-t6)/(t6_grid(1)-t6_grid(2))
      else if (log10_density >= ramp_start_density) then
!        If we get here, ramping in density is needed.
         needs_ramp = .true.
         ramp_factor = (table_edge_density-log10_density)/ &
              (table_edge_density-ramp_start_density)
      else
         needs_ramp = .false.
      end if

end subroutine find_opal_ramp

!----------------------------------------------------------------------
! quad
!----------------------------------------------------------------------
! Quadratic interpolation helper shared by the OPAL 1995/2001/2006 EOS
! table interpolators (esac*.f90/t6rint*.f90). Given three (x,y) pairs,
! returns the value of the interpolating quadratic at x. When
! recompute_flag is 0 the quadratic's coefficients are (re)computed and
! cached in slot cache_slot of the caller's opal_quad_cache (each
! vintage instance carries its own, opal95%quad etc.); otherwise the
! cached coefficients for that slot are reused as-is (the caller is
! responsible for knowing the cache is still valid for the grid points
! in question -- this is a straight transliteration of the original
! caching scheme, not a redesign of it).
!
! Readability W3 (2026): the former quadeos01.f90 and quadeos06.f90
! clones were identical to this routine apart from the cache they
! addressed (opal_eos%quad01_* / quad06_*; now the cache argument) and
! a redundant pair of parentheses around the single variable
! x1_plus_x2(cache_slot) in the 2006 copy; all three vintages now call
! this one function with their own opalXX%quad cache.
real(dp) function quad(cache, recompute_flag, cache_slot, eval_point, &
     y1, y2, y3, x1, x2, x3)


      type(opal_quad_cache), intent(inout) :: cache
      integer, intent(in) :: recompute_flag, cache_slot
      real(dp), intent(in) :: eval_point, y1, y2, y3, x1, x2, x3

      real(dp) :: grid_x(3), grid_y(3)
      real(dp) :: coef2, coef1, coef0
      grid_x(1) = x1
      grid_x(2) = x2
      grid_x(3) = x3
      grid_y(1) = y1
      grid_y(2) = y2
      grid_y(3) = y3
      if (recompute_flag == 0) then
         cache%x12_inv(cache_slot) = 1.0d0/(grid_x(1) - grid_x(2))
         cache%x13_inv(cache_slot) = 1.0d0/(grid_x(1) - grid_x(3))
         cache%x23_inv(cache_slot) = 1.0d0/(grid_x(2) - grid_x(3))
         cache%x1_squared(cache_slot) = grid_x(1)*grid_x(1)
         cache%x1_plus_x2(cache_slot) = grid_x(1) + grid_x(2)
      end if
      coef2 = (grid_y(1) - grid_y(2))*cache%x12_inv(cache_slot)
      coef2 = coef2 - (grid_y(2) - grid_y(3))*cache%x23_inv(cache_slot)
      coef2 = coef2*cache%x13_inv(cache_slot)
      coef1 = (grid_y(1) - grid_y(2))*cache%x12_inv(cache_slot) - &
           cache%x1_plus_x2(cache_slot)*coef2
      coef0 = grid_y(1) - grid_x(1)*coef1 - cache%x1_squared(cache_slot)*coef2
      quad = coef0 + eval_point*(coef1 + eval_point*coef2)
end function quad

!----------------------------------------------------------------------
! gmass01
!----------------------------------------------------------------------
! OPAL 2001 EOS analogue of gmass (see there for the general
! description); called from esac01 and, since the 2026 W3
! readability pass removed the identical gmass06.f90 copy, from
! esac06 as well. Note the difference from
! gmass: here the function accumulator starts at 0, so the
! electron-mass term (electron_mole_excess*atomic_weight(1) in gmass)
! is not included. Preserved verbatim rather than "fixed".
real(dp) function gmass01(hydrogen_fraction, metal_fraction, &
     total_moles, ground_state_energy, metal_mole_fraction, &
     species_number_fraction)


      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(out) :: total_moles, ground_state_energy, &
           metal_mole_fraction
! species_number_fraction(1:6): number (mole) fractions of Ne, O, N, C,
! He, H (same order as atomic_number); slot 7 is never assigned.
      real(dp), intent(out) :: species_number_fraction(:)

! atomic_number(i): atomic number (= electrons per fully-ionized atom)
! of species i, for i=1,6 = Ne, O, N, C, He, H.
      real(dp) :: atomic_number(6)
! ionization_energy(i): total ground-state ionization energy (eV) of
! species i, same ordering as atomic_number.
      real(dp) :: ionization_energy(6)
! atomic_weight(i): atomic weight (amu), i=1 electron, 2 Ne, 3 O, 4 N,
! 5 C, 6 He, 7 H.
      real(dp) :: atomic_weight(7)
! solar-relative mole-fraction-per-unit-Z of C, N, O, Ne in the metal
! mix (fixed composition assumed for all Z).
      real(dp) :: carbon_mix_fraction, nitrogen_mix_fraction, &
           oxygen_mix_fraction, neon_mix_fraction
      real(dp) :: carbon_moles, nitrogen_moles, oxygen_moles, &
           neon_moles
      real(dp) :: hydrogen_moles, helium_moles, total_moles_raw
      integer :: species_idx

      data (ionization_energy(species_idx), species_idx=1,6) &
           /-3394.873554d0, -1974.86545d0, -1433.92718d0, -993.326315d0, &
           -76.1959403d0, -15.29409d0/
      data (atomic_number(species_idx), species_idx=1,6) &
           /10.0d0, 8.0d0, 7.0d0, 6.0d0, 2.0d0, 1.0d0/
      carbon_mix_fraction = 0.247137766d0
      nitrogen_mix_fraction = 0.0620782d0
      oxygen_mix_fraction = 0.52837118d0
      neon_mix_fraction = 0.1624188d0
      atomic_weight(7) = 1.0079d0
      atomic_weight(6) = 4.0026d0
      atomic_weight(5) = 12.011d0
      atomic_weight(4) = 14.0067d0
      atomic_weight(3) = 15.9994d0
      atomic_weight(2) = 20.179d0
      atomic_weight(1) = 0.00054858d0
      metal_mole_fraction = metal_fraction/(carbon_mix_fraction* &
           atomic_weight(5) + nitrogen_mix_fraction*atomic_weight(4) + &
           oxygen_mix_fraction*atomic_weight(3) + neon_mix_fraction* &
           atomic_weight(2))
      carbon_moles = metal_mole_fraction*carbon_mix_fraction
      nitrogen_moles = metal_mole_fraction*nitrogen_mix_fraction
      oxygen_moles = metal_mole_fraction*oxygen_mix_fraction
      neon_moles = metal_mole_fraction*neon_mix_fraction
      hydrogen_moles = hydrogen_fraction/atomic_weight(7)
      helium_moles = (1.0d0 - hydrogen_fraction - metal_fraction)/ &
           atomic_weight(6)
      total_moles_raw = hydrogen_moles + helium_moles + carbon_moles + &
           nitrogen_moles + oxygen_moles + neon_moles
      species_number_fraction(6) = hydrogen_moles/total_moles_raw
      species_number_fraction(5) = helium_moles/total_moles_raw
      species_number_fraction(4) = carbon_moles/total_moles_raw
      species_number_fraction(3) = nitrogen_moles/total_moles_raw
      species_number_fraction(2) = oxygen_moles/total_moles_raw
      species_number_fraction(1) = neon_moles/total_moles_raw
      ground_state_energy = 0.0d0
      total_moles = 0.0d0
      do species_idx = 1, 6
         ground_state_energy = ground_state_energy + &
              ionization_energy(species_idx)*species_number_fraction(species_idx)
         total_moles = total_moles + (1.0d0 + atomic_number(species_idx))* &
              species_number_fraction(species_idx)
      end do
      gmass01 = 0.0d0
      do species_idx = 2, 7
         gmass01 = gmass01 + atomic_weight(species_idx)* &
              species_number_fraction(species_idx - 1)
      end do

end function gmass01

end module opal_eos_tables

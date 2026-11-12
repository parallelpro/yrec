!----------------------------------------------------------------------
! opal_eos_tables
!----------------------------------------------------------------------
! The OPAL EOS tables and interpolation state, one instance per
! vintage (opal_tables(iv_opal95), (iv_opal01), (iv_opal06)), plus
! opal_shared for what the vintages
! share: the table unit, the loaded flags and each table's Z. The
! interpolation routines take the instance as their first argument,
! "v". Also here: each vintage's read-only grids (x_grid,
! eos_var_order, t6_index_lo, density_index_edge) and its row of
! opal_layout (everything the code tells the vintages apart by), the
! slots of eos_output (i_opal_*, i_opal01_*), and the helpers the
! vintages share:
!     interpolate_opal_t6_rho  the T6-rho interpolation of the 2001 and
!                              2006 tables
!     find_opal_ramp           table coverage and edge-blend weight for
!                              the 1995 and 2001 tables
!     quad                     OPAL's quadratic interpolation
!
! The three vintages share one state type, opal_eos_vintage, dimensioned
! by the largest (2006) table; each addresses its first nt rows and nr
! columns. The 2001 and 2006 interpolation clamps every stencil to the
! table; the 1995 one does not (M94), and reads past its table what the
! 1995 code's own storage held there (opal_eos's stencil_entry).
!
! The 2001 and 2006 x_grid constants are written without the d0 suffix
! (0.0, 0.2, ...), as in OPAL's DATA statements: single-precision
! literals widened to double, which differ from the correctly rounded
! doubles in the last bits.
module opal_eos_tables
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_configuration, only: opal_vintage_1995, opal_vintage_2001, opal_vintage_2006
      use eos_messages, only: eos_report, message_length, unwritten_line
      use math_lib
      implicit none
      private
      public :: n_eos_nr_max, n_opal_vintages, find_opal_ramp, i_opal01_chi_rho, i_opal01_chi_t, &
           i_opal01_cv, i_opal01_e, i_opal01_gamma1, i_opal01_gamma2_ratio, &
           i_opal01_gamma3m1, i_opal01_p, i_opal01_s, i_opal_chi_rho, &
           i_opal_chi_t, i_opal_cv, i_opal_dedrho, i_opal_e, i_opal_gamma1, &
           i_opal_gamma2_ratio, i_opal_gamma3m1, i_opal_p, i_opal_s, &
           interpolate_opal_t6_rho, iv_opal01, iv_opal06, iv_opal95, &
           n_eos01_nr, n_eos01_nt, n_eos06_nr, n_eos06_nt, n_eos95_nr, &
           n_eos95_nt, n_eos_mv, n_eos_mx, opal_tables, opal01_density_index_edge, &
           opal01_t6_index_lo, opal01_x_grid, &
           opal06_density_index_edge, opal06_t6_index_lo, opal06_x_grid, &
           opal95_density_index_edge, opal95_t6_index_lo, &
           opal95_x_grid, opal_shared, opal_eos_vintage, &
           opal_flag_set, opal_rho_not_found, quad, &
           n_eos_nt_max, opal_empty_entry, opal_empty_sum, opal_usable_corner_sum, &
           opal_beyond_stencil, opal_x_interp_overflow, opal_range_overflow
      integer, parameter :: n_eos_mx = 5, n_eos_mv = 10
      integer, parameter :: n_eos95_nr = 77, n_eos95_nt = 56
      integer, parameter :: n_eos01_nr = 169, n_eos01_nt = 191
      integer, parameter :: n_eos06_nr = 169, n_eos06_nt = 197
! Largest 2001/2006 table dimensions: the physical size of every
! opal_eos_vintage table array (the 2001 loops run to n_eos01_*).
      integer, parameter :: n_eos_nr_max = n_eos06_nr, n_eos_nt_max = n_eos06_nt
! Vintage ids (eos_configuration's opal_vintage_*): the rows of
! opal_layout and the indices of opal_shared's per-vintage flags.
      integer, parameter :: n_opal_vintages = 3
      integer, parameter :: iv_opal95 = opal_vintage_1995, iv_opal01 = opal_vintage_2001, &
           iv_opal06 = opal_vintage_2006
! Variable numbers of the OPAL EOS results (the slots of
! opal95%eos_output and opal_tables(iv_opal06)%eos_output): the 1995 and 2006 tables
! carry the same variables in the same order (2006 without gamma3 - 1).
! Units as documented in esac (pressure in 1e12 dyn/cm**2, energy in
! 1e12 erg/g, ...).
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
! The 2001 tables (opal_tables(iv_opal01)%eos_output) carry no dE/dRho variable, so
! every slot after entropy is one lower.
      integer, parameter :: i_opal01_p = 1
      integer, parameter :: i_opal01_e = 2
      integer, parameter :: i_opal01_s = 3
      integer, parameter :: i_opal01_cv = 4
      integer, parameter :: i_opal01_chi_rho = 5
      integer, parameter :: i_opal01_chi_t = 6
      integer, parameter :: i_opal01_gamma1 = 7
      integer, parameter :: i_opal01_gamma2_ratio = 8
      integer, parameter :: i_opal01_gamma3m1 = 9
! The density opal_density_search returns when the (P, T6) point is not
! in the table, or its search does not converge.
      real(dp), parameter :: opal_rho_not_found = -999.0d0
! The readers fill a table entry the file does not give with
! opal_empty_entry, and the interpolation tells such entries by sums and
! values far above any real one: a coverage sum above opal_empty_sum
! holds one, a 2 x 2 corner sum below opal_usable_corner_sum holds none
! (the stencil can then shift onto it), opal_beyond_stencil stands for
! the sum of a cell the stencil lacks (above opal_empty_sum, so that it
! is not used), and an X-interpolated value beyond opal_x_interp_overflow
! either way (a huge entry with a negative Lagrange weight comes back
! negative, M292), or a T6-rho interpolated one above opal_range_overflow,
! came from one. The T6-rho check is one-sided: the 2006 tables' dE/dRho
! reaches -3.7e14, and every node reaching it has passed the X check.
      real(dp), parameter :: opal_empty_entry = 1.0d+35, opal_empty_sum = 1.0d+30, &
           opal_usable_corner_sum = 1.0d+25, opal_beyond_stencil = 2.0d+30, &
           opal_x_interp_overflow = 1.0d+20, opal_range_overflow = 1.0d+15
! a/3 in OPAL's units, Mb per T6^4 (radiation_constant_over_3 times
! 1e24/1e12), to the ten digits compute_opal_eos has always written:
! the radiation pressure's share of P is this T6^4/P12
      real(dp), parameter, public :: opal_radiation_constant_over_3 = 2.521971383d-3
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
! Hydrogen-fraction grid of each vintage's table set.
      real(dp), parameter :: opal95_x_grid(n_eos_mx) = &
           [0.0d0,0.2d0,0.4d0,0.6d0,0.8d0]
! Single-precision literals, as in OPAL (see the module header).
      real(dp), parameter :: opal01_x_grid(n_eos_mx) = [0.0, 0.2, 0.4, 0.6, 0.8]
      real(dp), parameter :: opal06_x_grid(n_eos_mx) = [0.0, 0.2, 0.4, 0.6, 0.8]
! Row of the lowest tabulated T6 at each density column.
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
! Column of the highest tabulated density at each T6 row.
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
! The 2001 density search's highest density column at each T6 row: the
! values of opal01_density_index_edge except at rows 78-80 (113 here,
! 115 there), as in OPAL's rhoofp01.
      integer, parameter :: opal01_search_density_index_edge(n_eos01_nt) = &
           [(169,i=1,16), 168, 167, 166, 165, (164,i=1,2), 163, (162,i=1,2), 161, 160, &
            (159,i=1,2), (143,i=1,4), (137,i=1,5), (134,i=1,6), (125,i=1,2), (123,i=1,5), &
            (122,i=1,2), (121,i=1,6), (119,i=1,4), (116,i=1,8), (115,i=1,6), (113,i=1,8), &
            (111,i=1,7), (110,i=1,6), (109,i=1,34), 107, 104, (100,i=1,40), (99,i=1,10), &
            98, 97, 96, 95, 94, 93, 92]
      integer, parameter :: opal06_density_index_edge(n_eos06_nt) = &
           [(169,i=1,26),168,(167,i=1,2),166,165,164,(163,i=1,2),162, &
            161,160,(159,i=1,2),158,(130,i=1,2),129,(128,i=1,2), &
            (126,i=1,2),(125,i=1,2),124,122,121,(120,i=1,2),(119,i=1,2), &
            (118,i=1,6),(117,i=1,2),(116,i=1,2),(115,i=1,2),(114,i=1,4), &
            (113,i=1,8),(111,i=1,22),(110,i=1,5),(109,i=1,6),(108,i=1,2), &
            (107,i=1,5),(106,i=1,2),105,(104,i=1,2),(103,i=1,8), &
            (102,i=1,16),(100,i=1,21),(99,i=1,9),(98,i=1,6),(97,i=1,4), &
            95,94,(87,i=1,6)]

! the ground-state ionization energies [eV] of Ne, O, N, C, He and H in
! the 2001 and 2006 codes' mass per particle (gmass01)
      real(dp), parameter :: gmass01_ionization_energy(6) = [-3394.873554d0, -1974.86545d0, &
           -1433.92718d0, -993.326315d0, -76.1959403d0, -15.29409d0]

! ==================== what tells the vintages apart ====================
! One row per vintage (opal_layout(iv_opal95), (iv_opal01),
! (iv_opal06)); the code reads its vintage's row instead of branching
! on the vintage (only the messages keep each vintage's own words).
      type, public :: opal_vintage_layout
! the table's density columns and T6 rows, and its hydrogen-fraction grid
         integer :: nr, nt
         real(dp) :: x_grid(n_eos_mx)
! the row of the lowest tabulated T6 at each density column, and the
! column of the highest density at each T6 row (zero past nr, nt)
         integer :: t6_index_lo(n_eos_nr_max), density_index_edge(n_eos_nt_max)
! how many variables an entry carries (esac's deriv_order for a full
! evaluation), and their eos_output slots (0 where the table lacks the
! variable: the 2001 table has no dE/dRho, so its later slots are one
! lower, and only the 1995 and 2001 ones carry gamma3 - 1)
         integer :: n_vars
         integer :: i_p, i_e, i_s, i_dedrho, i_cv, i_chi_rho, i_chi_t, i_gamma1, &
              i_gamma2_ratio, i_gamma3m1
! the radiation correction's a/c [Mb/T6^4] and molar gas constant
! [Mb cm^3/(mol T6)], as each vintage's code writes them (the 2001
! code's are single-precision literals, widened); whether it computes
! the gammas from gas plus radiation directly and rescales P and S by
! energy_unit_ratio, the ratio of N_A k in erg/eV (0.9648575e12) to the
! older value the table was computed with (1995), or adds to the table's
! gammas the difference that radiation makes to the same formulas
! (2001, 2006: add_radiation)
         real(dp) :: rad_const_over_c, molar_gas_constant_mbcc
         logical :: radiation_gammas_direct
         real(dp) :: energy_unit_ratio
! the mean mass per particle of the radiation correction (opal_eos's
! opal_mass_per_mole): the ground-state ionization energies [eV] of Ne,
! O, N, C, He and H, and whether the electrons' mass counts (1995)
         real(dp) :: ionization_energy(6)
         logical :: counts_electron_mass
! whether the interpolation clamps its stencils to the table (2001,
! 2006; not 1995, M94)
         logical :: clamps_stencil
! the reader's formats for an X block's header and for a T6 row, and the
! columns it skips before the variables (1995 and 2001: log10 R; 2006:
! amu and log10 n_e)
         character(len=48) :: header_format, row_format
         integer :: n_skipped_columns
! the prefix of the interpolator's messages
         character(len=15) :: routine_id
! compute_opal_eos's T6 range, and whether it checks that the
! interpolated pressure reproduces P (1995)
         real(dp) :: t6_min, t6_max
         logical :: checks_pressure
! the density search (opal_density_search; 1995: rhoofp): its highest
! density column at each T6 row, its a/c [Mb/T6^4] (2001 and 2006:
! OPAL's single-precision literal, widened), whether it subtracts the
! radiation pressure from the target and asks the interpolator for the
! gas alone (2006) or adds it to the table's pressure bounds and asks
! for gas plus radiation (1995, 2001), and its name in messages
         integer :: search_density_index_edge(n_eos_nt_max)
         real(dp) :: search_rad_const_over_c
         logical :: search_subtracts_radiation
         character(len=8) :: search_name
! the table's edge (find_opal_edge_ramp): the 2006 table's ragged edge
! (by column and row, find_ragged_edge_ramp), else one table edge in
! rho per T6 row (find_opal_ramp) below the lowest log10 rho the ramp
! test admits
         logical :: ragged_edge
         real(dp) :: ramp_log10_rho_min
! blend_eos's window, inside which the vintage is consulted at all: T in
! K, and the 1995 test of log10 T (the original literals, kept as
! written: 1995 bounds log10 T by 8, the others T by 1e8 and 2e8 K)
         real(dp) :: window_t_min, window_t_max, window_log10_t_max, window_log10_rho_max
      end type opal_vintage_layout
      type(opal_vintage_layout), parameter, public :: opal_layout(n_opal_vintages) = [ &
           opal_vintage_layout(nr=n_eos95_nr, nt=n_eos95_nt, x_grid=opal95_x_grid, &
                t6_index_lo=[opal95_t6_index_lo, (0, i=n_eos95_nr+1, n_eos_nr_max)], &
                density_index_edge=[opal95_density_index_edge, (0, i=n_eos95_nt+1, n_eos_nt_max)], &
                n_vars=10, i_p=1, i_e=2, i_s=3, i_dedrho=4, i_cv=5, i_chi_rho=6, i_chi_t=7, &
                i_gamma1=8, i_gamma2_ratio=9, i_gamma3m1=10, &
                rad_const_over_c=1.8914785d-3, molar_gas_constant_mbcc=83.1446304d0, &
                radiation_gammas_direct=.true., energy_unit_ratio=0.9648575d0/0.9652d0, &
                ionization_energy=[-3388.637d0, -1970.918d0, -1431.717d0, -993.2303d0, &
                -76.2315d0, -15.29409d0], counts_electron_mass=.true., clamps_stencil=.false., &
                header_format="(3X,F6.4,3X,F6.4,11X,F10.7,17X,F10.7)", &
                row_format="(F8.4,1X,F6.2,3E13.5,E11.3,6F8.4)", n_skipped_columns=1, routine_id="", &
                t6_min=0.0050d0, t6_max=100.0d0, checks_pressure=.true., &
                search_density_index_edge=[opal95_density_index_edge, (0, i=n_eos95_nt+1, n_eos_nt_max)], &
                search_rad_const_over_c=1.8914785d-3, search_subtracts_radiation=.false., &
                search_name="RHOOFP", ragged_edge=.false., ramp_log10_rho_min=-14.0d0, &
                window_t_min=5.0d3, window_t_max=huge(1.0d0), window_log10_t_max=8.0d0, &
                window_log10_rho_max=5.0d0), &
           opal_vintage_layout(nr=n_eos01_nr, nt=n_eos01_nt, x_grid=opal01_x_grid, &
                t6_index_lo=opal01_t6_index_lo, &
                density_index_edge=[opal01_density_index_edge, (0, i=n_eos01_nt+1, n_eos_nt_max)], &
                n_vars=9, i_p=1, i_e=2, i_s=3, i_dedrho=0, i_cv=4, i_chi_rho=5, i_chi_t=6, &
                i_gamma1=7, i_gamma2_ratio=8, i_gamma3m1=9, &
                rad_const_over_c=real(1.8914785e-3, dp), molar_gas_constant_mbcc=real(83.14510, dp), &
                radiation_gammas_direct=.false., energy_unit_ratio=1.0d0, &
                ionization_energy=gmass01_ionization_energy, counts_electron_mass=.false., &
                clamps_stencil=.true., header_format="(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)", &
                row_format="(F9.5,1X,F6.2,3E13.5,6F8.4)", n_skipped_columns=1, &
                routine_id="OPALEOS/ESAC01:", &
                t6_min=0.0020d0, t6_max=100.0d0, checks_pressure=.false., &
                search_density_index_edge=[opal01_search_density_index_edge, &
                (0, i=n_eos01_nt+1, n_eos_nt_max)], &
                search_rad_const_over_c=real(1.8914785e-3, dp), search_subtracts_radiation=.false., &
                search_name="RHOOFP01", ragged_edge=.false., ramp_log10_rho_min=-14.0d0, &
                window_t_min=2.0d3, window_t_max=100d6, window_log10_t_max=huge(1.0d0), &
                window_log10_rho_max=7.0d0), &
           opal_vintage_layout(nr=n_eos06_nr, nt=n_eos06_nt, x_grid=opal06_x_grid, &
                t6_index_lo=opal06_t6_index_lo, density_index_edge=opal06_density_index_edge, &
                n_vars=9, i_p=1, i_e=2, i_s=3, i_dedrho=4, i_cv=5, i_chi_rho=6, i_chi_t=7, &
                i_gamma1=8, i_gamma2_ratio=9, i_gamma3m1=0, &
                rad_const_over_c=1.8914785d-3, molar_gas_constant_mbcc=83.14510d0, &
                radiation_gammas_direct=.false., energy_unit_ratio=1.0d0, &
                ionization_energy=gmass01_ionization_energy, counts_electron_mass=.false., &
                clamps_stencil=.true., header_format="(3X,F6.4,3X,F12.9,11X,F10.7,17X,F10.7)", &
                row_format="(F11.6,1X,F6.4,E11.4,2E13.6,2E11.3,5F10.6)", n_skipped_columns=2, &
                routine_id="OPALEOS/ESAC06:", &
                t6_min=0.001870d0, t6_max=200.0d0, checks_pressure=.false., &
                search_density_index_edge=opal06_density_index_edge, &
                search_rad_const_over_c=real(1.8914785e-3, dp), search_subtracts_radiation=.true., &
                search_name="RHOOFP06", ragged_edge=.true., ramp_log10_rho_min=-14.0d0, &
                window_t_min=1.870d3, window_t_max=200d6, window_log10_t_max=huge(1.0d0), &
                window_log10_rho_max=7.0d0)]

! ==================== mutable state ====================
! quad's per-slot Lagrange-coefficient cache (recomputed when
! recompute_flag == 0, read on every other call); one per vintage.
      type :: opal_quad_cache
            real(dp) :: x12_inv(30), x13_inv(30), x23_inv(30), &
                 x1_squared(30), x1_plus_x2(30)
      end type opal_quad_cache

! State of one OPAL table set. Table arrays are dimensioned by the 2006
! table; the 1995 and 2001 code address only their first nr/nt
! columns/rows. Table arrays are dimensioned by the 2006
! table; the 2001 code addresses only its first n_eos01_nr/n_eos01_nt
! rows/columns.
      type :: opal_eos_vintage
            type(opal_quad_cache) :: quad
! the tables, their grids and inverse spacings, the X-interpolated slices
           real(dp) :: eos_table(n_eos_mx,n_eos_mv,n_eos_nt_max,n_eos_nr_max), &
                t6_list(n_eos_nr_max,n_eos_nt_max), density_grid(n_eos_nr_max), &
                t6_grid(n_eos_nt_max), x_interp_result(n_eos_nt_max,n_eos_nr_max), &
                x_interp_result_alt(n_eos_nt_max,n_eos_nr_max), &
                x_grid_spacing_inv(n_eos_mx), t6_grid_spacing_inv(n_eos_nt_max), &
                density_grid_spacing_inv(n_eos_nr_max)
           integer :: x_index_lo
! density-interpolation weights
           real(dp) :: rho_interp_hi(4), rho_interp_lo(4)
! each table's Z
           real(dp) :: z_table(n_eos_mx)
! the interpolation stencil and orders
           integer :: density_index_1, density_index_2, density_index_3, &
                density_index_4, t6_index_1, t6_index_2, t6_index_3, &
                t6_index_4, t6_interp_order, density_interp_order
! the last interpolated value (esact) and the results (eos_output)
           real(dp) :: esact, eos_output(n_eos_mv)
! the X-interpolation workspace
           real(dp) :: x_interp_workspace(n_eos_mx,n_eos_nt_max,n_eos_nr_max), &
                x_grid_copy(n_eos_mx)
! each block's X (for the reader's messages), each X table's densities
! and the number of T6 rows each density column has
           real(dp) :: hydrogen_fraction_header(n_eos_mx), &
                density_grid_table(n_eos_mx,n_eos_nr_max)
           integer :: temperature_count_used(n_eos_mx,n_eos_nr_max)
! the table edge at each T6 row and the cached row (find_opal_ramp;
! find_ragged_edge_ramp works from density_index_edge directly)
           real(dp) :: density_edge_at_t(n_eos_nt_max)
           integer :: t_row_index
      end type opal_eos_vintage


! Flags shared by, or indexed by, vintage. Kept apart from the vintage
! instances (15.5 MB apiece) so that those stay initializer-free and
! zero-filled at load rather than becoming initialized data.
      type :: opal_eos_state
! whether each vintage's table is loaded, indexed by iv_opal95/01/06
           integer :: table_loaded_flag(n_opal_vintages) = 0
! the unit open_opal_eos_table opened the run's OPAL table on
           integer :: table_unit = 0
! the readers' first-call guards and each loaded table's Z, indexed
! by iv_opal95/01/06
            integer :: readco_init_flag(n_opal_vintages) = 0
            real(dp) :: table_metal_fraction(n_opal_vintages) = 0.0d0
      end type opal_eos_state

      type(opal_eos_state) :: opal_shared
! the tables' instances, indexed by vintage
      type(opal_eos_vintage) :: opal_tables(iv_opal95:iv_opal06)


contains

!----------------------------------------------------------------------
! interpolate_opal_t6_rho (formerly t6rint_core)
!----------------------------------------------------------------------
! The T6/density interpolation of the OPAL EOS (opal_eos's t6rinteos
! adds each vintage's out-of-range message). (The 1995 code computed the
! upper-right quadratic esactq2 also where the T6 order is 2, with its
! weight dix unset; esactq is used only where that order is 3, so the
! one code gives its results.)
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
! The table-edge search and ramp factor of the OPAL 1995 and 2001
! tables, called by find_opal95_ramp and find_table_edge_ramp after
! their own density checks. find_ragged_edge_ramp searches the 2006
! table's ragged edge differently.
!
! Arguments: the vintage's row count nt, t6_grid (nt
! descending T6 rows), the caller's cached row index t_row_index
! (updated here), density_edge_at_t (table-edge density per row),
! density_index_edge (column of that edge per row) and density_grid
! (the density columns), then the point (t6, log10_density) and the
! find_opal95_ramp outputs.
subroutine find_opal_ramp(nt, t6_grid, t_row_index, density_edge_at_t, &
     density_index_edge, density_grid, t6, log10_density, ramp_factor, &
     in_opal_table, needs_ramp, ierr)


      integer, intent(in) :: nt
      real(dp), intent(in) :: t6_grid(:), density_edge_at_t(:)
      integer, intent(in) :: density_index_edge(:)
      real(dp), intent(in) :: density_grid(:)
      integer, intent(inout) :: t_row_index

      real(dp), intent(in) :: t6, log10_density
      real(dp), intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp

      integer, intent(out) :: ierr
! --- locals ---
      real(dp) :: table_edge_density, ramp_start_density
      real(dp) :: t6_top_of_table
      integer :: t6_scan_idx
      character(len=message_length) :: eos_lines(4)

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
         eos_lines = unwritten_line
         write(eos_lines,5) t6, t6_top_of_table, t_row_index
         call eos_report(eos_lines)
    5    format(' ERROR IN OPAL EOS: OUTSIDE TABLE IN T6',2F10.6,I5)
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
         eos_lines = unwritten_line
         write(eos_lines,5) t6, t6_top_of_table, t_row_index
         call eos_report(eos_lines)
         ierr = 1
         return
         end if
      end if
!     the table's edge in rho: the edge of the row below (in T6) the point
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
! Each vintage passes its own cache (opal95%quad,
! opal_tables(vintage)%quad).
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



!======================================================================
! Sources: the OPAL EOS interpolation codes of Rogers, Swenson &
! Iglesias (1996) and Rogers & Nayfonov (2002). This state was the
! COMMON blocks a, aa, b, bb, e, ee, eee and rmpopeos (1995); aeos,
! aaeos, beos, bbeos, eeos, eeeos, eeeeos and rmpopeos01 (2001); aeos06,
! aaeos06, beos06, bbeos06, eeeos06 and eeeeos06 (2006); and lreadco.

end module opal_eos_tables

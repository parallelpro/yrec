!----------------------------------------------------------------------
! opal_eos_lib
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Bundles the OPAL 1995/2001/2006 equation-of-
! state interpolation tables -- former common/a/, aa/, b/, bb/, e/,
! ee/, eee/, rmpopeos/ (1995), aeos/, aaeos/, beos/, bbeos/, eeos/,
! eeeos/, eeeeos/, rmpopeos01/ (2001), aeos06/, aaeos06/, beos06/,
! bbeos06/, eeeos06/, eeeeos06/ (2006), and the single common/lreadco/
! flag shared by name across all three years' readers -- into one
! derived type. Each EOS year's declaring files used byte-identical
! member names to the OTHER years' (e.g. esac.f90/esac01.f90/
! esac06.f90 all call their table array "eos_table"), even though each
! year's table is a physically distinct block/storage -- disambiguated
! here with the "_01"/"_06" suffixes the codebase already uses to
! distinguish these years elsewhere (eqbound01.f90 vs eqbound06.f90,
! etc.); the unsuffixed 1995-family names are the base case.
!
! DATA-statement defaults (x_grid/eos_var_order/t6_index_lo/
! density_index_edge_at_t/density_index_edge for all three years) are
! preserved as declaration-time initializers, transcribed verbatim
! from eos/blkdta000.f90 (1995 -- since deleted, its DATA statements
! were the only thing it did and became dead once these members moved
! here), eos/opal/esac01.f90 + eos/opal/readcoeos01.f90 (2001),
! eos/opal/esac06.f90 (2006), and eos/opal/readco.f90 (1995's rmpopeos
! ramp table) -- every repeat-count list's element count was
! independently verified against its table's nr/nt dimension before
! use. The 2001/2006 x_grid arrays are intentionally initialized
! WITHOUT a d0 suffix (0.0, 0.2, 0.4, 0.6, 0.8), exactly matching the
! original esac01.f/esac06.f DATA statements: those literals are
! parsed as single precision and then widened to double, which is NOT
! bit-identical to the correctly-rounded double value -- preserved
! verbatim rather than "fixed", per this project's policy on
! pre-existing behavior.
module opal_eos_lib
      implicit none
      integer, parameter :: n_eos_mx = 5, n_eos_mv = 10
      integer, parameter :: n_eos95_nr = 77, n_eos95_nt = 56
      integer, parameter :: n_eos01_nr = 169, n_eos01_nt = 191
      integer, parameter :: n_eos06_nr = 169, n_eos06_nt = 197
! Loop variable for the repeat-count array constructors used below to
! transcribe the original DATA statements' "N*value" runs -- must be a
! separately-declared integer in scope for the implied-do to resolve
! under `implicit none`, even though it is otherwise unused. Private
! so it can't collide with any `use`ing file's own local `i`.
      integer, private :: i

      type, public :: opal_eos_state
! former common/lreadco/: shared by COMMON block name across all three
! years' readers/interpolators in the original.
           integer :: table_loaded_flag

! ==================== 1995 (unsuffixed) ====================
! former common/a/
           double precision :: eos_table(n_eos_mx,n_eos_mv,n_eos95_nt,n_eos95_nr), &
                t6_list(n_eos95_nr,n_eos95_nt), density_grid(n_eos95_nr), &
                t6_grid(n_eos95_nt), x_interp_result(n_eos95_nt,n_eos95_nr), &
                x_interp_result_alt(n_eos95_nt,n_eos95_nr), &
                x_grid_spacing_inv(n_eos_mx), t6_grid_spacing_inv(n_eos95_nt), &
                density_grid_spacing_inv(n_eos95_nr)
           double precision :: x_grid(n_eos_mx) = &
                [0.0d0,0.2d0,0.4d0,0.6d0,0.8d0]
! x_loop_index is a plain local (not a module field): it is used as a
! DO-loop control variable in esac.f90/readco.f90, which the Fortran
! standard requires to be a named scalar variable, not a derived-type
! component -- see this file's header comment for the general note.
           integer :: x_index_lo
! former common/aa/
           double precision :: rho_interp_hi(4), rho_interp_lo(4), xxh
! former common/b/
           double precision :: z_table(n_eos_mx)
           integer :: eos_index_inverse(10)
           integer :: eos_var_order(10) = [1,2,3,4,5,6,7,8,9,10]
           integer :: t6_index_lo(n_eos95_nr) = &
                [(56,i=1,37),(54,i=1,7), &
                 39,37,36,34,33,31,31,30,29,28,27,26,26,25,24,23,22,21,21,20, &
                 19,18,17,16,15,15,13,13,11,11,9,9,7]
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
           integer :: density_index_edge_at_t(n_eos95_nt) = &
                [(77,i=1,7),(76,i=1,2),(74,i=1,2),(72,i=1,2),(70,i=1,2), &
                 68,67,66,65,64,63,61,60,59,58,57,55,54,53,52,51, &
                 (49,i=1,2),48,(47,i=1,2),46,(45,i=1,2),(44,i=1,15),(37,i=1,2)]
           integer :: t_row_index

! ==================== 2001 (_01 suffix) ====================
! former common/aeos/
           double precision :: eos_table_01(n_eos_mx,n_eos_mv,n_eos01_nt,n_eos01_nr), &
                t6_list_01(n_eos01_nr,n_eos01_nt), density_grid_01(n_eos01_nr), &
                t6_grid_01(n_eos01_nt), x_interp_result_01(n_eos01_nt,n_eos01_nr), &
                x_interp_result_alt_01(n_eos01_nt,n_eos01_nr), &
                x_grid_spacing_inv_01(n_eos_mx), t6_grid_spacing_inv_01(n_eos01_nt), &
                density_grid_spacing_inv_01(n_eos01_nr)
! x_loop_index_01 is a plain local, not a module field -- see the
! x_loop_index note in the 1995 section above.
           integer :: x_index_lo_01
           double precision :: x_grid_01(n_eos_mx) = [0.0, 0.2, 0.4, 0.6, 0.8]
! former common/aaeos/
           double precision :: rho_interp_hi_01(4), rho_interp_lo_01(4), xxh_01
! former common/beos/
           double precision :: z_table_01(n_eos_mx)
           integer :: eos_index_inverse_01(10)
           integer :: eos_var_order_01(10) = [1,2,3,4,5,6,7,8,9,10]
           integer :: t6_index_lo_01(n_eos01_nr) = &
                [(191,i=1,92),190,189,188,187,186,185,184,174,(134,i=1,4), &
                 (133,i=1,3),(132,i=1,2),98,92,(85,i=1,2),(77,i=1,2),71, &
                 (63,i=1,3),(59,i=1,2),53,51,(46,i=1,2),(44,i=1,9),(38,i=1,3), &
                 (33,i=1,6),(29,i=1,16),27,26,25,23,22,20,19,18,17,16]
! former common/bbeos/
           integer :: density_index_1_01, density_index_2_01, density_index_3_01, &
                density_index_4_01, t6_index_1_01, t6_index_2_01, t6_index_3_01, &
                t6_index_4_01, t6_interp_order_01, density_interp_order_01
! former common/eeos/
           double precision :: esact_01, eos_output_01(n_eos_mv)
! former common/eeeos/
           double precision :: x_interp_workspace_01(n_eos_mx,n_eos01_nt,n_eos01_nr), &
                x_grid_copy_01(n_eos_mx)
! former common/eeeeos/
           double precision :: moles_per_gram_table_01(n_eos_mx), &
                hydrogen_fraction_header_01(n_eos_mx), &
                mean_molecular_weight_header_01(n_eos_mx), &
                density_grid_table_01(n_eos_mx,n_eos01_nr), &
                species_fraction_header_01(n_eos_mx,6), &
                log10_r_value_01(n_eos01_nr,n_eos01_nt)
           integer :: temperature_count_used_01(n_eos_mx,n_eos01_nr)
! former common/rmpopeos01/
           double precision :: density_edge_at_t_01(n_eos01_nt)
           integer :: density_index_edge_at_t_01(n_eos01_nt) = &
                [(169,i=1,16),168,167,166,165,(164,i=1,2),163,(162,i=1,2), &
                 161,160,(159,i=1,2),(143,i=1,4),(137,i=1,5),(134,i=1,6), &
                 (125,i=1,2),(123,i=1,5),(122,i=1,2),(121,i=1,6),(119,i=1,4), &
                 (116,i=1,8),(115,i=1,9),(113,i=1,5),(111,i=1,7),(110,i=1,6), &
                 (109,i=1,34),107,104,(100,i=1,40),(99,i=1,10),98,97,96,95,94,93,92]
           integer :: t_row_index_01

! ==================== 2006 (_06 suffix) ====================
! former common/aeos06/
           double precision :: eos_table_06(n_eos_mx,n_eos_mv,n_eos06_nt,n_eos06_nr), &
                t6_list_06(n_eos06_nr,n_eos06_nt), density_grid_06(n_eos06_nr), &
                t6_grid_06(n_eos06_nt), x_interp_result_06(n_eos06_nt,n_eos06_nr), &
                x_interp_result_alt_06(n_eos06_nt,n_eos06_nr), &
                x_grid_spacing_inv_06(n_eos_mx), t6_grid_spacing_inv_06(n_eos06_nt), &
                density_grid_spacing_inv_06(n_eos06_nr)
! x_loop_index_06 is a plain local, not a module field -- see the
! x_loop_index note in the 1995 section above.
           integer :: x_index_lo_06
           double precision :: x_grid_06(n_eos_mx) = [0.0, 0.2, 0.4, 0.6, 0.8]
! former common/aaeos06/
           double precision :: rho_interp_hi_06(4), rho_interp_lo_06(4), xxh_06
! former common/beos06/
           double precision :: z_table_06(n_eos_mx)
           integer :: eos_index_inverse_06(10)
           integer :: eos_var_order_06(10) = [1,2,3,4,5,6,7,8,9,10]
           integer :: t6_index_lo_06(n_eos06_nr) = &
                [(197,i=1,87),(191,i=1,7),190,(189,i=1,2),185,179,170, &
                 (149,i=1,2),133,125,123,122,120,115,113,107,102,(80,i=1,2), &
                 72,68,66,64,62,56,54,52,51,(50,i=1,2),49,47,(45,i=1,2),43,42, &
                 (40,i=1,28),39,37,36,35,34,32,31,30,29,27,26]
           integer :: density_index_edge_06(n_eos06_nt) = &
                [(169,i=1,26),168,(167,i=1,2),166,165,164,(163,i=1,2),162, &
                 161,160,(159,i=1,2),158,(130,i=1,2),129,(128,i=1,2), &
                 (126,i=1,2),(125,i=1,2),124,122,121,(120,i=1,2),(119,i=1,2), &
                 (118,i=1,6),(117,i=1,2),(116,i=1,2),(115,i=1,2),(114,i=1,4), &
                 (113,i=1,8),(111,i=1,22),(110,i=1,5),(109,i=1,6),(108,i=1,2), &
                 (107,i=1,5),(106,i=1,2),105,(104,i=1,2),(103,i=1,8), &
                 (102,i=1,16),(100,i=1,21),(99,i=1,9),(98,i=1,6),(97,i=1,4), &
                 95,94,(87,i=1,6)]
! former common/bbeos06/
           integer :: density_index_1_06, density_index_2_06, density_index_3_06, &
                density_index_4_06, t6_index_1_06, t6_index_2_06, t6_index_3_06, &
                t6_index_4_06, t6_interp_order_06, density_interp_order_06
! former common/eeeos06/
           double precision :: x_interp_workspace_06(n_eos_mx,n_eos06_nt,n_eos06_nr), &
                x_grid_copy_06(n_eos_mx)
! former common/eeeeos06/
           double precision :: moles_per_gram_table_06(n_eos_mx), &
                hydrogen_fraction_header_06(n_eos_mx), &
                mean_molecular_weight_header_06(n_eos_mx), &
                amu_grid_06(n_eos06_nr,n_eos06_nt), &
                log10_ne_grid_06(n_eos06_nr,n_eos06_nt), &
                density_grid_table_06(n_eos_mx,n_eos06_nr), &
                species_fraction_header_06(n_eos_mx,6), &
                log10_r_value_06(n_eos06_nr,n_eos06_nt)
           integer :: temperature_count_used_06(n_eos_mx,n_eos06_nr)
      end type opal_eos_state

      type(opal_eos_state), public :: opal_eos

end module opal_eos_lib

!----------------------------------------------------------------------
! atm_table_lib
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Bundles atmosphere/opacity/EOS lookup-table
! state -- former common/ccr/, fac/, atmprt/, eeos06/, alatm01/,
! alatm02/, alatm05/, jtest/ -- into one derived type, following the
! pulse_diag_lib/run_diag_lib/rotdiff_lib precedent for grouping many
! small, unrelated, per-model or table-load-time blocks into a single
! module. None of these touch core/read_input.f90. Every declaring file
! for a given block used byte-identical member names/order.
module atm_table_lib
      implicit none
! Kurucz (nt) and Kurucz-Castelli (ntc) atmosphere tables: rows in
! log Teff, atm_table_ng columns in log g.
      integer, parameter :: atm_table_nt = 57, atm_table_ntc = 76
      integer, parameter :: atm_table_ng = 11
! Allard table maximum extents (rows in log Teff, columns in log g).
      integer, parameter :: atm_table_nta = 250, atm_table_nga = 25

! surface_p_table: one tabulated surface-pressure table -- log P at
! (log Teff row, log g column) -- with the log Teff rows, the
! atm_table_ng log g columns and the per-row lower/upper log-g edge
! indices (former common/fac/ gmin/gmax). 2026 wave 3 (R5): the former
! kurucz_teff_table/kurucz_logg_table/kurucz_log10_pressure_table +
! atm_table%kurucz_gmin_index/kurucz_gmax_index (atm_table_nt rows)
! and their kurucz_castelli_*/castelli_* twins (atm_table_ntc rows)
! became atm_table%kurucz_surface_p / castelli_surface_p so that one
! interpolator (tables/surfp.f90) serves both. The two tables have
! different row counts, so the components are allocatable: atm_init
! allocates each instance once, with exactly the extents the fixed
! arrays had (num_teff = atm_table_nt or atm_table_ntc, atm_table_ng
! columns), and zero-fills it as the static arrays were.
      type, public :: surface_p_table
           integer :: num_teff = 0
           double precision, allocatable :: teff(:), logg(:)
           double precision, allocatable :: log10_pressure(:,:)
           integer, allocatable :: gmin_index(:), gmax_index(:)
      end type surface_p_table

      type, public :: atm_table_state
! (former common/ccr/ -- the degenerate-electron Fermi-Dirac EOS
! table, fermi_table_* -- lived here through phases one and two
! purely as an accident of the original COMMON grouping; relocated
! to yale_eos_lib in 2026, ROADMAP.md stage 1, alongside moving its
! load from an inline setup/setups.f90 block into eos_lib's
! eos_init. Consumed only by eos/yale/fully_ionized_eos.f90.)
! The Kurucz (atm_table_nt rows) and Kurucz/Castelli (atm_table_ntc
! rows) surface-pressure tables, each with its former common/fac/
! log-g edge indices. See surface_p_table.
           type(surface_p_table) :: kurucz_surface_p, castelli_surface_p
! former common/atmprt/: current T-tau atmosphere integration point.
           double precision :: atm_tau, atm_log10_pressure, &
                atm_log10_temperature, atm_log10_density, atm_opacity
           double precision :: atm_ion_fraction(3)
! (former common/eeos06/ -- the OPAL 2006 EOS interpolator result,
! esact + eos_output(10) -- lived here through phase one purely as an
! accident of the original COMMON grouping; relocated to
! opal_eos_lib as esact_06/eos_output_06 in 2026, ROADMAP.md stage 1,
! so eos state lives in eos's state module.)
! former common/alatm01/: the Allard NextGen/BT-Settl model-atmosphere
! surface-pressure table.
           double precision :: allard_teffl_grid(atm_table_nta), &
                allard_gl_grid(atm_table_nga), allard_feh_grid(atm_table_nga), &
                allard_alpha_grid(atm_table_nga)
           double precision :: allard_log10_pressure(atm_table_nta,atm_table_nga), &
                allard_log10_pressure_tau100(atm_table_nta,atm_table_nga), &
                allard_log10_temp_tau100(atm_table_nta,atm_table_nga)
           logical :: allard_is_old_nextgen
           integer :: allard_num_teff, allard_num_gl, allard_num_feh, &
                allard_num_alpha
! former common/alatm02/: per-row Allard table interpolation bounds.
           double precision :: allard_gl_row_min(atm_table_nta), &
                allard_gl_row_max(atm_table_nta)
           integer :: allard_gl_index_min(atm_table_nta), &
                allard_gl_index_max(atm_table_nta)
           double precision :: allard_teffl_min, allard_teffl_max, &
                allard_gl_min, allard_gl_max
! former common/alatm05/: Allard alpha-enhanced table Teff range.
           double precision :: allard_al_teffl_min, allard_al_teffl_max
      end type atm_table_state

      type(atm_table_state), public :: atm_table


! 2026 (phase six, step 3 -- ROADMAP.md): evicted here from
! const_lib, where this table/working data had landed during the
! phase-one COMMON conversion; it belongs with this domain's state.
! former common/atmos2/: the table abundance kurucz_table_z
! (originally atmz; written by atm_init for the Kurucz and the
! Kurucz/Castelli table alike) and atm_table_file_unit (originally
! ioatm) are spelled identically to their canonical names --
! use-associated directly. The tables themselves (former
! common/atmos2/ atmpl/atmtl/atmgl and common/atmos2c/ atmplc/atmtlc/
! atmglc) are atm_table%kurucz_surface_p / castelli_surface_p above.
      double precision :: kurucz_table_z
      integer :: atm_table_file_unit


end module atm_table_lib

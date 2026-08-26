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
      integer, parameter :: atm_table_nt = 57, atm_table_ntc = 76
      integer, parameter :: atm_table_nta = 250, atm_table_nga = 25

      type, public :: atm_table_state
! (former common/ccr/ -- the degenerate-electron Fermi-Dirac EOS
! table, fermi_table_* -- lived here through phases one and two
! purely as an accident of the original COMMON grouping; relocated
! to yale_eos_lib in 2026, ROADMAP.md stage 1, alongside moving its
! load from an inline setup/setups.f90 block into eos_lib's
! eos_init. Consumed only by eos/yale/fully_ionized_eos.f90.)
! former common/fac/: lower-edge-of-table-in-log-g indices for the
! Kurucz (atm_table_nt) and Kurucz/Castelli (atm_table_ntc) surface-
! pressure tables.
           integer :: kurucz_gmin_index(atm_table_nt), &
                kurucz_gmax_index(atm_table_nt)
           integer :: teff_interp_start_index
           integer :: gravity_interp_indices(4)
           integer :: castelli_gmin_index(atm_table_ntc), &
                castelli_gmax_index(atm_table_ntc)
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
! former common/jtest/: dead-everywhere placeholder in both its
! declaring files.
           integer :: imax1_placeholder, imax2_placeholder
           logical :: ljvs_placeholder
      end type atm_table_state

      type(atm_table_state), public :: atm_table


! 2026 (phase six, step 3 -- ROADMAP.md): evicted here from
! const_lib, where this table/working data had landed during the
! phase-one COMMON conversion; it belongs with this domain's state.
! former common/atmos2/: the Kurucz surface-pressure table
! (kurucz_log10_pressure_table/kurucz_teff_table/kurucz_logg_table/
! kurucz_table_z, originally atmpl/atmtl/atmgl/atmz) and
! atm_table_file_unit (originally ioatm) are spelled identically to
! their canonical names -- use-associated directly.
      double precision :: kurucz_log10_pressure_table(57,11), &
           kurucz_teff_table(57), kurucz_logg_table(11), kurucz_table_z
      integer :: atm_table_file_unit

! former common/atmos2c/: the Kurucz/Castelli surface-pressure table
! (kurucz_castelli_log10_pressure_table/kurucz_castelli_teff_table/
! kurucz_castelli_logg_table, originally atmplc/atmtlc/atmglc) is
! spelled identically to its canonical name everywhere -- use-associated
! directly. Unused in core/read_input.f90.
      double precision :: kurucz_castelli_log10_pressure_table(76,11), &
           kurucz_castelli_teff_table(76), kurucz_castelli_logg_table(11)


end module atm_table_lib

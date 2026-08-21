!----------------------------------------------------------------------
! atm_table_lib
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Bundles atmosphere/opacity/EOS lookup-table
! state -- former common/ccr/, fac/, atmprt/, eeos06/, alatm01/,
! alatm02/, alatm05/, jtest/ -- into one derived type, following the
! pulse_diag_lib/run_diag_lib/rotdiff_lib precedent for grouping many
! small, unrelated, per-model or table-load-time blocks into a single
! module. None of these touch core/parmin.f90. Every declaring file
! for a given block used byte-identical member names/order.
module atm_table_lib
      implicit none
      integer, parameter :: atm_table_nt = 57, atm_table_ntc = 76
      integer, parameter :: atm_table_nta = 250, atm_table_nga = 25

      type, public :: atm_table_state
! former common/ccr/: the degenerate-electron (Fermi-Dirac) EOS table
! used by eos/yale/eqrelv.f90, loaded once at startup by setup/setups.f90.
           double precision :: fermi_table_x_grid(43), fermi_table_eta(43), &
                fermi_table_data(5,43,20)
           integer :: fermi_table_x_lookup(261)
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
! former common/eeos06/: OPAL 2006 EOS interpolator result.
           double precision :: esact
           double precision :: eos_output(10)
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

end module atm_table_lib

!----------------------------------------------------------------------
! controls_buffer
!----------------------------------------------------------------------
! The controls under their code names, each with its default as its
! initializer, the one place a control's default is written: the buffer
! the decks are read into. read_controls seeds it from the pristine
! star%control and star%physics (seed_controls_buffer) and reads the deck
! into it: an inlist's groups bind these names (or a local copied into
! one, for an inlist name that differs), a legacy pair's names are copied
! to them (tools/gen_inlist.py, from defaults/controls_registry.tsv, which
! gives each control its names). map_user_inputs and derive_controls then
! compute the controls derived from others, and store_controls_to_star
! copies the buffer into star%control and star%physics, where the rest of
! the code reads its controls; arm_run_copies sets the run's working
! copies of the controls the run changes. tools/gen_controls_state.py
! generates the two records' types and the copies (io/controls_sync.f90)
! from these declarations: regenerate after changing a member, and give
! it its registry row. Each member is declared on its own line, in the
! section of its group (&control, &physics) and the topic it belongs to,
! with its default as its initializer and its description as the comment
! directly above it (gen_inlist copies the description into the
! documents). Only the controls read's files
! (tools/check_boundaries.py's CONTROLS_BUFFER_IMPORTERS) use this module,
! each naming what it takes (use ..., only:); the sizes and slots it uses
! are state/control_slots.f90's and state/reaction_layout.f90's.

module controls_buffer
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use control_slots, only: max_runs
      use reaction_layout, only: num_reactions, num_extra_scales, num_sstandard
      implicit none
! The buffer is the read's shared state: every member, and the named
! constants, are public by name (the generated controls_buffer_public.inc);
! the read's files import what they use, each with an only: list
! (tools/check_boundaries.py's CONTROLS_BUFFER_IMPORTERS).
      private
      include 'controls_buffer_public.inc'

! The light-element reference cross sections (keV b), against which
! derive_controls turns a cross section given for a light-element rate
! (li6_rate_value etc.) into its scale, and which are the cross sections'
! defaults: Li6 5500 (Fowler et al. 1967), Li7 52 (Rolfs & Kavanagh 1986),
! Be9(p,g) 1.1 (approximate, from Fowler 1988), Be9(p,d) and Be9(p,a)
! 15000 (Fowler et al. 1967).
      real(dp), parameter :: li6_reference_xs = 5.5d3, li7_reference_xs = 5.2d1, &
           be9_pg_reference_xs = 1.1d0, be9_pd_reference_xs = 1.5d4, &
           be9_palpha_reference_xs = 1.5d4

! Each member below is a control, under its code name, with its default
! as its initializer and its description in the comment directly above
! it (gen_inlist copies the description into control.defaults,
! physics.defaults and controls_registry.md, in defaults/). The registry
! (defaults/controls_registry.tsv) gives its inlist and legacy names, its
! group, home, run copy and use.

! ==== &control =======================================================

! ---- The run cards (max_runs entries each, state/control_slots.f90; card n's
!      in element n): what each card does, from which model, its rescale
!      targets, its composition, mixing length and fitting point

! The number of runs (kind cards): each rescales a model, evolves it, or both
      integer :: num_runs = 1
! Per-run kind card: 1 = evolve, 2 = rescale (zero timestep, MS/CHeB), 3 =
! rescale and evolve (pre-MS)
      integer :: run_kind(max_runs) = 1
! Per run: the number of models to compute (the run stops earlier if a stop
! condition is met)
      integer :: num_models(max_runs) = 0
! Per card: start from the starting model (first_model_file) rather than the
! previous card's model; .true. for every card unless the deck says otherwise,
! so a continuation card sets it .false. (card 1 always starts from the
! starting model)
      logical :: first_call_flag(max_runs) = .true.
! Model number given to the first converged model after a fresh start-model
! load (first_call_flag card). <= 0 keeps the counter stored in the model file
! (restart continuation, the pre-2026 behavior).
      integer :: set_initial_model_number = 1
! Per-run description strings echoed in output headers
      character(len=256) :: run_description(2) = repeat(achar(0), 256)
! Per-run target mass (Msun) for rescale runs
      real(dp) :: rescale_mass(max_runs) = 0.0d0
! Per-run target hydrogen fraction for rescale runs
      real(dp) :: rescale_x(max_runs) = 0.0d0
! Per-run target metallicity for rescale runs
      real(dp) :: rescale_z(max_runs) = 0.0d0
! Per-run rescale: core mass bound
      real(dp) :: rescale_core_mass(max_runs) = 0.0d0
! Per card: evolve, rescale or both (run_kind for the cards the deck runs, 0
! beyond); derived by the run plan
      integer :: rescale_kind(max_runs)
! Per card: a rescaling card's targets M, X, Z and core mass (rescale_mass,
! rescale_x, rescale_z, rescale_core_mass); derived by the run plan
      real(dp) :: rescale_params(4,max_runs)
! Starting hydrogen mass fraction X, per run card
      real(dp) :: initial_x_array(max_runs) = 0
! Starting metal mass fraction Z, per run card
      real(dp) :: initial_z_array(max_runs) = 0
! Per card: the mixing-length parameter alpha
      real(dp) :: mixing_length_array(max_runs) = 0
! Per run: log10 of the fractional envelope mass of the outer fitting point
! (with has_senv0_array)
      real(dp) :: senv0_array(max_runs) = 1.26d-4
! Per run: move the outer fitting point to the fractional envelope mass
! senv0_array
      logical :: has_senv0_array(max_runs) = .false.
! On a rescaling run with rescale_core_z, rescale_z_mix1 and rescale_z_mix2
! positive: Z = rescale_core_z inside the mass fraction rescale_z_mix1, then
! linear to the surface Z at rescale_z_mix2
      logical :: use_z_ramp = .false.
! Per rescaling card, with use_z_ramp: the Z the ramp sets inside the mass
! fraction rescale_z_mix1 (positive to ramp; the default -1 is off)
      real(dp) :: rescale_core_z(max_runs) = -1.0d0
! Per rescaling card, with use_z_ramp: the mass fraction m/M inside which Z
! is rescale_core_z (positive to ramp)
      real(dp) :: rescale_z_mix1(max_runs) = -1.0d0
! Per rescaling card, with use_z_ramp: the mass fraction m/M at which the
! ramp reaches the surface Z (positive to ramp)
      real(dp) :: rescale_z_mix2(max_runs) = -1.0d0

! ---- Each card's stops and its fixed timestep (a structure limit at its
!      +-1d99 sentinel is off)

! Stop the run card at this age (yr); 0 = no age limit
      real(dp) :: target_end_age(max_runs) = 0.0d0
! Per run: stop when the central deuterium mass fraction falls below this (0:
! off)
      real(dp) :: central_deuterium_stop(max_runs) = 0.0d0
! Per run: stop when the central hydrogen mass fraction falls below this (0:
! off)
      real(dp) :: central_hydrogen_stop(max_runs) = 0.0d0
! Per run: stop when the central helium mass fraction falls below this (0:
! off)
      real(dp) :: central_helium_stop(max_runs) = 0.0d0
! Per card: an evolving card has a stop armed (max_age or a central-abundance
! stop positive); derived by the run plan (interpret_kind_cards)
      logical :: end_age_stop_active(max_runs)
! Per run: if positive, a fixed timestep (yr) for every model of the run,
! overriding the timestep limits
      real(dp) :: timestep_override(max_runs) = 0.0d0
! Per card: an evolving card takes the fixed timestep timestep_override
! (positive); derived by the run plan
      logical :: timestep_override_active(max_runs)
! Per kind card: stop the card when log(L/Lsun) exceeds this (1d99 = off)
      real(dp) :: log_L_upper_limit(max_runs) = 1.0d99
! Per kind card: stop the card when log(L/Lsun) drops below this (-1d99 = off)
      real(dp) :: log_L_lower_limit(max_runs) = -1.0d99
! Per kind card: stop the card when Teff [K] exceeds this (1d99 = off)
      real(dp) :: Teff_upper_limit(max_runs) = 1.0d99
! Per kind card: stop the card when Teff [K] drops below this (-1d99 = off)
      real(dp) :: Teff_lower_limit(max_runs) = -1.0d99
! Per kind card: stop the card when surface log g [cgs] exceeds this (1d99 =
! off)
      real(dp) :: log_g_upper_limit(max_runs) = 1.0d99
! Per kind card: stop the card when surface log g [cgs] drops below this
! (-1d99 = off)
      real(dp) :: log_g_lower_limit(max_runs) = -1.0d99
! Per kind card: stop the card when nu_max [uHz] (scaling relation) exceeds
! this (1d99 = off)
      real(dp) :: nu_max_upper_limit(max_runs) = 1.0d99
! Per kind card: stop the card when nu_max [uHz] (scaling relation) drops
! below this (-1d99 = off)
      real(dp) :: nu_max_lower_limit(max_runs) = -1.0d99

! ---- The calibrations: the solar model, and a star to a target luminosity,
!      Teff or radius

! Calibrate Y and the mixing length so that the model matches the solar L and
! R at target_solar_age
      logical :: calibrate_solar_model = .false.
! Solar calibration: also calibrate Z to match target_solar_zx (within
! zx_tolerance)
      logical :: calibrate_solar_zx = .false.
! Solar calibration: the target age of the Sun (yr)
      real(dp) :: target_solar_age = 4.57d9
! Solar calibration: the target surface Z/X (with calibrate_solar_zx)
      real(dp) :: target_solar_zx = 0.02292d0
! Solar calibration: tolerance on |log L/Lsun|
      real(dp) :: luminosity_tolerance = 1.0d-5
! Solar calibration: tolerance on |log R/Rsun|
      real(dp) :: radius_tolerance = 1.0d-4
! Solar calibration: the tolerance on log10(Z/X)
      real(dp) :: zx_tolerance = 1.0d-3
! Calibrate the mixing length and Y so that the star reaches
! target_luminosity_lsun at the target radius or Teff
      logical :: calibrate_star_flag = .false.
! Star calibration: match the target Teff (target_teff) instead of the target
! radius (target_radius_rsun)
      logical :: specify_teff_flag = .false.
! Star calibration: the target luminosity (Lsun)
      real(dp) :: target_luminosity_lsun = 0
! Star calibration: the tolerance on L/Lsun at the target radius or Teff
      real(dp) :: target_star_luminosity_tolerance = 0
! Star calibration: the target Teff (K), with specify_teff_flag
      real(dp) :: target_teff = 0
! Star calibration: the target radius (Rsun)
      real(dp) :: target_radius_rsun = 0

! ---- The starting model's CNO mixture, isotope ratios and light elements

! 1: set the starting model's C, N and O fractions of Z from cno_mixture_name
      integer :: mixture_change_mode = 0
! Heavy-element mixture selector (GS98/CUS/...) applied to the start model
      character(len=8) :: cno_mixture_name = 'GS98'
! Carbon fraction of the CNO metals for cno_mixture_name 'CUS' (with
! mixture_change_mode = 1)
      real(dp) :: target_carbon_cno_fraction = 0.172148d0
! Nitrogen fraction of the CNO metals for cno_mixture_name 'CUS'
      real(dp) :: target_nitrogen_cno_fraction = 0.050426d0
! Oxygen fraction of the CNO metals for cno_mixture_name 'CUS'
      real(dp) :: target_oxygen_cno_fraction = 0.468195d0
! Z/X for cno_mixture_name 'CUS' (inert: stored and used nowhere; the starting
! model keeps its Z)
      real(dp) :: target_metal_fraction = 0.02292d0
! The starting model's CNO mixture is changed (mixture_change_mode = 1 and a
! valid mixture)
      logical :: change_cno_mixture_active
! 1: set the starting model's CNO isotope ratios and light-element abundances
! from isotope_mixture_name
      integer :: isotope_change_mode = 0
! Isotope-ratio set selector (L21 default or CUS custom)
      character(len=8) :: isotope_mixture_name = 'L21'
! C12/C13 number ratio of the starting model for isotope_mixture_name 'CUS'
! (with isotope_change_mode = 1)
      real(dp) :: c12_to_c13_ratio = 88.27d0
! N14/N15 number ratio for isotope_mixture_name 'CUS'
      real(dp) :: n14_to_n15_ratio = 411.9d0
! O16/O17 number ratio for isotope_mixture_name 'CUS'
      real(dp) :: o16_to_o17_ratio = 471.4d0
! O16/O18 number ratio for isotope_mixture_name 'CUS'
      real(dp) :: o16_to_o18_ratio = 2693.0d0
! Initial H2 mass fraction for isotope_mixture_name 'CUS'
      real(dp) :: initial_h2_fraction = 2.781d-5
! Initial He3 mass fraction for isotope_mixture_name 'CUS'
      real(dp) :: initial_he3_fraction = 3.461d-5
! Initial Li6 mass fraction for isotope_mixture_name 'CUS'
      real(dp) :: initial_li6_fraction = 7.187d-10
! Initial Li7 mass fraction for isotope_mixture_name 'CUS'
      real(dp) :: initial_li7_fraction = 1.025d-8
! Initial Be9 mass fraction for isotope_mixture_name 'CUS'
      real(dp) :: initial_be9_fraction = 1.595d-10
! Initial B10 mass fraction for isotope_mixture_name 'CUS'
      real(dp) :: initial_b10_fraction = 1.002d-9
! Initial B11 mass fraction for isotope_mixture_name 'CUS'
      real(dp) :: initial_b11_fraction = 4.405d-9
! The starting model's isotope ratios are changed (isotope_change_mode = 1
! and a valid mixture)
      logical :: change_isotope_ratios_active

! ---- The deck and the run's own files (a path the deck leaves unset holds
!      NUL bytes)

! The deck's file: the single-file inlist, or the legacy pair's CONTROL
! namelist (.nml1); read_input finds it (the first command-line argument,
! default yrec8.nml1, or an embedding host's override)
      character(len=256) :: control_nml_file = repeat(achar(0), 256)
! The legacy pair's PHYSICS namelist (.nml2; the second command-line argument,
! default yrec8.nml2), or the inlist again
      character(len=256) :: physics_nml_file = repeat(achar(0), 256)
! Verbatim copy of the input namelist file(s), written at startup as run
! provenance; blank (the default): inlist_used in the output directory, the
! run log's
      character(len=256) :: inlist_used_file = ' '
! The starting model, which card 1 and every card with first_call_flag set
! read
      character(len=256) :: first_model_file = repeat(achar(0), 256)
! Path for the final model (.mod) written every model; reloadable
      character(len=256) :: last_model_file = '{YREC_OUTPUT}/final.mod'
! Monte Carlo runs: the file of each run's samples (the rate, composition,
! luminosity and age factors), read at the start
      character(len=256) :: dynamics_file = repeat(achar(0), 256)
! Monte Carlo runs: the file for each run's summary
      character(len=256) :: monte_carlo_file1_path = repeat(achar(0), 256)
! Monte Carlo runs: the file for each run's final model
      character(len=256) :: monte_carlo_file2_path = repeat(achar(0), 256)
! Write an isochrone record every model to isochrone_file
      logical :: isochrone_output_active = .false.
! Isochrone output file
      character(len=256) :: isochrone_file = repeat(achar(0), 256)
! Rewind the run log every model, so it holds only the latest model's solver
! diagnostics
      logical :: rewind_short_file = .false.
! Accepted and ignored: no surface-composition file is written any more (the
! light elements are history columns); kept so that decks setting it still
! read
      character(len=256) :: surface_composition_file = repeat(achar(0), 256)
! Unused: read by the legacy reader and ignored
      character(len=256) :: stored_change_file = repeat(achar(0), 256)

! ---- The output's columns and its interior GYRE files

! File listing the history columns (one name per line, ! comments); blank: the
! default selection of defaults/history_columns.list
      character(len=256) :: history_columns_file = ' '
! File listing the profile columns; blank: the default selection of
! defaults/profile_columns.list
      character(len=256) :: profile_columns_file = ' '
! Write the interior-only GYRE file every this many models (0 = off)
      integer :: pulse_gyre_interval = 0
! Inert: read and used nowhere (legacy NPENV, an envelope print count, whose
! one use in the F77 was commented out); the read makes 0 or less into 9999
      integer :: envelope_print_points = 0

! ---- The opacity tables: which families, and their metallicities

! Use the LAOL 1989 opacity tables (laol_opacity_table_file)
      logical :: use_laol89_tables = .false.
! Use the OPAL 1992 opacity tables (opal92_opacity_table_file)
      logical :: use_opal92_tables = .false.
! Use the OPAL 1995 / OP atomic opacity tables (opal95_opacity_table_file)
      logical :: use_opal95_tables = .false.
! Use the Kurucz 1990 low-temperature opacity tables
      logical :: use_kurucz90_tables = .false.
! Use the Alexander & Ferguson 1994 low-temperature opacity tables
! (alex95_table_z_list)
      logical :: use_alex95_tables = .false.
! Use the Ferguson et al. 2005 molecular opacity tables
! (alex06_opacity_table_file)
      logical :: use_alex06_tables = .false.
! Use the pure C and N opacity tables (pure_z_opacity_table_file) for
! metal-enriched mixtures
      logical :: use_pure_z_table = .false.
! Metallicity of the first LAOL89 opacity table
      real(dp) :: laol_table_z1 = 0
! Metallicity of the second LAOL89 opacity table (two-Z interpolation)
      real(dp) :: laol_table_z2 = 0
! Metallicity of the first OPAL92 opacity table
      real(dp) :: opal_table_z1 = 0
! Z of the second OPAL 1992 opacity table
      real(dp) :: opal_table_z2 = 0
! Z of the surface-mixture OPAL 1995 / OP table generated at the start
      real(dp) :: opal95_single_table_z = 0
! Inert: read and used nowhere (the OPAL 1995 tables take one Z,
! opal95_single_table_z); the controls read warns when it is set
      real(dp) :: opal95_z2 = 0
! Z of the Alexander & Ferguson 1994 low-temperature table
      real(dp) :: alex_table_z1 = 0
! Inert: read and used nowhere (the Alexander tables take one Z,
! alex_table_z1); the controls read warns when it is set
      real(dp) :: alex_z2 = 0
! Z of the first Kurucz 1990 opacity table
      real(dp) :: kurucz_table_z1 = 0
! Z of the second Kurucz 1990 opacity table
      real(dp) :: kurucz_table_z2 = 0

! ---- The table files, by the domain that reads them (the builders in
!      setup/<domain>_settings_builder.f90 copy them into each domain's paths
!      record); a path the deck leaves unset holds NUL bytes

! LAOL 1989 opacity table
      character(len=256) :: laol_table_path = repeat(achar(0), 256)
! LAOL 1989 opacity table, second Z
      character(len=256) :: laol_table2_path = repeat(achar(0), 256)
! OPAL 1995-format opacity table (use_opal95_tables)
      character(len=256) :: opal95_table_path = repeat(achar(0), 256)
! OPAL 1992 opacity table (use_opal92_tables)
      character(len=256) :: opal92_table_path = repeat(achar(0), 256)
! OPAL 1992 opacity table, second Z
      character(len=256) :: opal92_table2_path = repeat(achar(0), 256)
! Pure-Z opacity table (metal diffusion support)
      character(len=256) :: pure_z_table_path = repeat(achar(0), 256)
! Ferguson/Alexander 2006 low-T opacity tables
      character(len=256) :: alex06_table_path = repeat(achar(0), 256)
! The seven Alexander & Ferguson 1994 low-temperature opacity table files
      character(len=256) :: alex95_table_paths(7) = repeat(achar(0), 256)
! Kurucz 1990 low-temperature opacity table (use_kurucz90_tables)
      character(len=256) :: kurucz_table_path = repeat(achar(0), 256)
! Kurucz 1990 low-temperature opacity table, second Z
      character(len=256) :: kurucz_table2_path = repeat(achar(0), 256)
! Potekhin conductive opacity table
      character(len=256) :: conductive_table_path = repeat(achar(0), 256)
! Fermi-Dirac integral table (Yale/Prather degenerate-electron EOS)
      character(len=256) :: fermi_table_path = repeat(achar(0), 256)
! Saumon-Chabrier-Van Horn EOS: hydrogen table
      character(len=256) :: scv_h_table_path = repeat(achar(0), 256)
! SCV EOS: helium table
      character(len=256) :: scv_he_table_path = repeat(achar(0), 256)
! SCV EOS: metals table
      character(len=256) :: scv_z_table_path = repeat(achar(0), 256)
! MHD EOS table 1 of 8: the first ZAMS-type composition
      character(len=256) :: zams_a_table_path = repeat(achar(0), 256)
! MHD EOS table 2 of 8: the second ZAMS-type composition
      character(len=256) :: zams_b_table_path = repeat(achar(0), 256)
! MHD EOS table 3 of 8: the third ZAMS-type composition
      character(len=256) :: zams_c_table_path = repeat(achar(0), 256)
! MHD EOS table 4 of 8: the first centre-type composition
      character(len=256) :: centre1_table_path = repeat(achar(0), 256)
! MHD EOS table 5 of 8: the second centre-type composition
      character(len=256) :: centre2_table_path = repeat(achar(0), 256)
! MHD EOS table 6 of 8: the third centre-type composition
      character(len=256) :: centre3_table_path = repeat(achar(0), 256)
! MHD EOS table 7 of 8: the fourth centre-type composition
      character(len=256) :: centre4_table_path = repeat(achar(0), 256)
! MHD EOS table 8 of 8: the fifth centre-type composition
      character(len=256) :: centre5_table_path = repeat(achar(0), 256)
! OPAL 1995 EOS table (use_opal95_eos)
      character(len=256) :: opal95_eos_table_path = repeat(achar(0), 256)
! OPAL 2001 EOS table (use_opal2001_eos)
      character(len=256) :: opal2001_eos_table_path = repeat(achar(0), 256)
! OPAL 2006 EOS table (use_opal2006_eos)
      character(len=256) :: opal2006_eos_table_path = repeat(achar(0), 256)
! Kurucz (atm_choice 3 and 4) or Castelli/Kurucz (5) surface-pressure table
      character(len=256) :: atm_table_path = repeat(achar(0), 256)
! Allard NextGen model-atmosphere table
      character(len=256) :: allard_table_path = repeat(achar(0), 256)

! ==== &physics =======================================================

! ---- The output and the run log

! The history file's name, in the output directory
      character(len=256) :: star_history_name = 'history.data'
! Write profile{N}.data every N converged models (0 = off); the history's
! profile_number column maps models to profiles
      integer :: profile_interval = 50
! Profile/pulse file name prefix: <prefix><N>.data[.GYRE/.FGONG/.GSM] in the
! output directory
      character(len=256) :: profile_data_prefix = 'profile'
! Write the profile{N}.data files (every profile_interval models)
      logical :: write_profile_flag = .false.
! Write a pulsation file (pulse_format) with each profile
      logical :: write_pulse_flag = .false.
! The pulsation files' format (GYRE text, FGONG, or GSM = GYRE HDF5 [needs
! make USE_HDF5=1]), written beside each profile as profile{N}.data.{format}
      character(len=8) :: pulse_format = 'GYRE'
! Print the compact per-model progress line (terminal + run log) every N
! models; 0 = off; the final model always prints
      integer :: terminal_interval = 50
! Run log path; its directory part defines the output directory (legacy
! FSHORT)
      character(len=256) :: log_output_file = '{YREC_OUTPUT}/run.log'
! Write the solver forensics to the run log (Henyey iteration trace, envelope
! and atmosphere integration statistics, triangle and rezoning bookkeeping);
! the read installs it as run_messages' gate, solver_diagnostics(), which the
! run asks
      logical :: report_solver_diagnostics = .false.
! Diagnostic: write the per-zone solar neutrino production table (fluxes per
! shell, normalized profiles, Bahcall-Pinsonneault layout) to the run log at
! the start of each run
      logical :: compute_neutrino_fluxes = .false.

! ---- The solar units (prepare_run copies them to the star, where a
!      Monte-Carlo run rescales the luminosity), Newton's G and the
!      asteroseismic references

! Solar mass (g), which converts Msun to grams everywhere (the rescaling, log
! g, the seismic scalings)
      real(dp) :: solar_mass_cgs = 1.9891d33
! Solar luminosity (erg/s)
      real(dp) :: solar_luminosity_cgs = 3.8515d33
! The solar radius (cm) used for the tracks and the output
      real(dp) :: solar_radius_cgs = 6.9598d10
! Newton's G (cgs). Negative (the default) keeps YREC's log10 G = -7.17571
! (G = 6.6726e-8) bit for bit; a positive value (e.g. 6.6743d-8, CODATA
! 2018) replaces it
      real(dp) :: G_cgs = -1.0d0
! Solar nu_max [uHz] for the scaling relation (Huber et al. 2011)
      real(dp) :: nu_max_sun = 3090.0d0
! Solar delta_nu [uHz] for the mean-density scaling relation (Huber et al.
! 2011)
      real(dp) :: delta_nu_sun = 135.1d0
! Solar effective temperature [K] for the nu_max scaling relation (IAU nominal
! value)
      real(dp) :: Teff_sun = 5772.0d0

! ---- The timestep: its factors, and the readable tolerances map_user_inputs
!      maps onto them

! The timestep limits, slot by slot as the time_* controls (itime_* in
! state/control_slots.f90): time_core_min, time_dx_core_tot,
! time_dx_core_frac, time_dy_core_tot, time_dy_core_frac, time_dx_total,
! time_dx_shell, time_dt, time_dp, time_dr, time_dl, time_dy_shell,
! time_max_dt_frac, time_dy_total; used unless use_new_rate_variables, which
! takes the named controls instead
      real(dp) :: timestep_factors(14) = (/1.0d-3,2.0d-2,5.0d-1,2.0d-2, 3.0d-1, &
           1.5d-3,1.0d-1,2.0d-2,4.0d-2,2.0d-2,2.0d-2,0.25d0, 1.5d0,0.25d0/)
! Timestep limit: the largest fractional change of omega per step
      real(dp) :: max_domega_global = 1.0d1
! Limit the timestep by the structure changes (time_dt, time_dp, time_dr,
! time_dl)
      logical :: use_structure_dt_limits = .true.
! Limit the timestep by the envelope-fit triangle
      logical :: use_envelope_triangle_dt = .false.
! Take the zoning and timestep tolerances from the named tol_*, flag_* and
! time_* controls (and use_structure_time_limits) instead of
! zoning_tolerances and timestep_factors; nothing to do with reaction rates
      logical :: use_new_rate_variables = .false.
! With use_new_rate_variables: limit the timestep by the structure changes
! time_dt..time_dl (switched off once nuclear burning dominates)
      logical :: use_structure_time_limits = .false.
! With use_new_rate_variables: the smallest central X or Y for the
! core-burning timestep criteria
      real(dp) :: time_core_min = 1.0d-3
! With use_new_rate_variables: the largest change of L at any point per step
! (pre-main sequence)
      real(dp) :: time_dl = 2.0d-2
! With use_new_rate_variables: the largest change of P at any point per step
! (pre-main sequence)
      real(dp) :: time_dp = 4.0d-2
! With use_new_rate_variables: the largest change of r at any point per step
! (pre-main sequence)
      real(dp) :: time_dr = 2.0d-2
! With use_new_rate_variables: the largest change of T at any point per step
! (pre-main sequence)
      real(dp) :: time_dt = 2.0d-2
! With use_new_rate_variables: the largest change of omega per step
      real(dp) :: time_dw_global = 8.0d-2
! With use_new_rate_variables: the largest change of omega per
! rotational-mixing sub-step
      real(dp) :: time_dw_mix = 8.0d-2
! With use_new_rate_variables: the largest fraction of the central X burned
! per step
      real(dp) :: time_dx_core_frac = 0.5d0
! With use_new_rate_variables: the largest change of the central X per step
      real(dp) :: time_dx_core_tot = 2.0d-2
! With use_new_rate_variables: the largest fraction of X burned in the shell
! per step
      real(dp) :: time_dx_shell = 0.1d0
! With use_new_rate_variables: the largest mass of hydrogen burned per step
! (Msun)
      real(dp) :: time_dx_total = 1.5d-3
! With use_new_rate_variables: the largest fraction of the central Y burned
! per step
      real(dp) :: time_dy_core_frac = 0.5d0
! With use_new_rate_variables: the largest change of the central Y per step
      real(dp) :: time_dy_core_tot = 2.0d-2
! With use_new_rate_variables: the largest fraction of Y burned in the shell
! per step
      real(dp) :: time_dy_shell = 0.1d0
! With use_new_rate_variables: the largest mass of helium burned per step
! (Msun)
      real(dp) :: time_dy_total = 1.5d-3
! With use_new_rate_variables: the largest factor by which the timestep
! changes from one model to the next
      real(dp) :: time_max_dt_frac = 1.5d0

! ---- The zoning, and the readable tolerances map_user_inputs maps onto it

! The zoning tolerances, slot by slot as the tol_*/flag_* controls (ichi_* in
! state/control_slots.f90): tol_dm_min, tol_dm_max, flag_dx, flag_dz,
! tol_dx_max, tol_dz_max, tol_czbase_fine_width, tol_dp_env_max, tol_dl_max,
! tol_dp_czbase_max, tol_dp_core_max, flag_dw; used unless
! use_new_rate_variables, which takes the named controls instead
      real(dp) :: chi_grid_scale(12) = (/ 1.0d-8, 8.0d-2, 5.0d-2, 5.0d-2, 1.0d0, &
            1.0d0, 0.0d0, 5.0d-2, 2.0d-2, 5.0d-2, 5.0d-2, 1.0d-1 /)
! Zoning: an X jump between shells larger than this flags a shell
      real(dp) :: flag_dx = 0.05d0
! Zoning: the largest log-omega spacing (rotating models)
      real(dp) :: flag_dw = 0.10d0
! Zoning: a Z jump between shells larger than this flags a shell
      real(dp) :: flag_dz = 0.05d0
! With use_new_rate_variables: the width in log P of the finely zoned region
! at the base of the convection zone
      real(dp) :: tol_czbase_fine_width = 0.0d0
! With use_new_rate_variables: the largest spacing in L/L_surface between
! points
      real(dp) :: tol_dl_max = 0.02d0
! With use_new_rate_variables: the largest spacing in log m between points
      real(dp) :: tol_dm_max = 0.08d0
! With use_new_rate_variables: the smallest spacing in log m between points
      real(dp) :: tol_dm_min = 1.0d-8
! With use_new_rate_variables: the largest spacing in log P below the surface
! convection zone
      real(dp) :: tol_dp_core_max = 0.05d0
! With use_new_rate_variables: the largest spacing in log P near the base of
! the surface convection zone
      real(dp) :: tol_dp_czbase_max = 0.05d0
! With use_new_rate_variables: the largest spacing in log P in the surface
! convection zone
      real(dp) :: tol_dp_env_max = 0.05d0
! With use_new_rate_variables: the largest spacing in X between points
! (usually off)
      real(dp) :: tol_dx_max = 1.0d0
! With use_new_rate_variables: the largest spacing in Z between points
! (usually off)
      real(dp) :: tol_dz_max = 1.0d0

! ---- The Newton iteration of the structure

! The Henyey thresholds: rows 1-4 the per-variable convergence (column 1) and
! correction-too-large (column 2) limits, row 5 the level-5 residual test
      real(dp) :: solver_tolerances(5,2) = reshape((/6.0d-5,4.5d-5,3.0d-5, 9.0d-5, &
           3.0d-5,9.0d-1,5.0d-1,5.0d-1,2.0d0,2.5d-6/), (/5,2/))
! Henyey under-relaxation: each step starts with the correction factor
! |solver_underrelax_start| - solver_underrelax_increment; a positive value
! lets it grow each iteration
      real(dp) :: solver_underrelax_start = 0.8d0
! Henyey under-relaxation: the correction factor grows by this each corrector
! pass, up to 1
      real(dp) :: solver_underrelax_increment = 0.1d0
! Henyey iteration limit, level 1 (fresh surface BC)
      integer :: max_iter_level1 = 2
! Henyey iteration limit, level 2
      integer :: max_iter_level2 = 20
! Henyey iteration limit, level 3
      integer :: max_iter_level3 = 2
! Henyey iteration limit, level 4 (high precision)
      integer :: max_iter_level4 = 0
! Use stored structure changes for the next model's first guess
      logical :: improved_first_guess_flag = .false.

! ---- The log T cutoffs

! log T cutoffs (5) below which burning groups are switched off
      real(dp) :: nuclear_logT_cutoffs(5) = (/6.5d0,6.5d0,6.82d0,7.7d0,7.5d0/)
! log10 T above which the default EOS assumes full ionization
      real(dp) :: saha_log10t_cutoff = 6.0d0
! Above this log10 T a convective point takes the adiabatic gradient without
! the mixing-length theory (deep convective cores)
      real(dp) :: gradient_logt_cutoff = 6.9d0
! The rezone keeps the outermost interior point hotter than this log10 T:
! cooler outer points are deleted
      real(dp) :: rezone_env_logt0 = 3.0d0
! The rezone adds a point outside the outermost interior point when that
! point is hotter than this log10 T (at the stored envelope state)
      real(dp) :: rezone_env_logt1 = 9.0d0
! The log T at which an envelope integration stores its state for the surface
! boundary: the mean of rezone_env_logT0 and rezone_env_logT1
      real(dp) :: envelope_store_log_temperature

! ---- Nuclear burning: the network, the neutrinos and the reaction rates'
!      scales

! Take the cross sections from the s0_* controls (Adelberger et al. 2011
! defaults) instead of the ratios to Bahcall & Ulrich (1988)
      logical :: use_new_nuclear_rates = .false.
! Threshold of the screening parameter below which nuclear rates use weak
! screening
      real(dp) :: weak_screening_threshold = 0.03d0
! Neutrino losses from Itoh et al. (1996)
      logical :: use_itoh_neutrino_loss = .false.
! Compute solar neutrino fluxes
      logical :: calc_neutrinos = .false.
! Burning network: abundances below this are treated as zero
      real(dp) :: min_abundance = 1.0d-20
! Burning network: absolute tolerance
      real(dp) :: absolute_tolerance = 1.0d-5
! The burning network solver's relative tolerance
      real(dp) :: relative_tolerance = 1.0d-4
! The burning network solver's largest number of iterations
      integer :: max_burn_iterations = 50
! The cross sections as ratios to Bahcall & Ulrich (1988), one per slot of
! the legacy SSTANDARD array, used without use_new_nuclear_rates
      real(dp) :: cross_section_ratios(num_sstandard) = (/ 0.9828d0, 1.0485d0, &
            0.9815d0, 0.9241d0, 1.3818d0, 1.0542d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, &
            1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0108d0, 0.7819d0, 0.2875d0 /)
! Each reaction's cross-section (S0) scale against Bahcall & Ulrich (1988), by
! reaction slot r_*: from the S-factor inputs (use_new_nuclear_rates) or from
! cross_section_ratios
      real(dp) :: cross_section_scale(num_reactions) = 1.0d0
! The scales of the rates outside the reaction table, by extra slot xs_* (pep,
! Be7 + e, Be7 + p, hep): from the S-factor inputs (use_new_nuclear_rates) or
! from cross_section_ratios
      real(dp) :: extra_cross_section_scale(num_extra_scales) = 1.0d0
! Each reaction's S'/S scale against Solar Fusion I, by reaction slot (1
! without use_new_nuclear_rates)
      real(dp) :: s_prime_scale(num_reactions) = 1.0d0
! Each reaction's S''/S scale against Solar Fusion I, by reaction slot (1
! without use_new_nuclear_rates)
      real(dp) :: s_double_prime_scale(num_reactions) = 1.0d0
! Be7 + p's S'/S scale against Solar Fusion I
      real(dp) :: be7p_s_prime_scale = 1.0d0
! Be7 + p's S''/S scale against Solar Fusion I
      real(dp) :: be7p_s_double_prime_scale = 1.0d0
! The N15(p,gamma)O16 branch's scale (s0_n15_p_o16_branch over net's
! reference)
      real(dp) :: n15_pg_scale = 1.0d0
! The N15(p,alpha)C12 branch's scale (s0_n15_p_c12_branch over net's
! reference)
      real(dp) :: n15_pa_scale = 1.0d0
! Li6 burning rate scale (scale_li6_rate with li6_rate_value overrides it)
      real(dp) :: li6_rate_scale = 1.0d0
! Li7 burning rate scale (scale_li7_rate with li7_rate_value overrides it)
      real(dp) :: li7_rate_scale = 1.0d0
! Scale of the Be9(p,g) rate (scale_be9_rate1 with be9_rate_value1 overrides
! it)
      real(dp) :: be9_pg_rate_scale = 1.0d0
! Scale of the Be9(p,d) rate (scale_be9_rate2 with be9_rate_value2 overrides
! it)
      real(dp) :: be9_pd_rate_scale = 1.0d0
! Scale of the Be9(p,a) rate (scale_be9_rate3 with be9_rate_value3 overrides
! it)
      real(dp) :: be9_palpha_rate_scale = 1.0d0
! Take li6_rate_scale from the cross section li6_rate_value
      logical :: scale_li6_rate = .false.
! Take li7_rate_scale from the cross section li7_rate_value
      logical :: scale_li7_rate = .false.
! Take be9_pg_rate_scale from the cross section be9_rate_value1
      logical :: scale_be9_rate1 = .false.
! Take be9_pd_rate_scale from the cross section be9_rate_value2
      logical :: scale_be9_rate2 = .false.
! Take be9_palpha_rate_scale from the cross section be9_rate_value3
      logical :: scale_be9_rate3 = .false.
! The Li6 burning cross section (keV b): with scale_li6_rate, li6_rate_scale
! is this over li6_reference_xs (5.5e3)
      real(dp) :: li6_rate_value = li6_reference_xs
! The Li7 burning cross section (keV b): with scale_li7_rate, li7_rate_scale
! is this over li7_reference_xs (52)
      real(dp) :: li7_rate_value = li7_reference_xs
! The Be9(p,g) cross section (keV b): with scale_be9_rate1,
! be9_pg_rate_scale is this over be9_pg_reference_xs (1.1)
      real(dp) :: be9_rate_value1 = be9_pg_reference_xs
! The Be9(p,d) cross section (keV b): with scale_be9_rate2,
! be9_pd_rate_scale is this over be9_pd_reference_xs (1.5e4)
      real(dp) :: be9_rate_value2 = be9_pd_reference_xs
! The Be9(p,a) cross section (keV b): with scale_be9_rate3,
! be9_palpha_rate_scale is this over be9_palpha_reference_xs (1.5e4)
      real(dp) :: be9_rate_value3 = be9_palpha_reference_xs

! ---- The nuclear S-factors (keV b) and their first and second derivatives
!      (b, b/keV), used with use_new_nuclear_rates: Solar Fusion II
!      (Adelberger et al. 2011). pep is relative to pp, and Be7 + e- is the
!      constant of the linear term; the code uses T9, so a coefficient for
!      T/1e6 K is divided by sqrt(1000). map_user_inputs scales them against
!      the Bahcall & Ulrich (1988) S0 and the Solar Fusion I derivatives

! S-factor of p+p (keV b), used with use_new_nuclear_rates
      real(dp) :: s0_pp = 4.01d-22
! S-factor of He3+He3 (keV b), used with use_new_nuclear_rates
      real(dp) :: s0_he3he3 = 5.21d3
! S-factor of He3+He4 (keV b), used with use_new_nuclear_rates
      real(dp) :: s0_he3he4 = 5.6d-1
! S-factor of p+C12 (keV b), used with use_new_nuclear_rates
      real(dp) :: s0_p_c12 = 1.34d0
! S-factor of p+C13 (keV b), used with use_new_nuclear_rates
      real(dp) :: s0_p_c13 = 7.6d0
! S0 of N14(p,g)O15 (keV b), with use_new_nuclear_rates
      real(dp) :: s0_p_n14 = 1.66d0
! S0 of O16(p,g)F17 (keV b), with use_new_nuclear_rates
      real(dp) :: s0_p_o16 = 1.06d1
! The constant of the linear term of Be7 electron capture, with
! use_new_nuclear_rates
      real(dp) :: s0_be7_electron = 1.7709d-10
! S0 of Be7(p,g)B8 (keV b), with use_new_nuclear_rates
      real(dp) :: s0_be7_p = 0.0208d0
! S0 of N15(p,a)C12 (keV b), which sets the N15 + p branching
      real(dp) :: s0_n15_p_c12_branch = 7.3d4
! S0 of N15(p,g)O16 (keV b), which sets the N15 + p branching
      real(dp) :: s0_n15_p_o16_branch = 3.6d1
! dS/dE of H1(p,e+ nu)H2 (b), with use_new_nuclear_rates
      real(dp) :: s0p_pp = 4.49d-24
! dS/dE of He3(He3,2p)He4 (b), with use_new_nuclear_rates
      real(dp) :: s0p_he3he3 = -4.9d0
! dS/dE of He3(a,g)Be7 (b), with use_new_nuclear_rates
      real(dp) :: s0p_he3he4 = -3.6d-4
! dS/dE of C12(p,g)N13 (b), with use_new_nuclear_rates
      real(dp) :: s0p_p_c12 = 2.6d-3
! dS/dE of C13(p,g)N14 (b), with use_new_nuclear_rates
      real(dp) :: s0p_p_c13 = -7.83d-3
! dS/dE of N14(p,g)O15 (b), with use_new_nuclear_rates
      real(dp) :: s0p_p_n14 = -3.3d-3
! dS/dE of O16(p,g)F17 (b), with use_new_nuclear_rates
      real(dp) :: s0p_p_o16 = -5.4d-2
! d2S/dE2 of C12(p,g)N13 (b/keV), with use_new_nuclear_rates
      real(dp) :: s0pp_p_c12 = 8.3d-5
! d2S/dE2 of C13(p,g)N14 (b/keV), with use_new_nuclear_rates
      real(dp) :: s0pp_p_c13 = 7.29d-4
! d2S/dE2 of O16(p,g)F17 (b/keV), with use_new_nuclear_rates
      real(dp) :: s0pp_p_o16 = 0.0d0
! dS/dE of Be7(p,g)B8 (b), with use_new_nuclear_rates
      real(dp) :: s0p_be7_p = -3.12d-5
! d2S/dE2 of Be7(p,g)B8 (b/keV), with use_new_nuclear_rates
      real(dp) :: s0pp_be7_p = -2.288d-7
! The pep rate relative to pp, with use_new_nuclear_rates
      real(dp) :: s0_pep = 3.5734d-6
! S0 of He3(p,e+ nu)He4 (keV b), with use_new_nuclear_rates
      real(dp) :: s0_hep = 8.6d-20

! ---- The equation of state

! Use the OPAL 1995 EOS (a newer OPAL EOS selected wins)
      logical :: use_opal95_eos = .false.
! Use the OPAL 2001 EOS (the 2006 EOS, if also selected, wins)
      logical :: use_opal2001_eos = .false.
! Use the OPAL 2006 EOS
      logical :: use_opal2006_eos = .false.
! The OPAL EOS vintage the run uses: 0 none, 1 1995, 2 2001, 3 2006 (eos_lib's
! opal_vintage_*), the newest selected (set by derive_controls)
      integer :: opal_eos_vintage = 0
! The OPAL EOS table the run reads: the newest vintage selected
! (use_opal2006_eos, use_opal2001_eos, use_opal95_eos), blank for none (set by
! derive_controls)
      character(len=256) :: opal_eos_table_path = repeat(achar(0), 256)
! Use the MHD EOS
      logical :: use_mhd_eos = .false.
! Use the SCV EOS where it is valid and the OPAL EOS is not
      logical :: use_scv_eos = .false.
! Apply the Debye-Hueckel correction to the simplified EOS
      logical :: use_debye_huckel_correction = .false.
! Debye-Hueckel correction: the lower degeneracy (eta) bound of its domain
      real(dp) :: debye_huckel_eta_min = -1.0d0
! Debye-Hueckel correction: the upper degeneracy (eta) bound of its domain
      real(dp) :: debye_huckel_eta_max = 1.0d0
! Use numerical derivatives of the EOS instead of analytic ones (much slower)
      logical :: use_numerical_derivatives = .false.

! ---- The opacities

! Opacity tables: tolerance for matching the model's metallicity to a table's
      real(dp) :: metal_fraction_match_tolerance = 1.0d-8
! log10 T below which only the molecular opacities are used
      real(dp) :: molecular_opacity_logt_min = 4.0d0
! log10 T above which only the atomic opacities are used (a ramp joins the
! two)
      real(dp) :: molecular_opacity_logt_max = 4.1d0
! The LAOL 1989 interpolators extrapolate linearly within
! log10(laol_extrapolation_tolerance) of a table edge
      real(dp) :: laol_extrapolation_tolerance = 10.0d0
! Use the conductive opacity tables (conductive_opacity_table_file)
      logical :: use_conductive_opacity = .false.

! ---- Convection: overshoot and semiconvection

! Deep mixing: below 1, the composition step mixes every shell from the centre
! out to the mass envelope_logp_step times the star's (a mass fraction,
! whatever the name says); 1 or more is off
      real(dp) :: envelope_logp_step = 1.0d0
! Enable convective-core overshoot
      logical :: core_overshoot_active = .false.
! Core overshoot extent in pressure scale heights
      real(dp) :: overshoot_alpha_core = 0.0d0
! Overshoot below the surface convection zone, by overshoot_alpha_envelope
! pressure scale heights
      logical :: envelope_overshoot_active = .false.
! Envelope overshoot extent in pressure scale heights
      real(dp) :: overshoot_alpha_envelope = 0.0d0
! Overshoot at intermediate convection zones (neither the core nor the surface
! convection zone), by overshoot_alpha_min
      logical :: intermediate_overshoot_active = .false.
! Overshoot of intermediate convection zones (intermediate_overshoot_active),
! in pressure scale heights
      real(dp) :: overshoot_alpha_min = 0.0d0
! The overshoot region below the surface convection zone is adiabatic: its
! radiative points take the adiabatic gradient
      logical :: adiabatic_overshoot_active = .false.
! Cap the core overshoot extent at overshoot_beta_core times the core radius
      logical :: overshoot_max_extent_flag = .false.
! With overshoot_max_extent_flag: the core overshoot extent is capped at this
! fraction of the core radius
      real(dp) :: overshoot_beta_core = 0.15d0
! Semiconvection treatment (conflicts with core overshoot)
      logical :: use_semiconvection = .false.

! ---- The atmosphere and the envelope: the surface boundary, the integrations
!      and their Bulirsch-Stoer integrator (numerics' limits: 1 to 11 stages,
!      orders 2 to 7; check_controls checks them), the envelope-fit triangle,
!      the spots and the Allard tables

! Surface boundary condition: 0 Eddington gray T(tau), 1 Krishna-Swamy T(tau),
! 2 Harvard-Smithsonian atmosphere, 3 Kurucz tables, 4 Allard tables, 5
! Kurucz/Castelli tables
      integer :: atm_choice = 0
! Integrate the full envelope (false: the atmosphere to the photosphere only)
      logical :: calc_envelope_flag = .true.
! Atmosphere integration: error tolerance
      real(dp) :: atm_error_tol = 3.0d-4
! Atmosphere integration: the starting density (g/cm^3) at optical depth near
! zero; divided by 10 on each retry
      real(dp) :: atm_step_initial = 1.0d-10
! Atmosphere integration: the first step
      real(dp) :: atm_step_begin = 1.0d-1
! The atmosphere integration's smallest step (in log tau)
      real(dp) :: atm_step_min = 1.0d-1
! The atmosphere integration's largest step (in log tau)
      real(dp) :: atm_step_max = 5.0d-1
! Envelope integration: error tolerance
      real(dp) :: env_error_tol = 3.0d-4
! Envelope integration: first step
      real(dp) :: env_step_begin = 1.0d-1
! Envelope integration: minimum step
      real(dp) :: env_step_min = 1.0d-1
! The envelope integration's largest step (in log P)
      real(dp) :: env_step_max = 5.0d-1
! The envelope integration stops when log(m/M) is within this fraction of senv
! (also the atmosphere integration's step tolerance)
      real(dp) :: tolerance_fraction = 1.0d-3
! The Bulirsch-Stoer integrator (envelope and atmosphere): the number of
! stages, 1 to 11
      integer :: max_stage_index = 11
! The Bulirsch-Stoer integrator (envelope and atmosphere): the extrapolation
! order, 2 to 7 and at most max_stage_index
      integer :: extrap_order = 7
! Envelope-fit triangle: size in log Teff
      real(dp) :: tri_delta_teffl = 1.0d-2
! Envelope-fit triangle: size in log L
      real(dp) :: tri_delta_logl = 8.0d-2
! Force a fresh envelope-fit triangle at run start
      logical :: new_triangle_init_flag = .false.
! Envelope generation: the atmosphere step size
      real(dp) :: atm_step_size = 0.5d0
! Envelope generation: the envelope step size
      real(dp) :: envelope_step_size = 0.5d0
! Inert: read and stored, used nowhere (the outer layers are built by
! stitched_model)
      logical :: envelope_generation_flag = .false.
! Allard atmosphere tables: [Fe/H]
      real(dp) :: allard_target_feh = 0d0
! Allard atmosphere tables: [alpha/Fe]
      real(dp) :: allard_target_alpha = 0d0
! Allard atmosphere tables: the boundary at tau = 100 instead of the
! photosphere
      logical :: allard_use_tau100 = .false.
! Starspots: the filling factor
      real(dp) :: spot_filling_factor = 0.00d0
! Starspots: the spot temperature contrast (with spot_filling_factor > 0; 0.85
! is typical of sunspots)
      real(dp) :: spot_temp_contrast = 1.00d0
! Starspots: decrease the spot contrast with depth, at constant B
      logical :: spot_depth_varies = .false.

! ---- Rotation and the transport of angular momentum (a rotating run's
!      settings clamp cz_rotation_law_exponent to [-2, 0],
!      setup/rotation_settings_builder.f90)

! Rotation: evolve the rotation curve and the angular momentum transport
      logical :: rotation_active = .false.
! Initialize solid-body rotation at omega = initial_omega
      logical :: set_initial_omega = .false.
! Initial rotation rate when set_initial_omega
      real(dp) :: initial_omega = 0.0d0
! Exponent alpha of the surface convection zone's rotation law omega ~ r^alpha
! (0 or above: solid body; -2 or below: constant specific angular momentum);
! also gives the omega of shells added at the surface and the specific angular
! momentum the wind removes
      real(dp) :: cz_rotation_law_exponent = 0.0d0
! Rotationally induced mixing and hydrodynamic angular momentum transport
      logical :: instability_transport_active = .false.
! The rotating star's shape: the relative tolerance on r0^3 of the
! equipotential radius iteration (its cube root that on r0 and eta2); not a
! weight
      real(dp) :: fpft_average_weight = 1.0d-36
! The rotating star's shape: the most iterations of eta2 and r0 in turn
      integer :: fpft_iter1 = 5
! The rotating star's shape: the most Newton steps of each r0 solve
      integer :: fpft_iter2 = 20
! The passes of the structure-rotation iteration per step (each pass
! converges the structure, mixes the convection zones and updates the angular
! momentum and the shape); 1 without rotation
      integer :: num_rotation_structure_iters = 1
! The most iterations of each angular momentum transport sub-step (1 to 16;
! the controls read rejects a rotating run outside that)
      integer :: max_diffusion_iters = 10
! Sub-timestep factor for rotation diffusion
      real(dp) :: diffusion_timestep_factor = 1.0d-2
! Angular momentum transport: the convergence tolerance of the transport
! iterations
      real(dp) :: convergence_tolerance = 1.0d-4
! Rotational mixing: scale the composition diffusion coefficient by the Zahn
! (1992) coupling factor (D_mix = D_omega f_Zahn) instead of the fixed
! efficiencies (mixing_velocity_scale and the per-mechanism mixing scales)
      logical :: velocity_form_correction = .false.
! Treat meridional circulation as diffusion plus advection (Zahn) instead of
! diffusion (stops a rotating run: the path is commented out)
      logical :: use_diffusion_advection_transport = .false.
! No angular momentum transport in the convective core
      logical :: no_am_transport_in_core = .false.
! Force solid-body rotation (sets solid_body_mode_flag = 1)
      logical :: force_solid_body_rotation = .false.
! Angular momentum coupling: 0 none, 1 solid body everywhere, 2 solid-body
! core with a decoupled convection zone, 3 solid-body core coupled at the base
! of the convection zone
      integer :: solid_body_mode_flag = 0
! Inert: read and used nowhere (the F77 code only echoed it, as COUPLING, in
! the run log's header)
      real(dp) :: coupling_factor = 1.0d0
! Scale of the angular momentum diffusion coefficients in the
! diffusion-advection transport
      real(dp) :: difad_velocity_scale = 1.0d0
! Ratio of the chemical mixing to the angular momentum diffusion coefficients
! (calibrated to Li depletion)
      real(dp) :: mixing_velocity_scale = 1.0d0
! Scale of the Eddington-Sweet circulation velocity in the angular momentum
! transport
      real(dp) :: es_velocity_scale = 1.0d0
! Scale of the GSF and ABCD instability velocities in the angular momentum
! transport
      real(dp) :: gsf_velocity_scale = 1.0d0
! Scale of the secular shear instability velocity in the angular momentum
! transport
      real(dp) :: secular_shear_velocity_scale = 1.0d0
! Critical Reynolds number for the onset of the secular shear instability
      real(dp) :: critical_reynolds = 1.0d3
! Unused: read by the legacy reader and ignored
      integer :: enable_eddington_sweet = 1
! Unused: read by the legacy reader and ignored
      integer :: enable_mu_gradient = 1
! Unused: read from the legacy decks and ignored (the mu-gradient scale had no
! effect, docs/bug_sweep.md O15; dropped from the inlist format)
      real(dp) :: mu_gradient_scale = 1.0d0
! Scale of the Eddington-Sweet circulation's chemical mixing relative to its
! angular momentum transport
      real(dp) :: es_mixing_scale = 1.0d0
! Scale of the secular shear's chemical mixing relative to its angular
! momentum transport
      real(dp) :: secular_shear_mixing_scale = 1.0d0
! Scale of the GSF instability's chemical mixing relative to its angular
! momentum transport
      real(dp) :: gsf_mixing_scale = 1.0d0
! Rotational transport: the Solberg-Hoiland instability's velocity scale (0 =
! off; it can act only where the dynamical shear does; MESA D_SH_factor)
      real(dp) :: solberg_hoiland_velocity_scale = 0.0d0
! Rotational mixing: the Solberg-Hoiland instability's mixing efficiency
      real(dp) :: solberg_hoiland_mixing_scale = 1.0d0
! The stability condition of the GSF and ABCD instabilities (1 is recommended)
      integer :: gsf_inhibition_mode = 1
! Add a constant background diffusion (constant_background_diffusion_coeff) to
! the angular momentum transport
      logical :: use_constant_background_diffusion = .false.
! Constant background diffusion coefficient for angular momentum (cm^2/s)
      real(dp) :: constant_background_diffusion_coeff = 2.5d4

! ---- Magnetic braking (rotation/magnetic_braking.f90's derive_wind_law) and
!      disk locking

! Apply magnetic wind angular-momentum loss
      logical :: use_wind_torque = .true.
! Matt wind law: V13 (van Saders & Pinsonneault 2013), K97 (Krishnamurthi et
! al. 1997), or CUS (the pmm_exponent_* as given)
      character(len=3) :: wind_law_name = 'K97'
! Kawaler wind law: the index n (Jdot ~ Omega**(1+4n/3) ...); K97 takes m =
! n/3
      real(dp) :: wind_law_index = 1.5d0
! The constant factor of both wind laws (1 reproduces the solar angular
! momentum loss)
      real(dp) :: wind_law_constant = 1.0d0
! Saturation angular velocity of the wind law (rad/s); solar_omega_sat below
! 1000 overrides it in solar units
      real(dp) :: wind_saturation_omega = 3.0d-4
! The wind's saturation omega in units of the solar omega (pmm_solar_omega);
! below 1000 it overrides wind_saturation_omega
      real(dp) :: solar_omega_sat = 1000.0d0
! Use the Matt et al. (2012) magnetized wind law (wind_law_name) instead of
! the Kawaler law
      logical :: use_pmm_wind_law = .false.
! PMM wind law: exponent a (Lx/Lbol ~ omega**a); K97 and V13 set it
      real(dp) :: pmm_exponent_a = 2.0d0
! PMM wind law: exponent b, with Bmag ~ (omega R**c)**b; K97 and V13 set it
      real(dp) :: pmm_exponent_b = 1.0d0
! PMM wind law: exponent c, with Bmag ~ (omega R**c)**b; K97 and V13 set it
      real(dp) :: pmm_exponent_c = 2.0d0
! Matt wind law: exponent d of the convective turnover time in Bmag (CUS only;
! K97 and V13 set it)
      real(dp) :: pmm_exponent_d = 0.0d0
! Matt wind law: exponent m of the torque (CUS only; K97 and V13 set it)
      real(dp) :: pmm_exponent_m = 0.5d0
! Matt wind law: the solar wind torque (cgs), the zero point of the loss rate
! (adjust with wind_law_constant)
      real(dp) :: pmm_norm_jdot = 1.32d30
! Matt wind law: the solar mass-loss rate (g/s), its zero point
      real(dp) :: pmm_norm_mdot = 1.27d12
! Matt wind law: log10 of the solar photospheric pressure (dyn/cm^2), the zero
! point of the field strength
      real(dp) :: pmm_solar_pressure = 4.9304d0
! Matt wind law: the solar angular velocity (rad/s), the zero point of the
! torque
      real(dp) :: pmm_solar_omega = 2.836d-6
! Matt wind law: the solar convective turnover time (s), the zero point of the
! Rossby number
      real(dp) :: pmm_solar_turnover_timescale = 9.4805d5
! Matt wind law: scale the loss rate with the Rossby number (use with a
! saturation threshold)
      logical :: scale_by_rossby_number = .true.
! Matt wind law: scale the field with the equipartition B (set by
! wind_law_name for K97 and V13)
      logical :: scale_by_b_field = .false.
! Inert: the rotation-scaled solar wind is commented out
! (core/step/mass_change.f90); the controls read warns when it is set
      logical :: use_rotation_scaled_solar_wind = .false.
! Inert: read and used nowhere (the Matt law is use_pmm_wind_law); the
! controls read warns when it is set
      logical :: use_pmm_wind = .false.
! Inert: read and used nowhere (the saturation is wind_saturation_omega or
! solar_omega_sat); the controls read warns when it is set
      real(dp) :: pmm_omega_sat = 0
! Disk locking: hold the surface rotation at disk_omega_rad_s until
! disk_locking_age_gyr
      logical :: disk_locking_active = .false.
! Disk locking: the age (Gyr) up to which the surface rotation is held at
! disk_omega_rad_s (some older deck comments give it in yr, which the code
! never read)
      real(dp) :: disk_locking_age_gyr = 0.0d0
! Disk locking: the surface angular velocity (rad/s) held until
! disk_locking_age_gyr (some older deck comments call it a period in days,
! which the code never read)
      real(dp) :: disk_omega_rad_s = 7.2722d-6

! ---- Mass loss and accretion

! Mass change: a constant rate (mass_accretion_rate) or the Reimers wind;
! switched off for the run when the disk is exhausted
      logical :: use_mass_accretion = .false.
! Mass change rate (Msun/yr; negative for loss) when use_mass_change is set
      real(dp) :: mass_accretion_rate = -1.0d-14
! Mass change: the timestep changes the surface convection zone's mass by at
! most this fraction
      real(dp) :: accretion_cz_fraction = 0.1d0
! Mass change: the timestep changes the star's mass by at most this fraction
      real(dp) :: accretion_total_fraction = 1.0d-2
! Mass fractions of the accreted material (star%xa layout)
      real(dp) :: accreted_composition(15) = (/0.71668d0, 0.265721d0,0.01757d0, &
           2.9d-5,3.013d-3,3.385d-5,9.346d-4, 0.0d0,8.462d-3,0.0d0,1.696d-5,0.0d0, &
           2.0d-9,2.0d-9,3.0d-11/)
! Enable Reimers red-giant mass loss
      logical :: use_reimers_wind = .false.
! Reimers mass-loss efficiency eta
      real(dp) :: reimers_scaling_factor = -4.0d-13

! ---- Gravitational settling and element diffusion

! Gravitational settling of helium (and hydrogen)
      logical :: diffuse_helium_active = .false.
! Gravitational settling of the metals (implies helium settling and the Thoul
! et al. routines)
      logical :: use_diffusion_z = .false.
! Include Li/Be in gravitational settling
      logical :: diffuse_lithium = .false.
! Thoul et al. diffusion coefficients
      logical :: use_thoul_diffusion = .false.
! Use the published fit to the Thoul et al. settling coefficients instead of
! the full solution (not recommended)
      logical :: use_thoul_fit = .false.
! Use the new diffusion routines
      logical :: use_new_diffusion_routines = .false.
! Gravitational settling: the factor on helium's settling coefficients
      real(dp) :: settling_y_factor = 1.0d0
! Gravitational settling: the factor on the metals' settling coefficients
! (a Monte-Carlo run samples it)
      real(dp) :: settling_z_factor = 1.0d0
! Gravitational settling: solver tolerance
      real(dp) :: settling_tolerance = 1.0d-8
! Settling: the largest number of iterations of the abundance-change solver
      integer :: settling_num_iterations = 10
! Gravitational settling: sub-timestep as a fraction of the timestep
      real(dp) :: settling_timestep_fraction = 0.1d0
! The Coulomb logarithm of the settling coefficients (4, a constant 2.2, is
! recommended)
      integer :: coulomb_log_choice = 1
! Settling stops (for the rest of the run card) once the central X falls below
! this
      real(dp) :: hydrogen_diffusion_floor = 1.0d-3
! Settling's outer boundary moves inward to the first zone where Y exceeds
! this
      real(dp) :: helium_diffusion_min = 1.0d-3
! Constant added diffusion coefficient (mixing)
      real(dp) :: constant_mixing_coeff = 1.0d0
! Constant coefficient mimicking mixing to reduce settling
      real(dp) :: constant_settling_reduction = 1.0d0

! ---- The starting model: its composition array, a species rescaled, the core
!      extended inward, the helium flash

! Extended composition: also follow H2, Li6, Li7 and Be9 (15 species)
      logical :: use_extended_composition = .false.
! The surface mixture's relative mass fractions, one per slot of
! state/composition_layout.f90 (Na, Al, Mg, Fe, Si, C, H, O, N, Ar, Ne, He;
! H and He are set from the envelope X and Z): the starting model's envelope
! number fractions and the EOS mixture are made from the ten metals'
      real(dp) :: mixture_weights_seed(12) = (/0.001999d0, 0.003238d0, 0.037573d0, &
            0.071794d0, 0.040520d0, 0.173285d0, 0.000000d0, 0.482273d0, 0.053152d0, &
            0.005379d0, 0.098668d0, 0.000000d0/)
! Set one species of the starting model (change_species_name) to
! new_species_value
      logical :: rescale_species_active = .false.
! Element id whose abundance to override in the start model (with
! change_species_value)
      character(len=3) :: change_species_name = '   '
! With rescale_species_active: 'ABS' makes new_species_value an absolute mass
! fraction; anything else, relative to hydrogen
      character(len=3) :: new_species_value_mode = repeat(achar(0), 3)
! New abundance for change_species_name
      real(dp) :: new_species_value = 1.3d1
! With rescale_species_active: the composition slot of change_species_name
! (set by the reader)
      integer :: new_species_index
! With rescale_species_active: new_species_value is relative to hydrogen (from
! new_species_value_mode)
      logical :: value_relative_to_h
! Add shells inward of the starting model's innermost shell
      logical :: extend_core_inward = .false.
! With extend_core_inward: the number of points added inward (recomputed from
! core_mass_reduction_factor)
      integer :: num_core_shells_added = 0
! With extend_core_inward: the factor by which the innermost mass is reduced
      real(dp) :: core_mass_reduction_factor = 1.0d0
! Helium flash run: no rezoning, the starting model written to the .mod file
! for timestep cuts, and flash-specific solver handling
      logical :: helium_flash_active = .false.

! ---- Monte Carlo runs

! Monte-Carlo multi-realization mode
      logical :: monte_carlo_active = .false.
! First Monte-Carlo realization index
      integer :: mc_run_start = 1
! Last Monte-Carlo realization index
      integer :: mc_run_end = 1

end module controls_buffer

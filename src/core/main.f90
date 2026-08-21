!----------------------------------------------------------------------
! main
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original main.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model). This is the LAST file converted
! in the project (256 of 257 files were already converted before this
! one); every callee below was cross-checked against its own current
! (already-converted) signature, argument by argument.
!
! This is the top-level PROGRAM unit (not a SUBROUTINE) -- the original
! main.f has no PROGRAM statement (a legal, unnamed F77 main program);
! it is given the name "main" here, matching the file name, since
! free-form Fortran expects "end program <name>" to pair with
! "program <name>".
!
! MAIN drives the whole simulation: it reads user parameters (parmin)
! and physics tables (setups), then for each Monte-Carlo realization
! (if any) and each "kind card" run (the DO 200 NK loop) reads or
! rescales the starting model (starin), evolves it for the requested
! number of models (the DO 100 MODELN loop) via the Newton-Raphson
! relaxation driver (crrect, 4 iteration levels), mass loss (massloss),
! mixing (mix/mixcz), light-element burning (lirate88/liburn), the
! rotation curve (getw/fpft), and re-meshing (hpoint), writes output
! after every converged model (wrtout, plus wrtlst/putstore for stored
! models), and (for automatic solar/stellar calibration) checks and
! iterates the calibration parameters (chkcal/chkscal/setcal/setscal).
!
! CROSS-CALLEE NAMING NOTE: several locals here are threaded into more
! than one already-converted callee, and those callees do not always
! agree on a name for the same physical data (this file, being the
! PROGRAM unit, is free to choose its own local names; the callees'
! names were fixed by earlier batches). Judgment calls made below, all
! verified against the actual physics/usage in this file and cross-
! checked against each callee's current signature:
!   - BL (log10 L/Lsun, surface value) is named log10_luminosity,
!     matching crrect.f90/starin.f90/chkcal.f90/chkscal.f90's own
!     slot names for this exact scalar; wrtlst.f90/massloss.f90/
!     getnewenv.f90/putstore.f90/wrtout.f90 instead call the same
!     scalar "log_luminosity_lsun".
!   - HL (linear L/Lsun, per-shell array) is named luminosity_lsun,
!     matching crrect.f90/starin.f90; findsh.f90/htimer.f90/mix.f90/
!     hpoint.f90/getnewenv.f90/wrtlst.f90/wrtout.f90/putstore.f90/
!     massloss.f90/getw.f90 instead call the same array
!     "log_luminosity" (or plain "luminosity"/"luminosity_lsun" in
!     findsh.f90's own case) -- a misnomer inherited from those files'
!     own earlier conversions in most cases; out of scope to fix here.
!   - M (number of mass zones) is named num_zones, the majority
!     spelling among callees (massloss/mix/hpoint/getw/convec/
!     lirate88/liburn/getnewenv/mixcz); starin.f90/wrtlst.f90/
!     wrtout.f90/putstore.f90 instead call it num_shells,
!     crrect.f90/fpft.f90/findsh.f90/htimer.f90 call it num_points,
!     wrtmonte.f90 keeps it cryptic as "m".
!   - JENV (envelope CZ base zone index) is named
!     envelope_cz_bottom_index, matching wrtlst.f90/wrtout.f90/
!     putstore.f90/wrtmonte.f90 (majority); findsh.f90 calls it
!     envelope_edge, convec.f90/mix.f90 call it envelope_cz_edge,
!     hpoint.f90 calls it convective_envelope_edge_zone,
!     massloss.f90 calls it envelope_boundary_zone, starin.f90 calls
!     it envelope_zone_index.
!   - JCORE (core CZ top zone index) is named core_cz_top_index,
!     matching wrtlst.f90/wrtout.f90/putstore.f90/wrtmonte.f90
!     (majority); findsh.f90 calls it core_edge, convec.f90/mix.f90
!     call it core_cz_edge, hpoint.f90 calls it
!     convective_core_edge_zone, htimer.f90 calls it
!     convective_core_edge_zone as well.
!   - LARGE (T if the correction routine has diverged / the model has
!     failed) is named model_diverged_flag; crrect.f90's own slot name
!     for this same local is corrections_too_large, starin.f90's is
!     model_failed_flag.
!   - LNEW (T if a fresh envelope-fit triangle must be generated) is
!     named recompute_envelope_triangle; crrect.f90's own slot name is
!     start_new_triangle, starin.f90's is envelope_recomputed_flag,
!     getnewenv.f90's is new_points_added_flag.
!   - LRESET is named reset_triangle, matching crrect.f90 (used in 4 of
!     its 5 call sites here); hpoint.f90's own slot name is
!     point_reset_flag.
!   - DELTSH is a dual-purpose scratch slot: starin.f90 treats it as
!     delta_time_abs (= |delta_time|, set by starin itself) while
!     htimer.f90 treats it as hydrogen_dt (the H-burning-shell
!     timestep governor, computed fresh by htimer). It is named
!     hydrogen_dt here since that is its dominant role (read back by
!     the next model's htimer.f90 call); starin's transient use of the
!     same storage as |delta_time| is noted here for the record.
!   - ISTORE is named istore_flag, matching starin.f90; hpoint.f90's
!     own slot name for the same scalar is envelope_store_index.
!   - SENVOLD (set by massloss, read by getnewenv) is named
!     target_envelope_mass, matching getnewenv.f90 (how it is actually
!     consumed); massloss.f90's own slot name for the same scalar is
!     old_log_envelope_mass_fraction.
!   - V(12) is a dual-purpose work array: passed to setups.f90 as
!     laol_work_array (a scratch/table-loading array) and to
!     starin.f90 as species_mix_weights (the VNEW mixture-weight
!     copy). Named mixture_weights here as the more physically
!     meaningful of the two uses; flagged here since the two callees'
!     own names disagree on its purpose.
!   - FPATM/FPENV/FPMOD are named pulse_atm_path/pulse_env_path/
!     pulse_mod_path, matching pdist.f90 (their only consumer besides
!     parmin.f90, which keeps its own dummy spelling fpatm/fpenv/fpmod
!     per that file's NAMELIST-driven naming constraint -- COMMON/
!     argument passing is positional, so this does not need to match).
!   - The remaining PARMIN/SETUPS file-path locals (FALEX06, FALLARD,
!     FATM, FFERMI, FKUR, FKUR2, FLAOL, FLAOL2, FLIV95, FLLDAT, FMHD1-
!     FMHD8, FOPAL2, FPUREZ, FSCVH, FSCVHE, FSCVZ, OPECALEX) are named
!     to match setups.f90's own (descriptive) dummy-argument spelling
!     for the same slot, since parmin.f90 keeps its own dummy names at
!     cryptic lowercase spelling only (per that file's NAMELIST
!     constraint, documented in parmin.f90's own header) -- again,
!     positional COMMON/argument passing means this file's own local
!     names are free to differ from parmin.f90's.
!
! COMMON BLOCK NOTE: every COMMON block below was grepped against all
! other already-converted .f90 files and, where established, reuses
! the exact member names found there (per the project's COMMON-block
! rule). Two conflicts were found and resolved, both documented again
! at the point of declaration below:
!   - common/difus/'s 4th member (originally ITDIF2) is established as
!     "max_iterations" in checkc.f90, but this file already uses
!     max_iterations as a local name (matching crrect.f90's own slot
!     name for the NITER argument, threaded through all 4 of this
!     file's crrect calls) -- reusing it for the (here-unused) common
!     member would alias two different entities to the same name, so
!     the common member is instead kept at its lowercased-cryptic
!     spelling, itdif2, here.
!   - common/cals2/'s 1st member and common/calstar/'s 2nd member are
!     both independently established as "luminosity_tolerance" (in
!     chkcal.f90 and chkscal.f90 respectively) -- distinct physical
!     quantities (solar-calibration L tolerance vs. target-star L
!     tolerance) that never previously co-occurred in the same file.
!     common/cals2/'s member keeps the established "luminosity_
!     tolerance"; common/calstar/'s member is disambiguated here as
!     "target_star_luminosity_tolerance".
!   - Several common blocks establish "_placeholder"-suffixed names
!     for members that are unused in the file(s) that first declared
!     them, but that this file actually reads/writes (common/acdpth/'s
!     ageout_placeholder/lclcd_placeholder/ljlast_placeholder/
!     ljwrt_placeholder; common/rotprt/'s lprt0_placeholder;
!     common/chrone/'s lrwsh_placeholder; common/cenv/'s
!     lnew0_placeholder). Per the precedent set by envint.f90 (keeps
!     lclcd_placeholder despite noting its own active use) and
!     getw.f90 (keeps lprt0_placeholder despite noting its own active
!     use), these established names are reused verbatim here too
!     rather than renamed, even though this file actively assigns/
!     reads them; each active use is called out in a comment at its
!     point of use below.
program main

! the array size, i.e. max # of shells is specified in the parameter
! statement. it defines JSON. to change the array size do a global
! change on "JSON=2000" or whatever.
      use fluxes_lib
      use engeb_diag_lib
      use light_burn_lib
      use turnover_lib
      use oldmod_lib
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000
      integer, parameter :: numtt = 70
      integer, parameter :: numd = 19
      integer, parameter :: numx = 10
      integer, parameter :: numz = 13
      integer, parameter :: numxz = 126

! DBGLAOL - to save space make tables single precision
! MHP 8/25 Removed unused variables and added pass-through variables
! See the CROSS-CALLEE NAMING NOTE above for the file-path locals
! (declared further below, in the --- locals --- section, using
! setups.f90's/pdist.f90's own descriptive dummy-argument spelling).



! common/lufnm/ is commented out in the original (superseded by
! passing file names directly, per the MHP 8/25 note below) and is not
! declared.
! MHP 8/25 Removed character file names from common block
! DBGLAOL
! MHP 8/25 Removed character file names from common block
! common/nwlaol/: not used in this file's logic; declared only to
! preserve layout. Naming matches envint.f90/getopac.f90.
      double precision :: olaol(12,104,52), oxa(12), ot(52), orho(104), &
           tollaol
      integer :: iolaol, numofxyz, numrho, numt, iopurez
      logical :: llaol, use_pure_z_table
      common/nwlaol/ olaol, oxa, ot, orho, tollaol, iolaol, numofxyz, &
           numrho, numt, llaol, use_pure_z_table, iopurez

! common/track/: track_file_version (ITRVER); not used in this file's
! logic. Naming matches wrthead.f90.
      integer :: track_file_version
      common/track/ track_file_version

! common/label/: initial_envelope_x/initial_envelope_z (XENV0/ZENV0),
! both used here. Naming matches wrthead.f90.
      double precision :: initial_envelope_x, initial_envelope_z
      common/label/ initial_envelope_x, initial_envelope_z



! common/cenv/: requested_envelope_mass/change_envelope_mass_flag are
! used here, matching starin.f90's own active naming for the same two
! slots (that file explains why they were given descriptive names
! rather than kept as placeholders). lnew0_placeholder (5th member) is
! also actively used in this file (see the header COMMON BLOCK NOTE
! above for why it nonetheless keeps its established placeholder
! name); tri_delta_teffl/tri_delta_logl are unused placeholders here.
      double precision :: tri_delta_teffl, tri_delta_logl, &
           requested_envelope_mass
      logical :: change_envelope_mass_flag, lnew0_placeholder
      common/cenv/ tri_delta_teffl, tri_delta_logl, requested_envelope_mass, &
           change_envelope_mass_flag, lnew0_placeholder

! common/ckind/: all used here. Naming matches chkcal.f90/rscale.f90/
! wrthead.f90/starin.f90.
      double precision :: rescale_params(4,50)
      integer :: num_models(50), rescale_kind(50)
      logical :: first_call_flag(50)
      integer :: num_runs
      common/ckind/ rescale_params, num_models, rescale_kind, &
           first_call_flag, num_runs









! MHP 8/17 ADDED EXCEN, C_2 TO COMMON BLOCK FOR MATT ET AL. 2012 CENT. TERM
! common/cwind/: only ljdot0 used here; the rest are unused
! placeholders. Naming matches cowind.f90/getw.f90.
      double precision :: wind_saturation_omega, exmd, &
           wind_law_omega_exponent, extau, exr, exm, exl, expr, &
           constfactor, structfactor, excen, c_2
      logical :: ljdot0
      common/cwind/ wind_saturation_omega, exmd, wind_law_omega_exponent, &
           extau, exr, exm, exl, expr, constfactor, structfactor, excen, &
           c_2, ljdot0


      logical :: helium_flash_active
      common/heflsh/ helium_flash_active




! 7/91 ENTROPY TERM COMMON BLOCK ADDED.
      double precision :: temperature_entropy_term(json), &
           pressure_entropy_term(json), luminosity_entropy_term(json), &
           radius_entropy_term(json)
      common/entrop/ temperature_entropy_term, pressure_entropy_term, &
           luminosity_entropy_term, radius_entropy_term


! MHP 9/94
! common/rotprt/: lprt0_placeholder actively used here (see the header
! COMMON BLOCK NOTE for why the established placeholder name is kept).
! Naming matches wrtout.f90/getw.f90.
      logical :: lprt0_placeholder
      common/rotprt/ lprt0_placeholder




! MHP 8/25 Removed character file names from common block
! common/chrone/: lrwsh_placeholder actively used here (see the header
! COMMON BLOCK NOTE for why the established placeholder name is kept);
! the rest are unused placeholders. Naming matches wrthead.f90.
      logical :: lrwsh_placeholder, isochrone_output_active
      integer :: isochrone_file_unit
      common/chrone/ lrwsh_placeholder, isochrone_output_active, &
           isochrone_file_unit

! DBG 1/92 let XENV0, ZENV0, and CMIXL be arrays so can change during
! a set of runs.
! common/newxym/: all used here. Naming matches chkcal.f90/chkscal.f90.
      double precision :: initial_x_array(50), initial_z_array(50), &
           mixing_length_array(50)
      logical :: has_senv0_array(50)
      double precision :: senv0_array(50)
      common /newxym/ initial_x_array, initial_z_array, &
           mixing_length_array, has_senv0_array, senv0_array

! MHP 8/25 Removed character file names from common block
! common/zramp/: nk used here (the outer kind-card run index -- this
! IS the same storage as the DO 200 loop variable below, not a
! separate local); the rest unused placeholders. Naming matches
! wrthead.f90/setllo.f90.
      double precision :: rsclzc(50), rsclzm1(50), rsclzm2(50)
      integer :: iolaol2, ioopal2, nk
      logical :: use_z_ramp
      common/zramp/ rsclzc, rsclzm1, rsclzm2, iolaol2, ioopal2, nk, &
           use_z_ramp

! DBG 4/26/94 Tired of not have access to current age of model so...
! common/theage/: dage used here. Naming matches pdist.f90/mix.f90/
! microdiff_cod.f90/microdiff_run.f90 (established there as an unused
! placeholder; actively used here).
      double precision :: dage
      common/theage/ dage

      double precision :: luminosity_breakdown(8), trial_log_luminosity(3), &
           trial_log_temperature(3), fit_point_pressure(3), &
           fit_point_temperature(3), fit_point_radius(3), &
           envelope_fit_coeffs(9)
      double precision :: log_mass(json), luminosity_lsun(json), &
           log_radius(json), log_pressure(json), log_temperature(json), &
           log_density(json)
      logical :: convective_flag(json)
      double precision :: composition(15,json), enclosed_mass(json), &
           shell_mass(json), gravitational_luminosity(json)
      double precision :: elim_coeff(4,2,json), elim_rhs(4,json)
      double precision :: log_pressure_delta(json), log_temperature_delta(json)
      double precision :: pressure_rotation_factor(json), &
           temperature_rotation_factor(json)
      double precision :: max_residual(4)
      integer :: max_correction_index(4)
      double precision :: surface_bc(6), stored_envelope_state(4)
      double precision :: trial_sign_flag
      double precision :: omega(json), moment_of_inertia(json), &
           specific_angular_momentum(json), kinetic_energy_rot(json), &
           mean_radius(json), eta_squared(json)
      double precision :: qiw(json), mean_gravity(json), &
           mixture_weights(12), kinetic_energy_rot_old(json)
!     MHP 10/24 FLAG FOR END OF RUN
      logical :: end_kind_flag

! JVS 02/12
! common/stch/: not used in this file's logic; declared only to
! preserve layout. First appearance of this common block in the
! converted sources.
      double precision :: composition_final(15,json), log_radius_final(json), &
           log_pressure_final(json), log_density_final(json), &
           log_mass_final(json), log_temperature_final(json)
      common/stch/ composition_final, log_radius_final, log_pressure_final, &
           log_density_final, log_mass_final, log_temperature_final

! DBG PULSE OUT 7/92
! variables used to contral output of pulsation models. model
! output after has traveled pomax in HR diagram
! LPOUT and POMAX added to control common block, rest in physics
! common/po/: only po_output_enabled used here; the rest unused
! placeholders. Naming matches pdist.f90.
      double precision :: po_weight_l, po_weight_teff, po_weight_age, &
           po_max_len_sq
      logical :: po_output_enabled
      common /po/ po_weight_l, po_weight_teff, po_weight_age, &
           po_max_len_sq, po_output_enabled


! common/calsun/: dlum_dx/drad_dx/dlum_dalpha/drad_dalpha/
! solar_calibration_active used here; log_l_prev/log_r_prev/delta_x/
! delta_alpha unused placeholders. Naming matches chkcal.f90.
      double precision :: dlum_dx, drad_dx, dlum_dalpha, drad_dalpha, &
           log_l_prev, log_r_prev, delta_x, delta_alpha
      logical :: solar_calibration_active
      common/calsun/ dlum_dx, drad_dx, dlum_dalpha, drad_dalpha, &
           log_l_prev, log_r_prev, delta_x, delta_alpha, &
           solar_calibration_active

! MHP 6/13 ADD OPTION TO CALIBRATE SOLAR Z/X, SOLAR Z/X, SOLAR AGE
! common/cals2/: all used here. Naming matches chkcal.f90; see the
! header COMMON BLOCK NOTE above regarding the collision between this
! block's own established luminosity_tolerance member and
! common/calstar/'s.
      double precision :: luminosity_tolerance, radius_tolerance, &
           zx_tolerance
      logical :: calibrate_solar_model, calibrate_solar_zx
      double precision :: target_solar_zx, target_solar_age
      common/cals2/ luminosity_tolerance, radius_tolerance, zx_tolerance, &
           calibrate_solar_model, calibrate_solar_zx, target_solar_zx, &
           target_solar_age


! DBG 12/94 added calibrate stellar model
! common/calstar/: all used here. Naming matches chkscal.f90, except
! the 2nd member -- see the header COMMON BLOCK NOTE above regarding
! the collision with common/cals2/'s own established
! luminosity_tolerance member.
      double precision :: target_luminosity_lsun, &
           target_star_luminosity_tolerance, target_teff, &
           target_radius_rsun, log_l_prev_model, log_r_prev_model, &
           age_at_target_radius, log_l_at_target_radius, &
           log_l_at_target_radius_prev_run, age_prev_model
      logical :: star_found_flag, specify_teff_flag, &
           just_passed_target_radius_flag, calibrate_star_flag
      common/calstar/ target_luminosity_lsun, &
           target_star_luminosity_tolerance, target_teff, &
           target_radius_rsun, log_l_prev_model, log_r_prev_model, &
           age_at_target_radius, log_l_at_target_radius, &
           log_l_at_target_radius_prev_run, age_prev_model, &
           star_found_flag, specify_teff_flag, &
           just_passed_target_radius_flag, calibrate_star_flag

! MHP 7/96 common block added for sound speed
! common/sound/: only sound_speed_output_active used here; gam1
! unused placeholder. Naming matches coefft.f90.
      double precision :: adiabatic_index_gamma1(json)
      logical :: sound_speed_output_active
      common/sound/ adiabatic_index_gamma1, sound_speed_output_active




! MHP 8/96
! added monte carlo parameters for metal diffusion, solar L, solar age.
! common/monte2/: all used here. Naming matches wrtmonte.f90.
      double precision :: s11_rate(1000), s33_rate(1000), s34_rate(1000), &
           s17_rate(1000), metal_to_h_ratio(1000), &
           helium_fraction_param(1000), diffusion_factor(1000), &
           luminosity_target(1000), age_target(1000)
      common/monte2/ s11_rate, s33_rate, s34_rate, s17_rate, &
           metal_to_h_ratio, helium_fraction_param, diffusion_factor, &
           luminosity_target, age_target

! MHP 3/96 central T,P,RHO added
! common/cent/: all used here. Naming matches wrtmonte.f90.
      double precision :: central_log10_temperature, central_log10_pressure, &
           central_log10_density, envelope_mass, envelope_radius
      common/cent/ central_log10_temperature, central_log10_pressure, &
           central_log10_density, envelope_mass, envelope_radius

! common/llot95a/: not used in this file's logic; declared only to
! preserve layout. Naming matches op95ztab.f90.
      double precision :: opal95_grid_logt(numtt), opal95_grid_x(numx), &
           opal95_grid_logr(numd), opal95_grid_z(numz), &
           opal95_full_opacity(numxz,numtt,numd)
      integer :: opal95_num_x_at_z(numz), opal95_table_start_index(numz)
      common /llot95a/ opal95_grid_logt, opal95_grid_x, opal95_grid_logr, &
           opal95_grid_z, opal95_full_opacity, opal95_num_x_at_z, &
           opal95_table_start_index

! as above for the model Z.
! common/llot95/: not used in this file's logic; declared only to
! preserve layout. Naming matches op95xtab.f90.
      double precision :: opal95_fixed_z_opacity(numx,numtt,numd), &
           opal95_fixed_z
      common /llot95/ opal95_fixed_z_opacity, opal95_fixed_z

! as above for desired surface value of X.
! common/llot95e/: not used in this file's logic; declared only to
! preserve layout. Naming matches op95xtab.f90.
      double precision :: opal95_surface_opacity(numtt,numd), opal95_surface_x
      common /llot95e/ opal95_surface_opacity, opal95_surface_x

! indices for interpolation in Z,X,T, and R
! common/op95indx/: not used in this file's logic; declared only to
! preserve layout. Naming matches getopal95.f90.
      integer :: opal95_index_z, opal95_index_x(4,4), opal95_index_t, &
           opal95_index_rho(4)
      common/op95indx/ opal95_index_z, opal95_index_x, opal95_index_t, &
           opal95_index_rho

! interpolation factors for Z,X,T, and R, as well as derivative
! factors for T and RHO.
! common/op95fact/: not used in this file's logic; declared only to
! preserve layout. Naming matches getopal95.f90.
      double precision :: opal95_weight_z(4), opal95_weight_x(4,4), &
           opal95_weight_t(4), opal95_dweight_t(4), opal95_weight_rho(4,4), &
           opal95_dweight_rho(4,4)
      common/op95fact/ opal95_weight_z, opal95_weight_x, opal95_weight_t, &
           opal95_dweight_t, opal95_weight_rho, opal95_dweight_rho

! MHP 05/02 ADDED FOR ITERATION BETWEEN ROTATION AND STRUCTURE CALCULATION
! common/origstart/: both used here. First appearance of this common
! block in the converted sources.
      double precision :: orig_specific_angular_momentum(json), &
           orig_composition(15,json)
      common/origstart/ orig_specific_angular_momentum, orig_composition


!***MHP 1/04 inserted for test
! common/optab/ is commented out in the original and is not declared.
!*** end test
! LLP 8/07 Make LPTIME available for calibration
! common/ct3/: use_structure_dt_limits used here. Naming matches
! htimer.f90.
      logical :: use_structure_dt_limits
      common/ct3/ use_structure_dt_limits


! JVS 02/11
! KC 2025-05-30 reordered common block elements
! common/acdpth/: ageout_placeholder/lclcd_placeholder/
! ljlast_placeholder/ljwrt_placeholder/compute_acoustic_depth/
! acoustic_depth_output are actively used in this file (see the header
! COMMON BLOCK NOTE for why several nonetheless keep their established
! placeholder names); the rest are unused placeholders. Naming matches
! getopal95.f90/envint.f90.
      double precision :: taucz_placeholder, deladj_placeholder(json), &
           tauhe_placeholder, tnorm_placeholder, tcz_placeholder, &
           whe_placeholder
      double precision :: acatmr_placeholder(json), acatmd_placeholder(json), &
           acatmp_placeholder(json), acatmt_placeholder(json), &
           tatmos_placeholder
      double precision :: ageout_placeholder(5)
      logical :: lclcd_placeholder
      integer :: iclcd_placeholder, iacat_placeholder, ijlast_placeholder
      logical :: ljlast_placeholder, ljwrt_placeholder, &
           compute_acoustic_depth, laoly_placeholder
      integer :: ijvs_placeholder, ijent_placeholder, ijdel_placeholder
      logical :: acoustic_depth_output
      common/acdpth/ taucz_placeholder, deladj_placeholder, &
           tauhe_placeholder, tnorm_placeholder, tcz_placeholder, &
           whe_placeholder, acatmr_placeholder, acatmd_placeholder, &
           acatmp_placeholder, acatmt_placeholder, tatmos_placeholder, &
           ageout_placeholder, lclcd_placeholder, iclcd_placeholder, &
           iacat_placeholder, ijlast_placeholder, ljlast_placeholder, &
           ljwrt_placeholder, compute_acoustic_depth, laoly_placeholder, &
           ijvs_placeholder, ijent_placeholder, ijdel_placeholder, &
           acoustic_depth_output
      integer :: nao
      data nao/1/


! JVS 08/13 IF THE CZ IS BEYOND THE FITTING POINT, STORE ITS LOCATION
! common/envcz/: not used in this file's logic; declared only to
! preserve layout. Naming matches wrtout.f90/envint.f90.
      double precision :: convection_zone_radius_placeholder, rint_placeholder
      common/envcz/ convection_zone_radius_placeholder, rint_placeholder

! JVS 04/14 Common block for additional timestep governors
! common/govs/: not used in this file's logic; declared only to
! preserve layout. Naming matches htimer.f90.
      logical :: use_envelope_triangle_dt
      common/govs/ use_envelope_triangle_dt

! MHP 10/24 NEW VARIABLES FOR STOP CRITERIA ON CENTRAL ABUNDANCE are
! carried in common/sett/ above.

      double precision :: be7_mass_fraction_zone(json), &
           neutrino_flux_zone(10,json)
      double precision :: reaction_rate_1(json), reaction_rate_2(json), &
           reaction_rate_3(json), reaction_rate_4(json), &
           reaction_rate_5(json), reaction_rate_6(json), &
           reaction_rate_7(json), reaction_rate_8(json), &
           reaction_rate_9(json), reaction_rate_10(json), &
           reaction_rate_11(json), reaction_rate_12(json), &
           reaction_rate_13(json), n15_alpha_branch_fraction(json), &
           be7_electron_capture_fraction(json)
      integer :: mixed_zone_bounds(12,2), mixed_zone_bounds_no_overshoot(12,2)
      integer :: radiative_zone_bounds(13,2)
      logical :: am_transport_convective_flag(json)
! latest values (Bahcall and Pinsonneault 1996)-actual values set in
! subroutine PARMIN
      double precision :: bp96_scale_factor(17)
      data bp96_scale_factor/0.9558,0.9690,0.9712,1.0,1.0,0.992,1.0,1.0, &
           1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.92088,0.1625/
! MHP 3/96 added data for base solar age, L
!       DATA SUNAGE,SUNL/4.57D09,3.844D33/  ! KC 2025-05-31
      double precision :: reference_solar_luminosity
      data reference_solar_luminosity/3.844D33/

! --- locals ---
      integer :: mc_run_start, mc_run_end, monte_carlo_run_number
      double precision :: age_scale_factor
      logical :: saved_pulse_output_flag
      integer :: convergence_iterations
      double precision :: initial_x_guess, initial_alpha_guess
      logical :: saved_use_structure_dt_limits
      integer :: saved_atm_choice
      integer :: num_zones
      logical :: reset_triangle, model_diverged_flag
      double precision :: total_angular_momentum, total_rotational_ke
      double precision :: log10_luminosity, log_teff
      integer :: ikut_flag, istore_flag, envelope_cz_bottom_index, &
           model_number
      double precision :: delta_time, hydrogen_dt, timestep_yr
      double precision :: total_mass_msun, log_total_mass
      double precision :: convective_velocity
      logical :: compute_neutrino_fluxes
      double precision :: prev_mass_bound, curr_mass_bound, next_mass_bound
      double precision :: dlnrho_dlnt_unused, dlnrho_dlnp_unused
      integer :: i, j, k, ii
      double precision :: shell_log_density, shell_log_temperature, &
           hydrogen_fraction, helium_fraction, metal_fraction, he3_fraction, &
           c12_fraction, c13_fraction, n14_fraction, n15_fraction, &
           o16_fraction, o17_fraction, o18_fraction, deuterium_fraction
      double precision :: pp_chain_energy_gen, he3he4_be7_electron_energy_gen, &
           he3he4_be7_proton_energy_gen, cno_cycle_energy_gen, &
           triple_alpha_energy_gen, dlnepsilon_dlnrho, dlnepsilon_dlnt, &
           total_energy_gen_rate
      integer :: shell_index
      double precision :: t6_million_k, log_electron_density, &
           zone_mass_fraction, zone_radius_fraction
      double precision :: prev_log_l, prev_log_teff, prev_age, path_length_sq
      integer :: wrtlst_unit, model_iteration
      integer :: core_cz_top_index, h_shell_zone_begin, h_shell_midpoint_zone, &
           h_shell_end_index
      logical :: has_h_shell
      double precision :: delta_time_saved
      logical :: evolve_model_flag
      integer :: num_species
      double precision :: target_envelope_mass
      logical :: new_atmosphere_fit_needed
      integer :: envelope_cz_zone_prev, envelope_cz_zone_end
      logical :: converged, in_atmosphere, want_derivatives, mixing_active, &
           conductive_opacity_flag
      integer :: iteration_level, max_iterations, iterations_done
      double precision :: dlnrho_dlnt, dlnrho_dlnp
      logical :: recompute_surface_bc, recompute_envelope_triangle
      logical :: use_correct_gradients
      integer :: num_radiative_zones, num_mixed_zones, &
           num_mixed_zones_no_overshoot
      double precision :: log_gravity
      double precision :: teff_kelvin_unused
      double precision :: log_r_rsun, current_zx, surface_z_over_x
      double precision :: initial_helium_fraction, initial_metal_fraction
      double precision :: central_temperature_mk, central_pressure_scaled, &
           central_density_linear
      logical :: wind_loss_active
      double precision :: max_domega_frac
      integer :: itrot
      double precision :: delta_temp_step, delta_pressure_step, &
           delta_lum_step, delta_radius_step
      character(len=256) :: alex06_table_path, allard_table_path, &
           atm_table_path, fermi_table_path, kurucz_table_path, &
           kurucz_table2_path, laol_table_path, laol_table2_path, &
           opal95_table_path, opal92_table_path
      character(len=256) :: zams_a_table_path, zams_b_table_path, &
           zams_c_table_path, centre1_table_path, centre2_table_path, &
           centre3_table_path, centre4_table_path, centre5_table_path
      character(len=256) :: opal92_table2_path, pure_z_table_path, &
           scv_h_table_path, scv_he_table_path, scv_z_table_path
      character(len=256) :: alex95_table_paths(7)
      character(len=256) :: pulse_atm_path, pulse_env_path, pulse_mod_path
      double precision :: monte_helium_diffusion_fraction
      logical :: punch_pending_flag

      save

!*******
! START
!*******
      call setversion()

      iowr = 9
! LPUNCH is TRUE once first model is calculated
      punch_pending_flag = .false.
! read in user parameters
      call parmin(alex06_table_path,allard_table_path,atm_table_path,fermi_table_path,kurucz_table_path,kurucz_table2_path,laol_table_path, &
           laol_table2_path,opal95_table_path,opal92_table_path,zams_a_table_path,zams_b_table_path,zams_c_table_path,centre1_table_path,centre2_table_path,centre3_table_path,centre4_table_path, &
           centre5_table_path,opal92_table2_path,pulse_atm_path,pulse_env_path,pulse_mod_path,pure_z_table_path,scv_h_table_path,scv_he_table_path,scv_z_table_path,alex95_table_paths)
! set up constants and read in tabular data
! MHP 8/25 directly pass file names instead of using common blocks
      call setups(mixture_weights,alex06_table_path,allard_table_path,atm_table_path,fermi_table_path,kurucz_table_path,kurucz_table2_path, &
           laol_table_path,laol_table2_path,opal95_table_path,opal92_table_path,zams_a_table_path,zams_b_table_path,zams_c_table_path,centre1_table_path,centre2_table_path,centre3_table_path, &
           centre4_table_path,centre5_table_path,opal92_table2_path,pure_z_table_path,scv_h_table_path,scv_he_table_path,scv_z_table_path,alex95_table_paths)
! MHP 3/96 changed I/O to read in only up to max run needed.
      if (lmonte) then
!c MHP 8/25 moved file open to parmin
!     OPEN(UNIT=IDYN,FILE=FDYN,FORM='FORMATTED',STATUS='OLD')
         mc_run_start = imbeg
         imend = min(imend,1000)
         mc_run_end = imend
! read in monte carlo data
         do i = 1,imend
            read(dynamics_unit,1511)s11_rate(i),s33_rate(i),s34_rate(i), &
                 s17_rate(i),metal_to_h_ratio(i),helium_fraction_param(i), &
                 luminosity_target(i),age_target(i)
 1511       format(7X,1P7E10.3/E9.3)
            write(iowr,*)i,s11_rate(i),s33_rate(i),s34_rate(i),s17_rate(i), &
                 metal_to_h_ratio(i),helium_fraction_param(i), &
                 luminosity_target(i),age_target(i)
            diffusion_factor(i) = helium_fraction_param(i)
         end do
      else
         mc_run_start = 1
         mc_run_end = 1
      endif
      do 500 monte_carlo_run_number = mc_run_start,mc_run_end
! for monte carlo run, input values of parameters being changed.
      if (lmonte) then
         cross_section_scale(1) = s11_rate(monte_carlo_run_number)*bp96_scale_factor(1)
         cross_section_scale(2) = s33_rate(monte_carlo_run_number)*bp96_scale_factor(2)
         cross_section_scale(3) = s34_rate(monte_carlo_run_number)*bp96_scale_factor(3)
         cross_section_scale(16) = s17_rate(monte_carlo_run_number)*bp96_scale_factor(16)
         monte_helium_diffusion_fraction = helium_fraction_param(monte_carlo_run_number)
         fgrz = diffusion_factor(monte_carlo_run_number)
         solar_luminosity_cgs = reference_solar_luminosity*luminosity_target(monte_carlo_run_number)
         log10_solar_luminosity = dlog10(solar_luminosity_cgs)
         ln_solar_luminosity = ln10/solar_luminosity_cgs
         age_scale_factor = age_target(monte_carlo_run_number)
! timestep and final age are altered in SR SETCAL; input #s should be
! scaled for a solar age of 4.57 Gyr
         target_end_age(2)=1.0D8
         target_end_age(3)=4.57D9
      else
         age_scale_factor = 1.0D0
      endif
! DBG PULSE: save LPULSE flag, set LPULSE to F except on last model of
! last run, then set LPULSE to saved value of LPULSE.
      saved_pulse_output_flag = pulsation_output_active
! 02/11 JVS uncommented LPULSE=.FALSE.
!      LPULSE = .FALSE.
! MHP 1/93 add option to automatically calibrate solar model.
! MHP 3/96 added counter for # of iterations per converged model and
! starting estimate of ALPHA and X
      if (calibrate_solar_model) then
         call setcal(age_scale_factor)
         convergence_iterations = 1
         initial_x_guess = rescale_params(2,1)
         initial_alpha_guess = mixing_length_array(1)
         saved_use_structure_dt_limits = use_structure_dt_limits   ! save LPTIME for reuse during calibration
         saved_atm_choice  = atm_choice    ! save KTTAU for reuse during calibration
      else
         convergence_iterations = 0
      endif
! DBG 12/94 add option to automatically calculate a stellar model
! of specified Teff and L
      if (calibrate_star_flag) then
         call setscal
      endif

!**********
!     RUN THROUGH THE KIND CARDS IN ORDER
!**********
      do 200 nk = 1, num_runs
         sound_speed_output_active = .false.
!         LPULSE=.FALSE.
         initial_envelope_x = initial_x_array(nk)
         initial_envelope_z = initial_z_array(nk)
         cmixl = mixing_length_array(nk)
       change_envelope_mass_flag = has_senv0_array(nk)
       requested_envelope_mass = senv0_array(nk)
       reset_triangle = .false.
       model_diverged_flag = .false.
! MHP 10/02 ZERO OUT INITIAL ANGULAR MOMENTUM
         total_angular_momentum = 0.0D0
         total_rotational_ke = 0.0D0
! read in the initial model here
! STARIN also calls RSCALE to perform rescaling if requested
!        CALL STARIN(BL,CFENV,DAGE,DDAGE,DELTS,DELTSH,DELTS0,ETA2,  ! KC 2025-05-31
       call starin(log10_luminosity,envelope_fit_coeffs,dage,timestep_yr, &
            delta_time,hydrogen_dt,eta_squared, &
            pressure_rotation_factor,temperature_rotation_factor,trial_sign_flag, &
            composition,log_density,moment_of_inertia,specific_angular_momentum, &
            kinetic_energy_rot,luminosity_lsun,log_pressure,log_radius, &
            log_mass,enclosed_mass,shell_mass,log_total_mass,log_temperature, &
!      *    HSTOT,HT,IKUT,ISTORE,JCORE,JENV,LARGE,LC,LNEW,M,MODEL,  ! KC 2025-05-31
            ikut_flag,istore_flag,envelope_cz_bottom_index,model_diverged_flag, &
            convective_flag,recompute_envelope_triangle,num_zones,model_number, &
            nk,omega,fit_point_pressure,dlnrho_dlnp,dlnrho_dlnt,qiw, &
            mean_radius,fit_point_radius,total_angular_momentum,total_rotational_ke, &
            total_mass_msun,log_teff,luminosity_breakdown,trial_log_luminosity, &
            trial_log_temperature,fit_point_temperature,convective_velocity, &
            mean_gravity,mixture_weights)

      if ((omega(1) .eq. 0) .and. (rotation_active)) then

1611      format('LROT set to TRUE, but OMEGA(1) = 0. Stopping.', &
                 ' Initialize rotation rates or set LROT to', &
                 ' FALSE.')
          print 1611
          stop
      endif
!     MHP 10/24 CHECK STOP CONDITIONS AND DISABLE THEM IF THE STARTING VALUES ARE BELOW THE TARGET THRESHOLD
         if (end_age_stop_active(nk)) then
            if (central_deuterium_stop(nk).gt.0.0D0 .and. &
                 composition(12,1).lt.central_deuterium_stop(nk)) then
               central_deuterium_stop(nk)=-central_deuterium_stop(nk)
               write(*,101)composition(12,1),central_deuterium_stop(nk)
               write(short_file_unit,101)composition(12,1),central_deuterium_stop(nk)
 101           format('STARTING D ',E12.4,' BELOW STOP VALUE ', &
                      E12.4,' STOP DISABLED.')
            endif
            if (central_hydrogen_stop(nk).gt.0.0D0 .and. &
                 composition(1,1).lt.central_hydrogen_stop(nk)) then
               central_deuterium_stop(nk)=-central_hydrogen_stop(nk)
               write(*,102)composition(12,1),central_deuterium_stop(nk)
               write(short_file_unit,102)composition(12,1),central_deuterium_stop(nk)
 102           format('STARTING X ',E12.4,' BELOW STOP VALUE ', &
                      E12.4,' STOP DISABLED.')
            endif
            if (central_helium_stop(nk).gt.0.0D0 .and. &
                 composition(2,1).lt.central_helium_stop(nk)) then
               central_helium_stop(nk)=-central_helium_stop(nk)
               write(*,103)composition(12,1),central_deuterium_stop(nk)
               write(short_file_unit,103)composition(12,1),central_deuterium_stop(nk)
 103           format('STARTING Y ',E12.4,' BELOW STOP VALUE ', &
                      E12.4,' STOP DISABLED.')
            endif
         endif
!     MHP 2/04 NEUTRINO TABLE
!      LNUTAB = .TRUE.
      compute_neutrino_fluxes = .false.
      if (compute_neutrino_fluxes) then
! SET UP WEIGHTS AND MASSES.
! HS1 = LOCATION IN GM (UNLOGGED) OF SHELL CENTERS.
! HS2 = MASS IN GM OF EACH SHELL.
      curr_mass_bound = exp(ln10*log_mass(1))
      prev_mass_bound = - curr_mass_bound
      do i = 2,num_zones
         next_mass_bound = prev_mass_bound
         prev_mass_bound = curr_mass_bound
         curr_mass_bound = exp(ln10*log_mass(i))
         enclosed_mass(i-1) = prev_mass_bound
         shell_mass(i-1) = 0.5D0*(curr_mass_bound-next_mass_bound)
      end do
      enclosed_mass(num_zones) = curr_mass_bound
      shell_mass(num_zones) = exp(ln10*log_total_mass) - 0.5D0*(prev_mass_bound+curr_mass_bound)
      dlnrho_dlnt_unused = -1.0D0
      dlnrho_dlnp_unused = 1.0D0
      do j = 1,10
         flux_diag%neutrino_flux_total(j) = 0.0D0
         do k = 1,num_zones
            neutrino_flux_zone(j,k) = 0.0D0
         end do
      end do
! ASSIGN LOCAL VARIABLES FOR SR CALL FROM GLOBAL VECTORS.
      do i = 1,num_zones
         shell_log_density = log_density(i)
         shell_log_temperature = log_temperature(i)
! SKIP CALCULATIONS FOR LOW TEMPERATURES.
         if (shell_log_temperature.lt.6.0D0) goto 666
         hydrogen_fraction = composition(1,i)
         helium_fraction = composition(2,i)
         metal_fraction = composition(3,i)
         he3_fraction = composition(4,i)
         c12_fraction = composition(5,i)
         c13_fraction = composition(6,i)
         n14_fraction = composition(7,i)
         n15_fraction = composition(8,i)
         o16_fraction = composition(9,i)
         o17_fraction = composition(10,i)
         o18_fraction = composition(11,i)
         call engeb(pp_chain_energy_gen,he3he4_be7_electron_energy_gen, &
              he3he4_be7_proton_energy_gen,cno_cycle_energy_gen, &
              triple_alpha_energy_gen,dlnepsilon_dlnrho,dlnepsilon_dlnt, &
              total_energy_gen_rate,shell_log_density, &
!      *TL,PDT,PDP,X,Y,Z,XHE3,XC12,XC13,XN14,XN15,XO16,XO17,
!      *XO18,XH2,XLI6,XLI7,XBE9,I,HR1,HR2,HR3,HR4,HR5,HR6,HR7,  ! KC 2025-05-31
              shell_log_temperature,hydrogen_fraction,helium_fraction, &
              he3_fraction,c12_fraction,c13_fraction,n14_fraction,o16_fraction, &
              o18_fraction,deuterium_fraction,shell_index,reaction_rate_1, &
              reaction_rate_2,reaction_rate_3,reaction_rate_4,reaction_rate_5, &
              reaction_rate_6,reaction_rate_7,reaction_rate_8,reaction_rate_9, &
              reaction_rate_10,reaction_rate_11,reaction_rate_12, &
              reaction_rate_13,n15_alpha_branch_fraction, &
              be7_electron_capture_fraction)
! BE7 MASS FRACTION.
         be7_mass_fraction_zone(i) = engeb_diag%be7_mass_fraction
! CONVERT FROM ERG/GM/S TO ERG/S FOR EACH SHELL BY MULTIPLYING
! BY THE MASS OF EACH SHELL IN GM (HS2).
         do j = 1,10
            neutrino_flux_zone(j,i) = flux_diag%neutrino_flux(j)*shell_mass(i)
            flux_diag%neutrino_flux_total(j) = flux_diag%neutrino_flux_total(j) + neutrino_flux_zone(j,i)
         end do
         write(*,911)i,shell_mass(i),(neutrino_flux_zone(j,i),j=1,10)
 911     format(I5,1P11E10.3)
      end do
  666 continue
! WRITE OUT TOTAL NEUTRINO FLUXES.
! ***NOTE THAT THESE ARE IN UNITS OF 10**10. ***
      write(76,222)(flux_diag%neutrino_flux_total(i),i=1,10)
! NORMALIZE FLUXES.
      do j = 1,10
         do i = 1,num_zones
            neutrino_flux_zone(j,i) = neutrino_flux_zone(j,i)/flux_diag%neutrino_flux_total(j)
         end do
      end do
      do i = 1,num_zones
! TEMPERATURE IN UNITS OF 10**6 K.
         t6_million_k = exp(ln10*(log_temperature(i)-6.0D0))
         if (t6_million_k.lt.5.0D0) goto 141
! ELECTRON DENSITY.
         log_electron_density = log_density(i)+log10((1.0D0+composition(1,i))/2.0D0)
! MASS FRACTION.
         zone_mass_fraction = shell_mass(i)/1.9891D33
! RADIUS FRACTION.
         zone_radius_fraction = exp(ln10*log_radius(i))/solar_radius_cgs
! FLUXES ARE PRINTED IN THE SAME ORDER AS BAHCALL AND PINSONNEAULT.
         write(76,145)zone_radius_fraction,t6_million_k,log_electron_density, &
         zone_mass_fraction,be7_mass_fraction_zone(i),neutrino_flux_zone(1,i), &
         neutrino_flux_zone(5,i), &
         neutrino_flux_zone(6,i), &
         neutrino_flux_zone(7,i),neutrino_flux_zone(8,i),neutrino_flux_zone(4,i), &
         neutrino_flux_zone(2,i),neutrino_flux_zone(3,i)
  145    format(F9.5,F7.3,F6.3,1P10E10.3)
      end do
  141 continue
  222    format(1P10E10.3)
!         IF(M.GT.1)STOP999
      endif
! save mass in solar units
         pulsation_mass_msun=total_mass_msun
! MHP 08/02 STORE STARTING CZ PROPERTIES
         light_burn%jcz = envelope_cz_bottom_index
         turnover%convective_turnover_timescale = 0.0D0
! write out headers of the appropriate output files
      call wrthead(total_mass_msun)
! DBG PULSE OUT 7/92
! initialize variables for calculating when to dump pulse output
         prev_log_l = log10_luminosity
         prev_log_teff = log_teff
         prev_age = dage
         path_length_sq = 0.0D0

       if (helium_flash_active) then
! timestep cutting requires a model stored in logical unit ILAST
! or it will crash - so copy initial model to unit ILAST
          if (punch_pending_flag) then
             wrtlst_unit = ilast
             call wrtlst(wrtlst_unit,composition,log_density,luminosity_lsun, &
                  log_pressure,log_radius,log_mass,log_temperature,convective_flag, &
                  trial_log_temperature,trial_log_luminosity,fit_point_pressure, &
                  fit_point_temperature,fit_point_radius,envelope_fit_coeffs, &
                  trial_sign_flag,luminosity_breakdown,core_cz_top_index, &
                  envelope_cz_bottom_index,model_number,num_zones, &
                  total_mass_msun,log_teff,log10_luminosity,log_total_mass,dage, &
                  timestep_yr,omega)
          endif
       endif

! locate the hydrogen-burning shell and the boundaries of the central
! and surface convection zones (if applicable).
         call findsh(composition,luminosity_lsun,convective_flag,num_zones, &
              core_cz_top_index,envelope_cz_bottom_index,h_shell_zone_begin, &
              h_shell_end_index,h_shell_midpoint_zone,has_h_shell)
! determine timestep for model
! JVS 04/14 Added Teffl to passed variables
!        CALL HTIMER(DELTS,DELTSH,M,HD,HL,HS1,HS2,HT,LC,HCOMP,JCORE,
!      *               JXMID,TLUMX,DAGE,DDAGE,QDT,QDP,NK,HP,HR,OMEGA,  ! KC 2025-05-31
       call htimer(delta_time,hydrogen_dt,num_zones,log_density,luminosity_lsun, &
            enclosed_mass,shell_mass,log_temperature,composition,core_cz_top_index, &
            h_shell_midpoint_zone,luminosity_breakdown,dage,timestep_yr,nk, &
            log_pressure,log_radius,omega,max_domega_frac,h_shell_zone_begin, &
            log_teff)

       delta_time_saved = delta_time
! zero out entropy terms.
         do 99 i = 1,num_zones
            temperature_entropy_term(i) = 0.0D0
            pressure_entropy_term(i) = 0.0D0
            luminosity_entropy_term(i) = 0.0D0
            radius_entropy_term(i) = 0.0D0
   99    continue

! zero out light element burning rates in the surface CZ.
         if (use_extended_composition) then
            light_burn%log_rate_li6_prev = 0.0D0
            light_burn%log_rate_li7_prev = 0.0D0
            light_burn%log_rate_be9_prev = 0.0D0
         endif

! for a given kind card, evolve NMODLS(NK) times
! if rescaling is being performed, NMODLS(NK) is the number of times
! the new model is being relaxed
       do 100 model_iteration = 1,num_models(nk)
! rewind ISHORT if LRWSH is true (keeps ISHORT small)
          if (lrwsh_placeholder) then
             rewind(short_file_unit)
          endif
! DBG PULSE:  if last model of last run then set LPULSE to LSAVPU
            if (model_iteration.eq.num_models(nk) .and. nk .eq. num_runs) then
                 pulsation_output_active = saved_pulse_output_flag
            end if

! JVS 02/11: Also allow pulse output at particular ages along the way
!
!  If the step is bracketing an age of interest, turn on output. This will
! for the step before and step after the age in AGEOUT. Once the info has
! been printed out, AGEOUT is set to the next age.
!
! Turn on calcad:
      if (acoustic_depth_output) then
            compute_acoustic_depth=.true.
      else
            compute_acoustic_depth = .false.
      endif
! If output has been turned on for a previous step, keep it on for the next
! step, but then turn it off.
      if (acoustic_depth_output) then
            if (ljwrt_placeholder) then
                  print*, 'LJWRT on'
                  pulsation_output_active = saved_pulse_output_flag
                  nao=nao+1
                  lclcd_placeholder =.true.
                  ljlast_placeholder =.false.
                  ljwrt_placeholder=.false.
            else if (.not.ljwrt_placeholder) then
                  lclcd_placeholder=.false.
            endif
! If this is the step before one of the ages of interest, print everything out.
! Also, save model structure.
            if (nao.lt.6) then
                  if (dage+timestep_yr/1.0D9-ageout_placeholder(nao) .le. 0.0D0 .and. &
                  dage+2.0D0*timestep_yr/1.0D9-ageout_placeholder(nao) .ge. 0.0D0 .and. .not. ljwrt_placeholder) then
                        print*, 'AGEOUT reached'
                        pulsation_output_active = saved_pulse_output_flag
                        lclcd_placeholder = .true.
                        ljlast_placeholder = .true.
                        ljwrt_placeholder=.true.
                  endif
            endif
       endif


! JVS end
!
! DBG PULSE:  if endage reached then set LPULSE to LSAVPU
! MHP 10/24 GENERALIZE CHECK
         if (end_age_stop_active(nk).and.target_end_age(nk).gt.0.0D0 .and. &
         (abs(target_end_age(nk)-dage*1.0D9-timestep_yr) .le. 1.0D0)) then
                 pulsation_output_active = saved_pulse_output_flag
! MHP 7/96 compute sound speed for solar model
                 sound_speed_output_active = .true.
            end if


!FD echo LSOUND
!        print*,'MAIN LSOUND = ',LSOUND
!FD end
            if (po_output_enabled) then
! MHP 8/25 changed to add file names as declared variables
             call pdist(prev_log_l,prev_log_teff,prev_age,path_length_sq,log10_luminosity,log_teff,model_iteration,pulse_atm_path, &
             pulse_env_path,pulse_mod_path)
          endif

! STARIN called here for timestep cutting
   15       if (model_diverged_flag) then
!              CALL STARIN(BL,CFENV,DAGE,DDAGE,DELTS,DELTSH,DELTS0,ETA2,  ! KC 2025-05-31
             call starin(log10_luminosity,envelope_fit_coeffs,dage,timestep_yr, &
                  delta_time,hydrogen_dt,eta_squared, &
                  pressure_rotation_factor,temperature_rotation_factor,trial_sign_flag, &
                  composition,log_density,moment_of_inertia,specific_angular_momentum, &
                  kinetic_energy_rot,luminosity_lsun,log_pressure,log_radius, &
                  log_mass,enclosed_mass,shell_mass,log_total_mass,log_temperature, &
!      *          HSTOT,HT,IKUT,ISTORE,JCORE,JENV,LARGE,LC,LNEW,M,MODEL,  ! KC 2025-05-31
                  ikut_flag,istore_flag,envelope_cz_bottom_index,model_diverged_flag, &
                  convective_flag,recompute_envelope_triangle,num_zones,model_number, &
                  nk,omega,fit_point_pressure,dlnrho_dlnp,dlnrho_dlnt,qiw, &
                  mean_radius,fit_point_radius,total_angular_momentum,total_rotational_ke, &
                  total_mass_msun,log_teff,luminosity_breakdown,trial_log_luminosity, &
                  trial_log_temperature,fit_point_temperature,convective_velocity, &
                  mean_gravity,mixture_weights)
             if ((omega(1) .eq. 0) .and. (rotation_active)) then
18               format('LROT set to TRUE, but OMEGA(1) = 0. Stopping.', &
                        ' Initialize rotation rates or set LROT to', &
                        ' FALSE.')
                 print 18
                 stop
             endif
          endif
          punch_pending_flag = .true.

! skip this section if model not to be aged
! MHP 7/98
! need to add logic to permit resacling + time evolution for
! pre-main sequence models
            if (rescale_kind(nk).ne.2 .and. model_number.ge.0) then
               evolve_model_flag = .true.
            else if (model_number.ge.0 .and. log_temperature(1).lt.6.6D0) then
               evolve_model_flag = .true.
            else
               evolve_model_flag = .false.
            endif
            new_atmosphere_fit_needed = .false.
            if (evolve_model_flag) then
! ADD MASS LOSS CALCULATION
               call massloss(log10_luminosity,dage,delta_time,composition,log_density,specific_angular_momentum,log_pressure,log_radius, &
                             log_mass,enclosed_mass,shell_mass,log_total_mass,log_temperature,envelope_cz_bottom_index,recompute_envelope_triangle, &
                             num_zones,omega,total_mass_msun,log_teff,target_envelope_mass,new_atmosphere_fit_needed)
! STORE COMPOSITION MATRIX AT THE BEGINNING OF THE TIMESTEP.
               num_species = 11
               if (use_extended_composition) num_species=15
               do 33 i = 1,num_zones
                  do 32 j = 1,num_species
                     prev_model%old_composition(j,i) = composition(j,i)
   32             continue
   33          continue
               iteration_level=1
               call mix(delta_time,composition,log_density,luminosity_lsun,log_pressure,log_radius,log_mass,enclosed_mass,shell_mass,log_total_mass, &
!      *                     HT,ITLVL,LC,M,QDT,QDP,DDAGE,JCORE,JENV,  ! KC 2025-05-31
                            log_temperature,iteration_level,convective_flag,num_zones,timestep_yr,core_cz_top_index,envelope_cz_bottom_index, &
                            mixed_zone_bounds,mixed_zone_bounds_no_overshoot,log_teff)
             timestep_yr = delta_time/seconds_per_year
             dage = dage + 1.0D-9*timestep_yr
            endif
!***MHP 1/04 OPACITY TEST
!      IDT = 15
!      DO JJJ = 1,4
!         IDD(JJJ) = 5
!      END DO
!      XXX = 0.7D0
!      ZZZ = 0.02D0
!      Do JJJ = 1,2000
!         READ(75,*)TL,DL,XX,ZZ,OO
!         IF(TT.GT.1.0D9)STOP911
!         TL = LOG10(TT)
!         DL = LOG10(DD)
! DBG 12/95 GET OPACITY
!         CALL GETOPAC(DL, TL,XXX,ZZZ, O, OL, QOD, QOT)
!         DIFF = (O-OO)/O
!         RL = DL - 3*TL +18.0D0
!         WRITE(76,1554)TL,RL,XXX,ZZZ,O,OO,DIFF
! 1554    format(4f11.6,3e20.10)
!      END DO
!*** END TEST
! rezone new model, except rezoning not performed for He flash calculations
          if (.not.helium_flash_active) then
             call hpoint(num_zones,log_total_mass,log_mass,enclosed_mass, &
                  shell_mass,log_temperature,log_pressure,log_radius,convective_flag, &
                  luminosity_lsun,log_density,composition,fit_point_pressure, &
                  fit_point_radius,stored_envelope_state,istore_flag,reset_triangle,h_shell_zone_begin, &
                  model_number,has_h_shell,core_cz_top_index,envelope_cz_bottom_index,omega,eta_squared,mean_radius,moment_of_inertia, &
!      *              HJM,HKEROT,SJTOT,SKEROT,BL,DELTS,FP,FT,
!      *              HG,QIW,SMASS,TEFFL)  ! KC 2025-05-31
                  specific_angular_momentum,kinetic_energy_rot,total_angular_momentum,total_rotational_ke, &
                  pressure_rotation_factor,temperature_rotation_factor,mean_gravity,qiw,log_teff)
! STORE NEW CZ BASE
               light_burn%jcz = envelope_cz_bottom_index
            else
! save old model for PTIME
               do i=1, num_zones
                  prev_model%old_pressure(i) = log_pressure(i)
                  prev_model%old_temperature(i) = log_temperature(i)
                  prev_model%old_radius(i) = log_radius(i)
                  prev_model%old_luminosity(i) = luminosity_lsun(i)
                  prev_model%old_density(i) = log_density(i)
               end do
! JVS 04/14 Save Teff as well
               prev_model%old_teff = log_teff
!  JVS 05/25 Added model number to list of saved values
           prev_model%old_num_zones = num_zones

          endif
! store starting distribution of rotational kinetic energy.
            if (rotation_active) then
               do i = 1,num_zones
                  kinetic_energy_rot_old(i) = kinetic_energy_rot(i)
               end do
            endif
! changed for lithium burning with overshoot.
! store starting depth of C.Z. for light element burning.
            if (use_extended_composition) then
               light_burn%cz_base_radius_prev = 0.0D0
               envelope_cz_zone_prev = envelope_cz_bottom_index
               if (envelope_overshoot_active) then
                  light_burn%pressure_scale_height_start = alphae*exp(clndp*(log_pressure(envelope_cz_bottom_index)+2.0D0*log_radius(envelope_cz_bottom_index) &
                           -log_density(envelope_cz_bottom_index)-cgl-log_mass(envelope_cz_bottom_index)))
               else
                  light_burn%pressure_scale_height_start = 0.0D0
               endif
! find burning rates at the beginning of the time step.
               call lirate88(composition,log_density,log_temperature,num_zones,1)
            endif
! begin correction routines
! set flags for CRRECT
! CRRECT checks surface boundary conditions in every iteration
! if LNEW0 = T, new envelope triangle calculated the 1st iteration
! (i.e. old triangle ignored)
! LFINI = T if model has converged
! LARGE = T if model has diverged
          if (lnew0_placeholder) recompute_envelope_triangle = .true.
            if (.not.evolve_model_flag) delta_time = -dabs(delta_time)
            fcorr = dabs(fcorr0) - fcorri
            iterations_done = 0
            model_diverged_flag = .false.
            converged = .false.
            if (.not.lnews .or. delta_time.le.0.0D0) then
               do 20 i = 1,num_zones
! zero entropy terms
                  log_temperature_delta(i) = 0.0D0
                  log_pressure_delta(i) = 0.0D0
                  temperature_entropy_term(i) = 0.0D0
                  pressure_entropy_term(i) = 0.0D0
                  luminosity_entropy_term(i) = 0.0D0
                  radius_entropy_term(i) = 0.0D0
! zero gravitational energy terms.
                  gravitational_luminosity(i) = 0.0D0
 20            continue
            else
! use the rate of change in the previous model to estimate the new
! run of structure variables.
               do 30 i = 1,num_zones
                  delta_temp_step = temperature_entropy_term(i)*delta_time
                  delta_pressure_step = pressure_entropy_term(i)*delta_time
                  delta_lum_step = luminosity_lsun(i)*luminosity_entropy_term(i)*delta_time
                  delta_radius_step = radius_entropy_term(i)*delta_time
                  log_temperature_delta(i) = delta_temp_step
                  log_pressure_delta(i) = delta_pressure_step
                  log_temperature(i) = log_temperature(i) + delta_temp_step
                  log_pressure(i) = log_pressure(i) + delta_pressure_step
                  luminosity_lsun(i) = luminosity_lsun(i) + delta_lum_step
                  log_radius(i) = log_radius(i) + delta_radius_step
! zero gravitational energy terms.
                  gravitational_luminosity(i) = 0.0D0
 30            continue
            endif

! FIRST LEVEL OF ITERATIONS
! USE ENVELOPE TRIANGLE OF THE PREVIOUS MODEL;
! FOR THE FIRST MODEL OF A RUN,THE TRIANGLE IS GENERATED HERE.
            max_iterations = niter1
            recompute_surface_bc = .false.
! CALL TO CRRECT - ADDED ITERATION LEVEL
            iteration_level = 1
            call crrect(delta_time,num_zones,max_iterations,converged,model_diverged_flag,recompute_envelope_triangle,reset_triangle,recompute_surface_bc, &
                 trial_log_temperature,trial_log_luminosity,envelope_fit_coeffs,fit_point_pressure,fit_point_temperature,fit_point_radius,trial_sign_flag,istore_flag,stored_envelope_state,log_total_mass, &
                 log10_luminosity,log_teff,log_density,elim_coeff,elim_rhs,gravitational_luminosity,luminosity_lsun,max_residual,log_pressure,log_pressure_delta,log_radius,log_mass,enclosed_mass, &
                 shell_mass,surface_bc,log_temperature,log_temperature_delta,composition,convective_flag,max_correction_index,luminosity_breakdown,in_atmosphere,want_derivatives, &
!      *           LMIX,LOCOND,QDT,QDP,MODEL,FP,FT,ETA2,OMEGA,R0,ITDONE,
!      *           HG,HI,HJM, ITLVL,LCZ,MRZONE,MXZONE,NRZONE,NZONE,  ! KC 2025-05-31
                 mixing_active,conductive_opacity_flag,dlnrho_dlnt,dlnrho_dlnp,pressure_rotation_factor,temperature_rotation_factor,eta_squared,omega,mean_radius,iterations_done, &
                 mean_gravity,moment_of_inertia,specific_angular_momentum,iteration_level,am_transport_convective_flag,mixed_zone_bounds,qiw,kinetic_energy_rot,kinetic_energy_rot_old)
! SECOND LEVEL OF ITERATIONS
! CHECK ENVELOPE TRIANGLE BEFORE ITERATING FOR SOLUTION
            if (model_diverged_flag) goto 15
            recompute_surface_bc = .true.
            max_iterations = niter2
            iteration_level = 2
            call crrect(delta_time,num_zones,max_iterations,converged,model_diverged_flag,recompute_envelope_triangle,reset_triangle,recompute_surface_bc, &
                 trial_log_temperature,trial_log_luminosity,envelope_fit_coeffs,fit_point_pressure,fit_point_temperature,fit_point_radius,trial_sign_flag,istore_flag,stored_envelope_state,log_total_mass, &
                 log10_luminosity,log_teff,log_density,elim_coeff,elim_rhs,gravitational_luminosity,luminosity_lsun,max_residual,log_pressure,log_pressure_delta,log_radius,log_mass,enclosed_mass, &
                 shell_mass,surface_bc,log_temperature,log_temperature_delta,composition,convective_flag,max_correction_index,luminosity_breakdown,in_atmosphere,want_derivatives, &
!      *           LMIX,LOCOND,QDT,QDP,MODEL,FP,FT,ETA2,OMEGA,R0,ITDONE,
!      *           HG,HI,HJM, ITLVL,LCZ,MRZONE,MXZONE,NRZONE,NZONE,  ! KC 2025-05-31
                 mixing_active,conductive_opacity_flag,dlnrho_dlnt,dlnrho_dlnp,pressure_rotation_factor,temperature_rotation_factor,eta_squared,omega,mean_radius,iterations_done, &
                 mean_gravity,moment_of_inertia,specific_angular_momentum,iteration_level,am_transport_convective_flag,mixed_zone_bounds,qiw,kinetic_energy_rot,kinetic_energy_rot_old)
            if (model_diverged_flag) goto 15
! 7/91 STORE CHANGES IN THE STRUCTURE. THESE CHANGES ARE USED TO GET AN
! IMPROVED FIRST GUESS AT THE STRUCTURE FOR THE NEXT MODEL IF LNEWS=T.
            if (delta_time.gt.0.0D0) then
               do 27 ii = 1,num_zones
                  temperature_entropy_term(ii)=log_temperature_delta(ii)/delta_time
                  pressure_entropy_term(ii)=log_pressure_delta(ii)/delta_time
                  luminosity_entropy_term(ii)=2.0D0*(luminosity_lsun(ii)-prev_model%old_luminosity(ii))/(luminosity_lsun(ii)+prev_model%old_luminosity(ii))/delta_time
                  radius_entropy_term(ii)=(log_radius(ii)-prev_model%old_radius(ii))/delta_time
 27            continue
            endif
! THIRD LEVEL OF ITERATIONS
            recompute_surface_bc = .false.
            max_iterations = niter3
            iteration_level = 3
            call crrect(delta_time,num_zones,max_iterations,converged,model_diverged_flag,recompute_envelope_triangle,reset_triangle,recompute_surface_bc, &
                 trial_log_temperature,trial_log_luminosity,envelope_fit_coeffs,fit_point_pressure,fit_point_temperature,fit_point_radius,trial_sign_flag,istore_flag,stored_envelope_state,log_total_mass, &
                 log10_luminosity,log_teff,log_density,elim_coeff,elim_rhs,gravitational_luminosity,luminosity_lsun,max_residual,log_pressure,log_pressure_delta,log_radius,log_mass,enclosed_mass, &
                 shell_mass,surface_bc,log_temperature,log_temperature_delta,composition,convective_flag,max_correction_index,luminosity_breakdown,in_atmosphere,want_derivatives, &
!      *           LMIX,LOCOND,QDT,QDP,MODEL,FP,FT,ETA2,OMEGA,R0,ITDONE,
!      *           HG,HI,HJM, ITLVL,LCZ,MRZONE,MXZONE,NRZONE,NZONE,  ! KC 2025-05-31
                 mixing_active,conductive_opacity_flag,dlnrho_dlnt,dlnrho_dlnp,pressure_rotation_factor,temperature_rotation_factor,eta_squared,omega,mean_radius,iterations_done, &
                 mean_gravity,moment_of_inertia,specific_angular_momentum,iteration_level,am_transport_convective_flag,mixed_zone_bounds,qiw,kinetic_energy_rot,kinetic_energy_rot_old)
            if (model_diverged_flag) goto 15
            if (.not.rotation_active) then
               itdif1 = 1
            endif
! MHP 05/02
! IF THE CODE IS ITERATING BETWEEN THE STRUCTURE AND ROTATION
! SOLUTIONS, ENSURE THAT THE START-OF-TIMESTEP QUANTITIES
! HCOMPP (COMPOSITION) AND HJMSAV (ANGULAR MOMENTUM) ARE ONLY
! OVERWRITTEN ON THE LAST RUN THROUGH.
            if (itdif1.gt.1) then
               do i = 1,num_zones
                  orig_specific_angular_momentum(i) = specific_angular_momentum(i)
                  do j = 1,15
                     orig_composition(j,i) = prev_model%old_composition(j,i)
                  end do
               end do
            endif
            do itrot = 1, itdif1
! MHP 05/02 RESTORE ORIGINAL "START OF TIMESTEP"
! VALUES FOR THE COMPOSITION MATRIX
               if (itrot.gt.1) then
                  do i = 1,num_zones
                     do j = 1,15
                        prev_model%old_composition(j,i) = orig_composition(j,i)
                     end do
                  end do
               endif
! 7/91 THE FOURTH LEVEL OF ITERATION REPEATS THE ITERATION BETWEEN THE
! MIXING AND THE STRUCTURE VARIABLES.  IT SHOULD NOT BE USED FOR MODELS
! WHERE SEMI-CONVECTION IS IMPORTANT (ITERATING BETWEEN THE BURNING AND
! STRUCTURE GENERATES OSCILLATIONS). IT SHOULD BE USED FOR HIGH-PRECISION
! WORK (E.G. SOLAR MODELS).
! Surface boundary conditions checked again since we've changed the
! composition (and hence the structure) of the model in ITLVL=3
! (to be implemented when I know the rest of it works!)
            max_iterations = niter4
            recompute_surface_bc=.false.
            iteration_level = 4
            call crrect(delta_time,num_zones,max_iterations,converged,model_diverged_flag,recompute_envelope_triangle,reset_triangle,recompute_surface_bc, &
                 trial_log_temperature,trial_log_luminosity,envelope_fit_coeffs,fit_point_pressure,fit_point_temperature,fit_point_radius,trial_sign_flag,istore_flag,stored_envelope_state,log_total_mass, &
                 log10_luminosity,log_teff,log_density,elim_coeff,elim_rhs,gravitational_luminosity,luminosity_lsun,max_residual,log_pressure,log_pressure_delta,log_radius,log_mass,enclosed_mass, &
                 shell_mass,surface_bc,log_temperature,log_temperature_delta,composition,convective_flag,max_correction_index,luminosity_breakdown,in_atmosphere,want_derivatives, &
!      *           LMIX,LOCOND,QDT,QDP,MODEL,FP,FT,ETA2,OMEGA,R0,ITDONE,
!      *           HG,HI,HJM,ITLVL,LCZ,MRZONE,MXZONE,NRZONE,NZONE,  ! KC 2025-05-31
                 mixing_active,conductive_opacity_flag,dlnrho_dlnt,dlnrho_dlnp,pressure_rotation_factor,temperature_rotation_factor,eta_squared,omega,mean_radius,iterations_done, &
                 mean_gravity,moment_of_inertia,specific_angular_momentum,iteration_level,am_transport_convective_flag,mixed_zone_bounds,qiw,kinetic_energy_rot,kinetic_energy_rot_old)
!  25         CONTINUE
            if (.not.converged) then
! MODEL FAILED TO CONVERGE WITHIN(NITER1+NITER2+NITER3+NITER4)ITERATIONS
               model_diverged_flag = .true.
               goto 15
            endif

! MODEL HAS CONVERGED
! ENSURE CONVECTION ZONES ARE FULLY MIXED.
! MHP 02/12 INFER CONVECTIVE OVERTURN TIMESCALE (USED IN MDOT)
! JVS 02/12 CALL MIXCZ(HCOMP,HS2,LC,M)
! KC 2025-05-30 addressed warning messages from Makefile.legacy
! C G Somers 6/14, SET IMIX = .TRUE. SO THE CORRECT GRADS ARE USED.
!       IMIX = .TRUE.
!       CALL MIXCZ(HCOMP,HS2,HS1,LC,HR,HP,HD,HG,M,IMIX)
! G Somers 6/14, SET LIMIX = .TRUE. SO THE CORRECT GRADS ARE USED.
      use_correct_gradients = .true.
!       CALL MIXCZ(HCOMP,HS2,HS1,LC,HR,HP,HD,HG,M,LIMIX)  ! KC 2025-05-31
      call mixcz(composition,shell_mass,convective_flag,num_zones)
! G Somers END

! MHP 9/94 STORE TOTAL AGE IN SAGE
            disk_lifetime = dage
            if (rotation_active) then
! RESTORE ORIGINAL START OF TIMESTEP VALUES
! TO THE ANGULAR MOMENTUM DISTRIBUTION
               if (itrot.gt.1) then
                  do i = 1,num_zones
                     specific_angular_momentum(i) = orig_specific_angular_momentum(i)
                  end do
               endif
! MHP 9/94 ADDED FLAG TO TURN ON ROTATION OUTPUT WHEN END OF KIND
! CARD REACHED.
! MHP 10/24 GENERALIZE CHECK
         if (end_age_stop_active(nk).and.target_end_age(nk).gt.0.0D0 .and. &
         (abs(target_end_age(nk)-dage*1.0D9-timestep_yr) .le. 1.0D0)) then
!               IF(LENDAG(NK).AND.ENDAGE(NK)-DAGE*1.0D9.LE.1.0D0)THEN
                  lprt0_placeholder = .true.
               else
                  lprt0_placeholder = .false.
               endif
! FIND THE NEW RUN OF OMEGA
! JENV0 ADDED TO SR CALL.
               wind_loss_active = ljdot0
               call getw(log10_luminosity,delta_time,max_domega_frac,pressure_rotation_factor,temperature_rotation_factor,composition,log_density,specific_angular_momentum,luminosity_lsun,log_pressure,log_radius, &
!      *              HS,HS1,HS2,HSTOT,HT,LC,LJDOT,M,MODEL,SJTOT,SMASS,  ! KC 2025-05-31
                    log_mass,enclosed_mass,shell_mass,log_total_mass,log_temperature,convective_flag,wind_loss_active,num_zones,total_mass_msun, &
                    log_teff,eta_squared,mean_gravity,moment_of_inertia,omega,qiw,mean_radius,envelope_cz_zone_prev)
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
               call fpft(log_density,log_radius,log_mass,num_zones,omega,eta_squared,pressure_rotation_factor,temperature_rotation_factor,mean_gravity,mean_radius)
            endif
            end do
! LOCATE THE HYDROGEN-BURNING SHELL AND THE BOUNDARIES OF THE CENTRAL
! AND SURFACE CONVECTION ZONES (IF APPLICABLE).
       call findsh(composition,luminosity_lsun,convective_flag,num_zones, &
              core_cz_top_index,envelope_cz_bottom_index,h_shell_zone_begin,h_shell_end_index,h_shell_midpoint_zone, &
              has_h_shell)
! PERFORM LIGHT ELEMENT BURNING
         if (use_extended_composition .and. model_number.ge.0 .and. delta_time.gt.0.0D0) then
! ONLY FOR MODELS WITHOUT ROTATION, OR WITHOUT ROTATIONAL MIXING.
            if (.not.rotation_active .or. .not.instability_transport_active) then
! FIND CONVECTION ZONE DEPTH AT THE END OF THE TIME STEP.
               call convec(composition,log_density,log_pressure,log_radius,log_mass,log_temperature,convective_flag,num_zones,radiative_zone_bounds,mixed_zone_bounds, &
                            mixed_zone_bounds_no_overshoot,core_cz_top_index,envelope_cz_bottom_index,num_radiative_zones,num_mixed_zones,num_mixed_zones_no_overshoot)
! CHANGED FOR LITHIUM BURNING WITH OVERSHOOT.
               envelope_cz_zone_end = envelope_cz_bottom_index
               if (envelope_overshoot_active) then
                  light_burn%pressure_scale_height_end = alphae*exp(clndp*(log_pressure(envelope_cz_bottom_index)+2.0D0*log_radius(envelope_cz_bottom_index) &
                           -log_density(envelope_cz_bottom_index)-cgl-log_mass(envelope_cz_bottom_index)))
               else
                  light_burn%pressure_scale_height_end = 0.0D0
               endif
! FIND BURNING RATES AT THE END OF THE TIME STEP.
               call lirate88(composition,log_density,log_temperature,num_zones,2)
!                CALL LIBURN(DELTS,HCOMP,HD,HR,HS1,HS2,HT,JENV1,JENV0,M)  ! KC 2025-05-31
               call liburn(delta_time,composition,log_radius,enclosed_mass,shell_mass,log_temperature,envelope_cz_zone_end,envelope_cz_zone_prev,num_zones)
            endif
         endif
! MHP 07/02 RESTORE PRIOR FITTING POINT IF MASS ACCRETION IS BEING
! INCLUDED
         if (new_atmosphere_fit_needed) then
            call getnewenv(target_envelope_mass,composition,log_density,luminosity_lsun,log_pressure,log_radius,log_mass,enclosed_mass,shell_mass, &
!     *                     HSTOT,HT,LC,ETA2,HG,HI,HJM,QIW,R0,  ! KC 2025-05-31
                            log_total_mass,log_temperature,convective_flag,eta_squared,moment_of_inertia,specific_angular_momentum,qiw,mean_radius, &
                            kinetic_energy_rot,log10_luminosity,total_angular_momentum,total_rotational_ke,log_teff,num_zones,recompute_envelope_triangle)
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
            call fpft(log_density,log_radius,log_mass,num_zones,omega,eta_squared,pressure_rotation_factor,temperature_rotation_factor,mean_gravity,mean_radius)
            new_atmosphere_fit_needed = .false.
         endif
! DETERMINE TIMESTEP FOR NEXT MODEL
! HTIMER ALSO LOCATES THE H-BURNING SHELL
! JVS 04/14 added teffl to passed htimer variables
       delta_time = dabs(delta_time)
       delta_time_saved = delta_time
!        CALL HTIMER(DELTS,DELTSH,M,HD,HL,HS1,HS2,HT,LC,HCOMP,JCORE,
!      *        JXMID,TLUMX,DAGE,DDAGE,QDT,QDP,NK,HP,HR,OMEGA,  ! KC 2025-05-31
       call htimer(delta_time,hydrogen_dt,num_zones,log_density,luminosity_lsun,enclosed_mass,shell_mass,log_temperature,composition,core_cz_top_index, &
              h_shell_midpoint_zone,luminosity_breakdown,dage,timestep_yr,nk,log_pressure,log_radius,omega, &
              max_domega_frac,h_shell_zone_begin,log_teff)
! IF EVOLVING TO A GIVEN AGE AND KIND CARD IS DONE, AVOID ZEROING OUT
! TIMESTEP WRITTEN TO MODEL (AS THIS MAKES CONTINUING A SEQUENCE AWKWARD.)
!     INSTEAD WRITE THE PREVIOUS MODEL TIMESTEP TO MODEL.
! ONLY IF A FIXED END AGE IS USED, NOT FOR OTHER STOPS
       if (end_age_stop_active(nk) .and. target_end_age(nk).gt.0.0D0) then
          if (target_end_age(nk)-dage*1.0D9.le.1.0D0) then
             delta_time = max(delta_time_saved,1.0D-3*dage*seconds_per_year)
             timestep_yr = delta_time/seconds_per_year
          else
             delta_time_saved = delta_time
          endif
       else
          delta_time_saved = delta_time
       endif
       if (rescale_kind(nk).ne.2) model_number = model_number+1
! WRTOUT IS THE OUTPUT DRIVER ROUTINE
       call wrtout(composition,log_density,luminosity_lsun,log_pressure,log_radius,log_mass,enclosed_mass,log_temperature,convective_flag,num_zones,model_number,dage, &
              timestep_yr,total_mass_msun,log_teff,log10_luminosity,log_gravity,has_h_shell,h_shell_zone_begin,h_shell_midpoint_zone,h_shell_end_index,core_cz_top_index, &
              envelope_cz_bottom_index,luminosity_breakdown,trial_log_temperature,trial_log_luminosity,fit_point_pressure,fit_point_temperature,fit_point_radius,envelope_fit_coeffs,trial_sign_flag,log_total_mass,omega, &
!      *        LPUNCH,FP,FT,ETA2,R0,HJM,HI,SJTOT,SKEROT,HS2,NK)  ! KC 2025-05-31
              punch_pending_flag,pressure_rotation_factor,temperature_rotation_factor,eta_squared,mean_radius,specific_angular_momentum,moment_of_inertia,total_angular_momentum,total_rotational_ke,shell_mass)

! MHP 10/24 GENERALIZED STOP CONDITIONS
!     IF EVOLVING TO A GIVEN AGE AND AGE IS REACHED, KIND CARD IS DONE
!       IF(LENDAG(NK).AND.ENDAGE(NK)-DAGE*1.0D9.LE.1.0D0)GOTO 110
       if (end_age_stop_active(nk).and.target_end_age(nk).gt.0.0D0 .and. &
         (target_end_age(nk)-dage*1.0D9).le.1.0D0) goto 110
! MHP 10/24 CHECK ALL STOP CONDITIONS, EXIT IF ANY SATISFIED
         end_kind_flag = .false.
         if (end_age_stop_active(nk).and.central_deuterium_stop(nk).gt.0.0D0 .and. &
              composition(12,1).lt.central_deuterium_stop(nk)) then
            write(*,104)composition(12,1),central_deuterium_stop(nk)
 104        format('CENTRAL D ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
         if (end_age_stop_active(nk).and.central_hydrogen_stop(nk).gt.0.0D0 .and. &
              composition(1,1).lt.central_hydrogen_stop(nk)) then
            write(*,105)composition(1,1),central_hydrogen_stop(nk)
 105        format('CENTRAL X ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
         if (end_age_stop_active(nk).and.central_helium_stop(nk).gt.0.0D0 .and. &
              composition(2,1).lt.central_helium_stop(nk)) then
            write(*,106)composition(2,1),central_helium_stop(nk)
 106        format('CENTRAL Y ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
! IF EXITING, SET I/O FLAGS PROPERLY AND EXIT LOOP
         if (end_kind_flag) then
            pulsation_output_active = saved_pulse_output_flag
            sound_speed_output_active = .true.
            lprt0_placeholder = .true.
            goto 110
         endif
! TEST IF MODEL IS NEAR DESIRED Teff AND L. IF NOT RESCALE AND TRY AGAIN.
         if (calibrate_star_flag .and. .not. star_found_flag) then
            if (mod(nk,2).eq.0) then
             if (model_iteration.eq.1) then
                teff_kelvin_unused = 10.0D0**log_teff
             else
                call chkscal(log10_luminosity, log_teff, dage, nk)
                if (just_passed_target_radius_flag) goto 200
             end if
          endif
       endif

! END OF RUN
  100    continue

! G Somers 11/14, CHANGE CALL TO PUTSTORE INSTEAD OF WRTLST.
! STORE LAST MODEL IN ISTOR IF LSTORE, LSTPCH, AND LPUNCH ARE .TRUE.
  110    if (lstore.and.lstpch.and.punch_pending_flag) then
          call putstore(composition,log_density,luminosity_lsun,log_pressure,log_radius,log_mass,log_temperature,convective_flag,trial_log_temperature,trial_log_luminosity,fit_point_pressure,fit_point_temperature,fit_point_radius, &
                 envelope_fit_coeffs,trial_sign_flag,luminosity_breakdown,core_cz_top_index,envelope_cz_bottom_index,model_number,num_zones,total_mass_msun,log_teff,log10_luminosity,log_total_mass, &
                 dage,timestep_yr,omega,enclosed_mass,eta_squared,mean_radius,pressure_rotation_factor,temperature_rotation_factor,specific_angular_momentum,moment_of_inertia)
            punch_pending_flag = .false.
       endif
! 110  CONTINUE
! G Somers END

! MHP 1/93 CHECK AUTOMATIC CALIBRATATION OF SOLAR MODEL.
!c MHP 5/96 changed solar calibration to perform solar models in 3 kind cards
         if (calibrate_solar_model) then
! JVS Turn off calcad - speeds things up
            compute_acoustic_depth=.false.
            if (mod(nk,3).eq.0) then
               log_r_rsun = 0.5D0*(log10_luminosity+log10_solar_luminosity-c4pil-csigl-4.0D0*log_teff)-log10_solar_radius
! MHP 06/13 Add solar Z/X to observables
               current_zx = composition(3,num_zones)/composition(1,num_zones)
               call chkcal(log10_luminosity,log_r_rsun,nk,current_zx)
!               CALL CHKCAL(BL,RLL,NK)
               use_structure_dt_limits = saved_use_structure_dt_limits  ! Restore LPTIME to original value for next cycle
               atm_choice  = saved_atm_choice    ! Restore KTTAU to original value for next cycle
               if (solar_calibration_active) then
                  go to 250
               else
!c MHP 8/96 added counter for # of runs needed for calibration
                  convergence_iterations = convergence_iterations + 1
! MHP 6/97 STOP AFTER 10 ATTEMPTS AT CALIBRATION
!                  IF(ICONV.GE.11) GOTO 250
                  if (convergence_iterations.ge.15) goto 250
                  if (pulsation_output_active) then
! DBG 6/93 Need to delete pulse output because have not got ultimate
! model yet.
! MHP 8/25 Replaced delete file with rewind file. This is functionally the same and avoids the need to pass the character string for the file name from parmin.
                     rewind(opal_model_unit)
                     rewind(opal_envelope_unit)
                     rewind(opal_atm_unit)
!                     CLOSE(IOPMOD, STATUS='DELETE')
!                     CLOSE(IOPENV, STATUS='DELETE')
!                     CLOSE(IOPATM, STATUS='DELETE')
!                     OPEN(IOPMOD, FILE=FPMOD,STATUS='UNKNOWN',
!     *                    FORM='FORMATTED')
!                     OPEN(IOPENV, FILE=FPENV,STATUS='UNKNOWN',
!     *                    FORM='FORMATTED')
!                     OPEN(IOPATM, FILE=FPATM,STATUS='UNKNOWN',
!     *                    FORM='FORMATTED')
                  end if
               end if
            endif
         endif

! DBG 12/94 NO MORE RUNS NEEDED. HAVE CALIBRATED STELLAR MODEL
         if (calibrate_star_flag .and. star_found_flag.and.(mod(nk,2).eq.0)) goto 250

! END RUN LOOP
 200  continue
! EXIT RUN LOOP
 250  continue

! FOR MONTE CARLO, REWIND OUTPUT FILES AND WRITE OUT SNU FLUXES AND
! MODEL PARAMTERS TO AN OUTPUT FILE.
! RUN FAILED TO CONVERGE.  WRITE FINAL INFO WITH WARNING NOTE.
      if (lmonte .and. convergence_iterations.ge.11 .and. .not.solar_calibration_active) then
         rewind(ilast)
         rewind(first_unit)
         rewind(idebug)
         rewind(itrack)
         rewind(short_file_unit)
         rewind(imodpt)
         rewind(istor)
         write(neutrino_unit,1525)log10_luminosity,log_r_rsun
 1525    format(5X,'DID NOT CONVERGE WITHIN 10 ATTEMPTS L,R',2F10.6)
! MONTE CARLO #, CONVERGED MIXING LENGTH AND INITIAL H, SURFACE X,
! SURFACE Z, Z/X, CENTRAL X, CENTRAL Z
         write(neutrino_unit,1519) monte_carlo_run_number,mixing_length_array(nk),rescale_params(2,nk-2),composition(1,num_zones), &
              composition(3,num_zones),surface_z_over_x,composition(1,1),composition(3,1)
 1519    format(1X,I5,3F10.6,4E10.3)
! NUMERICAL DATA : #OF RUNS NEEDED FOR A CONVERGED MODEL, INITIAL X
! AND ALPHA, FINAL DL/DX,DR/DX,DL/D ALPHA, DR/D ALPHA
         write(neutrino_unit,1518)convergence_iterations,initial_x_guess,initial_alpha_guess,dlum_dx,drad_dx,dlum_dalpha,drad_dalpha
! SUMMARY OF STRUCTURE : TC, RHOC, PC
         write(neutrino_unit, 1517)central_log10_temperature,central_log10_pressure,central_log10_density, &
              composition(1,1),composition(3,1)
! NEUTRINO FLUXES (SEE ENGEB FOR DETAILS)
         write(neutrino_unit, 1516) flux_diag%cl37_snu_rate,flux_diag%ga71_snu_rate,(flux_diag%neutrino_flux_total(i),i=1,8)
!          CALL WRTMONTE(HCOMP,HD,HL,HP,HR,HS,HT,LC,M,MODEL,DAGE,
!      *        DDAGE,SMASS,TEFFL,BL,GL,LSHELL,JXBEG,JXMID,
!      *        JXEND,JCORE,JENV,TLUMX,TRIT,TRIL,PS,TS,RS,
!      *        CFENV,FTRI,HSTOT,OMEGA,RLL,ICONV,NK,NN)  ! KC 2025-05-31
         call wrtmonte(composition,log_density,luminosity_lsun,log_pressure,log_radius,log_mass,log_temperature,convective_flag,num_zones,dage, &
              timestep_yr,total_mass_msun,log_teff,log10_luminosity, &
              core_cz_top_index,envelope_cz_bottom_index,luminosity_breakdown,trial_log_temperature,trial_log_luminosity,fit_point_pressure,fit_point_temperature,fit_point_radius, &
              envelope_fit_coeffs,trial_sign_flag,log_total_mass,omega,log_r_rsun,convergence_iterations,nk,monte_carlo_run_number)
      else if (calibrate_solar_model .and. lsnu .and. solar_calibration_active) then
         rewind(ilast)
         rewind(first_unit)
         rewind(idebug)
         rewind(itrack)
         rewind(short_file_unit)
         rewind(imodpt)
         rewind(istor)

         surface_z_over_x = composition(3,num_zones)/composition(1,num_zones)
! HEADER FILE:  MONTE CARLO PARAMETERS
         if (lmonte) then
            write(neutrino_unit,1520)monte_carlo_run_number,s11_rate(monte_carlo_run_number),s33_rate(monte_carlo_run_number),s34_rate(monte_carlo_run_number),s17_rate(monte_carlo_run_number), &
                 metal_to_h_ratio(monte_carlo_run_number),helium_fraction_param(monte_carlo_run_number),diffusion_factor(monte_carlo_run_number),luminosity_target(monte_carlo_run_number),age_target(monte_carlo_run_number)
         endif
 1520    format(I7,1P9E10.3)
! NUMERICAL DATA : #OF RUNS NEEDED FOR A CONVERGED MODEL, INITIAL X
! AND ALPHA, FINAL DL/DX,DR/DX,DL/D ALPHA, DR/D ALPHA
         write(neutrino_unit,1518)convergence_iterations,initial_x_guess,initial_alpha_guess,dlum_dx,drad_dx,dlum_dalpha,drad_dalpha
 1518    format(1X,I2,2F10.6,1P4E11.4)
! NEUTRINO FLUXES (SEE ENGEB FOR DETAILS)
         write(neutrino_unit, 1516) flux_diag%cl37_snu_rate,flux_diag%ga71_snu_rate,(flux_diag%neutrino_flux_total(i),i=1,10)
 1516    format(1X,2F8.3,1P10E10.3)
! SUMMARY OF STRUCTURE : TC, RHOC, PC, XC, ZC (ADD MU C)
         central_temperature_mk = 10.0D0**(central_log10_temperature-6.0D0)
         central_pressure_scaled = 10.0D0**(central_log10_pressure-17.0D0)
         central_density_linear = 10.0D0**central_log10_density
         write(neutrino_unit, 1517)central_temperature_mk,central_density_linear,central_pressure_scaled,composition(1,1),composition(3,1)
 1517    format(1X,F7.3,F7.2,F6.3,2F8.5)
! INITIAL ALPHA,Y,Z,ALPHA; FINAL R, L
         initial_helium_fraction = 1.0D0 - rescale_params(2,nk-2) - rescale_params(3,nk-2)
         initial_metal_fraction = rescale_params(3,nk-2)
         write(neutrino_unit,1521)mixing_length_array(nk),initial_helium_fraction,initial_metal_fraction,log10_luminosity,log_r_rsun
 1521    format(F7.4,2F8.5,1P2E10.3)
! CZ DEPTH (R,M), SURFACE Y, Z, Z/X (ADD T CZ BASE, RHO CZ BASE)
         write(neutrino_unit,1522)envelope_radius,envelope_mass,composition(2,num_zones),composition(3,num_zones),surface_z_over_x
 1522    format(F8.5,F9.6,2F8.5,F9.6)
! ENERGY GENERATION FRACTIONS PP I,II,III,CNO,EGRAV
         write(neutrino_unit,1523)(luminosity_breakdown(j),j=1,4),luminosity_breakdown(7)
 1523    format(1P5E10.3)
         if (lmonte) then
!             CALL WRTMONTE(HCOMP,HD,HL,HP,HR,HS,HT,LC,M,MODEL,DAGE,
!      *           DDAGE,SMASS,TEFFL,BL,GL,LSHELL,JXBEG,JXMID,
!      *           JXEND,JCORE,JENV,TLUMX,TRIT,TRIL,PS,TS,RS,
!      *           CFENV,FTRI,HSTOT,OMEGA,RLL,ICONV,NK,NN)  ! KC 2025-05-31
            call wrtmonte(composition,log_density,luminosity_lsun,log_pressure,log_radius,log_mass,log_temperature,convective_flag,num_zones,dage, &
                 timestep_yr,total_mass_msun,log_teff,log10_luminosity, &
                 core_cz_top_index,envelope_cz_bottom_index,luminosity_breakdown,trial_log_temperature,trial_log_luminosity,fit_point_pressure,fit_point_temperature,fit_point_radius, &
                 envelope_fit_coeffs,trial_sign_flag,log_total_mass,omega,log_r_rsun,convergence_iterations,nk,monte_carlo_run_number)
         endif
      endif
 500  end do

      stop
end program main

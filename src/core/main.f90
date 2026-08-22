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
!   - HL (linear L/Lsun, per-shell array) is named star%luminosity_lsun,
!     matching crrect.f90/starin.f90; findsh.f90/htimer.f90/mix.f90/
!     hpoint.f90/getnewenv.f90/wrtlst.f90/wrtout.f90/putstore.f90/
!     massloss.f90/getw.f90 instead call the same array
!     "log_luminosity" (or plain "luminosity"/"star%luminosity_lsun" in
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
!     ljwrt_placeholder; common/rotprt/'s run_diag%lprt0_placeholder;
!     common/chrone/'s lrwsh_placeholder; common/cenv/'s
!     lnew0). Per the precedent set by atm_lib.f90 (keeps
!     lclcd_placeholder despite noting its own active use) and
!     getw.f90 (keeps run_diag%lprt0_placeholder despite noting its own active
!     use), these established names are reused verbatim here too
!     rather than renamed, even though this file actively assigns/
!     reads them; each active use is called out in a comment at its
!     point of use below.
program main

! the array size, i.e. max # of shells is specified in the parameter
! statement. it defines JSON. to change the array size do a global
! change on "JSON=2000" or whatever.
      use nuclear_lib
      use run_diag_lib
      use fluxes_lib
      use engeb_diag_lib
      use light_burn_lib
      use turnover_lib
      use oldmod_lib
      use luout_lib
      use const_lib
      use star_info_lib, only: star
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


































      double precision :: trial_sign_flag
      double precision :: mixture_weights(12)
!     MHP 10/24 FLAG FOR END OF RUN
      logical :: end_kind_flag























      integer :: nao
      data nao/1/




! MHP 10/24 NEW VARIABLES FOR STOP CRITERIA ON CENTRAL ABUNDANCE are
! carried in common/sett/ above.

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
      ! 2026 (ROADMAP.md stage 3): library errors return here via ierr;
      ! this driver-side call site preserves the historical stop.
      integer :: jerr

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
            read(dynamics_unit,1511)run_diag%s11_rate(i),run_diag%s33_rate(i),run_diag%s34_rate(i), &
                 run_diag%s17_rate(i),run_diag%metal_to_h_ratio(i),run_diag%helium_fraction_param(i), &
                 run_diag%luminosity_target(i),run_diag%age_target(i)
 1511       format(7X,1P7E10.3/E9.3)
            write(iowr,*)i,run_diag%s11_rate(i),run_diag%s33_rate(i),run_diag%s34_rate(i),run_diag%s17_rate(i), &
                 run_diag%metal_to_h_ratio(i),run_diag%helium_fraction_param(i), &
                 run_diag%luminosity_target(i),run_diag%age_target(i)
            run_diag%diffusion_factor(i) = run_diag%helium_fraction_param(i)
         end do
      else
         mc_run_start = 1
         mc_run_end = 1
      endif
      do 500 monte_carlo_run_number = mc_run_start,mc_run_end
! for monte carlo run, input values of parameters being changed.
      if (lmonte) then
         cross_section_scale(1) = run_diag%s11_rate(monte_carlo_run_number)*bp96_scale_factor(1)
         cross_section_scale(2) = run_diag%s33_rate(monte_carlo_run_number)*bp96_scale_factor(2)
         cross_section_scale(3) = run_diag%s34_rate(monte_carlo_run_number)*bp96_scale_factor(3)
         cross_section_scale(16) = run_diag%s17_rate(monte_carlo_run_number)*bp96_scale_factor(16)
         monte_helium_diffusion_fraction = run_diag%helium_fraction_param(monte_carlo_run_number)
         fgrz = run_diag%diffusion_factor(monte_carlo_run_number)
         solar_luminosity_cgs = reference_solar_luminosity*run_diag%luminosity_target(monte_carlo_run_number)
         log10_solar_luminosity = dlog10(solar_luminosity_cgs)
         ln_solar_luminosity = ln10/solar_luminosity_cgs
         age_scale_factor = run_diag%age_target(monte_carlo_run_number)
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
         run_diag%sound_speed_output_active = .false.
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
       call starin(log10_luminosity, run_diag%dage, timestep_yr, &
            delta_time, hydrogen_dt, trial_sign_flag, log_total_mass, &
            ikut_flag, istore_flag, envelope_cz_bottom_index, &
            model_diverged_flag, recompute_envelope_triangle, num_zones, &
            model_number, nk, dlnrho_dlnp, dlnrho_dlnt, &
            total_angular_momentum, total_rotational_ke, total_mass_msun, &
            log_teff, convective_velocity, mixture_weights)

      if ((star%omega(1) .eq. 0) .and. (rotation_active)) then

1611      format('LROT set to TRUE, but OMEGA(1) = 0. Stopping.', &
                 ' Initialize rotation rates or set LROT to', &
                 ' FALSE.')
          print 1611
          stop
      endif
!     MHP 10/24 CHECK STOP CONDITIONS AND DISABLE THEM IF THE STARTING VALUES ARE BELOW THE TARGET THRESHOLD
         if (end_age_stop_active(nk)) then
            if (central_deuterium_stop(nk).gt.0.0D0 .and. &
                 star%composition(12,1).lt.central_deuterium_stop(nk)) then
               central_deuterium_stop(nk)=-central_deuterium_stop(nk)
               write(*,101)star%composition(12,1),central_deuterium_stop(nk)
               write(short_file_unit,101)star%composition(12,1),central_deuterium_stop(nk)
 101           format('STARTING D ',E12.4,' BELOW STOP VALUE ', &
                      E12.4,' STOP DISABLED.')
            endif
            if (central_hydrogen_stop(nk).gt.0.0D0 .and. &
                 star%composition(1,1).lt.central_hydrogen_stop(nk)) then
               central_deuterium_stop(nk)=-central_hydrogen_stop(nk)
               write(*,102)star%composition(12,1),central_deuterium_stop(nk)
               write(short_file_unit,102)star%composition(12,1),central_deuterium_stop(nk)
 102           format('STARTING X ',E12.4,' BELOW STOP VALUE ', &
                      E12.4,' STOP DISABLED.')
            endif
            if (central_helium_stop(nk).gt.0.0D0 .and. &
                 star%composition(2,1).lt.central_helium_stop(nk)) then
               central_helium_stop(nk)=-central_helium_stop(nk)
               write(*,103)star%composition(12,1),central_deuterium_stop(nk)
               write(short_file_unit,103)star%composition(12,1),central_deuterium_stop(nk)
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
      curr_mass_bound = exp(ln10*star%log_mass(1))
      prev_mass_bound = - curr_mass_bound
      do i = 2,num_zones
         next_mass_bound = prev_mass_bound
         prev_mass_bound = curr_mass_bound
         curr_mass_bound = exp(ln10*star%log_mass(i))
         star%enclosed_mass(i-1) = prev_mass_bound
         star%shell_mass(i-1) = 0.5D0*(curr_mass_bound-next_mass_bound)
      end do
      star%enclosed_mass(num_zones) = curr_mass_bound
      star%shell_mass(num_zones) = exp(ln10*log_total_mass) - 0.5D0*(prev_mass_bound+curr_mass_bound)
      dlnrho_dlnt_unused = -1.0D0
      dlnrho_dlnp_unused = 1.0D0
      do j = 1,10
         flux_diag%neutrino_flux_total(j) = 0.0D0
         do k = 1,num_zones
            star%neutrino_flux_zone(j,k) = 0.0D0
         end do
      end do
! ASSIGN LOCAL VARIABLES FOR SR CALL FROM GLOBAL VECTORS.
      do i = 1,num_zones
         shell_log_density = star%log_density(i)
         shell_log_temperature = star%log_temperature(i)
! SKIP CALCULATIONS FOR LOW TEMPERATURES.
         if (shell_log_temperature.lt.6.0D0) goto 666
         hydrogen_fraction = star%composition(1,i)
         helium_fraction = star%composition(2,i)
         metal_fraction = star%composition(3,i)
         he3_fraction = star%composition(4,i)
         c12_fraction = star%composition(5,i)
         c13_fraction = star%composition(6,i)
         n14_fraction = star%composition(7,i)
         n15_fraction = star%composition(8,i)
         o16_fraction = star%composition(9,i)
         o17_fraction = star%composition(10,i)
         o18_fraction = star%composition(11,i)
         call engeb(pp_chain_energy_gen,he3he4_be7_electron_energy_gen, &
              he3he4_be7_proton_energy_gen,cno_cycle_energy_gen, &
              triple_alpha_energy_gen,dlnepsilon_dlnrho,dlnepsilon_dlnt, &
              total_energy_gen_rate,shell_log_density, &
!      *TL,PDT,PDP,X,Y,Z,XHE3,XC12,XC13,XN14,XN15,XO16,XO17,
!      *XO18,XH2,XLI6,XLI7,XBE9,I,HR1,HR2,HR3,HR4,HR5,HR6,HR7,  ! KC 2025-05-31
              shell_log_temperature,hydrogen_fraction,helium_fraction, &
              he3_fraction,c12_fraction,c13_fraction,n14_fraction,o16_fraction, &
              o18_fraction,deuterium_fraction,shell_index,star%reaction_rate_1, &
              star%reaction_rate_2,star%reaction_rate_3,star%reaction_rate_4,star%reaction_rate_5, &
              star%reaction_rate_6,star%reaction_rate_7,star%reaction_rate_8,star%reaction_rate_9, &
              star%reaction_rate_10,star%reaction_rate_11,star%reaction_rate_12, &
              star%reaction_rate_13,star%n15_alpha_branch_fraction, &
              star%be7_electron_capture_fraction)
! BE7 MASS FRACTION.
         star%be7_mass_fraction_zone(i) = engeb_diag%be7_mass_fraction
! CONVERT FROM ERG/GM/S TO ERG/S FOR EACH SHELL BY MULTIPLYING
! BY THE MASS OF EACH SHELL IN GM (HS2).
         do j = 1,10
            star%neutrino_flux_zone(j,i) = flux_diag%neutrino_flux(j)*star%shell_mass(i)
            flux_diag%neutrino_flux_total(j) = flux_diag%neutrino_flux_total(j) + star%neutrino_flux_zone(j,i)
         end do
         write(*,911)i,star%shell_mass(i),(star%neutrino_flux_zone(j,i),j=1,10)
 911     format(I5,1P11E10.3)
      end do
  666 continue
! WRITE OUT TOTAL NEUTRINO FLUXES.
! ***NOTE THAT THESE ARE IN UNITS OF 10**10. ***
      write(76,222)(flux_diag%neutrino_flux_total(i),i=1,10)
! NORMALIZE FLUXES.
      do j = 1,10
         do i = 1,num_zones
            star%neutrino_flux_zone(j,i) = star%neutrino_flux_zone(j,i)/flux_diag%neutrino_flux_total(j)
         end do
      end do
      do i = 1,num_zones
! TEMPERATURE IN UNITS OF 10**6 K.
         t6_million_k = exp(ln10*(star%log_temperature(i)-6.0D0))
         if (t6_million_k.lt.5.0D0) goto 141
! ELECTRON DENSITY.
         log_electron_density = star%log_density(i)+log10((1.0D0+star%composition(1,i))/2.0D0)
! MASS FRACTION.
         zone_mass_fraction = star%shell_mass(i)/1.9891D33
! RADIUS FRACTION.
         zone_radius_fraction = exp(ln10*star%log_radius(i))/solar_radius_cgs
! FLUXES ARE PRINTED IN THE SAME ORDER AS BAHCALL AND PINSONNEAULT.
         write(76,145)zone_radius_fraction,t6_million_k,log_electron_density, &
         zone_mass_fraction,star%be7_mass_fraction_zone(i),star%neutrino_flux_zone(1,i), &
         star%neutrino_flux_zone(5,i), &
         star%neutrino_flux_zone(6,i), &
         star%neutrino_flux_zone(7,i),star%neutrino_flux_zone(8,i),star%neutrino_flux_zone(4,i), &
         star%neutrino_flux_zone(2,i),star%neutrino_flux_zone(3,i)
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
         prev_age = run_diag%dage
         path_length_sq = 0.0D0

       if (helium_flash_active) then
! timestep cutting requires a model stored in logical unit ILAST
! or it will crash - so copy initial model to unit ILAST
          if (punch_pending_flag) then
             wrtlst_unit = ilast
             call wrtlst(wrtlst_unit,star%composition,star%log_density,star%luminosity_lsun, &
                  star%log_pressure,star%log_radius,star%log_mass,star%log_temperature,star%convective_flag, &
                  star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure, &
                  star%fit_point_temperature,star%fit_point_radius,star%envelope_fit_coeffs, &
                  trial_sign_flag,star%luminosity_breakdown,core_cz_top_index, &
                  envelope_cz_bottom_index,model_number,num_zones, &
                  total_mass_msun,log_teff,log10_luminosity,log_total_mass,run_diag%dage, &
                  timestep_yr,star%omega)
          endif
       endif

! locate the hydrogen-burning shell and the boundaries of the central
! and surface convection zones (if applicable).
         call findsh(star%composition,star%luminosity_lsun,star%convective_flag,num_zones, &
              core_cz_top_index,envelope_cz_bottom_index,h_shell_zone_begin, &
              h_shell_end_index,h_shell_midpoint_zone,has_h_shell)
! determine timestep for model
! JVS 04/14 Added Teffl to passed variables
!        CALL HTIMER(DELTS,DELTSH,M,HD,HL,HS1,HS2,HT,LC,HCOMP,JCORE,
!      *               JXMID,TLUMX,DAGE,DDAGE,QDT,QDP,NK,HP,HR,OMEGA,  ! KC 2025-05-31
       call htimer(delta_time,hydrogen_dt,num_zones,star%log_density,star%luminosity_lsun, &
            star%enclosed_mass,star%shell_mass,star%log_temperature,star%composition,core_cz_top_index, &
            h_shell_midpoint_zone,star%luminosity_breakdown,run_diag%dage,timestep_yr,nk, &
            star%log_pressure,star%log_radius,star%omega,max_domega_frac,h_shell_zone_begin, &
            log_teff)

       delta_time_saved = delta_time
! zero out entropy terms.
         do 99 i = 1,num_zones
            run_diag%temperature_entropy_term(i) = 0.0D0
            run_diag%pressure_entropy_term(i) = 0.0D0
            run_diag%luminosity_entropy_term(i) = 0.0D0
            run_diag%radius_entropy_term(i) = 0.0D0
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
                  if (run_diag%dage+timestep_yr/1.0D9-ageout_placeholder(nao) .le. 0.0D0 .and. &
                  run_diag%dage+2.0D0*timestep_yr/1.0D9-ageout_placeholder(nao) .ge. 0.0D0 .and. .not. ljwrt_placeholder) then
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
         (abs(target_end_age(nk)-run_diag%dage*1.0D9-timestep_yr) .le. 1.0D0)) then
                 pulsation_output_active = saved_pulse_output_flag
! MHP 7/96 compute sound speed for solar model
                 run_diag%sound_speed_output_active = .true.
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
             call starin(log10_luminosity, run_diag%dage, timestep_yr, &
                  delta_time, hydrogen_dt, trial_sign_flag, &
                  log_total_mass, ikut_flag, istore_flag, &
                  envelope_cz_bottom_index, model_diverged_flag, &
                  recompute_envelope_triangle, num_zones, model_number, &
                  nk, dlnrho_dlnp, dlnrho_dlnt, total_angular_momentum, &
                  total_rotational_ke, total_mass_msun, log_teff, &
                  convective_velocity, mixture_weights)
             if ((star%omega(1) .eq. 0) .and. (rotation_active)) then
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
            else if (model_number.ge.0 .and. star%log_temperature(1).lt.6.6D0) then
               evolve_model_flag = .true.
            else
               evolve_model_flag = .false.
            endif
            new_atmosphere_fit_needed = .false.
            if (evolve_model_flag) then
! ADD MASS LOSS CALCULATION
               call massloss(log10_luminosity,run_diag%dage,delta_time,star%composition,star%log_density,star%specific_angular_momentum,star%log_pressure,star%log_radius, &
                             star%log_mass,star%enclosed_mass,star%shell_mass,log_total_mass,star%log_temperature,envelope_cz_bottom_index,recompute_envelope_triangle, &
                             num_zones,star%omega,total_mass_msun,log_teff,target_envelope_mass,new_atmosphere_fit_needed)
! STORE COMPOSITION MATRIX AT THE BEGINNING OF THE TIMESTEP.
               num_species = 11
               if (use_extended_composition) num_species=15
               do 33 i = 1,num_zones
                  do 32 j = 1,num_species
                     prev_model%old_composition(j,i) = star%composition(j,i)
   32             continue
   33          continue
               iteration_level=1
! mixed_zone_bounds_no_overshoot stays an ARGUMENT of mix (not read as
! star% inside it) because crrect passes its own local array there --
! storage deliberately separate from main's. main passes the star copy
! explicitly.
               call mix(delta_time, log_total_mass, iteration_level, &
                    num_zones, timestep_yr, core_cz_top_index, &
                    envelope_cz_bottom_index, &
                    star%mixed_zone_bounds_no_overshoot, log_teff, jerr)
               if (jerr /= 0) stop
             timestep_yr = delta_time/seconds_per_year
             run_diag%dage = run_diag%dage + 1.0D-9*timestep_yr
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
             call hpoint(num_zones, log_total_mass, istore_flag, &
                  reset_triangle, h_shell_zone_begin, model_number, &
                  has_h_shell, core_cz_top_index, &
                  envelope_cz_bottom_index, total_angular_momentum, &
                  total_rotational_ke, log_teff)
! STORE NEW CZ BASE
               light_burn%jcz = envelope_cz_bottom_index
            else
! save old model for PTIME
               do i=1, num_zones
                  prev_model%old_pressure(i) = star%log_pressure(i)
                  prev_model%old_temperature(i) = star%log_temperature(i)
                  prev_model%old_radius(i) = star%log_radius(i)
                  prev_model%old_luminosity(i) = star%luminosity_lsun(i)
                  prev_model%old_density(i) = star%log_density(i)
               end do
! JVS 04/14 Save Teff as well
               prev_model%old_teff = log_teff
!  JVS 05/25 Added model number to list of saved values
           prev_model%old_num_zones = num_zones

          endif
! store starting distribution of rotational kinetic energy.
            if (rotation_active) then
               do i = 1,num_zones
                  star%kinetic_energy_rot_old(i) = star%kinetic_energy_rot(i)
               end do
            endif
! changed for lithium burning with overshoot.
! store starting depth of C.Z. for light element burning.
            if (use_extended_composition) then
               light_burn%cz_base_radius_prev = 0.0D0
               envelope_cz_zone_prev = envelope_cz_bottom_index
               if (envelope_overshoot_active) then
                  light_burn%pressure_scale_height_start = alphae*exp(clndp*(star%log_pressure(envelope_cz_bottom_index)+2.0D0*star%log_radius(envelope_cz_bottom_index) &
                           -star%log_density(envelope_cz_bottom_index)-cgl-star%log_mass(envelope_cz_bottom_index)))
               else
                  light_burn%pressure_scale_height_start = 0.0D0
               endif
! find burning rates at the beginning of the time step.
               call lirate88(star%composition,star%log_density,star%log_temperature,num_zones,1)
            endif
! begin correction routines
! set flags for CRRECT
! CRRECT checks surface boundary conditions in every iteration
! if LNEW0 = T, new envelope triangle calculated the 1st iteration
! (i.e. old triangle ignored)
! LFINI = T if model has converged
! LARGE = T if model has diverged
          if (lnew0) recompute_envelope_triangle = .true.
            if (.not.evolve_model_flag) delta_time = -dabs(delta_time)
            fcorr = dabs(fcorr0) - fcorri
            iterations_done = 0
            model_diverged_flag = .false.
            converged = .false.
            if (.not.lnews .or. delta_time.le.0.0D0) then
               do 20 i = 1,num_zones
! zero entropy terms
                  star%log_temperature_delta(i) = 0.0D0
                  star%log_pressure_delta(i) = 0.0D0
                  run_diag%temperature_entropy_term(i) = 0.0D0
                  run_diag%pressure_entropy_term(i) = 0.0D0
                  run_diag%luminosity_entropy_term(i) = 0.0D0
                  run_diag%radius_entropy_term(i) = 0.0D0
! zero gravitational energy terms.
                  star%gravitational_luminosity(i) = 0.0D0
 20            continue
            else
! use the rate of change in the previous model to estimate the new
! run of structure variables.
               do 30 i = 1,num_zones
                  delta_temp_step = run_diag%temperature_entropy_term(i)*delta_time
                  delta_pressure_step = run_diag%pressure_entropy_term(i)*delta_time
                  delta_lum_step = star%luminosity_lsun(i)*run_diag%luminosity_entropy_term(i)*delta_time
                  delta_radius_step = run_diag%radius_entropy_term(i)*delta_time
                  star%log_temperature_delta(i) = delta_temp_step
                  star%log_pressure_delta(i) = delta_pressure_step
                  star%log_temperature(i) = star%log_temperature(i) + delta_temp_step
                  star%log_pressure(i) = star%log_pressure(i) + delta_pressure_step
                  star%luminosity_lsun(i) = star%luminosity_lsun(i) + delta_lum_step
                  star%log_radius(i) = star%log_radius(i) + delta_radius_step
! zero gravitational energy terms.
                  star%gravitational_luminosity(i) = 0.0D0
 30            continue
            endif

! FIRST LEVEL OF ITERATIONS
! USE ENVELOPE TRIANGLE OF THE PREVIOUS MODEL;
! FOR THE FIRST MODEL OF A RUN,THE TRIANGLE IS GENERATED HERE.
            max_iterations = niter1
            recompute_surface_bc = .false.
! CALL TO CRRECT - ADDED ITERATION LEVEL
            iteration_level = 1
            call crrect(delta_time, num_zones, max_iterations, converged, &
                 model_diverged_flag, recompute_envelope_triangle, &
                 reset_triangle, recompute_surface_bc, trial_sign_flag, &
                 istore_flag, log_total_mass, log10_luminosity, log_teff, &
                 in_atmosphere, want_derivatives, mixing_active, &
                 conductive_opacity_flag, dlnrho_dlnt, dlnrho_dlnp, &
                 iterations_done, iteration_level)
! SECOND LEVEL OF ITERATIONS
! CHECK ENVELOPE TRIANGLE BEFORE ITERATING FOR SOLUTION
            if (model_diverged_flag) goto 15
            recompute_surface_bc = .true.
            max_iterations = niter2
            iteration_level = 2
            call crrect(delta_time, num_zones, max_iterations, converged, &
                 model_diverged_flag, recompute_envelope_triangle, &
                 reset_triangle, recompute_surface_bc, trial_sign_flag, &
                 istore_flag, log_total_mass, log10_luminosity, log_teff, &
                 in_atmosphere, want_derivatives, mixing_active, &
                 conductive_opacity_flag, dlnrho_dlnt, dlnrho_dlnp, &
                 iterations_done, iteration_level)
            if (model_diverged_flag) goto 15
! 7/91 STORE CHANGES IN THE STRUCTURE. THESE CHANGES ARE USED TO GET AN
! IMPROVED FIRST GUESS AT THE STRUCTURE FOR THE NEXT MODEL IF LNEWS=T.
            if (delta_time.gt.0.0D0) then
               do 27 ii = 1,num_zones
                  run_diag%temperature_entropy_term(ii)=star%log_temperature_delta(ii)/delta_time
                  run_diag%pressure_entropy_term(ii)=star%log_pressure_delta(ii)/delta_time
                  run_diag%luminosity_entropy_term(ii)=2.0D0*(star%luminosity_lsun(ii)-prev_model%old_luminosity(ii))/(star%luminosity_lsun(ii)+prev_model%old_luminosity(ii))/delta_time
                  run_diag%radius_entropy_term(ii)=(star%log_radius(ii)-prev_model%old_radius(ii))/delta_time
 27            continue
            endif
! THIRD LEVEL OF ITERATIONS
            recompute_surface_bc = .false.
            max_iterations = niter3
            iteration_level = 3
            call crrect(delta_time, num_zones, max_iterations, converged, &
                 model_diverged_flag, recompute_envelope_triangle, &
                 reset_triangle, recompute_surface_bc, trial_sign_flag, &
                 istore_flag, log_total_mass, log10_luminosity, log_teff, &
                 in_atmosphere, want_derivatives, mixing_active, &
                 conductive_opacity_flag, dlnrho_dlnt, dlnrho_dlnp, &
                 iterations_done, iteration_level)
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
                  run_diag%orig_specific_angular_momentum(i) = star%specific_angular_momentum(i)
                  do j = 1,15
                     run_diag%orig_composition(j,i) = prev_model%old_composition(j,i)
                  end do
               end do
            endif
            do itrot = 1, itdif1
! MHP 05/02 RESTORE ORIGINAL "START OF TIMESTEP"
! VALUES FOR THE COMPOSITION MATRIX
               if (itrot.gt.1) then
                  do i = 1,num_zones
                     do j = 1,15
                        prev_model%old_composition(j,i) = run_diag%orig_composition(j,i)
                     end do
                  end do
               endif
! 7/91 THE FOURTH LEVEL OF ITERATION REPEATS THE ITERATION BETWEEN THE
! MIXING AND THE STRUCTURE VARIABLES.  IT SHOULD NOT BE USED FOR MODELS
! WHERE SEMI-CONVECTION IS IMPORTANT (ITERATING BETWEEN THE BURNING AND
! STRUCTURE GENERATES OSCILLATIONS). IT SHOULD BE USED FOR HIGH-PRECISION
! WORK (E.G. SOLAR MODELS).
! Surface boundary conditions checked again since we've changed the
! star%composition (and hence the structure) of the model in ITLVL=3
! (to be implemented when I know the rest of it works!)
            max_iterations = niter4
            recompute_surface_bc=.false.
            iteration_level = 4
            call crrect(delta_time, num_zones, max_iterations, converged, &
                 model_diverged_flag, recompute_envelope_triangle, &
                 reset_triangle, recompute_surface_bc, trial_sign_flag, &
                 istore_flag, log_total_mass, log10_luminosity, log_teff, &
                 in_atmosphere, want_derivatives, mixing_active, &
                 conductive_opacity_flag, dlnrho_dlnt, dlnrho_dlnp, &
                 iterations_done, iteration_level)
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
      call mixcz(star%composition,star%shell_mass,star%convective_flag,num_zones)
! G Somers END

! MHP 9/94 STORE TOTAL AGE IN SAGE
            disk_lifetime = run_diag%dage
            if (rotation_active) then
! RESTORE ORIGINAL START OF TIMESTEP VALUES
! TO THE ANGULAR MOMENTUM DISTRIBUTION
               if (itrot.gt.1) then
                  do i = 1,num_zones
                     star%specific_angular_momentum(i) = run_diag%orig_specific_angular_momentum(i)
                  end do
               endif
! MHP 9/94 ADDED FLAG TO TURN ON ROTATION OUTPUT WHEN END OF KIND
! CARD REACHED.
! MHP 10/24 GENERALIZE CHECK
         if (end_age_stop_active(nk).and.target_end_age(nk).gt.0.0D0 .and. &
         (abs(target_end_age(nk)-run_diag%dage*1.0D9-timestep_yr) .le. 1.0D0)) then
!               IF(LENDAG(NK).AND.ENDAGE(NK)-DAGE*1.0D9.LE.1.0D0)THEN
                  run_diag%lprt0_placeholder = .true.
               else
                  run_diag%lprt0_placeholder = .false.
               endif
! FIND THE NEW RUN OF OMEGA
! JENV0 ADDED TO SR CALL.
               wind_loss_active = ljdot0
               call getw(log10_luminosity, delta_time, max_domega_frac, &
                    log_total_mass, wind_loss_active, num_zones, &
                    total_mass_msun, log_teff, envelope_cz_zone_prev, &
                    jerr)
               if (jerr /= 0) stop
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
               call fpft(star%log_density,star%log_radius,star%log_mass,num_zones,star%omega,star%eta_squared,star%pressure_rotation_factor,star%temperature_rotation_factor,star%mean_gravity,star%mean_radius)
            endif
            end do
! LOCATE THE HYDROGEN-BURNING SHELL AND THE BOUNDARIES OF THE CENTRAL
! AND SURFACE CONVECTION ZONES (IF APPLICABLE).
       call findsh(star%composition,star%luminosity_lsun,star%convective_flag,num_zones, &
              core_cz_top_index,envelope_cz_bottom_index,h_shell_zone_begin,h_shell_end_index,h_shell_midpoint_zone, &
              has_h_shell)
! PERFORM LIGHT ELEMENT BURNING
         if (use_extended_composition .and. model_number.ge.0 .and. delta_time.gt.0.0D0) then
! ONLY FOR MODELS WITHOUT ROTATION, OR WITHOUT ROTATIONAL MIXING.
            if (.not.rotation_active .or. .not.instability_transport_active) then
! FIND CONVECTION ZONE DEPTH AT THE END OF THE TIME STEP.
               call convec(star%composition,star%log_density,star%log_pressure,star%log_radius,star%log_mass,star%log_temperature,star%convective_flag,num_zones,star%radiative_zone_bounds,star%mixed_zone_bounds, &
                            star%mixed_zone_bounds_no_overshoot,core_cz_top_index,envelope_cz_bottom_index,num_radiative_zones,num_mixed_zones,num_mixed_zones_no_overshoot)
! CHANGED FOR LITHIUM BURNING WITH OVERSHOOT.
               envelope_cz_zone_end = envelope_cz_bottom_index
               if (envelope_overshoot_active) then
                  light_burn%pressure_scale_height_end = alphae*exp(clndp*(star%log_pressure(envelope_cz_bottom_index)+2.0D0*star%log_radius(envelope_cz_bottom_index) &
                           -star%log_density(envelope_cz_bottom_index)-cgl-star%log_mass(envelope_cz_bottom_index)))
               else
                  light_burn%pressure_scale_height_end = 0.0D0
               endif
! FIND BURNING RATES AT THE END OF THE TIME STEP.
               call lirate88(star%composition,star%log_density,star%log_temperature,num_zones,2)
!                CALL LIBURN(DELTS,HCOMP,HD,HR,HS1,HS2,HT,JENV1,JENV0,M)  ! KC 2025-05-31
               call liburn(delta_time,star%composition,star%log_radius,star%enclosed_mass,star%shell_mass,star%log_temperature,envelope_cz_zone_end,envelope_cz_zone_prev,num_zones)
            endif
         endif
! MHP 07/02 RESTORE PRIOR FITTING POINT IF MASS ACCRETION IS BEING
! INCLUDED
         if (new_atmosphere_fit_needed) then
            call getnewenv(target_envelope_mass,star%composition,star%log_density,star%luminosity_lsun,star%log_pressure,star%log_radius,star%log_mass,star%enclosed_mass,star%shell_mass, &
!     *                     HSTOT,HT,LC,ETA2,HG,HI,HJM,QIW,R0,  ! KC 2025-05-31
                            log_total_mass,star%log_temperature,star%convective_flag,star%eta_squared,star%moment_of_inertia,star%specific_angular_momentum,star%qiw,star%mean_radius, &
                            star%kinetic_energy_rot,log10_luminosity,total_angular_momentum,total_rotational_ke,log_teff,num_zones,recompute_envelope_triangle)
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
            call fpft(star%log_density,star%log_radius,star%log_mass,num_zones,star%omega,star%eta_squared,star%pressure_rotation_factor,star%temperature_rotation_factor,star%mean_gravity,star%mean_radius)
            new_atmosphere_fit_needed = .false.
         endif
! DETERMINE TIMESTEP FOR NEXT MODEL
! HTIMER ALSO LOCATES THE H-BURNING SHELL
! JVS 04/14 added teffl to passed htimer variables
       delta_time = dabs(delta_time)
       delta_time_saved = delta_time
!        CALL HTIMER(DELTS,DELTSH,M,HD,HL,HS1,HS2,HT,LC,HCOMP,JCORE,
!      *        JXMID,TLUMX,DAGE,DDAGE,QDT,QDP,NK,HP,HR,OMEGA,  ! KC 2025-05-31
       call htimer(delta_time,hydrogen_dt,num_zones,star%log_density,star%luminosity_lsun,star%enclosed_mass,star%shell_mass,star%log_temperature,star%composition,core_cz_top_index, &
              h_shell_midpoint_zone,star%luminosity_breakdown,run_diag%dage,timestep_yr,nk,star%log_pressure,star%log_radius,star%omega, &
              max_domega_frac,h_shell_zone_begin,log_teff)
! IF EVOLVING TO A GIVEN AGE AND KIND CARD IS DONE, AVOID ZEROING OUT
! TIMESTEP WRITTEN TO MODEL (AS THIS MAKES CONTINUING A SEQUENCE AWKWARD.)
!     INSTEAD WRITE THE PREVIOUS MODEL TIMESTEP TO MODEL.
! ONLY IF A FIXED END AGE IS USED, NOT FOR OTHER STOPS
       if (end_age_stop_active(nk) .and. target_end_age(nk).gt.0.0D0) then
          if (target_end_age(nk)-run_diag%dage*1.0D9.le.1.0D0) then
             delta_time = max(delta_time_saved,1.0D-3*run_diag%dage*seconds_per_year)
             timestep_yr = delta_time/seconds_per_year
          else
             delta_time_saved = delta_time
          endif
       else
          delta_time_saved = delta_time
       endif
       if (rescale_kind(nk).ne.2) model_number = model_number+1
! WRTOUT IS THE OUTPUT DRIVER ROUTINE
       call wrtout(num_zones, model_number, run_diag%dage, timestep_yr, &
            total_mass_msun, log_teff, log10_luminosity, log_gravity, &
            has_h_shell, h_shell_zone_begin, h_shell_midpoint_zone, &
            h_shell_end_index, core_cz_top_index, &
            envelope_cz_bottom_index, trial_sign_flag, log_total_mass, &
            punch_pending_flag, total_angular_momentum, &
            total_rotational_ke)

! MHP 10/24 GENERALIZED STOP CONDITIONS
!     IF EVOLVING TO A GIVEN AGE AND AGE IS REACHED, KIND CARD IS DONE
!       IF(LENDAG(NK).AND.ENDAGE(NK)-DAGE*1.0D9.LE.1.0D0)GOTO 110
       if (end_age_stop_active(nk).and.target_end_age(nk).gt.0.0D0 .and. &
         (target_end_age(nk)-run_diag%dage*1.0D9).le.1.0D0) goto 110
! MHP 10/24 CHECK ALL STOP CONDITIONS, EXIT IF ANY SATISFIED
         end_kind_flag = .false.
         if (end_age_stop_active(nk).and.central_deuterium_stop(nk).gt.0.0D0 .and. &
              star%composition(12,1).lt.central_deuterium_stop(nk)) then
            write(*,104)star%composition(12,1),central_deuterium_stop(nk)
 104        format('CENTRAL D ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
         if (end_age_stop_active(nk).and.central_hydrogen_stop(nk).gt.0.0D0 .and. &
              star%composition(1,1).lt.central_hydrogen_stop(nk)) then
            write(*,105)star%composition(1,1),central_hydrogen_stop(nk)
 105        format('CENTRAL X ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
         if (end_age_stop_active(nk).and.central_helium_stop(nk).gt.0.0D0 .and. &
              star%composition(2,1).lt.central_helium_stop(nk)) then
            write(*,106)star%composition(2,1),central_helium_stop(nk)
 106        format('CENTRAL Y ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
! IF EXITING, SET I/O FLAGS PROPERLY AND EXIT LOOP
         if (end_kind_flag) then
            pulsation_output_active = saved_pulse_output_flag
            run_diag%sound_speed_output_active = .true.
            run_diag%lprt0_placeholder = .true.
            goto 110
         endif
! TEST IF MODEL IS NEAR DESIRED Teff AND L. IF NOT RESCALE AND TRY AGAIN.
         if (calibrate_star_flag .and. .not. star_found_flag) then
            if (mod(nk,2).eq.0) then
             if (model_iteration.eq.1) then
                teff_kelvin_unused = 10.0D0**log_teff
             else
                call chkscal(log10_luminosity, log_teff, run_diag%dage, nk)
                if (just_passed_target_radius_flag) goto 200
             end if
          endif
       endif

! END OF RUN
  100    continue

! G Somers 11/14, CHANGE CALL TO PUTSTORE INSTEAD OF WRTLST.
! STORE LAST MODEL IN ISTOR IF LSTORE, LSTPCH, AND LPUNCH ARE .TRUE.
  110    if (lstore.and.lstpch.and.punch_pending_flag) then
          call putstore(star%composition,star%log_density,star%luminosity_lsun,star%log_pressure,star%log_radius,star%log_mass,star%log_temperature,star%convective_flag,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
                 star%envelope_fit_coeffs,trial_sign_flag,star%luminosity_breakdown,core_cz_top_index,envelope_cz_bottom_index,model_number,num_zones,total_mass_msun,log_teff,log10_luminosity,log_total_mass, &
                 run_diag%dage,timestep_yr,star%omega,star%enclosed_mass,star%eta_squared,star%mean_radius,star%pressure_rotation_factor,star%temperature_rotation_factor,star%specific_angular_momentum,star%moment_of_inertia)
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
               current_zx = star%composition(3,num_zones)/star%composition(1,num_zones)
               call chkcal(log10_luminosity,log_r_rsun,nk,current_zx)
!               CALL CHKCAL(BL,RLL,NK)
               use_structure_dt_limits = saved_use_structure_dt_limits  ! Restore LPTIME to original value for next cycle
               atm_choice  = saved_atm_choice    ! Restore KTTAU to original value for next cycle
               if (run_diag%solar_calibration_active) then
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
      if (lmonte .and. convergence_iterations.ge.11 .and. .not.run_diag%solar_calibration_active) then
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
         write(neutrino_unit,1519) monte_carlo_run_number,mixing_length_array(nk),rescale_params(2,nk-2),star%composition(1,num_zones), &
              star%composition(3,num_zones),surface_z_over_x,star%composition(1,1),star%composition(3,1)
 1519    format(1X,I5,3F10.6,4E10.3)
! NUMERICAL DATA : #OF RUNS NEEDED FOR A CONVERGED MODEL, INITIAL X
! AND ALPHA, FINAL DL/DX,DR/DX,DL/D ALPHA, DR/D ALPHA
         write(neutrino_unit,1518)convergence_iterations,initial_x_guess,initial_alpha_guess,run_diag%dlum_dx,run_diag%drad_dx,run_diag%dlum_dalpha,run_diag%drad_dalpha
! SUMMARY OF STRUCTURE : TC, RHOC, PC
         write(neutrino_unit, 1517)run_diag%central_log10_temperature,run_diag%central_log10_pressure,run_diag%central_log10_density, &
              star%composition(1,1),star%composition(3,1)
! NEUTRINO FLUXES (SEE ENGEB FOR DETAILS)
         write(neutrino_unit, 1516) flux_diag%cl37_snu_rate,flux_diag%ga71_snu_rate,(flux_diag%neutrino_flux_total(i),i=1,8)
!          CALL WRTMONTE(HCOMP,HD,HL,HP,HR,HS,HT,LC,M,MODEL,DAGE,
!      *        DDAGE,SMASS,TEFFL,BL,GL,LSHELL,JXBEG,JXMID,
!      *        JXEND,JCORE,JENV,TLUMX,TRIT,TRIL,PS,TS,RS,
!      *        CFENV,FTRI,HSTOT,OMEGA,RLL,ICONV,NK,NN)  ! KC 2025-05-31
         call wrtmonte(star%composition,star%log_density,star%luminosity_lsun,star%log_pressure,star%log_radius,star%log_mass,star%log_temperature,star%convective_flag,num_zones,run_diag%dage, &
              timestep_yr,total_mass_msun,log_teff,log10_luminosity, &
              core_cz_top_index,envelope_cz_bottom_index,star%luminosity_breakdown,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
              star%envelope_fit_coeffs,trial_sign_flag,log_total_mass,star%omega,log_r_rsun,convergence_iterations,nk,monte_carlo_run_number)
      else if (calibrate_solar_model .and. lsnu .and. run_diag%solar_calibration_active) then
         rewind(ilast)
         rewind(first_unit)
         rewind(idebug)
         rewind(itrack)
         rewind(short_file_unit)
         rewind(imodpt)
         rewind(istor)

         surface_z_over_x = star%composition(3,num_zones)/star%composition(1,num_zones)
! HEADER FILE:  MONTE CARLO PARAMETERS
         if (lmonte) then
            write(neutrino_unit,1520)monte_carlo_run_number,run_diag%s11_rate(monte_carlo_run_number),run_diag%s33_rate(monte_carlo_run_number),run_diag%s34_rate(monte_carlo_run_number),run_diag%s17_rate(monte_carlo_run_number), &
                 run_diag%metal_to_h_ratio(monte_carlo_run_number),run_diag%helium_fraction_param(monte_carlo_run_number),run_diag%diffusion_factor(monte_carlo_run_number),run_diag%luminosity_target(monte_carlo_run_number),run_diag%age_target(monte_carlo_run_number)
         endif
 1520    format(I7,1P9E10.3)
! NUMERICAL DATA : #OF RUNS NEEDED FOR A CONVERGED MODEL, INITIAL X
! AND ALPHA, FINAL DL/DX,DR/DX,DL/D ALPHA, DR/D ALPHA
         write(neutrino_unit,1518)convergence_iterations,initial_x_guess,initial_alpha_guess,run_diag%dlum_dx,run_diag%drad_dx,run_diag%dlum_dalpha,run_diag%drad_dalpha
 1518    format(1X,I2,2F10.6,1P4E11.4)
! NEUTRINO FLUXES (SEE ENGEB FOR DETAILS)
         write(neutrino_unit, 1516) flux_diag%cl37_snu_rate,flux_diag%ga71_snu_rate,(flux_diag%neutrino_flux_total(i),i=1,10)
 1516    format(1X,2F8.3,1P10E10.3)
! SUMMARY OF STRUCTURE : TC, RHOC, PC, XC, ZC (ADD MU C)
         central_temperature_mk = 10.0D0**(run_diag%central_log10_temperature-6.0D0)
         central_pressure_scaled = 10.0D0**(run_diag%central_log10_pressure-17.0D0)
         central_density_linear = 10.0D0**run_diag%central_log10_density
         write(neutrino_unit, 1517)central_temperature_mk,central_density_linear,central_pressure_scaled,star%composition(1,1),star%composition(3,1)
 1517    format(1X,F7.3,F7.2,F6.3,2F8.5)
! INITIAL ALPHA,Y,Z,ALPHA; FINAL R, L
         initial_helium_fraction = 1.0D0 - rescale_params(2,nk-2) - rescale_params(3,nk-2)
         initial_metal_fraction = rescale_params(3,nk-2)
         write(neutrino_unit,1521)mixing_length_array(nk),initial_helium_fraction,initial_metal_fraction,log10_luminosity,log_r_rsun
 1521    format(F7.4,2F8.5,1P2E10.3)
! CZ DEPTH (R,M), SURFACE Y, Z, Z/X (ADD T CZ BASE, RHO CZ BASE)
         write(neutrino_unit,1522)run_diag%envelope_radius,run_diag%envelope_mass,star%composition(2,num_zones),star%composition(3,num_zones),surface_z_over_x
 1522    format(F8.5,F9.6,2F8.5,F9.6)
! ENERGY GENERATION FRACTIONS PP I,II,III,CNO,EGRAV
         write(neutrino_unit,1523)(star%luminosity_breakdown(j),j=1,4),star%luminosity_breakdown(7)
 1523    format(1P5E10.3)
         if (lmonte) then
!             CALL WRTMONTE(HCOMP,HD,HL,HP,HR,HS,HT,LC,M,MODEL,DAGE,
!      *           DDAGE,SMASS,TEFFL,BL,GL,LSHELL,JXBEG,JXMID,
!      *           JXEND,JCORE,JENV,TLUMX,TRIT,TRIL,PS,TS,RS,
!      *           CFENV,FTRI,HSTOT,OMEGA,RLL,ICONV,NK,NN)  ! KC 2025-05-31
            call wrtmonte(star%composition,star%log_density,star%luminosity_lsun,star%log_pressure,star%log_radius,star%log_mass,star%log_temperature,star%convective_flag,num_zones,run_diag%dage, &
                 timestep_yr,total_mass_msun,log_teff,log10_luminosity, &
                 core_cz_top_index,envelope_cz_bottom_index,star%luminosity_breakdown,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
                 star%envelope_fit_coeffs,trial_sign_flag,log_total_mass,star%omega,log_r_rsun,convergence_iterations,nk,monte_carlo_run_number)
         endif
      endif
 500  end do

      stop
end program main

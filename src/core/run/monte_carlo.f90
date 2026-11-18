!----------------------------------------------------------------------
! monte_carlo
!----------------------------------------------------------------------
! The built-in Monte-Carlo runs. Each run takes sampled nuclear
! cross-section factors, a metal settling factor, and solar luminosity
! and age factors from the sample file (control%dynamics_file; the
! fields of star%mc_samples, state/README.md):
!
!  * first_mc_run, last_mc_run -- the run range, functions of the
!    controls (run_yrec's loop bounds), which io's check_controls keeps
!    within 1 to max_mc_runs;
!  * setup_monte_carlo_runs -- each run's samples from the sample file
!    (once, from prepare_run; a short or unreadable file is an error);
!  * begin_monte_carlo_run -- the start of every Monte-Carlo run, the
!    per-run hook: every run copy of a control re-armed from its control,
!    the other state found so far that a run leaves and the next would
!    read cleared (clear_previous_run_state), then the run's sampled
!    values (apply_monte_carlo_parameters); a value the run changes per
!    Monte-Carlo run is set here;
!  * write_run_summaries -- the end-of-run bookkeeping, including the
!    Monte-Carlo model summary (write_monte_carlo_model, legacy
!    wrtmonte), which reads star% directly.
!
! The files' layouts, reads, writes and rewinds are io's
! (io/monte_carlo_files.f90); this module holds no I/O statement.
! The Monte-Carlo loop itself is in run_yrec. The deck
! examples/run_config_matrix/cm_monte_carlo.inlist makes one Monte-Carlo
! run, with the sample file mc_input.dat beside it.
module monte_carlo
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: i_h1, i_he4, i_metals, i_lum_grav, star_info
      use control_slots, only: i_target_x, i_target_z
      use core_messages, only: core_report, core_terminal_values, &
           message_length, unwritten_line
      use run_files, only: mc_model_unit
      use model_files, only: write_mod_model
      use monte_carlo_files, only: num_mc_sample_values, read_monte_carlo_sample, &
           monte_carlo_summary, write_monte_carlo_summary, rewind_run_cycle_files
      use calibration, only: dlogl_dx_fixed, dlogr_dx_fixed, dlogl_dalpha_fixed, &
           dlogr_dalpha_fixed
      implicit none
      private
      public :: first_mc_run, last_mc_run, setup_monte_carlo_runs, &
           begin_monte_carlo_run, write_run_summaries

contains

! ---------------------------------------------------------------
! The run range: the first and the last Monte-Carlo run, the controls
! mc_run_start and mc_run_end (legacy IMBEG, IMEND), which io's
! check_controls requires to satisfy 1 <= mc_run_start <= mc_run_end <=
! max_mc_runs when Monte Carlo is on (M115; the F77 clamped the last to
! 1000 without a word); without Monte Carlo the single run 1.
!   star reads: physics%mc_run_start, physics%monte_carlo_active
!   star writes: none
pure integer function first_mc_run(s)
      type(star_info), intent(in) :: s
      first_mc_run = 1
      if (s%physics%monte_carlo_active) first_mc_run = s%physics%mc_run_start
end function first_mc_run

!   star reads: physics%mc_run_end, physics%monte_carlo_active
!   star writes: none
pure integer function last_mc_run(s)
      type(star_info), intent(in) :: s
      last_mc_run = 1
      if (s%physics%monte_carlo_active) last_mc_run = s%physics%mc_run_end
end function last_mc_run

! ---------------------------------------------------------------
! For a Monte-Carlo job, each run's samples from the sample file (io's
! read_monte_carlo_sample, one two-line record per run): the S11, S33,
! S34 and S17 factors, Z/X, the metal settling factor, and the
! luminosity and age factors. Only the records up to last_mc_run (at
! most max_mc_runs) are read (MHP 3/96). A file that ends before them,
! or a record that does not read, is an error (ierr = 1, reported on the
! terminal and in the run log; M115).
!   star reads: control%dynamics_file, physics%monte_carlo_active
!   star writes: mc_samples
subroutine setup_monte_carlo_runs(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      integer :: i, ios
      real(dp) :: samples(num_mc_sample_values)
      character(len=256) :: msg
      character(len=message_length) :: lines(1)
      ierr = 0
      if (.not. s%physics%monte_carlo_active) return
      do i = 1,last_mc_run(s)
         call read_monte_carlo_sample(samples, ios, msg)
         if (ios /= 0) then
            lines = unwritten_line
            write(lines,'(a,i0,3a)') ' ERROR: Monte Carlo: cannot read run ', i, &
                 '''s samples from the dynamics file ', &
                 trim(s%control%dynamics_file), ': '//trim(msg)//'; run stopped'
            call core_report(lines)
            ierr = 1
            return
         end if
         s%mc_samples%s11_factor(i) = samples(1)
         s%mc_samples%s33_factor(i) = samples(2)
         s%mc_samples%s34_factor(i) = samples(3)
         s%mc_samples%s17_factor(i) = samples(4)
         s%mc_samples%metal_to_h_ratio(i) = samples(5)
         s%mc_samples%settling_z_factor(i) = samples(6)
         s%mc_samples%luminosity_factor(i) = samples(7)
         s%mc_samples%age_factor(i) = samples(8)
         call core_terminal_values(samples, int_value=i)
      end do
end subroutine setup_monte_carlo_runs

! ---------------------------------------------------------------
! The start of every Monte-Carlo run (and of the one run without Monte
! Carlo): the per-run hook, the first call of run_yrec's Monte-Carlo
! loop. Each Monte-Carlo run is the deck's run again with other samples,
! so it starts from the deck's values: every run copy of a control is
! re-armed from its control (controls_sync's arm_run_copies, as the
! controls read armed them), whatever the run before changed -- the
! stops init_stop_conditions disarmed, the solar calibration's
! age-scaled fixed timesteps and laid-out cycles, the structure-timestep
! and mass-change latches, the atmosphere the surface boundary switched
! to, the core extension's point count. The other state found so far that a run
! leaves on the star and the next would read before writing it is
! cleared (clear_previous_run_state). Then the run's samples are applied
! (apply_monte_carlo_parameters). On the first run the re-arm and the
! clearing change nothing: no run copy is written between the controls
! read and the loop, and the cleared state is still the
! template star's. A value the run changes per Monte-Carlo run is a
! star member of the Monte-Carlo runs class, set here from the run's
! sample (as settling_z_factor_in_use); state/README.md, "Add a value
! the run changes per Monte-Carlo run". The contract below covers this
! body alone: what the run starts from is its three callees' contracts.
!   star reads: none
!   star writes: none
subroutine begin_monte_carlo_run(s, monte_carlo_run_number, &
     age_scale_factor)
      use controls_sync, only: arm_run_copies
      type(star_info), intent(inout) :: s
      integer, intent(in) :: monte_carlo_run_number
      real(dp), intent(out) :: age_scale_factor
      call arm_run_copies(s)
      call clear_previous_run_state(s)
      call apply_monte_carlo_parameters(s, monte_carlo_run_number, &
           age_scale_factor)
end subroutine begin_monte_carlo_run

! ---------------------------------------------------------------
! The members found so far that a run leaves on the star, other than
! the run copies, and that the next run would read before writing them,
! put back to their value at the start of the first run (the template
! star's zero). The list is what a bisection of the star's classes on
! four two-run decks found (M169), not one shown to be complete
! (disk_gate_age_gyr, the disk-locking tests' age, would carry over too:
! every model read from file sets it, starting_model's
! reset_state_after_read; M190). The members: the
! start-of-step deuterium rates. compute_burn_rates stores them only in
! the zones above the burning cutoff (and while there is deuterium to
! burn), and burn_deuterium averages them over a mixed region, the
! zones beyond the cutoff included, so the next run's first steps would
! take the previous run's rates there (docs/bug_sweep.md M169); and the
! settling routines' report flag, so that a run that ended with settling
! suspended does not make the next one print "SETTLING RESUMED" for a
! suspension it never reported.
! test_monte_carlo.py runs two Monte-Carlo runs of identical samples,
! which must be the same run: state that carries over fails it.
!   star reads: none
!   star writes: deuterium_burning_rate_start,
!     settling_suspended_reported
subroutine clear_previous_run_state(s)
      type(star_info), intent(inout) :: s
      s%deuterium_burning_rate_start = 0.0d0
      s%settling_suspended_reported = .false.
end subroutine clear_previous_run_state

! ---------------------------------------------------------------
! For a Monte-Carlo run, apply the current run's sampled values: the
! factors on the S-factors of p+p, He3+He3, He3+He4 and Be7+p (the
! cross-section scales r_pp, r_he3he3, r_he3he4 and the extra scale
! xs_be7p), each times the factor that takes the code's reference S0 to
! the Bahcall & Pinsonneault (1996) value; the metal settling factor;
! the solar luminosity; and the age scale factor. With the solar
! calibration on, it also sets cards 2 and 3's stop ages, which no sample
! gives (below). Outside Monte Carlo it only sets age_scale_factor = 1.
!   star reads: control%calibrate_solar_model, mc_samples,
!     physics%monte_carlo_active, physics%solar_luminosity_cgs
!   star writes: card_max_age_yr, cross_section_scale_in_use,
!     extra_cross_section_scale_in_use, ln10_over_lsun,
!     log10_solar_luminosity, settling_z_factor_in_use,
!     solar_luminosity_cgs
subroutine apply_monte_carlo_parameters(s, monte_carlo_run_number, &
     age_scale_factor)
      use math_lib
      use physical_constants, only: ln10
      use net_lib, only: r_pp, r_he3he3, r_he3he4, xs_be7p
      type(star_info), intent(inout) :: s
      integer, intent(in) :: monte_carlo_run_number
      real(dp), intent(out) :: age_scale_factor
! The ratios of the Bahcall & Pinsonneault (1996) S0 values to the
! code's reference ones, for the four sampled reactions. (The original
! data statement also had 0.992 for N14 + p and 0.1625 for hep, which
! are not sampled.) The literals are default (single-precision) reals on
! purpose: they reproduce the original data statement's values exactly,
! and appending d0 would shift them in the 8th decimal.
      real(dp), parameter :: bp96_pp_factor = 0.9558, &
           bp96_he3he3_factor = 0.9690, bp96_he3he4_factor = 0.9712, &
           bp96_be7p_factor = 0.92088
! Cards 2 and 3's stop ages [yr] for the solar calibration (below).
      real(dp), parameter :: mc_settling_age_yr = 1.0D8, mc_solar_age_yr = 4.57D9

      if (s%physics%monte_carlo_active) then
         s%cross_section_scale_in_use(r_pp) = s%mc_samples%s11_factor(monte_carlo_run_number)*bp96_pp_factor
         s%cross_section_scale_in_use(r_he3he3) = s%mc_samples%s33_factor(monte_carlo_run_number)*bp96_he3he3_factor
         s%cross_section_scale_in_use(r_he3he4) = s%mc_samples%s34_factor(monte_carlo_run_number)*bp96_he3he4_factor
         s%extra_cross_section_scale_in_use(xs_be7p) = s%mc_samples%s17_factor(monte_carlo_run_number)*bp96_be7p_factor
! The sampled metal settling factor (legacy FHE, which the F77 code also
! copied to FZDIF and FGRSET) sets element diffusion's metal factor fgrz
! only; the helium factor fgry keeps its control value. The sampled Z/X
! (metal_to_h_ratio) is only written to the summary.
         s%settling_z_factor_in_use = s%mc_samples%settling_z_factor(monte_carlo_run_number)
! The sampled solar luminosity: the deck's (physics%solar_luminosity_cgs)
! times the run's luminosity factor (M97; it was 3.844e33 erg/s, MHP 3/96,
! whatever the deck set).
         s%solar_luminosity_cgs = s%physics%solar_luminosity_cgs*s%mc_samples%luminosity_factor(monte_carlo_run_number)
         s%log10_solar_luminosity = log10(s%solar_luminosity_cgs)
         s%ln10_over_lsun = ln10/s%solar_luminosity_cgs
         age_scale_factor = s%mc_samples%age_factor(monte_carlo_run_number)
! With the solar calibration on, the stop ages of cards 2 and 3 for a
! solar age of 4.57 Gyr, the F77's (main.f:243-244), neither sampled nor
! scaled by the age factor; setup_solar_calibration then sets both
! again, with 1e8 yr and target_solar_age times age_scale_factor.
! Without the calibration the deck's max_age stands (M191).
         if (s%control%calibrate_solar_model) then
            s%card_max_age_yr(2) = mc_settling_age_yr
            s%card_max_age_yr(3) = mc_solar_age_yr
         end if
      else
         age_scale_factor = 1.0D0
      end if
end subroutine apply_monte_carlo_parameters

! ---------------------------------------------------------------
! End-of-run bookkeeping for the Monte-Carlo and solar-calibration
! runs: the final model, the starting model and the run log rewound
! (io's rewind_run_cycle_files), so that the next cycle rewrites them
! from the start, and the Monte-Carlo model summary written. Two cases
! do this:
!  - a Monte-Carlo run whose solar calibration did not converge by its
!    mc_unconverged_cycle-th cycle (the calibration gives up at 15),
!    whose model is written with its cycle count;
!  - a converged solar calibration with neutrinos on, whose model is
!    written when this is a Monte-Carlo run.
! The neutrino and structure summaries once written to the .snu file
! are in the history file: snu_cl37, snu_ga71 and the ten neutrino
! fluxes are history columns, and the structure, convection-zone and
! energy summary is the final model's history row; the calibration
! verdicts are in the .log file.
!
!   star reads: control%calibrate_solar_model, current_card,
!     physics%calc_neutrinos, physics%monte_carlo_active,
!     solar_calibration_converged
!   star writes: none
subroutine write_run_summaries(s, monte_carlo_run_number, &
     convergence_iterations, log_r_rsun, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(in) :: monte_carlo_run_number, convergence_iterations
      real(dp), intent(in) :: log_r_rsun
      integer, intent(inout) :: ierr
! A Monte-Carlo run still unconverged at this calibration cycle is
! written out as it is.
      integer, parameter :: mc_unconverged_cycle = 11

      if (s%physics%monte_carlo_active .and. convergence_iterations >= mc_unconverged_cycle .and. &
           .not.s%solar_calibration_converged) then
         call rewind_run_cycle_files()
         call write_monte_carlo_model(s, log_r_rsun, convergence_iterations, &
              s%current_card, monte_carlo_run_number, ierr)
      else if (s%control%calibrate_solar_model .and. s%physics%calc_neutrinos .and. &
           s%solar_calibration_converged) then
         call rewind_run_cycle_files()
         if (s%physics%monte_carlo_active) then
            call write_monte_carlo_model(s, log_r_rsun, convergence_iterations, &
                 s%current_card, monte_carlo_run_number, ierr)
         end if
      end if
end subroutine write_run_summaries

! ---------------------------------------------------------------
! write_monte_carlo_model (legacy wrtmonte)
! Appends the last converged model to the Monte-Carlo model file
! (write_mod_model on mc_model_unit, one model per run) and the run's
! summary to the summary file (io's write_monte_carlo_summary, whose
! monte_carlo_summary record documents the seven records): the run's
! sampled parameters; the calibration cycle count and partial
! derivatives; the neutrino fluxes; the central conditions; the shell
! count, alpha, initial composition and final L and R; the surface
! convection zone and surface composition; and the luminosities of the
! burning channels and of gravitational contraction. local_log_radius is
! log(R/Rsun), and run_index = current_card is the last card of the
! calibration triple, whose rescale card (run_index-2) holds the initial
! composition. After a failed write of the model or of the summary,
! ierr is 1.
!   star reads: card_mixing_length_alpha, card_rescale_targets,
!     central_log10_density, central_log10_pressure,
!     central_log10_temperature, cl37_snu_rate,
!     control%monte_carlo_file1_path, envelope_mass_msun,
!     envelope_radius_rsun, ga71_snu_rate, log_L_lsun,
!     luminosity_breakdown_lsun, mc_samples,
!     neutrino_flux_total_e10_per_cm2_s, nz, xa
!   star writes: none
subroutine write_monte_carlo_model(s, local_log_radius, &
     convergence_iterations, run_index, monte_carlo_run_number, ierr)
      use math_lib
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: local_log_radius
      integer, intent(in) :: convergence_iterations, run_index, &
           monte_carlo_run_number
      integer, intent(inout) :: ierr
      type(monte_carlo_summary) :: summary
      integer :: run

      call write_mod_model(s, mc_model_unit, replace=.false., ierr=ierr)
      if (ierr /= 0) return
      run = monte_carlo_run_number
      summary%run_number = run
      summary%samples = [s%mc_samples%s11_factor(run), s%mc_samples%s33_factor(run), &
           s%mc_samples%s34_factor(run), s%mc_samples%s17_factor(run), &
           s%mc_samples%metal_to_h_ratio(run), s%mc_samples%settling_z_factor(run), &
           s%mc_samples%luminosity_factor(run), s%mc_samples%age_factor(run)]
! The solar calibration's fixed partial derivatives (the F77 printed its
! COMMON copies of them, zero when the first run converged at its first
! check; M34).
      summary%calibration_cycles = convergence_iterations
      summary%calibration_derivatives = [dlogl_dx_fixed, dlogr_dx_fixed, &
           dlogl_dalpha_fixed, dlogr_dalpha_fixed]
      summary%cl37_snu = s%cl37_snu_rate
      summary%ga71_snu = s%ga71_snu_rate
      summary%neutrino_fluxes = s%neutrino_flux_total_e10_per_cm2_s(1:8)
      summary%central_t6 = exp10((s%central_log10_temperature-6.0d0))
      summary%central_p17 = exp10((s%central_log10_pressure-17.0d0))
      summary%central_density = exp10(s%central_log10_density)
      summary%central_x = s%xa(i_h1,1)
      summary%central_z = s%xa(i_metals,1)
      summary%num_shells = s%nz
      summary%alpha = s%card_mixing_length_alpha(run_index)
      summary%initial_y = 1.0d0 - s%card_rescale_targets(i_target_x,run_index-2) - &
           s%card_rescale_targets(i_target_z,run_index-2)
      summary%initial_z = s%card_rescale_targets(i_target_z,run_index-2)
      summary%log_l_lsun = s%log_L_lsun
      summary%log_r_rsun = local_log_radius
      summary%cz_base_radius_rsun = s%envelope_radius_rsun
      summary%cz_mass_msun = s%envelope_mass_msun
      summary%surface_y = s%xa(i_he4,s%nz)
      summary%surface_z = s%xa(i_metals,s%nz)
      summary%surface_z_over_x = s%xa(i_metals,s%nz)/s%xa(i_h1,s%nz)
      summary%luminosities_lsun = [s%luminosity_breakdown_lsun(1:4), &
           s%luminosity_breakdown_lsun(i_lum_grav)]
      call write_monte_carlo_summary(summary, s%control%monte_carlo_file1_path, ierr)
end subroutine write_monte_carlo_model

end module monte_carlo

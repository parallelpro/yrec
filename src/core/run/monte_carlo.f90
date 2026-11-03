!----------------------------------------------------------------------
! monte_carlo
!----------------------------------------------------------------------
! The built-in Monte-Carlo runs. Each run takes sampled nuclear
! cross-section factors, a metal settling factor, and solar luminosity
! and age factors from the sample file (control%dynamics_file; the
! fields of star%mc_samples, state/README.md):
!
!  * setup_monte_carlo_runs -- the run range, and each run's samples
!    read from the sample file (once, from prepare_run);
!  * begin_monte_carlo_run -- the start of every Monte-Carlo run, the
!    per-run hook: every run copy of a control re-armed from its control,
!    the state the previous run left that the next would read cleared
!    (clear_previous_run_state), then the run's sampled values
!    (apply_monte_carlo_parameters); a value the run changes per
!    Monte-Carlo run is set here;
!  * write_run_summaries -- the end-of-run bookkeeping, including the
!    Monte-Carlo model summary (write_monte_carlo_model, legacy
!    wrtmonte), which reads star% directly.
!
! The Monte-Carlo loop itself is in run_yrec. The deck
! examples/run_config_matrix/cm_monte_carlo.inlist makes one Monte-Carlo
! run, with the sample file mc_input.dat beside it.
module monte_carlo
      use star_info_lib, only: json, i_h1, i_he4, i_metals, i_lum_grav, &
           max_mc_runs, star_info
      use control_slots, only: i_target_x, i_target_z
      use run_messages, only: run_log_unit, terminal_unit
      use run_files, only: last_model_unit, starting_model_unit, &
           dynamics_unit, mc_summary_unit, mc_model_unit
      use model_files, only: write_mod_model
      use output_files, only: flush_output_file
      implicit none
      private
      public :: setup_monte_carlo_runs, begin_monte_carlo_run, &
           write_run_summaries

contains

! ---------------------------------------------------------------
! Set the run range and, for a Monte-Carlo job, read each run's samples
! from the sample file (opened on dynamics_unit with the run's files):
! the S11, S33, S34 and S17 factors, Z/X, the metal settling factor, and
! the luminosity and age factors. Only the records up to last_mc_run (at
! most max_mc_runs) are read (MHP 3/96).
!   star reads: last_mc_run, physics%monte_carlo_active
!   star writes: mc_samples
subroutine setup_monte_carlo_runs(s)
      type(star_info), intent(inout) :: s
      integer :: i
      call set_monte_carlo_range(s)
      if (s%physics%monte_carlo_active) then
! Read the samples, one two-line record per run.
         do i = 1,s%last_mc_run
            read(dynamics_unit,1511)s%mc_samples%s11_factor(i),s%mc_samples%s33_factor(i), &
                 s%mc_samples%s34_factor(i),s%mc_samples%s17_factor(i),s%mc_samples%metal_to_h_ratio(i), &
                 s%mc_samples%settling_z_factor(i),s%mc_samples%luminosity_factor(i), &
                 s%mc_samples%age_factor(i)
 1511       format(7X,1P7E10.3/E9.3)
            write(terminal_unit,*)i,s%mc_samples%s11_factor(i),s%mc_samples%s33_factor(i), &
                 s%mc_samples%s34_factor(i),s%mc_samples%s17_factor(i),s%mc_samples%metal_to_h_ratio(i), &
                 s%mc_samples%settling_z_factor(i),s%mc_samples%luminosity_factor(i), &
                 s%mc_samples%age_factor(i)
         end do
      endif
end subroutine setup_monte_carlo_runs

! ---------------------------------------------------------------
! The run range: first_mc_run and last_mc_run, the star's copies of the
! controls mc_run_start and mc_run_end (legacy IMBEG, IMEND), with the
! last clamped to max_mc_runs; without Monte Carlo the single run 1.
!   star reads: physics%monte_carlo_active
!   star writes: first_mc_run, last_mc_run
subroutine set_monte_carlo_range(s)
      type(star_info), intent(inout) :: s
      if (s%physics%monte_carlo_active) then
         s%last_mc_run = min(s%last_mc_run,max_mc_runs)
      else
         s%first_mc_run = 1
         s%last_mc_run = 1
      endif
end subroutine set_monte_carlo_range

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
! to, the core extension's factor. The run range, re-armed with them, is
! derived again as setup_monte_carlo_runs derived it (the loop's bounds
! were fixed when it began). The state the previous run left on the star
! that the next would read before writing it is cleared
! (clear_previous_run_state). Then the run's samples are applied
! (apply_monte_carlo_parameters). On the first run the re-arm and the
! clearing change nothing: no run copy but the range is written between
! the controls read and the loop, and the cleared state is still the
! template star's. A value the run changes per Monte-Carlo run is a
! star member of the Monte-Carlo runs class, set here from the run's
! sample (as settling_z_factor_in_use); state/README.md, "Add a value
! the run changes per Monte-Carlo run".
!   star reads: none
!   star writes: none
subroutine begin_monte_carlo_run(s, monte_carlo_run_number, &
     age_scale_factor)
      use controls_sync, only: arm_run_copies
      type(star_info), intent(inout) :: s
      integer, intent(in) :: monte_carlo_run_number
      double precision, intent(out) :: age_scale_factor
      call arm_run_copies(s)
      call set_monte_carlo_range(s)
      call clear_previous_run_state(s)
      call apply_monte_carlo_parameters(s, monte_carlo_run_number, &
           age_scale_factor)
end subroutine begin_monte_carlo_run

! ---------------------------------------------------------------
! The state a run leaves on the star, other than the run copies, that
! the next run would read before writing it, put back to its value at
! the start of the first run (the template star's zero): the
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
! the solar luminosity; and the age scale factor. Outside Monte Carlo it
! only sets age_scale_factor = 1.
!   star reads: mc_samples, physics%monte_carlo_active
!   star writes: card_max_age_yr, cross_section_scale_in_use,
!     extra_cross_section_scale_in_use, ln10_over_lsun,
!     log10_solar_luminosity, settling_z_factor_in_use,
!     solar_luminosity_cgs
subroutine apply_monte_carlo_parameters(s, monte_carlo_run_number, &
     age_scale_factor)
      use math_lib
      use physical_constants
      use net_lib, only: r_pp, r_he3he3, r_he3he4, xs_be7p
      type(star_info), intent(inout) :: s
      integer, intent(in) :: monte_carlo_run_number
      double precision, intent(out) :: age_scale_factor
! The ratios of the Bahcall & Pinsonneault (1996) S0 values to the
! code's reference ones, for the four sampled reactions. (The original
! data statement also had 0.992 for N14 + p and 0.1625 for hep, which
! are not sampled.) The literals are default (single-precision) reals on
! purpose: they reproduce the original data statement's values exactly,
! and appending d0 would shift them in the 8th decimal.
      double precision, parameter :: bp96_pp_factor = 0.9558, &
           bp96_he3he3_factor = 0.9690, bp96_he3he4_factor = 0.9712, &
           bp96_be7p_factor = 0.92088
! The base solar luminosity in erg/s (MHP 3/96), which
! luminosity_factor scales.
      double precision, parameter :: reference_solar_luminosity = 3.844D33

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
         s%solar_luminosity_cgs = reference_solar_luminosity*s%mc_samples%luminosity_factor(monte_carlo_run_number)
         s%log10_solar_luminosity = log10(s%solar_luminosity_cgs)
         s%ln10_over_lsun = ln10/s%solar_luminosity_cgs
         age_scale_factor = s%mc_samples%age_factor(monte_carlo_run_number)
! Stop ages of cards 2 and 3 for a solar age of 4.57 Gyr. With the
! solar calibration on, setup_solar_calibration overwrites both, with
! 1e8 yr and target_solar_age times age_scale_factor.
         s%card_max_age_yr(2)=1.0D8
         s%card_max_age_yr(3)=4.57D9
      else
         age_scale_factor = 1.0D0
      endif
end subroutine apply_monte_carlo_parameters

! ---------------------------------------------------------------
! End-of-run bookkeeping for the Monte-Carlo and solar-calibration
! runs: rewind last_model_unit, the starting-model unit (first_unit)
! and the run log, so that the next cycle rewrites them from the
! start, and write the Monte-Carlo model summary. Two cases do this:
!  - a Monte-Carlo run whose solar calibration did not converge
!    (convergence_iterations >= 11; the calibration gives up at 15),
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
      double precision, intent(in) :: log_r_rsun
      integer, intent(inout) :: ierr

! A Monte-Carlo run that failed to converge: rewind the output units
! and write the model with its cycle count.
      if (s%physics%monte_carlo_active .and. convergence_iterations.ge.11 .and. .not.s%solar_calibration_converged) then
         rewind(last_model_unit)
         rewind(starting_model_unit)
         rewind(run_log_unit)
         call write_monte_carlo_model(s, log_r_rsun, convergence_iterations, &
              s%current_card, monte_carlo_run_number, ierr)
      else if (s%control%calibrate_solar_model .and. s%physics%calc_neutrinos .and. s%solar_calibration_converged) then
         rewind(last_model_unit)
         rewind(starting_model_unit)
         rewind(run_log_unit)
         if (s%physics%monte_carlo_active) then
            call write_monte_carlo_model(s, log_r_rsun, convergence_iterations, &
                 s%current_card, monte_carlo_run_number, ierr)
         endif
      endif
end subroutine write_run_summaries

! ---------------------------------------------------------------
! write_monte_carlo_model (legacy wrtmonte)
! Appends the last converged model to monte_carlo_unit2
! (write_mod_model, one model per run) and a seven-record summary to
! monte_carlo_unit1: the run's sampled parameters; the calibration
! cycle count and partial derivatives; the neutrino fluxes; the central
! conditions; the shell count, alpha, initial composition and final
! L and R; the surface convection zone and surface composition; and
! the luminosities of the burning channels and of gravitational
! contraction. local_log_radius is log(R/Rsun), and run_index = current_card is
! the last card of the calibration triple.
!   star reads: card_mixing_length_alpha, card_rescale_targets,
!     central_log10_density, central_log10_pressure,
!     central_log10_temperature, cl37_snu_rate,
!     control%monte_carlo_file1_path, dlum_dalpha, dlum_dx, drad_dalpha,
!     drad_dx, envelope_mass_msun, envelope_radius_rsun, ga71_snu_rate,
!     log_L_lsun, luminosity_breakdown_lsun, mc_samples,
!     neutrino_flux_total_e10_per_cm2_s, nz, xa
!   star writes: none
subroutine write_monte_carlo_model(s, local_log_radius, &
     convergence_iterations, run_index, monte_carlo_run_number, ierr)
      use math_lib
      type(star_info), intent(inout) :: s
      double precision, intent(in) :: local_log_radius
      integer, intent(in) :: convergence_iterations, run_index, &
           monte_carlo_run_number
! a failed write of the model (write_mod_model) or of the summary
! returns ierr = 1
      integer, intent(inout) :: ierr
! --- locals ---
      integer :: iwrite, j, ios
      character(len=256) :: msg
      double precision :: surface_z_over_x, tcen, pcen, dcen, yini, zini

! Append the last converged model to monte_carlo_unit2 (legacy IMONTE2).
      iwrite = mc_model_unit
      call write_mod_model(s, iwrite, replace=.false., ierr=ierr)
      if (ierr /= 0) return
! The summary on monte_carlo_unit1, each record written unless an
! earlier one failed, then checked (output_files' flush_output_file: the
! first failed write, or a file cut short on disk). The surface Z/X:
      ios = 0
      surface_z_over_x = s%xa(i_metals,s%nz)/s%xa(i_h1,s%nz)
! Record 1: the run number and its samples. The metal settling factor
! is written twice, as WRTMONTE wrote FHE and its copy FZDIF.
      if (ios == 0) write(mc_summary_unit,10,iostat=ios,iomsg=msg)monte_carlo_run_number,s%mc_samples%s11_factor(monte_carlo_run_number), &
              s%mc_samples%s33_factor(monte_carlo_run_number),s%mc_samples%s34_factor(monte_carlo_run_number), &
              s%mc_samples%s17_factor(monte_carlo_run_number), &
              s%mc_samples%metal_to_h_ratio(monte_carlo_run_number),s%mc_samples%settling_z_factor(monte_carlo_run_number), &
              s%mc_samples%settling_z_factor(monte_carlo_run_number),s%mc_samples%luminosity_factor(monte_carlo_run_number), &
              s%mc_samples%age_factor(monte_carlo_run_number)
   10 format(I7,1P9E10.3)
! Record 2: the number of calibration cycles, and the partial
! derivatives d log L/dX, d log R/dX, d log L/d alpha, d log R/d alpha.
      if (ios == 0) write(mc_summary_unit,20,iostat=ios,iomsg=msg)convergence_iterations,s%dlum_dx,s%drad_dx, &
           s%dlum_dalpha,s%drad_dalpha
 20   format(1X,I2,1P4E11.4)
! Record 3: the Cl37 and Ga71 capture rates (SNU) and the first eight
! neutrino fluxes, pp, pep, hep, Be7, B8, N13, O15 and F17 (see
! compute_solar_neutrino_fluxes).
      if (ios == 0) write(mc_summary_unit,30,iostat=ios,iomsg=msg) s%cl37_snu_rate,s%ga71_snu_rate,(s%neutrino_flux_total_e10_per_cm2_s(j),j=1,8)
 30   format(1X,2F8.3,1P8E10.3)
! Record 4: the central T/1e6 K, density (g cm^-3),
! P/1e17 dyn cm^-2, X and Z.
      tcen = exp10((s%central_log10_temperature-6.0d0))
      pcen = exp10((s%central_log10_pressure-17.0d0))
      dcen = exp10(s%central_log10_density)
      if (ios == 0) write(mc_summary_unit,40,iostat=ios,iomsg=msg)tcen,dcen,pcen,s%xa(i_h1,1),s%xa(i_metals,1)
 40   format(1X,F7.3,F7.2,F6.3,2F8.5)
! Record 5: the number of shells, alpha, the initial Y and Z (from the
! triple's rescale card, current_card-2), and the final log(L/Lsun) and
! log(R/Rsun).
      yini = 1.0d0 - s%card_rescale_targets(i_target_x,run_index-2) - s%card_rescale_targets(i_target_z,run_index-2)
      zini = s%card_rescale_targets(i_target_z,run_index-2)
      if (ios == 0) write(mc_summary_unit,50,iostat=ios,iomsg=msg)s%nz,s%card_mixing_length_alpha(run_index),yini,zini, &
           s%log_L_lsun,local_log_radius
 50   format(I5,F7.4,2F8.5,1P2E10.3)
! Record 6: the radius (Rsun) of the base of the surface convection
! zone and the mass above it (Msun), and the surface Y, Z and Z/X.
      if (ios == 0) write(mc_summary_unit,60,iostat=ios,iomsg=msg)s%envelope_radius_rsun,s%envelope_mass_msun,s%xa(i_he4,s%nz), &
           s%xa(i_metals,s%nz),surface_z_over_x
 60   format(F8.5,F9.6,2F8.5,F9.6)
! Record 7: the luminosities (Lsun) of the pp I, pp II, pp III and CNO
! channels and of gravitational energy release.
      if (ios == 0) write(mc_summary_unit,70,iostat=ios,iomsg=msg)(s%luminosity_breakdown_lsun(j),j=1,4),s%luminosity_breakdown_lsun(i_lum_grav)
 70   format(1P5E10.3)
      call flush_output_file(mc_summary_unit, s%control%monte_carlo_file1_path, ios, msg, ierr)
end subroutine write_monte_carlo_model

end module monte_carlo

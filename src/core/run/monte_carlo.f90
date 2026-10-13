!----------------------------------------------------------------------
! monte_carlo
!----------------------------------------------------------------------
! The built-in Monte-Carlo runs. Each run takes sampled nuclear
! cross-section factors, a diffusion factor, and solar luminosity and
! age factors from the dynamics file:
!
!  * setup_monte_carlo_runs -- the run range, and each run's sampled
!    parameters read from the dynamics file;
!  * apply_monte_carlo_parameters -- applies the current run's values;
!  * write_run_summaries -- the end-of-run bookkeeping, including the
!    Monte-Carlo model summary (write_monte_carlo_model, legacy
!    wrtmonte), which reads star% directly.
!
! The Monte-Carlo loop itself is in run_yrec. No Stage-0 regression
! deck sets monte_carlo_active, so the regression suite checks only
! that the non-Monte-Carlo path is unchanged; the Monte-Carlo path is
! compile-verified only.
module monte_carlo
      use star_info_lib, only: json, i_h1, i_he4, i_metals, i_lum_grav, &
           max_mc_runs, star_info
      use io_units
      use model_files, only: write_mod_model
      implicit none
      private
      public :: setup_monte_carlo_runs, apply_monte_carlo_parameters, &
           write_run_summaries

contains

! ---------------------------------------------------------------
! Set the run range and, for a Monte-Carlo job, read each run's
! sampled parameters from the dynamics file (opened on dynamics_unit
! when the namelist is read): the S11, S33, S34 and S17 factors, Z/X,
! the diffusion factor, and the luminosity and age factors.
! mc_run_start and mc_run_end are star%job members set from the
! namelist (legacy IMBEG, IMEND). Only the records up to mc_run_end (at
! most max_mc_runs) are read (MHP 3/96). Without Monte Carlo the range
! is the single run 1.
subroutine setup_monte_carlo_runs(s)
      type(star_info), intent(inout) :: s
      integer :: i
      if (s%ctrl%monte_carlo_active) then
         s%job%mc_run_end = min(s%job%mc_run_end,max_mc_runs)
! Read the sampled parameters, one two-line record per run.
         do i = 1,s%job%mc_run_end
            read(s%ctrl%dynamics_unit,1511)s%job%s11_rate(i),s%job%s33_rate(i),s%job%s34_rate(i), &
                 s%job%s17_rate(i),s%job%metal_to_h_ratio(i),s%job%helium_fraction_param(i), &
                 s%job%luminosity_target(i),s%job%age_target(i)
 1511       format(7X,1P7E10.3/E9.3)
            write(terminal_unit,*)i,s%job%s11_rate(i),s%job%s33_rate(i),s%job%s34_rate(i),s%job%s17_rate(i), &
                 s%job%metal_to_h_ratio(i),s%job%helium_fraction_param(i), &
                 s%job%luminosity_target(i),s%job%age_target(i)
            s%job%diffusion_factor(i) = s%job%helium_fraction_param(i)
         end do
      else
         s%job%mc_run_start = 1
         s%job%mc_run_end = 1
      endif
end subroutine setup_monte_carlo_runs

! ---------------------------------------------------------------
! For a Monte-Carlo run, apply the current run's sampled values: the
! S-factor scales of p+p, He3+He3, He3+He4 and Be7+p
! (cross_section_scale slots 1, 2, 3 and 16), each times the factor
! that takes the code's reference S0 to the Bahcall & Pinsonneault
! (1996) value; the metal diffusion factor; the solar luminosity; and
! the age scale factor. Outside Monte Carlo it only sets
! age_scale_factor = 1.
subroutine apply_monte_carlo_parameters(s, monte_carlo_run_number, &
     age_scale_factor)
      use math_lib
      use physical_constants
      type(star_info), intent(inout) :: s
      integer, intent(in) :: monte_carlo_run_number
      double precision, intent(out) :: age_scale_factor
! The ratios of the Bahcall & Pinsonneault (1996) S0 values to the
! code's reference ones, by cross_section_scale slot (only slots 1, 2,
! 3 and 16 are used). The literals are default (single-precision)
! reals on purpose: they reproduce the original data statement's values
! exactly, and appending d0 would shift them in the 8th decimal.
      double precision, parameter :: bp96_scale_factor(17) = &
           [0.9558,0.9690,0.9712,1.0,1.0,0.992,1.0,1.0, &
           1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.92088,0.1625]
! The base solar luminosity in erg/s (MHP 3/96), which
! luminosity_target scales.
      double precision, parameter :: reference_solar_luminosity = 3.844D33

      if (s%ctrl%monte_carlo_active) then
         s%cross_section_scale(1) = s%job%s11_rate(monte_carlo_run_number)*bp96_scale_factor(1)
         s%cross_section_scale(2) = s%job%s33_rate(monte_carlo_run_number)*bp96_scale_factor(2)
         s%cross_section_scale(3) = s%job%s34_rate(monte_carlo_run_number)*bp96_scale_factor(3)
         s%cross_section_scale(16) = s%job%s17_rate(monte_carlo_run_number)*bp96_scale_factor(16)
! The sampled diffusion factor (helium_fraction_param, legacy FHE,
! copied to diffusion_factor) sets the metal diffusion factor fgrz
! only; the helium factor fgry keeps its control value. (The F77 code
! also copied it to FGRSET, which nothing read.) The sampled Z/X
! (metal_to_h_ratio) is only written to the summary.
         s%job%fgrz = s%job%diffusion_factor(monte_carlo_run_number)
         s%solar_luminosity_cgs = reference_solar_luminosity*s%job%luminosity_target(monte_carlo_run_number)
         s%log10_solar_luminosity = log10(s%solar_luminosity_cgs)
         s%ln10_over_lsun = ln10/s%solar_luminosity_cgs
         age_scale_factor = s%job%age_target(monte_carlo_run_number)
! Stop ages of cards 2 and 3 for a solar age of 4.57 Gyr. With the
! solar calibration on, setup_solar_calibration overwrites both, with
! 1e8 yr and target_solar_age times age_scale_factor.
         s%job%target_end_age(2)=1.0D8
         s%job%target_end_age(3)=4.57D9
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
subroutine write_run_summaries(s, monte_carlo_run_number, &
     convergence_iterations, log_r_rsun)
      type(star_info), intent(inout) :: s
      integer, intent(in) :: monte_carlo_run_number, convergence_iterations
      double precision, intent(in) :: log_r_rsun

! A Monte-Carlo run that failed to converge: rewind the output units
! and write the model with its cycle count.
      if (s%ctrl%monte_carlo_active .and. convergence_iterations.ge.11 .and. .not.s%solar_calibration_converged) then
         rewind(last_model_unit)
         rewind(s%ctrl%first_unit)
         rewind(run_log_unit)
         call write_monte_carlo_model(s, log_r_rsun, convergence_iterations, &
              s%job%nk, monte_carlo_run_number)
      else if (s%ctrl%calibrate_solar_model .and. s%ctrl%calc_neutrinos .and. s%solar_calibration_converged) then
         rewind(last_model_unit)
         rewind(s%ctrl%first_unit)
         rewind(run_log_unit)
         if (s%ctrl%monte_carlo_active) then
            call write_monte_carlo_model(s, log_r_rsun, convergence_iterations, &
                 s%job%nk, monte_carlo_run_number)
         endif
      endif
end subroutine write_run_summaries

! ---------------------------------------------------------------
! write_monte_carlo_model (legacy wrtmonte)
! Stores the last converged model on monte_carlo_unit2
! (write_mod_model) and appends a seven-record summary to
! monte_carlo_unit1: the run's sampled parameters; the calibration
! cycle count and partial derivatives; the neutrino fluxes; the central
! conditions; the shell count, alpha, initial composition and final
! L and R; the surface convection zone and surface composition; and
! the luminosities of the burning channels and of gravitational
! contraction. local_log_radius is log(R/Rsun), and run_index = nk is
! the last card of the calibration triple.
subroutine write_monte_carlo_model(s, local_log_radius, &
     convergence_iterations, run_index, monte_carlo_run_number)
      use math_lib
      type(star_info), intent(inout) :: s
      double precision, intent(in) :: local_log_radius
      integer, intent(in) :: convergence_iterations, run_index, &
           monte_carlo_run_number
! --- locals ---
      integer :: iwrite, j
      double precision :: surface_z_over_x, tcen, pcen, dcen, yini, zini

! Store the last converged model on monte_carlo_unit2 (legacy IMONTE2).
      iwrite = s%ctrl%monte_carlo_unit2
      call write_mod_model(iwrite)
! The summary on monte_carlo_unit1. The surface Z/X:
      surface_z_over_x = s%xa(i_metals,s%nz)/s%xa(i_h1,s%nz)
! Record 1: the run number and its sampled parameters.
      write(s%ctrl%monte_carlo_unit1,10)monte_carlo_run_number,s%job%s11_rate(monte_carlo_run_number), &
              s%job%s33_rate(monte_carlo_run_number),s%job%s34_rate(monte_carlo_run_number), &
              s%job%s17_rate(monte_carlo_run_number), &
              s%job%metal_to_h_ratio(monte_carlo_run_number),s%job%helium_fraction_param(monte_carlo_run_number), &
              s%job%diffusion_factor(monte_carlo_run_number),s%job%luminosity_target(monte_carlo_run_number), &
              s%job%age_target(monte_carlo_run_number)
   10 format(I7,1P9E10.3)
! Record 2: the number of calibration cycles, and the partial
! derivatives d log L/dX, d log R/dX, d log L/d alpha, d log R/d alpha.
      write(s%ctrl%monte_carlo_unit1,20)convergence_iterations,s%dlum_dx,s%drad_dx, &
           s%dlum_dalpha,s%drad_dalpha
 20   format(1X,I2,1P4E11.4)
! Record 3: the Cl37 and Ga71 capture rates (SNU) and the first eight
! neutrino fluxes, pp, pep, hep, Be7, B8, N13, O15 and F17 (see
! compute_solar_neutrino_fluxes).
      write(s%ctrl%monte_carlo_unit1,30) s%cl37_snu_rate,s%ga71_snu_rate,(s%neutrino_flux_total(j),j=1,8)
 30   format(1X,2F8.3,1P8E10.3)
! Record 4: the central T/1e6 K, density (g cm^-3),
! P/1e17 dyn cm^-2, X and Z.
      tcen = exp10((s%central_log10_temperature-6.0d0))
      pcen = exp10((s%central_log10_pressure-17.0d0))
      dcen = exp10(s%central_log10_density)
      write(s%ctrl%monte_carlo_unit1,40)tcen,dcen,pcen,s%xa(i_h1,1),s%xa(i_metals,1)
 40   format(1X,F7.3,F7.2,F6.3,2F8.5)
! Record 5: the number of shells, alpha, the initial Y and Z (from the
! triple's rescale card, nk-2), and the final log(L/Lsun) and
! log(R/Rsun).
      yini = 1.0d0 - s%job%rescale_params(2,run_index-2) - s%job%rescale_params(3,run_index-2)
      zini = s%job%rescale_params(3,run_index-2)
      write(s%ctrl%monte_carlo_unit1,50)s%nz,s%job%mixing_length_array(run_index),yini,zini, &
           s%log_L,local_log_radius
 50   format(I5,F7.4,2F8.5,1P2E10.3)
! Record 6: the radius (Rsun) of the base of the surface convection
! zone and the mass above it (Msun), and the surface Y, Z and Z/X.
      write(s%ctrl%monte_carlo_unit1,60)s%envelope_radius,s%envelope_mass,s%xa(i_he4,s%nz), &
           s%xa(i_metals,s%nz),surface_z_over_x
 60   format(F8.5,F9.6,2F8.5,F9.6)
! Record 7: the luminosities (Lsun) of the pp I, pp II, pp III and CNO
! channels and of gravitational energy release.
      write(s%ctrl%monte_carlo_unit1,70)(s%luminosity_breakdown(j),j=1,4),s%luminosity_breakdown(i_lum_grav)
 70   format(1P5E10.3)
end subroutine write_monte_carlo_model

end module monte_carlo

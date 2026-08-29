!----------------------------------------------------------------------
! monte_carlo
!----------------------------------------------------------------------
! New (2026): the built-in Monte-Carlo machinery, gathered into one
! standalone home from the three places it used to live, so run_yrec
! reads as a plain driver:
!
!  * setup_monte_carlo_runs -- the run range and the per-run sampled
!    parameter read from the dynamics file (from star_setup);
!  * apply_monte_carlo_parameters -- per-run application of the
!    sampled values (from run_yrec's contained subroutine);
!  * write_run_summaries -- the end-of-cycle bookkeeping (from io/),
!    hosting the former io/write_monte_carlo_model.f90 (its only
!    caller) as a private procedure that reads star% directly
!    instead of taking the historical 24-argument wrtlst tramp; the
!    unused initial_x_guess/initial_alpha_guess arguments (written
!    only by a format commented out in the original F77) are gone.
!
! The MC loop itself stays in run_yrec -- it is the run driver; this
! module is everything the loop body calls. No Stage-0 deck sets
! monte_carlo_active, so the regression suite verifies the non-MC path is
! byte-identical; the MC path is compile-verified only.
module monte_carlo_lib
      use star_info_lib, only: star, json, i_h1, i_metals, i_lum_grav
      use luout_lib
      implicit none
      private
      public :: setup_monte_carlo_runs, apply_monte_carlo_parameters, &
           write_run_summaries

contains

! ---------------------------------------------------------------
! Establish the run range (mc_run_start/mc_run_end) and, for a
! Monte-Carlo job, read each run's sampled parameters from the
! dynamics file (opened on dynamics_unit at namelist-read time).
! MHP 3/96 changed I/O to read in only up to max run needed.
subroutine setup_monte_carlo_runs
      integer :: i
      if (star%ctrl%monte_carlo_active) then
!c MHP 8/25 moved file open to parmin
!     OPEN(UNIT=IDYN,FILE=FDYN,FORM='FORMATTED',STATUS='OLD')
         star%job%mc_run_start = star%ctrl%imbeg
         star%job%imend = min(star%job%imend,1000)
         star%job%mc_run_end = star%job%imend
! read in monte carlo data
         do i = 1,star%job%imend
            read(star%ctrl%dynamics_unit,1511)star%job%s11_rate(i),star%job%s33_rate(i),star%job%s34_rate(i), &
                 star%job%s17_rate(i),star%job%metal_to_h_ratio(i),star%job%helium_fraction_param(i), &
                 star%job%luminosity_target(i),star%job%age_target(i)
 1511       format(7X,1P7E10.3/E9.3)
            write(terminal_unit,*)i,star%job%s11_rate(i),star%job%s33_rate(i),star%job%s34_rate(i),star%job%s17_rate(i), &
                 star%job%metal_to_h_ratio(i),star%job%helium_fraction_param(i), &
                 star%job%luminosity_target(i),star%job%age_target(i)
            star%job%diffusion_factor(i) = star%job%helium_fraction_param(i)
         end do
      else
         star%job%mc_run_start = 1
         star%job%mc_run_end = 1
      endif
end subroutine setup_monte_carlo_runs

! ---------------------------------------------------------------
! For a Monte-Carlo run, apply the current run's sampled parameters:
! nuclear cross-section scales (against the Bahcall & Pinsonneault
! 1996 reference values), the metal diffusion factor, and the solar
! luminosity/age targets. Outside Monte Carlo only the age scale
! factor (1.0) is set.
subroutine apply_monte_carlo_parameters(monte_carlo_run_number, &
     age_scale_factor)
      use phys_const_lib
      integer, intent(in) :: monte_carlo_run_number
      double precision, intent(out) :: age_scale_factor
! latest values (Bahcall and Pinsonneault 1996). NOTE: the literals
! are default-real on purpose -- the original data statement's
! single-precision constants, widened exactly as before; do not
! append d0 (it would shift the values in the 8th decimal).
      double precision, parameter :: bp96_scale_factor(17) = &
           [0.9558,0.9690,0.9712,1.0,1.0,0.992,1.0,1.0, &
           1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.92088,0.1625]
! MHP 3/96 added data for base solar age, L
      double precision, parameter :: reference_solar_luminosity = 3.844D33
      double precision :: monte_helium_diffusion_fraction

! for monte carlo run, input values of parameters being changed.
      if (star%ctrl%monte_carlo_active) then
         star%cross_section_scale(1) = star%job%s11_rate(monte_carlo_run_number)*bp96_scale_factor(1)
         star%cross_section_scale(2) = star%job%s33_rate(monte_carlo_run_number)*bp96_scale_factor(2)
         star%cross_section_scale(3) = star%job%s34_rate(monte_carlo_run_number)*bp96_scale_factor(3)
         star%cross_section_scale(16) = star%job%s17_rate(monte_carlo_run_number)*bp96_scale_factor(16)
! NOTE (2026): write-only since the original F77 (FGRSET = FHE(NN))
! -- the sampled helium diffusion factor never reaches the physics;
! only the metal factor (fgrz) is wired through. Preserved, not
! fixed; a candidate for an upstream report.
         monte_helium_diffusion_fraction = star%job%helium_fraction_param(monte_carlo_run_number)
         star%job%fgrz = star%job%diffusion_factor(monte_carlo_run_number)
         star%solar_luminosity_cgs = reference_solar_luminosity*star%job%luminosity_target(monte_carlo_run_number)
         star%log10_solar_luminosity = dlog10(star%solar_luminosity_cgs)
         star%ln_solar_luminosity = ln10/star%solar_luminosity_cgs
         age_scale_factor = star%job%age_target(monte_carlo_run_number)
! timestep and final age are altered in SR SETCAL; input #s should be
! scaled for a solar age of 4.57 Gyr
         star%job%target_end_age(2)=1.0D8
         star%job%target_end_age(3)=4.57D9
      else
         age_scale_factor = 1.0D0
      endif
end subroutine apply_monte_carlo_parameters

! ---------------------------------------------------------------
! The Monte-Carlo / solar-calibration end-of-cycle bookkeeping.
! Legacy-stream mode only: rewinds the legacy units between
! calibration/MC cycles (so each cycle rewrites them from scratch)
! and dumps the Monte-Carlo model. In MESA mode the history file
! simply accumulates every calibration cycle instead.
!
! 2026 retire-legacy: the .snu summary records that used to be
! written here are deleted with the .snu file itself. Everything
! they carried is in the history file (snu_cl37/snu_ga71 and the ten
! neutrino fluxes are history columns; the structure/CZ/energy
! summary is the final model's history row) or in the .log
! calibration verdicts.
!
! surface_z_over_x is intent(inout) deliberately: the historical
! SAVE semantics of the driver, preserved exactly (the converged
! branch refreshes it for the caller).
subroutine write_run_summaries(monte_carlo_run_number, &
     convergence_iterations, log_r_rsun, surface_z_over_x)
      integer, intent(in) :: monte_carlo_run_number, convergence_iterations
      double precision, intent(in) :: log_r_rsun
      double precision, intent(inout) :: surface_z_over_x

! FOR MONTE CARLO, REWIND OUTPUT FILES BETWEEN CYCLES.
! RUN FAILED TO CONVERGE: DUMP THE MODEL WITH THE FAILURE COUNT.
      if (star%ctrl%monte_carlo_active .and. convergence_iterations.ge.11 .and. .not.star%solar_calibration_active) then
         rewind(last_model_unit)
         rewind(star%ctrl%first_unit)
         rewind(run_log_unit)
         call write_monte_carlo_model(log_r_rsun, convergence_iterations, &
              star%job%nk, monte_carlo_run_number)
      else if (star%ctrl%calibrate_solar_model .and. star%ctrl%calc_neutrinos .and. star%solar_calibration_active) then
         rewind(last_model_unit)
         rewind(star%ctrl%first_unit)
         rewind(run_log_unit)

         surface_z_over_x = star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz)
         if (star%ctrl%monte_carlo_active) then
            call write_monte_carlo_model(log_r_rsun, convergence_iterations, &
                 star%job%nk, monte_carlo_run_number)
         endif
      endif
end subroutine write_run_summaries

! ---------------------------------------------------------------
! wrtmonte (formerly io/write_monte_carlo_model.f90; this is its
! only caller). Stores the last converged model (write_mod_model)
! and appends a summary of the run's sampled parameters, convergence
! diagnostics, central conditions, and energy generation fractions
! to the Monte-Carlo output files.
subroutine write_monte_carlo_model(local_log_radius, &
     convergence_iterations, run_index, monte_carlo_run_number)
      double precision, intent(in) :: local_log_radius
      integer, intent(in) :: convergence_iterations, run_index, &
           monte_carlo_run_number
! --- locals ---
      integer :: iwrite, j
      double precision :: surface_z_over_x, tcen, pcen, dcen, yini, zini

!  STORE LAST CONVERGED MODEL IN LOGICAL UNIT IMONTE2
      iwrite = star%ctrl%monte_carlo_unit2
      call write_mod_model(iwrite)
!  GLOBAL INFORMATION SENT TO FIRST MONTE CARLO OUTPUT FILE
!  SURFACE Z/X
      surface_z_over_x = star%xa(3,star%nz)/star%xa(1,star%nz)
!  HEADER FILE:  MONTE CARLO PARAMETERS
      write(star%ctrl%monte_carlo_unit1,10)monte_carlo_run_number,star%job%s11_rate(monte_carlo_run_number), &
              star%job%s33_rate(monte_carlo_run_number),star%job%s34_rate(monte_carlo_run_number), &
              star%job%s17_rate(monte_carlo_run_number), &
              star%job%metal_to_h_ratio(monte_carlo_run_number),star%job%helium_fraction_param(monte_carlo_run_number), &
              star%job%diffusion_factor(monte_carlo_run_number),star%job%luminosity_target(monte_carlo_run_number), &
              star%job%age_target(monte_carlo_run_number)
   10 format(I7,1P9E10.3)
!  #OF RUNS NEEDED FOR A CONVERGED MODEL, INITIAL X
!  AND ALPHA, FINAL DL/DX,DR/DX,DL/D ALPHA, DR/D ALPHA
      write(star%ctrl%monte_carlo_unit1,20)convergence_iterations,star%dlum_dx,star%drad_dx, &
           star%dlum_dalpha,star%drad_dalpha
!      WRITE(IMONTE1,20)ICONV,XGUESS,AGUESS,DLDX,DRDX,DLDA,DRDA
! 20   FORMAT(1X,I2,2F10.6,1P4E11.4)
 20   format(1X,I2,1P4E11.4)
!  NEUTRINO FLUXES (SEE ENGEB FOR DETAILS)
      write(star%ctrl%monte_carlo_unit1,30) star%cl37_snu_rate,star%ga71_snu_rate,(star%neutrino_flux_total(j),j=1,8)
 30   format(1X,2F8.3,1P8E10.3)
!  SUMMARY OF STRUCTURE : TC, RHOC, PC, XC, ZC (ADD MU C)
      tcen = 10.0d0**(star%central_log10_temperature-6.0d0)
      pcen = 10.0d0**(star%central_log10_pressure-17.0d0)
      dcen = 10.0d0**star%central_log10_density
      write(star%ctrl%monte_carlo_unit1,40)tcen,dcen,pcen,star%xa(1,1),star%xa(3,1)
 40   format(1X,F7.3,F7.2,F6.3,2F8.5)
!  #SHELLS, INITIAL ALPHA, Y, Z; FINAL R, L
      yini = 1.0d0 - star%job%rescale_params(2,run_index-2) - star%job%rescale_params(3,run_index-2)
      zini = star%job%rescale_params(3,run_index-2)
      write(star%ctrl%monte_carlo_unit1,50)star%nz,star%job%mixing_length_array(run_index),yini,zini, &
           star%log_L,local_log_radius
 50   format(I5,F7.4,2F8.5,1P2E10.3)
!  CZ DEPTH (R,M), SURFACE Y, Z, Z/X (ADD T CZ BASE, RHO CZ BASE)
      write(star%ctrl%monte_carlo_unit1,60)star%envelope_radius,star%envelope_mass,star%xa(2,star%nz), &
           star%xa(3,star%nz),surface_z_over_x
 60   format(F8.5,F9.6,2F8.5,F9.6)
!  ENERGY GENERATION FRACTIONS PP I,II,III,CNO,EGRAV
      write(star%ctrl%monte_carlo_unit1,70)(star%luminosity_breakdown(j),j=1,4),star%luminosity_breakdown(i_lum_grav)
 70   format(1P5E10.3)
end subroutine write_monte_carlo_model

end module monte_carlo_lib

!----------------------------------------------------------------------
! wrtmonte
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original wrtmonte.f; only variable names, source form, and comment
! style were updated.
!
! Writes Monte Carlo calibration output: stores the last converged
! model (via wrtlst) and appends a summary line of the run's
! parameters, convergence diagnostics, central conditions, and energy
! generation fractions to the Monte Carlo output files.
subroutine wrtmonte(hcomp, hd, hl, hp, hr, hs, ht, lc, m, age_gyr, &
     timestep_yr, total_mass_msun, log_teff, log_luminosity, core_cz_top_index, &
     envelope_cz_bottom_index, luminosity_breakdown, trial_log_temperature, &
     trial_log_luminosity, fit_point_pressure, fit_point_temperature, &
     fit_point_radius, envelope_fit_coeffs, trial_sign_flag, log_total_mass, &
     omega, local_log_radius, convergence_iterations, run_index, &
     monte_carlo_run_number)

      use star_info_lib, only: star, i_lum_grav, json
      use const_lib
      implicit none

      double precision, intent(in) :: hcomp(15,json), hd(json), hl(json), &
           hp(json), hr(json), ht(json), hs(json)
      logical, intent(in) :: lc(json)
      integer, intent(in) :: m
      double precision, intent(in) :: age_gyr, timestep_yr, total_mass_msun, &
           log_teff, log_luminosity
      integer, intent(in) :: core_cz_top_index, envelope_cz_bottom_index
! luminosity_breakdown (TLUMX) is passed through to wrtlst -> putmodel2,
! which may convert it from solar units to erg/s in place; declared
! inout to reflect that it can come back changed.
      double precision, intent(inout) :: luminosity_breakdown(8)
      double precision, intent(in) :: trial_log_temperature(3), &
           trial_log_luminosity(3), fit_point_pressure(3), &
           fit_point_temperature(3), fit_point_radius(3), &
           envelope_fit_coeffs(9), trial_sign_flag, log_total_mass
      double precision, intent(in) :: omega(json)
      double precision, intent(in) :: local_log_radius
      integer, intent(in) :: convergence_iterations, run_index, &
           monte_carlo_run_number







! former common/iomonte/: only monte_carlo_unit1/monte_carlo_unit2 are
! used here; the file-path members (monte_carlo_file1_path/
! monte_carlo_file2_path) are unused placeholders here, now
! use-associated from const_lib along with the two used members.
! --- locals ---
      integer :: iwrite, j
      double precision :: surface_z_over_x, tcen, pcen, dcen, yini, zini

!  STORE LAST CONVERGED MODEL IN LOGICAL UNIT IMONTE2
      iwrite = star%ctrl%monte_carlo_unit2
      call wrtlst(iwrite,hcomp,hd,hl,hp,hr,hs,ht,lc,trial_log_temperature, &
           trial_log_luminosity,fit_point_pressure,fit_point_temperature, &
           fit_point_radius,envelope_fit_coeffs,trial_sign_flag, &
           luminosity_breakdown,core_cz_top_index,envelope_cz_bottom_index, &
           monte_carlo_run_number,m,total_mass_msun,log_teff,log_luminosity, &
           log_total_mass,age_gyr,timestep_yr,omega)
!  GLOBAL INFORMATION SENT TO FIRST MONTE CARLO OUTPUT FILE
!  SURFACE Z/X
      surface_z_over_x = hcomp(3,m)/hcomp(1,m)
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
      write(star%ctrl%monte_carlo_unit1,40)tcen,dcen,pcen,hcomp(1,1),hcomp(3,1)
 40   format(1X,F7.3,F7.2,F6.3,2F8.5)
!  #SHELLS, INITIAL ALPHA, Y, Z; FINAL R, L
      yini = 1.0d0 - star%job%rescale_params(2,run_index-2) - star%job%rescale_params(3,run_index-2)
      zini = star%job%rescale_params(3,run_index-2)
      write(star%ctrl%monte_carlo_unit1,50)m,star%job%mixing_length_array(run_index),yini,zini, &
           log_luminosity,local_log_radius
 50   format(I5,F7.4,2F8.5,1P2E10.3)
!  CZ DEPTH (R,M), SURFACE Y, Z, Z/X (ADD T CZ BASE, RHO CZ BASE)
      write(star%ctrl%monte_carlo_unit1,60)star%envelope_radius,star%envelope_mass,hcomp(2,m), &
           hcomp(3,m),surface_z_over_x
 60   format(F8.5,F9.6,2F8.5,F9.6)
!  ENERGY GENERATION FRACTIONS PP I,II,III,CNO,EGRAV
      write(star%ctrl%monte_carlo_unit1,70)(luminosity_breakdown(j),j=1,4),luminosity_breakdown(i_lum_grav)
 70   format(1P5E10.3)
      return
end subroutine wrtmonte

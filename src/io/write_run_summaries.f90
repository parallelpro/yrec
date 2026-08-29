!----------------------------------------------------------------------
! write_run_summaries
!----------------------------------------------------------------------
! New location (2026, core/ readability phase 4): the Monte-Carlo /
! solar-calibration end-of-cycle bookkeeping, moved verbatim from
! run_yrec (it was the last output writer living in the driver).
! Legacy mode only: rewinds the legacy units between calibration/MC
! cycles (so each cycle rewrites them from scratch) and dumps the
! Monte-Carlo model. In MESA mode the history file simply accumulates
! every calibration cycle instead.
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
     convergence_iterations, initial_x_guess, initial_alpha_guess, &
     log_r_rsun, surface_z_over_x)
      use star_info_lib, only: star, i_h1, i_metals
      use luout_lib
      implicit none

      integer, intent(in) :: monte_carlo_run_number, convergence_iterations
      double precision, intent(in) :: initial_x_guess, initial_alpha_guess, &
           log_r_rsun
      double precision, intent(inout) :: surface_z_over_x

! FOR MONTE CARLO, REWIND OUTPUT FILES BETWEEN CYCLES.
! RUN FAILED TO CONVERGE: DUMP THE MODEL WITH THE FAILURE COUNT.
      if (star%ctrl%lmonte .and. convergence_iterations.ge.11 .and. .not.star%solar_calibration_active) then
         rewind(last_model_unit)
         rewind(star%ctrl%first_unit)
         rewind(debug_file_unit)
         rewind(run_log_unit)
         call write_monte_carlo_model(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%nz,star%dage, &
              star%timestep_yr,star%star_mass,star%log_Teff,star%log_L, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,star%luminosity_breakdown,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
              star%envelope_fit_coeffs,star%trial_sign_flag,star%log_total_mass,star%omega,log_r_rsun,convergence_iterations,star%job%nk,monte_carlo_run_number)
      else if (star%ctrl%calibrate_solar_model .and. star%ctrl%lsnu .and. star%solar_calibration_active) then
         rewind(last_model_unit)
         rewind(star%ctrl%first_unit)
         rewind(debug_file_unit)
         rewind(run_log_unit)

         surface_z_over_x = star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz)
         if (star%ctrl%lmonte) then
            call write_monte_carlo_model(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%nz,star%dage, &
                 star%timestep_yr,star%star_mass,star%log_Teff,star%log_L, &
                 star%core_cz_top_index,star%envelope_cz_bottom_index,star%luminosity_breakdown,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
                 star%envelope_fit_coeffs,star%trial_sign_flag,star%log_total_mass,star%omega,log_r_rsun,convergence_iterations,star%job%nk,monte_carlo_run_number)
         endif
      endif
      return
end subroutine write_run_summaries

!----------------------------------------------------------------------
! write_run_summaries
!----------------------------------------------------------------------
! New location (2026, core/ readability phase 4): the Monte-Carlo /
! solar-calibration end-of-cycle summary writer, moved verbatim from
! run_yrec (it was the last output writer living in the driver).
! Legacy mode only: rewinds the legacy units between calibration/MC
! cycles and writes the .snu summaries (+ wrtmonte); in MESA mode the
! history file simply accumulates every calibration cycle instead.
!
! surface_z_over_x is intent(inout) deliberately: the
! failed-convergence branch prints it without computing it first,
! reporting the value carried over from a previous cycle -- the
! historical SAVE semantics of the driver, preserved exactly.
subroutine write_run_summaries(monte_carlo_run_number, &
     convergence_iterations, initial_x_guess, initial_alpha_guess, &
     log_r_rsun, surface_z_over_x)
      use star_info_lib, only: star, i_h1, i_he4, i_lum_grav, i_metals
      use luout_lib
      use const_lib
      implicit none

      integer, intent(in) :: monte_carlo_run_number, convergence_iterations
      double precision, intent(in) :: initial_x_guess, initial_alpha_guess, &
           log_r_rsun
      double precision, intent(inout) :: surface_z_over_x

      double precision :: central_temperature_mk, central_pressure_scaled, &
           central_density_linear
      double precision :: initial_helium_fraction, initial_metal_fraction
      integer :: i, j


! FOR MONTE CARLO, REWIND OUTPUT FILES AND WRITE OUT SNU FLUXES AND
! MODEL PARAMTERS TO AN OUTPUT FILE.
! RUN FAILED TO CONVERGE.  WRITE FINAL INFO WITH WARNING NOTE.
! 2026 MESA-style output: this whole chain is legacy-file machinery
! (it rewinds the legacy units between calibration/MC cycles and
! writes the .snu summaries) -- legacy mode only. In MESA mode the
! history file simply accumulates every calibration cycle instead.
      if (use_legacy_output) then
      if (lmonte .and. convergence_iterations.ge.11 .and. .not.star%run%solar_calibration_active) then
         rewind(ilast)
         rewind(first_unit)
         rewind(idebug)
         rewind(itrack)
         rewind(short_file_unit)
         rewind(imodpt)
         rewind(istor)
         write(neutrino_unit,1525)star%log_L,log_r_rsun
 1525    format(5X,'DID NOT CONVERGE WITHIN 10 ATTEMPTS L,R',2F10.6)
! MONTE CARLO #, CONVERGED MIXING LENGTH AND INITIAL H, SURFACE X,
! SURFACE Z, Z/X, CENTRAL X, CENTRAL Z
         write(neutrino_unit,1519) monte_carlo_run_number,mixing_length_array(star%job%nk),rescale_params(2,star%job%nk-2),star%xa(i_h1,star%nz), &
              star%xa(i_metals,star%nz),surface_z_over_x,star%xa(i_h1,1),star%xa(i_metals,1)
 1519    format(1X,I5,3F10.6,4E10.3)
! NUMERICAL DATA : #OF RUNS NEEDED FOR A CONVERGED MODEL, INITIAL X
! AND ALPHA, FINAL DL/DX,DR/DX,DL/D ALPHA, DR/D ALPHA
         write(neutrino_unit,1518)convergence_iterations,initial_x_guess,initial_alpha_guess,star%run%dlum_dx,star%run%drad_dx,star%run%dlum_dalpha,star%run%drad_dalpha
! SUMMARY OF STRUCTURE : TC, RHOC, PC
         write(neutrino_unit, 1517)star%run%central_log10_temperature,star%run%central_log10_pressure,star%run%central_log10_density, &
              star%xa(i_h1,1),star%xa(i_metals,1)
! NEUTRINO FLUXES (SEE ENGEB FOR DETAILS)
         write(neutrino_unit, 1516) star%flux%cl37_snu_rate,star%flux%ga71_snu_rate,(star%flux%neutrino_flux_total(i),i=1,8)
!          CALL WRTMONTE(HCOMP,HD,HL,HP,HR,HS,HT,LC,M,MODEL,DAGE,
!      *        DDAGE,SMASS,TEFFL,BL,GL,LSHELL,JXBEG,JXMID,
!      *        JXEND,JCORE,JENV,TLUMX,TRIT,TRIL,PS,TS,RS,
!      *        CFENV,FTRI,HSTOT,OMEGA,RLL,ICONV,NK,NN)  ! KC 2025-05-31
         call wrtmonte(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%nz,star%run%dage, &
              star%evo%timestep_yr,star%star_mass,star%log_Teff,star%log_L, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,star%luminosity_breakdown,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
              star%envelope_fit_coeffs,star%evo%trial_sign_flag,star%log_total_mass,star%omega,log_r_rsun,convergence_iterations,star%job%nk,monte_carlo_run_number)
      else if (calibrate_solar_model .and. lsnu .and. star%run%solar_calibration_active) then
         rewind(ilast)
         rewind(first_unit)
         rewind(idebug)
         rewind(itrack)
         rewind(short_file_unit)
         rewind(imodpt)
         rewind(istor)

         surface_z_over_x = star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz)
! HEADER FILE:  MONTE CARLO PARAMETERS
         if (lmonte) then
            write(neutrino_unit,1520)monte_carlo_run_number,star%run%s11_rate(monte_carlo_run_number),star%run%s33_rate(monte_carlo_run_number),star%run%s34_rate(monte_carlo_run_number),star%run%s17_rate(monte_carlo_run_number), &
                 star%run%metal_to_h_ratio(monte_carlo_run_number),star%run%helium_fraction_param(monte_carlo_run_number),star%run%diffusion_factor(monte_carlo_run_number),star%run%luminosity_target(monte_carlo_run_number),star%run%age_target(monte_carlo_run_number)
         endif
 1520    format(I7,1P9E10.3)
! NUMERICAL DATA : #OF RUNS NEEDED FOR A CONVERGED MODEL, INITIAL X
! AND ALPHA, FINAL DL/DX,DR/DX,DL/D ALPHA, DR/D ALPHA
         write(neutrino_unit,1518)convergence_iterations,initial_x_guess,initial_alpha_guess,star%run%dlum_dx,star%run%drad_dx,star%run%dlum_dalpha,star%run%drad_dalpha
 1518    format(1X,I2,2F10.6,1P4E11.4)
! NEUTRINO FLUXES (SEE ENGEB FOR DETAILS)
         write(neutrino_unit, 1516) star%flux%cl37_snu_rate,star%flux%ga71_snu_rate,(star%flux%neutrino_flux_total(i),i=1,10)
 1516    format(1X,2F8.3,1P10E10.3)
! SUMMARY OF STRUCTURE : TC, RHOC, PC, XC, ZC (ADD MU C)
         central_temperature_mk = 10.0D0**(star%run%central_log10_temperature-6.0D0)
         central_pressure_scaled = 10.0D0**(star%run%central_log10_pressure-17.0D0)
         central_density_linear = 10.0D0**star%run%central_log10_density
         write(neutrino_unit, 1517)central_temperature_mk,central_density_linear,central_pressure_scaled,star%xa(i_h1,1),star%xa(i_metals,1)
 1517    format(1X,F7.3,F7.2,F6.3,2F8.5)
! INITIAL ALPHA,Y,Z,ALPHA; FINAL R, L
         initial_helium_fraction = 1.0D0 - rescale_params(2,star%job%nk-2) - rescale_params(3,star%job%nk-2)
         initial_metal_fraction = rescale_params(3,star%job%nk-2)
         write(neutrino_unit,1521)mixing_length_array(star%job%nk),initial_helium_fraction,initial_metal_fraction,star%log_L,log_r_rsun
 1521    format(F7.4,2F8.5,1P2E10.3)
! CZ DEPTH (R,M), SURFACE Y, Z, Z/X (ADD T CZ BASE, RHO CZ BASE)
         write(neutrino_unit,1522)star%run%envelope_radius,star%run%envelope_mass,star%xa(i_he4,star%nz),star%xa(i_metals,star%nz),surface_z_over_x
 1522    format(F8.5,F9.6,2F8.5,F9.6)
! ENERGY GENERATION FRACTIONS PP I,II,III,CNO,EGRAV
         write(neutrino_unit,1523)(star%luminosity_breakdown(j),j=1,4),star%luminosity_breakdown(i_lum_grav)
 1523    format(1P5E10.3)
         if (lmonte) then
!             CALL WRTMONTE(HCOMP,HD,HL,HP,HR,HS,HT,LC,M,MODEL,DAGE,
!      *           DDAGE,SMASS,TEFFL,BL,GL,LSHELL,JXBEG,JXMID,
!      *           JXEND,JCORE,JENV,TLUMX,TRIT,TRIL,PS,TS,RS,
!      *           CFENV,FTRI,HSTOT,OMEGA,RLL,ICONV,NK,NN)  ! KC 2025-05-31
            call wrtmonte(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%nz,star%run%dage, &
                 star%evo%timestep_yr,star%star_mass,star%log_Teff,star%log_L, &
                 star%core_cz_top_index,star%envelope_cz_bottom_index,star%luminosity_breakdown,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
                 star%envelope_fit_coeffs,star%evo%trial_sign_flag,star%log_total_mass,star%omega,log_r_rsun,convergence_iterations,star%job%nk,monte_carlo_run_number)
         endif
      endif
      end if
      return
end subroutine write_run_summaries

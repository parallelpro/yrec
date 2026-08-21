!----------------------------------------------------------------------
! eqstat2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eqstat2.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Core equation-of-state dispatcher. Given log10(T), log10(P), X, Z:
! decides whether the gas needs a Saha ionization solve (via eqsaha,
! optionally backed by the SCV table lookup in eqscve) or can be
! treated as fully ionized (via eqrelv), interpolates between the two
! regimes near the ionization-cutoff temperature, and finally checks
! whether an OPAL table (1995/2001/2006) covers the point and, if so,
! blends its result in (via oeqos/oeqos01/oeqos06 and the matching
! eqbound*/ramp logic).
subroutine eqstat2(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, ion_fraction, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, &
     electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
     specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
     dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
     specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives, &
     in_atmosphere, saha_state)

      use envelope_comp_lib
      use luout_lib
      use const_lib
      implicit none

      integer, parameter :: nts = 63, nps = 76

      double precision, intent(inout) :: log10_temperature
      double precision, intent(out) :: temperature
      double precision, intent(inout) :: log10_pressure
      double precision, intent(out) :: pressure
      double precision, intent(inout) :: log10_density
      double precision, intent(out) :: density
      double precision, intent(in) :: hydrogen_fraction
      double precision, intent(inout) :: metal_fraction
      double precision, intent(inout) :: beta
      double precision, intent(out) :: beta_inverse, beta14
      double precision, intent(inout) :: ion_fraction(3)
      double precision, intent(out) :: specific_gas_constant
      double precision, intent(inout) :: ion_mean_weight_inverse
      double precision, intent(out) :: electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      double precision, intent(inout) :: dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient
      double precision, intent(out) :: dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
           adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state









! common/scveos/: only use_scv_eos is used here (to select whether the
! Saumon-Chabrier-Van Horn table lookup, eqscve, backs the Saha
! solve). The table data itself is read/used elsewhere.
      double precision :: tlogx(nts), tablex(nts,nps,12), &
           tabley(nts,nps,12), smix(nts,nps), tablez(nts,nps,13), &
           tablenv(nts,nps,12)
      integer :: nptsx(nts), idtt, idp
      logical :: use_scv_eos
      common/scveos/ tlogx, tablex, tabley, smix, tablez, tablenv, nptsx, &
           use_scv_eos, idtt, idp

      save

! --- locals ---
      integer :: num_species, species_idx
      double precision :: atomic_weights(4), atomic_weights_full(12)
      double precision :: saha_ramp_width, saha_ramp_scale
      data num_species/12/
      data atomic_weights/0.9921d0, 0.24975d0, 0.08322d0, 0.4995d0/
      data saha_ramp_width, saha_ramp_scale/0.500d0, 2.000d0/
      data atomic_weights_full/23.0d0, 26.99d0, 24.32d0, 55.86d0, 28.1d0, &
           12.015d0, 1.008d0, 16.0d0, 14.01d0, 39.96d0, 20.19d0, 4.004d0/

      logical :: need_saha_solution, skip_relativistic_eos
      double precision :: saha_mass_fractions(12)
      double precision :: original_metal_fraction
      double precision :: metal_ratio, amu_correction, h_excess, y_excess
      double precision :: envelope_amu_over_amu
      double precision :: dfx1, dfx12, dfx4, amu_inverse, envelope_amu_frac
      logical :: valid_table_point
      double precision :: opal_log10_density, opal_density, opal_beta, &
           opal_beta_inverse, opal_beta14, opal_specific_gas_constant, &
           opal_ion_mean_weight_inverse, opal_electron_mean_weight_inverse, &
           opal_dlnrho_dlnt, opal_dlnrho_dlnp, opal_specific_heat_cp, &
           opal_adiabatic_gradient
      double precision :: ramp_factor
      logical :: in_opal_table, needs_ramp

! saved (across the dead-code numerical-derivative branch, preserved
! verbatim from the original -- see the "LSCVD" block below)
      double precision :: dtl, dtl2, ttl
      double precision :: dlnrho_dlnt_1, dlnrho_dlnp_1, specific_heat_cp_1, &
           adiabatic_gradient_1, dlnrho_dlnt_2, dlnrho_dlnp_2, &
           specific_heat_cp_2, adiabatic_gradient_2
      double precision :: dpl, dpl2, ppl
      double precision :: log10_density_1, density_1
      logical :: valid_table_point_1, valid_table_point_2
      logical :: do_scv_derivatives

! values saved across the eqsaha/eqrelv interpolation near the
! ionization cutoff
      double precision :: saha_log10_density, saha_dlnrho_dlnt, &
           saha_dlnrho_dlnp, saha_specific_heat_cp, saha_adiabatic_gradient
      double precision :: saha_dlnrho_dlnt_dt, saha_dlnrho_dlnp_dt, &
           saha_specific_heat_cp_dt, saha_specific_heat_cp_dp, &
           saha_adiabatic_gradient_dt, saha_adiabatic_gradient_dp
      double precision :: ramp_weight, ramp_weight_sq

!     LSAHA = T IF T LOW ENOUGH TO REQUIRE CALL TO EQSAHA TO SOLVE
!     SAHA IONIZATION EQUATION. IF LSAHA = F, FULLY IONIZED GAS ASSUMED.
!     LONLYS = T IF SAHA EQUATION MUST BE SOLVED BUT RELATIVISTIC
!     EQUATION OF STATE NOT NEEDED.
      need_saha_solution = log10_temperature.lt.saha_log10t_cutoff
      skip_relativistic_eos = log10_temperature.le. &
           (saha_log10t_cutoff - saha_ramp_width)

!     MHP 3/94 METAL DIFFUSION ADDED.  ASSUME ALL METALS SCALE EQUALLY.
      original_metal_fraction = metal_fraction
      if (use_diffusion_z) then
!        THE VECTOR env_comp%fxenv IS DEFINED IN STARIN AS
!        (MASS FRACTION OF SPECIES/ATOMIC WT/AMUENV) FOR THE COMPOSITION
!        XENV,ZENV.
!        AND AMU IS DEFINED AS THE SUM OF THE MASS FRACTIONS DIVIDED BY
!        THEIR ATOMIC WEIGHTS.  FOR METAL DIFFUSION, ALL THE METALS ARE
!        ASSUMED TO CHANGE EQUALLY.
         metal_ratio = metal_fraction/env_comp%envelope_metal_fraction
         amu_correction = (metal_ratio - 1.0d0)*env_comp%amuenv
         ion_mean_weight_inverse = env_comp%amuenv
         do species_idx = 1, 6
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 amu_correction*env_comp%fxenv(species_idx)
         end do
         h_excess = (hydrogen_fraction - env_comp%envelope_hydrogen_fraction)/ &
              atomic_weights_full(7)
         ion_mean_weight_inverse = ion_mean_weight_inverse + h_excess
         do species_idx = 8, 11
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 amu_correction*env_comp%fxenv(species_idx)
         end do
         y_excess = (env_comp%envelope_hydrogen_fraction + env_comp%envelope_metal_fraction - &
              hydrogen_fraction - metal_fraction)/atomic_weights_full(12)
         ion_mean_weight_inverse = ion_mean_weight_inverse + y_excess
         envelope_amu_over_amu = metal_ratio*env_comp%amuenv/ &
              ion_mean_weight_inverse
         if (need_saha_solution) then
            do species_idx = 1, 6
               saha_mass_fractions(species_idx) = envelope_amu_over_amu* &
                    env_comp%fxenv(species_idx)
            end do
            saha_mass_fractions(7) = (env_comp%fxenv(7)* &
                 env_comp%amuenv + h_excess)/ion_mean_weight_inverse
            do species_idx = 8, 11
               saha_mass_fractions(species_idx) = envelope_amu_over_amu* &
                    env_comp%fxenv(species_idx)
            end do
            saha_mass_fractions(12) = (env_comp%fxenv(12)* &
                 env_comp%amuenv + y_excess)/ion_mean_weight_inverse
         end if
      else
!        SET UP FRACTIONAL ABUNDANCES
         dfx1 = (hydrogen_fraction - env_comp%envelope_hydrogen_fraction)
         dfx12 = (metal_fraction - env_comp%envelope_metal_fraction)
         if (dabs(dfx1) + dabs(dfx12).lt.1.0d-5) then
!           USE ENVELOPE ABUNDANCES
            ion_mean_weight_inverse = env_comp%amuenv
            if (need_saha_solution) then
               do species_idx = 1, num_species
                  saha_mass_fractions(species_idx) = &
                       env_comp%fxenv(species_idx)
               end do
            end if
         else
            dfx1 = dfx1*atomic_weights(1)
            dfx12 = dfx12*atomic_weights(3)
            dfx4 = (env_comp%envelope_hydrogen_fraction + env_comp%envelope_metal_fraction - &
                 hydrogen_fraction - metal_fraction)*atomic_weights(2)
!           ASSUME EXCESS Z(METALS) IS IN THE FORM OF CARBON(12)
            ion_mean_weight_inverse = env_comp%amuenv + dfx1 + dfx4 + dfx12
            amu_inverse = 1.0d0/ion_mean_weight_inverse
            if (need_saha_solution) then
               envelope_amu_frac = env_comp%amuenv*amu_inverse
               do species_idx = 1, num_species
                  saha_mass_fractions(species_idx) = envelope_amu_frac* &
                       env_comp%fxenv(species_idx)
               end do
               saha_mass_fractions(6) = saha_mass_fractions(6) + &
                    dfx12*amu_inverse
               saha_mass_fractions(7) = saha_mass_fractions(7) + &
                    dfx1*amu_inverse
               saha_mass_fractions(12) = saha_mass_fractions(12) + &
                    dfx4*amu_inverse
            end if
         end if
      end if
!     COMPUTE RADIATION PRESSURE
      temperature = dexp(ln10*log10_temperature)
      pressure = dexp(ln10*log10_pressure)
      beta14 = radiation_constant_over_3*(temperature**2)**2/pressure
      beta = 1 - beta14

      if (beta.lt.0d0) then
         write(short_file_unit,*) &
              'EQSTAT2: BETA is negative: TL,PL,BETA,1-BETA', &
              log10_temperature, log10_pressure, beta, 1d0 - beta
      end if

      beta_inverse = 1.0d0/beta
      beta14 = 4.0d0*beta_inverse*dmin1(1.0d0, beta14)
      electron_degeneracy_parameter = 0.0d0
      if (.not.need_saha_solution) then
!        COMPUTE VALUES FOR FULLY IONIZED GAS AND RETURN
         electron_mean_weight_inverse = hydrogen_fraction*atomic_weights(1) &
              + (1.0d0 - hydrogen_fraction)*atomic_weights(4)
         specific_gas_constant = gas_constant*(ion_mean_weight_inverse + &
              electron_mean_weight_inverse)
         ion_fraction(1) = 1.0d0
         ion_fraction(2) = 0.0d0
         ion_fraction(3) = 1.0d0
         call eqrelv(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, beta, &
              ion_mean_weight_inverse, electron_mean_weight_inverse, &
              electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
              specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
              dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp)
         goto 200
      end if
!     CHECK IF SAUMON, CHABRIER, AND VAN HORN EQUATION OF STATE NEEDED.
!     THIS EOS REPLACES THE CALL TO EQSAHA, EXCEPT FOR DERIVATIVE PURPOSES.
      if (use_scv_eos) then
         call eqsaha(saha_mass_fractions, log10_temperature, temperature, &
              log10_pressure, pressure, log10_density, density, beta, &
              beta_inverse, beta14, ion_fraction, specific_gas_constant, &
              ion_mean_weight_inverse, electron_mean_weight_inverse, &
              want_derivatives, in_atmosphere, dlnrho_dlnt, dlnrho_dlnp, &
              specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
              dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp, saha_state)
         call eqscve(log10_temperature, temperature, pressure, &
              log10_density, density, hydrogen_fraction, metal_fraction, &
              beta, ion_fraction, dlnrho_dlnt, dlnrho_dlnp, &
              specific_heat_cp, adiabatic_gradient, valid_table_point)

         do_scv_derivatives = .false.   ! Do not do SCV derivatives
         if (do_scv_derivatives .and. valid_table_point) then
!           LLP 9/6/03 To get reasonable accuracy in numerical derivatives
!           in the 4-5 decimal place SCV tables, appropriate sizes for
!           the stepouts in the numerical derivatives must be obtained.
!           The row to row and column to column spacings are .20 in PL
!           and .08 in TL. Maximum row to row changes are of the order
!           of 0.2000 out of 10.0000 (in density). Maximum column to
!           column changes are of the order of .0500 out of 3.5000 (in
!           density. It appears that stepouts of plus and minus half a
!           row and column are needed to get appropriate accuracy for
!           the derivatives.
!
!           NOTE: preserved verbatim from the original, including its
!           pre-existing bug -- specific_heat_cp_2/adiabatic_gradient_2
!           below are read before ever being assigned (the "-dtl"/"-dpl"
!           calls both write into the "_1" locals, never "_2"), so this
!           whole branch, already unreachable since do_scv_derivatives
!           is hardcoded .false. just above, would also be broken if it
!           somehow ran. Not fixed here to keep this a pure
!           transliteration; flagged for anyone revisiting eqsaha/eqscve.

            dpl = .010d0
            dtl = .040d0
            ttl = log10_temperature + dtl
            temperature = 10.0d0**ttl
            call eqscve(ttl, temperature, pressure, log10_density_1, &
                 density_1, hydrogen_fraction, metal_fraction, beta, &
                 ion_fraction, dlnrho_dlnt_1, dlnrho_dlnp_1, &
                 specific_heat_cp_1, adiabatic_gradient_1, &
                 valid_table_point_1)
            ttl = log10_temperature - dtl
            temperature = 10.0d0**ttl
            call eqscve(ttl, temperature, pressure, log10_density_1, &
                 density_1, hydrogen_fraction, metal_fraction, beta, &
                 ion_fraction, dlnrho_dlnt_1, dlnrho_dlnp_1, &
                 specific_heat_cp_1, adiabatic_gradient_1, &
                 valid_table_point_1)
            dtl2 = 2d0*dtl
            dlnrho_dlnt_dt = (dlnrho_dlnt_1 - dlnrho_dlnt_2)/dtl2/ln10
            specific_heat_cp_dt = (dlog10(specific_heat_cp_1) - &
                 dlog10(specific_heat_cp_2))/dtl2
            adiabatic_gradient_dt = (dlog10(adiabatic_gradient_1) - &
                 dlog10(adiabatic_gradient_2))/dtl2

            temperature = 10.0d0**log10_temperature
            ppl = log10_pressure + dpl
            pressure = 10.0d0**ppl
            call eqscve(log10_temperature, temperature, pressure, &
                 log10_density_1, density_1, hydrogen_fraction, &
                 metal_fraction, beta, ion_fraction, dlnrho_dlnt_1, &
                 dlnrho_dlnp_1, specific_heat_cp_1, &
                 adiabatic_gradient_1, valid_table_point_1)
            ppl = log10_pressure - dpl
            pressure = 10.0d0**ppl
            call eqscve(log10_temperature, temperature, pressure, &
                 log10_density_1, density_1, hydrogen_fraction, &
                 metal_fraction, beta, ion_fraction, dlnrho_dlnt_2, &
                 dlnrho_dlnp_2, specific_heat_cp_2, &
                 adiabatic_gradient_2, valid_table_point_2)
            pressure = 10.0d0**log10_pressure
            dpl2 = 2d0*dpl
            dlnrho_dlnp_dt = (dlnrho_dlnt_1 - dlnrho_dlnt_2)/dpl2/ln10
            specific_heat_cp_dp = (dlog10(specific_heat_cp_1) - &
                 dlog10(specific_heat_cp_2))/dpl2
            adiabatic_gradient_dp = (dlog10(adiabatic_gradient_1) - &
                 dlog10(adiabatic_gradient_2))/dpl2
            dlnrho_dlnt_dt = dlnrho_dlnt_dt
            dlnrho_dlnp_dt = dlnrho_dlnp_dt
            adiabatic_gradient_dt = adiabatic_gradient_dt
            adiabatic_gradient_dp = adiabatic_gradient_dp
            specific_heat_cp_dp = specific_heat_cp_dp
            specific_heat_cp_dt = specific_heat_cp_dt
         end if
      else
!        CALL TO PRATHER EOS - Either because SCV was not requested or
!        it has failed.
         call eqsaha(saha_mass_fractions, log10_temperature, temperature, &
              log10_pressure, pressure, log10_density, density, beta, &
              beta_inverse, beta14, ion_fraction, specific_gas_constant, &
              ion_mean_weight_inverse, electron_mean_weight_inverse, &
              want_derivatives, in_atmosphere, dlnrho_dlnt, dlnrho_dlnp, &
              specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
              dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp, saha_state)
      end if
!     STORE EQSAHA VALUES FOR INTERPOLATION WITH EQRELV VALUES
      saha_log10_density = log10_density
      saha_dlnrho_dlnt = dlnrho_dlnt
      saha_dlnrho_dlnp = dlnrho_dlnp
      saha_specific_heat_cp = specific_heat_cp
      saha_adiabatic_gradient = adiabatic_gradient
      if (want_derivatives) then
         saha_dlnrho_dlnt_dt = dlnrho_dlnt_dt
         saha_dlnrho_dlnp_dt = dlnrho_dlnp_dt
         saha_specific_heat_cp_dt = specific_heat_cp_dt
         saha_specific_heat_cp_dp = specific_heat_cp_dp
         saha_adiabatic_gradient_dt = adiabatic_gradient_dt
         saha_adiabatic_gradient_dp = adiabatic_gradient_dp
      end if
!     Bypass relativistic EOS (eqrelv) if in low temperature region
      if (skip_relativistic_eos) goto 200
!     COMPUTE VALUES FOR FULLY IONIZED GAS
      electron_mean_weight_inverse = hydrogen_fraction*atomic_weights(1) + &
           (1.0d0 - hydrogen_fraction)*atomic_weights(4)
      specific_gas_constant = gas_constant*(ion_mean_weight_inverse + &
           electron_mean_weight_inverse)
      ion_fraction(1) = 1.0d0
      ion_fraction(2) = 0.0d0
      ion_fraction(3) = 1.0d0
      call eqrelv(log10_temperature, temperature, log10_pressure, &
           pressure, log10_density, density, beta, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp)
!     INTERPOLATE VALUES
      ramp_weight = saha_ramp_scale*(saha_log10t_cutoff - log10_temperature)
      ramp_weight_sq = ramp_weight*ramp_weight
      ramp_factor = ramp_weight_sq*ramp_weight* &
           (6.0d0*ramp_weight_sq - 15.0d0*ramp_weight + 10.0d0)
      saha_log10_density = saha_log10_density - log10_density
      saha_dlnrho_dlnt = saha_dlnrho_dlnt - dlnrho_dlnt
      saha_dlnrho_dlnp = saha_dlnrho_dlnp - dlnrho_dlnp
      saha_specific_heat_cp = saha_specific_heat_cp - specific_heat_cp
      saha_adiabatic_gradient = saha_adiabatic_gradient - adiabatic_gradient
      saha_specific_heat_cp_dt = saha_specific_heat_cp_dt - specific_heat_cp_dt
      saha_specific_heat_cp_dp = saha_specific_heat_cp_dp - specific_heat_cp_dp
      saha_adiabatic_gradient_dt = saha_adiabatic_gradient_dt - adiabatic_gradient_dt
      saha_adiabatic_gradient_dp = saha_adiabatic_gradient_dp - adiabatic_gradient_dp
      log10_density = log10_density + ramp_factor*saha_log10_density
      density = dexp(ln10*log10_density)
      specific_gas_constant = beta*pressure/(density*temperature)
      electron_mean_weight_inverse = specific_gas_constant/gas_constant - &
           ion_mean_weight_inverse
      ion_fraction(1) = 1.0d0 + ramp_factor*(ion_fraction(1) - 1.0d0)
      ion_fraction(2) = ramp_factor*ion_fraction(2)
      ion_fraction(3) = 1.0d0 + ramp_factor*(ion_fraction(3) - 1.0d0)
      dlnrho_dlnt = dlnrho_dlnt + ramp_factor*saha_dlnrho_dlnt
      dlnrho_dlnp = dlnrho_dlnp + ramp_factor*saha_dlnrho_dlnp
      specific_heat_cp = specific_heat_cp + ramp_factor*saha_specific_heat_cp
      adiabatic_gradient = adiabatic_gradient + &
           ramp_factor*saha_adiabatic_gradient
      if (want_derivatives) then
         dlnrho_dlnt_dt = dlnrho_dlnt_dt + ramp_factor*saha_dlnrho_dlnt_dt
         dlnrho_dlnp_dt = dlnrho_dlnp_dt + ramp_factor*saha_dlnrho_dlnp_dt
         specific_heat_cp_dt = specific_heat_cp_dt + &
              ramp_factor*saha_specific_heat_cp_dt
         specific_heat_cp_dp = specific_heat_cp_dp + &
              ramp_factor*saha_specific_heat_cp_dp
         adiabatic_gradient_dt = adiabatic_gradient_dt + &
              ramp_factor*saha_adiabatic_gradient_dt
         adiabatic_gradient_dp = adiabatic_gradient_dp + &
              ramp_factor*saha_adiabatic_gradient_dp
      end if
  200 continue

!     1995 OPAL eqos
      if (use_opal95_eos) then
      if (temperature.ge.5.0d3 .and. log10_temperature.le.8.0d0 .and. &
           log10_density.le.5.0d0) then

         call oeqos(log10_temperature, temperature, log10_pressure, &
              pressure, opal_log10_density, opal_density, &
              hydrogen_fraction, metal_fraction, opal_beta, &
              opal_beta_inverse, opal_beta14, opal_specific_gas_constant, &
              opal_ion_mean_weight_inverse, &
              opal_electron_mean_weight_inverse, opal_dlnrho_dlnt, &
              opal_dlnrho_dlnp, opal_specific_heat_cp, &
              opal_adiabatic_gradient, *998)

         call eqbound(temperature, opal_log10_density, ramp_factor, &
              in_opal_table, needs_ramp)

         if (.not.in_opal_table) goto 998  ! Point is not in OPAL 1995 EOS table, so exit.

         if (.not.needs_ramp) then
!           No ramping needed between OPAL 1995 EOS and Yale/SCV. Result
!           is fully in the OPAL 1995 table
            log10_density = opal_log10_density
            density = 10.0d0**log10_density
            beta = opal_beta
            beta_inverse = 1.0d0/beta
            beta14 = 1.0d0 - beta
            specific_gas_constant = opal_specific_gas_constant
            ion_mean_weight_inverse = opal_ion_mean_weight_inverse
            electron_mean_weight_inverse = opal_electron_mean_weight_inverse
            dlnrho_dlnt = opal_dlnrho_dlnt
            dlnrho_dlnp = opal_dlnrho_dlnp
            specific_heat_cp = opal_specific_heat_cp
            adiabatic_gradient = opal_adiabatic_gradient
         else
!           Ramping required. Result is on border between OPAL 1995 EOS
!           and Yale/SCV.
            log10_density = log10_density + &
                 ramp_factor*(opal_log10_density - log10_density)
            density = 10.0d0**log10_density
            beta = beta + ramp_factor*(opal_beta - beta)
            beta_inverse = 1.0d0/beta
            beta14 = 1.0d0 - beta
            specific_gas_constant = specific_gas_constant + &
                 ramp_factor*(opal_specific_gas_constant - specific_gas_constant)
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 ramp_factor*(opal_ion_mean_weight_inverse - ion_mean_weight_inverse)
            electron_mean_weight_inverse = electron_mean_weight_inverse + &
                 ramp_factor*(opal_electron_mean_weight_inverse - electron_mean_weight_inverse)
            dlnrho_dlnt = dlnrho_dlnt + &
                 ramp_factor*(opal_dlnrho_dlnt - dlnrho_dlnt)
            dlnrho_dlnp = dlnrho_dlnp + &
                 ramp_factor*(opal_dlnrho_dlnp - dlnrho_dlnp)
            specific_heat_cp = specific_heat_cp + &
                 ramp_factor*(opal_specific_heat_cp - specific_heat_cp)
            adiabatic_gradient = adiabatic_gradient + &
                 ramp_factor*(opal_adiabatic_gradient - adiabatic_gradient)
         end if
      end if
      end if

!     2001 OPAL eqos  LLP 6/17/03
      if (use_opal2001_eos) then
      if (temperature.ge.2.0d3 .and. temperature.le.100d6 .and. &
           log10_density.le.7.0d0) then
         call oeqos01(log10_temperature, temperature, log10_pressure, &
              pressure, opal_log10_density, opal_density, &
              hydrogen_fraction, metal_fraction, opal_beta, &
              opal_beta_inverse, opal_beta14, opal_specific_gas_constant, &
              opal_ion_mean_weight_inverse, &
              opal_electron_mean_weight_inverse, opal_dlnrho_dlnt, &
              opal_dlnrho_dlnp, opal_specific_heat_cp, &
              opal_adiabatic_gradient, *998)

         call eqbound01(temperature, opal_log10_density, ramp_factor, &
              in_opal_table, needs_ramp)
!        eqbound01 determines whether or not the point is in the OPAL
!        2001 EOS table

         if (.not.in_opal_table) goto 998  ! Point is not in OPAL 2001 table, so exit.
!        USE OPAL RESULTS IF NOT IN (RHO,T) REGIME WHERE RAMP NEEDED
         if (.not.needs_ramp) then
!           No ramping needed between OPAL 2001 EOS and Yale/SCV. Result
!           is fully in the OPAL 2001 table
            log10_density = opal_log10_density
            density = 10.0d0**log10_density
            beta = opal_beta
            beta_inverse = 1.0d0/beta
            beta14 = 1.0d0 - beta
            specific_gas_constant = opal_specific_gas_constant
            ion_mean_weight_inverse = opal_ion_mean_weight_inverse
            electron_mean_weight_inverse = opal_electron_mean_weight_inverse
            dlnrho_dlnt = opal_dlnrho_dlnt
            dlnrho_dlnp = opal_dlnrho_dlnp
            specific_heat_cp = opal_specific_heat_cp
            adiabatic_gradient = opal_adiabatic_gradient
         else
!           Ramping required. Result is on border between OPAL 2001 EOS
!           and Yale/SCV.
            log10_density = log10_density + &
                 ramp_factor*(opal_log10_density - log10_density)
            density = 10.0d0**log10_density
            beta = beta + ramp_factor*(opal_beta - beta)
            beta_inverse = 1.0d0/beta
            beta14 = 1.0d0 - beta
            specific_gas_constant = specific_gas_constant + &
                 ramp_factor*(opal_specific_gas_constant - specific_gas_constant)
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 ramp_factor*(opal_ion_mean_weight_inverse - ion_mean_weight_inverse)
            electron_mean_weight_inverse = electron_mean_weight_inverse + &
                 ramp_factor*(opal_electron_mean_weight_inverse - electron_mean_weight_inverse)
            dlnrho_dlnt = dlnrho_dlnt + &
                 ramp_factor*(opal_dlnrho_dlnt - dlnrho_dlnt)
            dlnrho_dlnp = dlnrho_dlnp + &
                 ramp_factor*(opal_dlnrho_dlnp - dlnrho_dlnp)
            specific_heat_cp = specific_heat_cp + &
                 ramp_factor*(opal_specific_heat_cp - specific_heat_cp)
            adiabatic_gradient = adiabatic_gradient + &
                 ramp_factor*(opal_adiabatic_gradient - adiabatic_gradient)
         end if
      end if
      end if

!     2006 OPAL eqos  LLP 10/13/2996
      if (use_opal2006_eos) then
      if (temperature.ge.1.870d3 .and. temperature.le.200d6 .and. &
           log10_density.le.7.0d0) then
         call oeqos06(log10_temperature, temperature, log10_pressure, &
              pressure, opal_log10_density, opal_density, &
              hydrogen_fraction, metal_fraction, opal_beta, &
              opal_beta_inverse, opal_beta14, opal_specific_gas_constant, &
              opal_ion_mean_weight_inverse, &
              opal_electron_mean_weight_inverse, opal_dlnrho_dlnt, &
              opal_dlnrho_dlnp, opal_specific_heat_cp, &
              opal_adiabatic_gradient, *998)

         call eqbound06(temperature, opal_log10_density, ramp_factor, &
              in_opal_table, needs_ramp)
!        eqbound06 determines whether or not the point is in the OPAL
!        2006 EOS table.
!        If in_opal_table is true, point is in Opal 2006 table.
!        If needs_ramp is true, point is in ramp area, and ramping is
!        required. ramp_factor is the weight of the point for ramping
!        purposes.
!        Also, to eliminate a point, one can set needs_ramp to true and
!        ramp_factor to zero.

         if (.not.in_opal_table) goto 998  ! Point is not in OPAL 2006 EOS table, so exit.

!        USE OPAL 2006 RESULTS ONLY IF NOT IN (RHO,T) REGIME WHERE
!        RAMPING is NEEDED
         if (.not.needs_ramp) then
!           No ramping needed between OPAL 2006 EOS and Yale/SCV. Result
!           is fully in the OPAL 2006 table
            log10_density = opal_log10_density
            density = 10.0d0**log10_density
            beta = opal_beta
            beta_inverse = 1.0d0/beta
            beta14 = 1.0d0 - beta
            specific_gas_constant = opal_specific_gas_constant
            ion_mean_weight_inverse = opal_ion_mean_weight_inverse
            electron_mean_weight_inverse = opal_electron_mean_weight_inverse
            dlnrho_dlnt = opal_dlnrho_dlnt
            dlnrho_dlnp = opal_dlnrho_dlnp
            specific_heat_cp = opal_specific_heat_cp
            adiabatic_gradient = opal_adiabatic_gradient
         else
!           Ramping required. Result is on border between OPAL 2006 EOS
!           and Yale/SCV.
            log10_density = log10_density + &
                 ramp_factor*(opal_log10_density - log10_density)
            density = 10.0d0**log10_density
            beta = beta + ramp_factor*(opal_beta - beta)
            beta_inverse = 1.0d0/beta
            beta14 = 1.0d0 - beta
            specific_gas_constant = specific_gas_constant + &
                 ramp_factor*(opal_specific_gas_constant - specific_gas_constant)
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 ramp_factor*(opal_ion_mean_weight_inverse - ion_mean_weight_inverse)
            electron_mean_weight_inverse = electron_mean_weight_inverse + &
                 ramp_factor*(opal_electron_mean_weight_inverse - electron_mean_weight_inverse)
            dlnrho_dlnt = dlnrho_dlnt + &
                 ramp_factor*(opal_dlnrho_dlnt - dlnrho_dlnt)
            dlnrho_dlnp = dlnrho_dlnp + &
                 ramp_factor*(opal_dlnrho_dlnp - dlnrho_dlnp)
            specific_heat_cp = specific_heat_cp + &
                 ramp_factor*(opal_specific_heat_cp - specific_heat_cp)
            adiabatic_gradient = adiabatic_gradient + &
                 ramp_factor*(opal_adiabatic_gradient - adiabatic_gradient)
         end if
      end if
      end if

  998 continue
      metal_fraction = original_metal_fraction
      return
end subroutine eqstat2

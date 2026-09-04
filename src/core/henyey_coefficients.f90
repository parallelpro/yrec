!----------------------------------------------------------------------
! henyey_coefficients (formerly coefft)
!----------------------------------------------------------------------
! Builds the Henyey structure-equation coefficients (star%elim_coeff/
! star%elim_rhs, via henyey_eliminate) for every mesh point: at each
! shell, calls the equation of state (eos_lib's eos_get), opacity
! (kap_lib's kap_get) and temperature-gradient
! (temperature_gradients) routines to get the local physics,
! optionally the nuclear energy generation (burn_lib's engeb) and
! gravitational/entropy ("Kelvin-Helmholtz") energy term, assembles
! the pressure/temperature/radius/luminosity equation residuals and
! derivatives (a henyey_shell_terms record, cur; prev carries the
! shell below), and forward-eliminates each shell pair via
! henyey_eliminate. Also stores the diagnostic per-zone physics used
! elsewhere for output (star%eps_*, star%pulse_*, star%gradr/gradT/
! grada, star%adiabatic_index_gamma1, rotation scratch, ...).
!
! Three of the shell terms (pt_scr%qt, pt_scr%qp, pt_scr%qtl -- the
! historical COMMON/PULSE2/ slots, now point_scratch_lib) are
! simultaneously pulse-output scratch: this file computes them as the
! live temperature/pressure Henyey-equation terms and mirrors them
! into cur%qt/cur%qp/cur%qt_dl for the elimination, while the pulse
! writers read them from pt_scr afterwards. They keep the historical
! pt_scr names here despite those names no longer being very
! descriptive of their role in this file.
subroutine henyey_coefficients(delta_time, in_atmosphere, &
     want_derivatives, mixing_active, conductive_opacity_flag, &
     dlnrho_dlnt, dlnrho_dlnp, saha_state, envelope_zone_index, &
     ierr)
      use temperature_gradients_lib
      use rotation_scratch_lib
      use henyey_eliminate_lib

      use net_lib
      use star_info_lib, only: star, i_eps_grav, i_eps_neu
      use point_scratch_lib
      use phys_const_lib
      use eos_lib
      use kap_lib
      use burn_lib
      use math_lib
      implicit none

      double precision, intent(in) :: delta_time
      logical, intent(out) :: in_atmosphere, want_derivatives, &
           mixing_active, conductive_opacity_flag
      double precision, intent(out) :: dlnrho_dlnt, dlnrho_dlnp
      integer, intent(inout) :: saha_state
      integer, intent(in) :: envelope_zone_index

! --- locals ---
      double precision :: energy_gen_component(6)
      double precision :: zone_energy_luminosity, zone_log_mass, &
           zone_log_temperature, zone_log_pressure, zone_log_radius, &
           zone_luminosity_lsun
      double precision :: hydrogen_fraction, helium_fraction, &
           metal_fraction, he3_fraction, c12_fraction, c13_fraction, &
           n14_fraction, n15_fraction, o16_fraction, o17_fraction, &
           o18_fraction, deuterium_fraction, li6_fraction, li7_fraction, &
           be9_fraction
      integer :: shell_index
      double precision :: zone_log10_density
      double precision :: pressure_rotation_factor, temperature_rotation_factor
! 2026 named-index results (same shape as shell_physics): the eos/kap
! relay soup is two arrays; the dlnrho_dlnt/p DUMMIES stay explicit
! (they are this routine's outputs to the driver) and are packed/
! unpacked around the eos call to keep the historical inout guess
! chain through the caller's storage.
      double precision :: eos_res(num_eos_results), kap_res(num_kap_results)
      double precision :: actual_gradient, radiative_gradient, &
           dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
           convective_velocity
      logical :: is_convective
      double precision :: qtemp
! Henyey structure-equation terms at the current shell (cur) and,
! carried from the previous loop iteration, at the shell below (prev);
! the qp/qt/qt_dl slots live in pt_scr (pulse scratch) and are
! mirrored into cur just before the elimination.
      type(henyey_shell_terms) :: prev, cur
      double precision :: pp_chain_gen, he3he4_be7_electron_gen, &
           he3he4_be7_proton_gen, cno_gen, triple_alpha_gen, &
           zone_dlnepsilon_dlnrho, zone_dlnepsilon_dlnt, total_energy_gen
      double precision :: energy_gen_rate, alpha_capture_energy_zone
      double precision :: zone_dt, entropy_term1, entropy_term2, &
           entropy_term, entropy_term3, egrav
      double precision :: zone_log_temperature_delta, zone_log_pressure_delta
      logical :: compute_entropy_term
      double precision :: delta_time_inv, one_year_sec, one_year_sec_inv
      double precision :: cccql
      integer :: im, im1, j
      double precision :: energy_sum_inverse
      double precision :: chi_rho, chi_t, specific_heat_cv
      double precision :: total_energy_sum, neutrino_and_grav_sum
! 7/91 MHP
! ZERO OUT SOLAR NEUTRINO FLUXES.
! FLUXTOT = TOTAL FLUX OF EACH OF THE NEUTRINOS PRODUCED IN THE SUN
      integer, intent(out) :: ierr

      ierr = 0

      if (star%ctrl%calc_neutrinos) then
         do j = 1,10
            star%neutrino_flux_total(j) = 0.0d0
         end do
      end if
      conductive_opacity_flag = .true.
      want_derivatives = .true.
      in_atmosphere = .false.
      mixing_active = .false.
      compute_entropy_term = delta_time.gt.0.0d0
      if (compute_entropy_term) then
       delta_time_inv = 1.0d0/delta_time
       one_year_sec = 3.1558d7
       one_year_sec_inv = 3.1688d-8
      end if
      do j = 1,8
       star%luminosity_breakdown(j) = 0.0d0
      end do
      do im = 1,star%nz
! SET UP LOCAL VARIABLES FOR CALLS TO BASIC PHYSICS ROUTINES
       zone_energy_luminosity = 0.0d0
       zone_log_mass = star%log_mass(im)
       zone_log_temperature = star%logT(im)
       zone_log_pressure = star%logP(im)
       zone_log_radius = star%logR(im)
       zone_luminosity_lsun = star%luminosity_lsun(im)
       hydrogen_fraction = star%xa(1,im)
       helium_fraction = star%xa(2,im)
       metal_fraction = star%xa(3,im)
       he3_fraction = star%xa(4,im)
       c12_fraction = star%xa(5,im)
       c13_fraction = star%xa(6,im)
       n14_fraction = star%xa(7,im)
       n15_fraction = star%xa(8,im)
       o16_fraction = star%xa(9,im)
       o17_fraction = star%xa(10,im)
       o18_fraction = star%xa(11,im)
! the extended-composition species are always defined; they are
! passed to engeb anyway.
       deuterium_fraction = star%xa(12,im)
       li6_fraction = star%xa(13,im)
       li7_fraction = star%xa(14,im)
       be9_fraction = star%xa(15,im)
       shell_index = im
       zone_log10_density = star%logRho(im)
       pressure_rotation_factor = star%fp_rot(im)
        temperature_rotation_factor = star%ft_rot(im)
         eos_res(i_log10_density) = zone_log10_density
         eos_res(i_dlnrho_dlnt) = dlnrho_dlnt
         eos_res(i_dlnrho_dlnp) = dlnrho_dlnp
         call eos_get(zone_log_temperature, zone_log_pressure, &
              hydrogen_fraction, metal_fraction, eos_res, &
              want_derivatives, in_atmosphere, saha_state, &
              composition_at_zone=star%xa(:,im), ierr=ierr)
         if (ierr /= 0) return
         dlnrho_dlnt = eos_res(i_dlnrho_dlnt)
         dlnrho_dlnp = eos_res(i_dlnrho_dlnp)
! eqstat historically updated zone_log10_density in place; the tail
! (cur%qr, energy generation, star%elim_rhs) reads the updated value.
         zone_log10_density = eos_res(i_log10_density)
! DBG 12/95 GET OPACITY (at eqstat's returned density -- the
! historical inout dataflow)
         call kap_get(eos_res(i_log10_density), zone_log_temperature, &
              hydrogen_fraction, metal_fraction, kap_res, &
              eos_res(i_fxion:i_fxion+2), ierr=ierr)
         if (ierr /= 0) return
         star%iovim = im
         call temperature_gradients(zone_log_temperature, zone_log_pressure, &
              eos_res, kap_res, zone_log_radius, zone_log_mass, &
              zone_luminosity_lsun, actual_gradient, radiative_gradient, &
              dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
              convective_velocity, want_derivatives, is_convective, &
              pressure_rotation_factor, temperature_rotation_factor, &
              star%log_Teff, ierr)
         if (ierr /= 0) return
       star%logRho(im) = eos_res(i_log10_density)
! COMPUTE DERIVATIVES
       qtemp = c4pil + zone_log_radius + zone_log_radius + zone_log_radius
       cur%qr =+exp(ln10*(zone_log_mass - zone_log10_density - qtemp))
       cur%qr_dr = - cur%qr - cur%qr - cur%qr
       cur%qr_dp = -cur%qr*dlnrho_dlnp
       cur%qr_dt = -cur%qr*dlnrho_dlnt
       pt_scr%qp =-exp(ln10*(cgl + zone_log_mass + zone_log_mass - &
            zone_log_pressure - qtemp - zone_log_radius ))*star%fp_rot(im)
!       QPR = -QP - QP - QP - QP*(1.0D0 - QFPR)
       cur%qp_dr = -pt_scr%qp - pt_scr%qp - pt_scr%qp - pt_scr%qp
       cur%qp_dp = -pt_scr%qp
       star%convective_flag(im) = is_convective
       pt_scr%qt = actual_gradient*pt_scr%qp
       cur%qt_dr = -pt_scr%qt - pt_scr%qt - pt_scr%qt - pt_scr%qt
!       QTR = -QT - QT - QT - QT*(1.0D0 - QFTR)
       if (.not.is_convective) then
! TEMPERATURE GRADIENT IS RADIATIVE
          pt_scr%qtl = clni*pt_scr%qt/zone_luminosity_lsun
          cur%qt_dp = pt_scr%qt*kap_res(i_dlnkap_dlnrho)*dlnrho_dlnp
          cur%qt_dt = pt_scr%qt*(-4.0d0 + kap_res(i_dlnkap_dlnt) + kap_res(i_dlnkap_dlnrho)*dlnrho_dlnt)
       else
! TEMPERATURE GRADIENT IS CONVECTIVE
          pt_scr%qtl = 0.0d0
          cur%qt_dp = pt_scr%qt*(-1.0d0 + dgrad_dp_component)
          cur%qt_dt = pt_scr%qt*dgrad_dt_component
          cur%qt_dr = cur%qt_dr + pt_scr%qt*dgrad_dr_component
       end if
       cur%ql = 0.0d0
       cur%ql_dt = 0.0d0
       cur%ql_dp = 0.0d0
       if (zone_log_temperature.gt.star%ctrl%nuclear_logT_cutoffs(1)) then
! SET UP NUCLEAR ENERGY TERMS
            call engeb(pp_chain_gen, he3he4_be7_electron_gen, &
                 he3he4_be7_proton_gen, cno_gen, triple_alpha_gen, &
                 zone_dlnepsilon_dlnrho, zone_dlnepsilon_dlnt, &
                 total_energy_gen, zone_log10_density, &
                 zone_log_temperature, hydrogen_fraction, helium_fraction, &
                 he3_fraction, c12_fraction, c13_fraction, n14_fraction, &
                 o16_fraction, o18_fraction, deuterium_fraction, &
                 shell_index)
            energy_gen_rate = total_energy_gen
            energy_gen_component(1) = pp_chain_gen
            energy_gen_component(2) = he3he4_be7_electron_gen
            energy_gen_component(3) = he3he4_be7_proton_gen
            energy_gen_component(4) = cno_gen
            energy_gen_component(5) = triple_alpha_gen
            energy_gen_component(6) = star%neutrino_loss_rate
            alpha_capture_energy_zone = star%alpha_capture_energy
! 7/91 MHP
! CONVERT NEUTRINO FLUX RATES (UNITS 10**10 ERGS PER GM)
! TO UNITS OF 10**10 ERGS BY MULTIPLYING BY THE MASS.
            if (star%ctrl%calc_neutrinos) then
               do j = 1,10
                  star%neutrino_flux_total(j) = star%neutrino_flux_total(j) + &
                       star%neutrino_flux(j)*star%dm(im)
               end do
            end if
            do j = 1,6
               star%luminosity_breakdown(j) = star%luminosity_breakdown(j) + &
                    (star%dm(im)/star%solar_luminosity_cgs)* &
                    energy_gen_component(j)
               zone_energy_luminosity = zone_energy_luminosity + &
                    (star%dm(im)/star%solar_luminosity_cgs)* &
                    energy_gen_component(j)
            end do
! JVS 10/11 Calculate the He3+He3 and sum of He3+He3 and He3+He4 luminosity
            star%he3_he3_luminosity_zone(im) = (star%dm(im)/ &
                 star%solar_luminosity_cgs)*star%he3_he3_energy_rate
            star%he3_burning_luminosity_zone(im) = (star%dm(im)/ &
                 star%solar_luminosity_cgs)*star%he3_burning_energy_rate
! JVS end
            star%luminosity_breakdown(8)=star%luminosity_breakdown(8)+(star%dm(im)/ &
                 star%solar_luminosity_cgs)*alpha_capture_energy_zone
            zone_energy_luminosity = zone_energy_luminosity + &
                 (star%dm(im)/star%solar_luminosity_cgs)* &
                 alpha_capture_energy_zone
            cur%ql = energy_gen_rate
            cur%ql_dt = cur%ql_dt + zone_dlnepsilon_dlnt + &
                 zone_dlnepsilon_dlnrho*dlnrho_dlnt
            cur%ql_dp = cur%ql_dp + zone_dlnepsilon_dlnrho*dlnrho_dlnp
         end if
         if (compute_entropy_term) then
! SET UP ENTROPY TERMS
            zone_dt = delta_time_inv
            if (star%job%use_mass_accretion.and.star%ctrl%mass_accretion_rate.gt.0.0d0) then
               if (im.ge.envelope_zone_index) then
                  zone_log_temperature_delta = star%log_temperature_delta(im)+ &
                       rot_scr%delta_log_temperature
                  zone_log_pressure_delta = star%log_pressure_delta(im)+ &
                       rot_scr%delta_log_pressure
               else
                  zone_log_temperature_delta = star%log_temperature_delta(im)
                  zone_log_pressure_delta = star%log_pressure_delta(im)
               end if
            else
               zone_log_temperature_delta = star%log_temperature_delta(im)
               zone_log_pressure_delta = star%log_pressure_delta(im)
            end if
            if (star%xa(1,im).gt.0.01d0 .and. delta_time.lt.one_year_sec) &
                 zone_dt = one_year_sec_inv
            entropy_term1 = eos_res(i_pressure)*dlnrho_dlnt/eos_res(i_density)
            entropy_term2 = entropy_term1/eos_res(i_grada)
            entropy_term = (entropy_term2*zone_log_temperature_delta - &
                 entropy_term1*zone_log_pressure_delta)*ln10
            entropy_term3 = entropy_term2*ln10*zone_log_temperature_delta
!            ENTR = (ENTR2*HTT(IM) - ENTR1*HPP(IM))*CLN
!            ENTR3 = ENTR2*CLN*HTT(IM)
            egrav = zone_dt*entropy_term
            star%gravitational_luminosity(im) = egrav
            star%luminosity_breakdown(7) = star%luminosity_breakdown(7) + (star%dm(im)/ &
                 star%solar_luminosity_cgs)*egrav
            cur%ql = cur%ql + egrav
            cur%ql_dp = cur%ql_dp + zone_dt*(entropy_term*(1.0d0-dlnrho_dlnp+ &
                 eos_res(i_dlnrho_dlnp_dt))-entropy_term1 - entropy_term3* &
                 eos_res(i_grada_dp))
            cur%ql_dt = cur%ql_dt + zone_dt*(entropy_term*(-dlnrho_dlnt+ &
                 eos_res(i_dlnrho_dlnt_dt)) + entropy_term2 - entropy_term3* &
                 eos_res(i_grada_dt))
! 7/92 INCLUDE CHANGE IN ROTATIONAL KINETIC ENERGY IN ENERGY EQUATION.
            if (star%job%rotation_active) then
               rot_scr%rotational_energy_term(im) = zone_dt*(star%kinetic_energy_rot(im)- &
                    star%kinetic_energy_rot_old(im))/star%dm(im)
               cur%ql = cur%ql - rot_scr%rotational_energy_term(im)
            end if
! ADD CHANGE IN ENTROPY FROM ACCRETED MATERIAL
         end if
         cccql = star%ln_solar_luminosity*star%m(im)
         cur%ql = cccql*cur%ql
         cur%ql_dp = cccql*cur%ql_dp
         cur%ql_dt = cccql*cur%ql_dt
! the qp/qt/qt_dl slots live in pt_scr (pulse scratch); mirror them
! into the current shell's term record before eliminating
         cur%qp = pt_scr%qp
         cur%qt = pt_scr%qt
         cur%qt_dl = pt_scr%qtl
         if (im.gt.1) then
! REDUCE MATRIX FOR PAIR OF POINTS (IM-1,IM)
            im1 = im
            call henyey_eliminate(im1, prev, cur)
         else
! SETUP CENTRAL BOUNDARY CONDITIONS
            star%elim_coeff(3,1,1) = cc13*dlnrho_dlnp
            star%elim_coeff(3,2,1) = cc13*dlnrho_dlnt
            star%elim_rhs(3,1) = -cc13*(c4pi3l + zone_log10_density - &
                 zone_log_mass) - zone_log_radius
            star%elim_coeff(4,1,1) = -cur%ql_dp
            star%elim_coeff(4,2,1) = -cur%ql_dt
            star%elim_rhs(4,1) = clni*cur%ql - zone_luminosity_lsun
         end if
! carry this shell's terms into the next pair's elimination
         prev = cur
!  STORE THE PER-ZONE PHYSICS FOR OUTPUT (always, for every model).
!  ZERO OUT NUCLEAR ENERGY TERMS IF T < NUCLEAR CUTOFF.
         if (star%logT(im).lt.star%ctrl%nuclear_logT_cutoffs(1)) then
            star%eps_total(im) = 0.0d0
            star%eps_channels(i_eps_grav,im) = star%gravitational_luminosity(im)
            do j = 1,6
               star%eps_channels(j,im) = 0.0d0
           end do
         else
            star%eps_total(im) = energy_gen_component(1)+energy_gen_component(2)+ &
                 energy_gen_component(3)+energy_gen_component(4)+ &
                 energy_gen_component(5)
            star%eps_channels(i_eps_neu,im) = energy_gen_component(6)
            star%eps_channels(i_eps_grav,im) = star%gravitational_luminosity(im)
            if (star%eps_total(im).gt.1.0d-22) then
               energy_sum_inverse = 1.0d0/star%eps_total(im)
            else
               energy_sum_inverse = 0.0d0
            end if
            do j = 1,5
               star%eps_channels(j,im) = energy_gen_component(j)*energy_sum_inverse
              end do
         end if
         star%beta(im) = eos_res(i_beta)
         star%eta(im) = eos_res(i_eta)
         star%converged_zone(im) = conductive_opacity_flag
         star%opacity_zone(im) = kap_res(i_kap)
         star%gradr(im) = radiative_gradient
         star%gradT(im) = actual_gradient
         star%grada(im) = eos_res(i_grada)
         do j = 1,3
            star%fxion_zone(j,im) = eos_res(i_fxion-1+j)
         end do
         star%conv_vel(im) = convective_velocity
         star%scp(im) = eos_res(i_cp)
! GAMMA1 FROM GUENTHER 1995 P.C. (always computed):
            chi_rho = 1.0d0/dlnrho_dlnp
            chi_t = -chi_rho*dlnrho_dlnt
            specific_heat_cv = eos_res(i_cp) - exp(ln10*(star%logP(im)- &
                 star%logRho(im)-star%logT(im)))*chi_t**2/chi_rho
            star%adiabatic_index_gamma1(im) = chi_rho*eos_res(i_cp)/ &
                 specific_heat_cv

         if (star%job%rotation_active) then
            rot_scr%dlnkappa_dlnrho(im) = kap_res(i_dlnkap_dlnrho)
            rot_scr%dlnkappa_dlnt(im) = kap_res(i_dlnkap_dlnt)
! MHP 10/02 variable index error
            if (star%eps_total(im).gt.0.0d0) then
               total_energy_sum = star%eps_total(im)
               neutrino_and_grav_sum = star%eps_channels(i_eps_neu,im)+star%eps_channels(i_eps_grav,im)
               rot_scr%neutrino_loss_fraction(im) = (total_energy_sum - &
                    neutrino_and_grav_sum)/total_energy_sum
            else
               rot_scr%neutrino_loss_fraction(im) = 0.0d0
            end if
            rot_scr%dlnepsilon_dlnrho(im) = zone_dlnepsilon_dlnrho
            rot_scr%dlnepsilon_dlnt(im) = zone_dlnepsilon_dlnt
         end if
! DBG PULSE
! MHP 8/25 unconditional: previously gated on pulsation_output_active
! (the legacy path-length-triggered OPAL pulsation writer's flag), but
! the GYRE-format pulsation writer (io/write_gyre_pulse.f90) is
! triggered independently (by model-number interval, checked in
! write_legacy_output.f90) and needs these populated for every converged model, not
! just ones flagged by the older mechanism. All source locals here
! (dlnrho_dlnp, dlnkap_dlnrho/dlnt, zone_dlnepsilon_dlnrho/dlnt, etc.)
! are already computed unconditionally above, so this is just always
! copying them into the pulse1/mixing-length output arrays -- no
! change to any existing output (.short/.track/.store/.pmod/.penv/
! .patm) values, since none of those read from these arrays.
         star%pulse_dlnrho_dlnp(im) = dlnrho_dlnp
         star%pulse_dlneps_dlnrho(im) = zone_dlnepsilon_dlnrho
         star%pulse_dlneps_dlnt(im) = zone_dlnepsilon_dlnt
         star%pulse_dlnkap_dlnrho(im) = kap_res(i_dlnkap_dlnrho)
         star%pulse_dlnkap_dlnt(im) = kap_res(i_dlnkap_dlnt)
         star%pulse_specific_heat(im) = eos_res(i_cp)
! 2026 (bugsweep sec-11): eos_res(i_gas_constant) is R/mu (the
! specific gas constant), not mu -- store the real mean molecular
! weight; the electron slot is 1/mu_e straight from the eos.
         if (eos_res(i_gas_constant) > 0.0d0) then
            star%pulse_mean_molecular_weight(im) = gas_constant/eos_res(i_gas_constant)
         else
            star%pulse_mean_molecular_weight(im) = 0.0d0
         end if
         star%pulse_electron_mean_weight_inverse(im) = eos_res(i_mu_e_inv)
         star%pulse_dlnrho_dlnt(im) = dlnrho_dlnt
! star%alfmlt/phmlt/cmxmlt are only ever set to zero (mixing/
! temperature_gradients.f90), so these output arrays are always zero;
! pending an author decision -- see audit SUMMARY.md 1.1 #10.
         star%valfmlt(im) = star%alfmlt
         star%vphmlt(im) = star%phmlt
         star%vcmxmlt(im) = star%cmxmlt
      end do

      return
end subroutine henyey_coefficients

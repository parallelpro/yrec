!----------------------------------------------------------------------
! coefft
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original coefft.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! 2/91 MHP FLAG TO TOGGLE BETWEEN OLD/NEW ENERGY GENERATION ROUTINES
! ADDED (COMMON BLOCK NEWENG).
!
! Builds the Henyey structure-equation coefficients (elim_coeff/
! elim_rhs, via reduce) for every mesh point: at each shell, calls the
! equation of state (meqos or eqstat), opacity (getopac), and
! temperature-gradient (tpgrad) routines to get the local physics,
! optionally the nuclear energy generation (engeb) and gravitational/
! entropy ("Kelvin-Helmholtz") energy term, assembles the pressure/
! temperature/radius/luminosity equation residuals and derivatives
! (eq_p_val/dqp_dr/.../eq_l_val/dql_dt and their previous-shell "0"
! counterparts), and forward-eliminates each shell pair via reduce.
! Also stores the diagnostic per-zone physics used elsewhere for
! output (common/scrtch/, common/pulse1/, common/sound/, common/
! rotder/, common/roten/).
!
! Three dummy arguments below (pulse_diag%qt, pulse_diag%qp, pulse_diag%qtl -- see common/pulse2/) are
! simultaneously common-block storage: COMMON/PULSE2/ was first
! declared (unused, as placeholders) in an earlier-converted file
! (wrtmod.f90) with generic member names kept close to the original;
! this file is the first to actually USE those three slots, as the
! live temperature/pressure Henyey-equation scratch terms (they also
! double as the eq_t_val/eq_p_val/dqt_dl arguments passed to reduce),
! so per the project's COMMON-block-reuse rule they keep the
! wrtmod.f90 names (pulse_diag%qt/pulse_diag%qp/pulse_diag%qtl) here too, despite those names no
! longer being very descriptive of their role in this file.
!
! KC 2025-05-31 removed the unused MODEL argument and reordered the
! trailing argument list slightly (see the commented-out original
! signature below).
!       SUBROUTINE COEFFT(DELTS,M,HD,HHA,HHB,HHC,HL,HMAX,HP,HPP,HR,HS,
!      *HS1,HS2,HT,HTT,HCOMP,LC,TLUMX,LATMO,LDERIV,LMIX,LOCOND,QDT,QDP,
!      *KSAHA,MODEL,FP,FT,HKEROT,HKEROTO,JENV,TEFFL)  ! KC 2025-05-31
subroutine coefft(delta_time, num_points, log10_density, elim_coeff, &
     elim_rhs, gravitational_luminosity, luminosity_lsun, max_residual, &
     log_pressure, log_pressure_delta, log_radius, log_mass, &
     mass_weight_ln, shell_mass, log_temperature, log_temperature_delta, &
     composition, convective_flag, luminosity_terms, in_atmosphere, &
     want_derivatives, mixing_active, conductive_opacity_flag, &
     dlnrho_dlnt, dlnrho_dlnp, saha_state, &
     rotation_p_factor, rotation_t_factor, kinetic_energy_rot, &
     kinetic_energy_rot_old, envelope_zone_index, log_teff)

      use pulse_diag_lib
      use fluxes_lib
      use engeb_diag_lib
      use scrtch_lib
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json=5000

      double precision, intent(in) :: delta_time
      integer, intent(in) :: num_points
      double precision, intent(inout) :: log10_density(json)
      double precision, intent(inout) :: elim_coeff(4,2,json), &
           elim_rhs(4,json)
      double precision, intent(inout) :: gravitational_luminosity(json)
      double precision, intent(in) :: luminosity_lsun(json)
      double precision, intent(inout) :: max_residual(4)
      double precision, intent(in) :: log_pressure(json), &
           log_pressure_delta(json), log_radius(json), log_mass(json), &
           mass_weight_ln(json), shell_mass(json), log_temperature(json), &
           log_temperature_delta(json)
      double precision, intent(in) :: composition(15,json)
      logical, intent(out) :: convective_flag(json)
      double precision, intent(out) :: luminosity_terms(8)
      logical, intent(out) :: in_atmosphere, want_derivatives, &
           mixing_active, conductive_opacity_flag
      double precision, intent(out) :: dlnrho_dlnt, dlnrho_dlnp
      integer, intent(inout) :: saha_state
      double precision, intent(in) :: rotation_p_factor(json), &
           rotation_t_factor(json), kinetic_energy_rot(json), &
           kinetic_energy_rot_old(json)
      integer, intent(in) :: envelope_zone_index
      double precision, intent(in) :: log_teff

! DBG 7/92 COMMON BLOCK ADDED TO COMPUTE DEBYE-HUCKEL CORRECTION.
      double precision :: cdh, etadh0, etadh1, zdh(18), xxdy, yydh, zzdh, &
           dhnue(18)
      logical :: ldh
      common/debhu/cdh,etadh0,etadh1,zdh,xxdy,yydh,zzdh,dhnue,ldh
! DBG 7/95 To store variables for pulse output
      double precision :: alfmlt, phmlt, cmxmlt
      double precision :: valfmlt(json), vphmlt(json), vcmxmlt(json)
      common/pualpha/alfmlt,phmlt,cmxmlt,valfmlt,vphmlt,vcmxmlt
! MHP 7/96 COMMON BLOCK ADDED FOR SOUND SPEED
      double precision :: adiabatic_index_gamma1(json)
      logical :: sound_speed_output_active
      common/sound/adiabatic_index_gamma1,sound_speed_output_active
      double precision :: rotational_energy_term(json)
      common/roten/rotational_energy_term
! MHP 06/02 COMMON BLOCK ADDED FOR DERIVATIVES OF
! CAPPA AND EPSILON
      double precision :: dlnkappa_dlnrho(json), dlnkappa_dlnt(json), &
           dlnepsilon_dlnrho(json), dlnepsilon_dlnt(json), &
           neutrino_loss_fraction(json)
      common/rotder/dlnkappa_dlnrho,dlnkappa_dlnt,dlnepsilon_dlnrho, &
           dlnepsilon_dlnt,neutrino_loss_fraction
      double precision :: accretion_specific_entropy, &
           envelope_specific_entropy, updated_mass_msun, delta_log_pressure, &
           delta_log_temperature
      common/masschg2/accretion_specific_entropy,envelope_specific_entropy, &
           updated_mass_msun,delta_log_pressure,delta_log_temperature

! JVS end

! --- locals ---
      double precision :: ion_fraction(3), energy_gen_component(6)
      double precision :: hf1(json), hf2(json), hr1(json), hr2(json), &
           hr3(json), hr4(json), hr5(json), hr6(json), hr7(json), &
           hr8(json), hr9(json), hr10(json), hr11(json), hr12(json), &
           hr13(json)
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
      double precision :: temperature, pressure, density
      double precision :: beta, beta_inverse, beta14, specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      double precision :: specific_heat_cp, adiabatic_gradient, &
           dlnrho_dlnt_dt, dlnrho_dlnp_dt, adiabatic_gradient_dt, &
           adiabatic_gradient_dp, specific_heat_cp_dt, specific_heat_cp_dp
      double precision :: opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
      double precision :: actual_gradient, radiative_gradient, &
           dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
           convective_velocity
      logical :: is_convective
      double precision :: qtemp
      double precision :: eq_p_val0, dqp_dr, dqp_dr0, dqp_dp, dqp_dp0
      double precision :: eq_t_val0, dqt_dr, dqt_dr0, dqt_dl0, dqt_dp, &
           dqt_dp0, dqt_dt, dqt_dt0
      double precision :: eq_r_val0, eq_r_val, dqr_dr, dqr_dr0, dqr_dp, &
           dqr_dp0, dqr_dt, dqr_dt0
      double precision :: eq_l_val0, eq_l_val, dql_dp, dql_dp0, dql_dt, &
           dql_dt0
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
      save

! 7/91 MHP
! ZERO OUT SOLAR NEUTRINO FLUXES.
! FLUXTOT = TOTAL FLUX OF EACH OF THE NEUTRINOS PRODUCED IN THE SUN
      if (lsnu) then
         do 5 j = 1,10
            flux_diag%neutrino_flux_total(j) = 0.0d0
   5     continue
      end if
! MHP 10/02 QFPR,QFTR NOT USED - OMIT
!      IF(.NOT.LROT) THEN
!       QFPR = 0.0D0
!       QFTR = 0.0D0
!      ENDIF
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
      do 10 j = 1,8
       luminosity_terms(j) = 0.0d0
   10 continue
      idt = 15
      do 15 j = 1,4
       idd(j) = 5
   15 continue
      do 30 im = 1,num_points
! SET UP LOCAL VARIABLES FOR CALLS TO BASIC PHYSICS ROUTINES
       zone_energy_luminosity = 0.0d0
       zone_log_mass = log_mass(im)
       zone_log_temperature = log_temperature(im)
       zone_log_pressure = log_pressure(im)
       zone_log_radius = log_radius(im)
       zone_luminosity_lsun = luminosity_lsun(im)
       hydrogen_fraction = composition(1,im)
       helium_fraction = composition(2,im)
       metal_fraction = composition(3,im)
       he3_fraction = composition(4,im)
       c12_fraction = composition(5,im)
       c13_fraction = composition(6,im)
       n14_fraction = composition(7,im)
       n15_fraction = composition(8,im)
       o16_fraction = composition(9,im)
       o17_fraction = composition(10,im)
       o18_fraction = composition(11,im)
! MHP 05/02 DEFINE THESE ALWAYS; THEY
! ARE PASSED TO THE SR ANYWAY.
!       IF(LEXCOM) THEN
          deuterium_fraction = composition(12,im)
          li6_fraction = composition(13,im)
          li7_fraction = composition(14,im)
          be9_fraction = composition(15,im)
!       ENDIF
       shell_index = im
       zone_log10_density = log10_density(im)
       pressure_rotation_factor = rotation_p_factor(im)
        temperature_rotation_factor = rotation_t_factor(im)
! YC   IF LMHD USE MHD EQUATION OF STATE.
         if (use_mhd_eos) then
            call meqos(zone_log_temperature, temperature, &
                 zone_log_pressure, pressure, zone_log10_density, density, &
                 hydrogen_fraction, metal_fraction, beta, beta_inverse, &
                 beta14, ion_fraction, specific_gas_constant, &
                 ion_mean_weight_inverse, electron_mean_weight_inverse, &
                 electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
                 specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
                 dlnrho_dlnp_dt, adiabatic_gradient_dt, &
                 adiabatic_gradient_dp, specific_heat_cp_dt, &
                 specific_heat_cp_dp)
         else
            if (ldh) then
               xxdy = composition(1,im)
               yydh = composition(2,im)+composition(4,im)
               zzdh = composition(3,im)
               zdh(1) = composition(5,im)+composition(6,im)
               zdh(2) = composition(7,im)+composition(8,im)
               zdh(3) = composition(9,im)+composition(10,im)+composition(11,im)
            end if
            call eqstat(zone_log_temperature, temperature, &
                 zone_log_pressure, pressure, zone_log10_density, density, &
                 hydrogen_fraction, metal_fraction, beta, beta_inverse, &
                 beta14, ion_fraction, specific_gas_constant, &
                 ion_mean_weight_inverse, electron_mean_weight_inverse, &
                 electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
                 specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
                 dlnrho_dlnp_dt, adiabatic_gradient_dt, &
                 adiabatic_gradient_dp, specific_heat_cp_dt, &
                 specific_heat_cp_dp, want_derivatives, in_atmosphere, &
                 saha_state)
         end if
! DBG 12/95 GET OPACITY
         call getopac(zone_log10_density, zone_log_temperature, &
              hydrogen_fraction, metal_fraction, opacity, log10_opacity, &
              dlnkap_dlnrho, dlnkap_dlnt, ion_fraction)
         iovim = im
         call tpgrad(zone_log_temperature, temperature, zone_log_pressure, &
              pressure, density, zone_log_radius, zone_log_mass, &
              zone_luminosity_lsun, opacity, dlnrho_dlnt, dlnrho_dlnp, &
              dlnkap_dlnt, dlnkap_dlnrho, specific_heat_cp, &
              actual_gradient, radiative_gradient, adiabatic_gradient, &
              dlnrho_dlnt_dt, dlnrho_dlnp_dt, adiabatic_gradient_dt, &
              adiabatic_gradient_dp, dgrad_dt_component, &
              dgrad_dp_component, dgrad_dr_component, specific_heat_cp_dt, &
              specific_heat_cp_dp, convective_velocity, want_derivatives, &
              is_convective, pressure_rotation_factor, &
              temperature_rotation_factor, log_teff)
       log10_density(im) = zone_log10_density
! COMPUTE DERIVATIVES
!       IF(LROT) THEN
!  CALCULATE D(LOG FP)/D(LOG R) AND D(LOG FT)/D(LOG R)
!            IF(IM.GT.1) THEN
!             IF(IM.LT.M) THEN
!              QFPR = (DLOG(FP(IM+1)) - DLOG(FP(IM-1)))/
!     *                 (CLN*(HR(IM+1) - HR(IM-1)))
!              QFTR = (DLOG(FT(IM+1)) - DLOG(FT(IM-1)))/
!     *                 (CLN*(HR(IM+1) - HR(IM-1)))
!             ELSE
!              QFPR = (DLOG(FP(M)) - DLOG(FP(M-1)))/
!     *                 (CLN*(HR(M) - HR(M-1)))
!              QFTR = (DLOG(FT(M)) - DLOG(FT(M-1)))/
!     *                 (CLN*(HR(M) - HR(M-1)))
!             ENDIF
!          ELSE
!             QFPR = (DLOG(FP(2)) - DLOG(FP(1)))/
!     *              (CLN*(HR(2) - HR(1)))
!             QFTR = (DLOG(FT(2)) - DLOG(FT(1)))/
!     *              (CLN*(HR(2) - HR(1)))
!          ENDIF
!       ENDIF
       qtemp = c4pil + zone_log_radius + zone_log_radius + zone_log_radius
       eq_r_val =+dexp(ln10*(zone_log_mass - zone_log10_density - qtemp))
       dqr_dr = - eq_r_val - eq_r_val - eq_r_val
       dqr_dp = -eq_r_val*dlnrho_dlnp
       dqr_dt = -eq_r_val*dlnrho_dlnt
       pulse_diag%qp =-dexp(ln10*(cgl + zone_log_mass + zone_log_mass - &
            zone_log_pressure - qtemp - zone_log_radius ))*rotation_p_factor(im)
!       QPR = -QP - QP - QP - QP*(1.0D0 - QFPR)
       dqp_dr = -pulse_diag%qp - pulse_diag%qp - pulse_diag%qp - pulse_diag%qp
       dqp_dp = -pulse_diag%qp
       convective_flag(im) = is_convective
       pulse_diag%qt = actual_gradient*pulse_diag%qp
       dqt_dr = -pulse_diag%qt - pulse_diag%qt - pulse_diag%qt - pulse_diag%qt
!       QTR = -QT - QT - QT - QT*(1.0D0 - QFTR)
       if (.not.is_convective) then
! TEMPERATURE GRADIENT IS RADIATIVE
          pulse_diag%qtl = clni*pulse_diag%qt/zone_luminosity_lsun
          dqt_dp = pulse_diag%qt*dlnkap_dlnrho*dlnrho_dlnp
          dqt_dt = pulse_diag%qt*(-4.0d0 + dlnkap_dlnt + dlnkap_dlnrho*dlnrho_dlnt)
       else
! TEMPERATURE GRADIENT IS CONVECTIVE
          pulse_diag%qtl = 0.0d0
          dqt_dp = pulse_diag%qt*(-1.0d0 + dgrad_dp_component)
          dqt_dt = pulse_diag%qt*dgrad_dt_component
          dqt_dr = dqt_dr + pulse_diag%qt*dgrad_dr_component
       end if
       eq_l_val = 0.0d0
       dql_dt = 0.0d0
       dql_dp = 0.0d0
       if (zone_log_temperature.gt.tcut(1)) then
! SET UP NUCLEAR ENERGY TERMS
            call engeb(pp_chain_gen, he3he4_be7_electron_gen, &
                 he3he4_be7_proton_gen, cno_gen, triple_alpha_gen, &
                 zone_dlnepsilon_dlnrho, zone_dlnepsilon_dlnt, &
                 total_energy_gen, zone_log10_density, &
                 zone_log_temperature, hydrogen_fraction, helium_fraction, &
                 he3_fraction, c12_fraction, c13_fraction, n14_fraction, &
                 o16_fraction, o18_fraction, deuterium_fraction, &
                 shell_index, hr1, hr2, hr3, hr4, hr5, hr6, hr7, &
                 hr8, hr9, hr10, hr11, hr12, hr13, hf1, hf2)
            energy_gen_rate = total_energy_gen
            energy_gen_component(1) = pp_chain_gen
            energy_gen_component(2) = he3he4_be7_electron_gen
            energy_gen_component(3) = he3he4_be7_proton_gen
            energy_gen_component(4) = cno_gen
            energy_gen_component(5) = triple_alpha_gen
            energy_gen_component(6) = engeb_diag%neutrino_loss_rate
            alpha_capture_energy_zone = engeb_diag%alpha_capture_energy
! 7/91 MHP
! CONVERT NEUTRINO FLUX RATES (UNITS 10**10 ERGS PER GM)
! TO UNITS OF 10**10 ERGS BY MULTIPLYING BY THE MASS.
            if (lsnu) then
               do 17 j = 1,10
                  flux_diag%neutrino_flux_total(j) = flux_diag%neutrino_flux_total(j) + &
                       flux_diag%neutrino_flux(j)*shell_mass(im)
 17            continue
            end if
            do 20 j = 1,6
               luminosity_terms(j) = luminosity_terms(j) + &
                    (shell_mass(im)/solar_luminosity_cgs)* &
                    energy_gen_component(j)
               zone_energy_luminosity = zone_energy_luminosity + &
                    (shell_mass(im)/solar_luminosity_cgs)* &
                    energy_gen_component(j)
 20         continue
! JVS 10/11 Calculate the He3+He3 and sum of He3+He3 and He3+He4 luminosity
            engeb_diag%he3_he3_rate_placeholder(im) = (shell_mass(im)/ &
                 solar_luminosity_cgs)*engeb_diag%he3_luminosity_placeholder
            engeb_diag%he3_he4_rate_placeholder(im) = (shell_mass(im)/ &
                 solar_luminosity_cgs)*engeb_diag%he3_total_placeholder
! JVS end
            luminosity_terms(8)=luminosity_terms(8)+(shell_mass(im)/ &
                 solar_luminosity_cgs)*alpha_capture_energy_zone
            zone_energy_luminosity = zone_energy_luminosity + &
                 (shell_mass(im)/solar_luminosity_cgs)* &
                 alpha_capture_energy_zone
            eq_l_val = energy_gen_rate
            dql_dt = dql_dt + zone_dlnepsilon_dlnt + &
                 zone_dlnepsilon_dlnrho*dlnrho_dlnt
            dql_dp = dql_dp + zone_dlnepsilon_dlnrho*dlnrho_dlnp
         end if
         if (compute_entropy_term) then
! SET UP ENTROPY TERMS
            zone_dt = delta_time_inv
            if (use_mass_accretion.and.mass_accretion_rate.gt.0.0d0) then
               if (im.ge.envelope_zone_index) then
                  zone_log_temperature_delta = log_temperature_delta(im)+ &
                       delta_log_temperature
                  zone_log_pressure_delta = log_pressure_delta(im)+ &
                       delta_log_pressure
               else
                  zone_log_temperature_delta = log_temperature_delta(im)
                  zone_log_pressure_delta = log_pressure_delta(im)
               end if
            else
               zone_log_temperature_delta = log_temperature_delta(im)
               zone_log_pressure_delta = log_pressure_delta(im)
            end if
            if (composition(1,im).gt.0.01d0 .and. delta_time.lt.one_year_sec) &
                 zone_dt = one_year_sec_inv
            entropy_term1 = pressure*dlnrho_dlnt/density
            entropy_term2 = entropy_term1/adiabatic_gradient
            entropy_term = (entropy_term2*zone_log_temperature_delta - &
                 entropy_term1*zone_log_pressure_delta)*ln10
            entropy_term3 = entropy_term2*ln10*zone_log_temperature_delta
!            ENTR = (ENTR2*HTT(IM) - ENTR1*HPP(IM))*CLN
!            ENTR3 = ENTR2*CLN*HTT(IM)
            egrav = zone_dt*entropy_term
            gravitational_luminosity(im) = egrav
            luminosity_terms(7) = luminosity_terms(7) + (shell_mass(im)/ &
                 solar_luminosity_cgs)*egrav
            eq_l_val = eq_l_val + egrav
            dql_dp = dql_dp + zone_dt*(entropy_term*(1.0d0-dlnrho_dlnp+ &
                 dlnrho_dlnp_dt)-entropy_term1 - entropy_term3* &
                 adiabatic_gradient_dp)
            dql_dt = dql_dt + zone_dt*(entropy_term*(-dlnrho_dlnt+ &
                 dlnrho_dlnt_dt) + entropy_term2 - entropy_term3* &
                 adiabatic_gradient_dt)
! 7/92 INCLUDE CHANGE IN ROTATIONAL KINETIC ENERGY IN ENERGY EQUATION.
            if (rotation_active) then
               rotational_energy_term(im) = zone_dt*(kinetic_energy_rot(im)- &
                    kinetic_energy_rot_old(im))/shell_mass(im)
               eq_l_val = eq_l_val - rotational_energy_term(im)
            end if
! ADD CHANGE IN ENTROPY FROM ACCRETED MATERIAL
!            IF(LMDOT.AND.DMDT0.GT.0.0D0)THEN
!               IF(IM.GE.JENV)THEN
!                  QACC = - T*SCEN*DMDT0/CSECYR/SMASS0
!                  WRITE(*,*)QL,QACC
!                  QL = QL + QACC
!               ENDIF
!            ENDIF
         end if
         cccql = ln_solar_luminosity*mass_weight_ln(im)
         eq_l_val = cccql*eq_l_val
         dql_dp = cccql*dql_dp
         dql_dt = cccql*dql_dt
         if (im.gt.1) then
! REDUCE MATRIX FOR PAIR OF POINTS (IM-1,IM)
            im1 = im
            call reduce(im1,elim_coeff,elim_rhs,luminosity_lsun,max_residual, &
                 log_pressure,log_radius,log_mass,log_temperature, &
                 eq_p_val0,pulse_diag%qp,dqp_dr0,dqp_dr,dqp_dp0,dqp_dp,eq_t_val0,pulse_diag%qt, &
                 dqt_dr0,dqt_dr,dqt_dl0,pulse_diag%qtl,dqt_dp0,dqt_dp,dqt_dt0,dqt_dt, &
                 eq_r_val0,eq_r_val,dqr_dr0,dqr_dr,dqr_dp0,dqr_dp,dqr_dt0, &
                 dqr_dt,eq_l_val0,eq_l_val,dql_dp0,dql_dp,dql_dt0,dql_dt)
         else
! SETUP CENTRAL BOUNDARY CONDITIONS
            elim_coeff(3,1,1) = cc13*dlnrho_dlnp
            elim_coeff(3,2,1) = cc13*dlnrho_dlnt
            elim_rhs(3,1) = -cc13*(c4pi3l + zone_log10_density - &
                 zone_log_mass) - zone_log_radius
            elim_coeff(4,1,1) = -dql_dp
            elim_coeff(4,2,1) = -dql_dt
            elim_rhs(4,1) = clni*eq_l_val - zone_luminosity_lsun
         end if
         eq_p_val0 = pulse_diag%qp
         dqp_dr0 = dqp_dr
         dqp_dp0 = dqp_dp
         eq_t_val0 = pulse_diag%qt
         dqt_dr0 = dqt_dr
         dqt_dl0 = pulse_diag%qtl
         dqt_dp0 = dqt_dp
         dqt_dt0 = dqt_dt
         eq_r_val0 = eq_r_val
         dqr_dr0 = dqr_dr
         dqr_dp0 = dqr_dp
         dqr_dt0 = dqr_dt
         eq_l_val0 = eq_l_val
         dql_dp0 = dql_dp
         dql_dt0 = dql_dt
! MHP 02/12 REMOVED RESTRICTIONS ON WHERE INTERMEDIATE VARIABLES
! SUCH AS OPACITY ARE SAVED; PRIOR RESTRICTIONS WERE BASED ON OBSOLETE
! MEMORY RESTRICTIONS IN LEGACY CODE
!         IF(LMDOT.AND.DMDT0.GT.0.0D0)THEN
         shell_diag%svel(im) = convective_velocity
!         ENDIF
!  STORE VARIABLES FOR OUTPUT IN SCRIB2 IF MODEL IS TO BE PRINTED OUT
! DBG PULSE STORE VARIABLES FOR PULSATION OUPUT
! DBG 3/91 CHANGED TO ALWAYS EXECUTE THIS STUFF
!         LONG = MOD(MODEL,NPRT2).EQ.0 .OR. LROT
! MHP 10/02 LSHORT NOT USED, OMIT
!         LSHORT = .NOT.LONG .AND. MOD(MODEL,NPRT1).EQ.0
!  ZERO OUT NUCLEAR ENERGY TERMS IF T < NUCLEAR CUTOFF.
         if (log_temperature(im).lt.tcut(1)) then
            shell_diag%sesum(im) = 0.0d0
            shell_diag%seg(7,im) = gravitational_luminosity(im)
            do j = 1,6
               shell_diag%seg(j,im) = 0.0d0
           end do
         else
!         ELSE IF(LONG) THEN
!  LONG OUTPUT NEEDED
            shell_diag%sesum(im) = energy_gen_component(1)+energy_gen_component(2)+ &
                 energy_gen_component(3)+energy_gen_component(4)+ &
                 energy_gen_component(5)
            shell_diag%seg(6,im) = energy_gen_component(6)
            shell_diag%seg(7,im) = gravitational_luminosity(im)
            if (shell_diag%sesum(im).gt.1.0d-22) then
               energy_sum_inverse = 1.0d0/shell_diag%sesum(im)
            else
               energy_sum_inverse = 0.0d0
            end if
            do j = 1,5
               shell_diag%seg(j,im) = energy_gen_component(j)*energy_sum_inverse
              end do
!  SHORT OUTPUT ONLY
!         ELSE
!            SESUM(IM) = EG(1)+EG(2)+EG(3)+EG(4)+EG(5)
!            SEG(6,IM) = EG(6)
!            SEG(7,IM) = HHC(IM)
         end if
         shell_diag%sbeta(im) = beta
         shell_diag%seta(im) = electron_degeneracy_parameter
         shell_diag%locons(im) = conductive_opacity_flag
         shell_diag%so(im) = opacity
         shell_diag%del_grad(1,im) = radiative_gradient
         shell_diag%del_grad(2,im) = actual_gradient
         shell_diag%del_grad(3,im) = adiabatic_gradient
         do j = 1,3
            shell_diag%sfxion(j,im) = ion_fraction(j)
         end do
         shell_diag%svel(im) = convective_velocity
         shell_diag%scp(im) = specific_heat_cp
! MHP 02/12 COMMENTED CODE OUT, AS REPLICATED BELOW
!         IF(LSOUND) THEN
! MHP 7/96 CALCULATION OF GAMMA1 FROM GUENTHER 1995 P.C.
!            CHRH = 1.0D0/QDP
!            CHT = -CHRH*QDT
!            CV = QCP - EXP(CLN*(HP(IM)-HD(IM)-HT(IM)))*CHT**2/CHRH
!            GAM1(IM) = CHRH*QCP/CV
!            PQDP(IM) = QDP
!         ENDIF

! JVS 01/11 always want gamma:
            chi_rho = 1.0d0/dlnrho_dlnp
            chi_t = -chi_rho*dlnrho_dlnt
            specific_heat_cv = specific_heat_cp - exp(ln10*(log_pressure(im)- &
                 log10_density(im)-log_temperature(im)))*chi_t**2/chi_rho
            adiabatic_index_gamma1(im) = chi_rho*specific_heat_cp/ &
                 specific_heat_cv
            pulse_diag%pulse_dlnrho_dlnp(im) = dlnrho_dlnp
            pulse_diag%pulse_dlnrho_dlnt(im) = dlnrho_dlnt
! JVS END


         if (rotation_active) then
            dlnkappa_dlnrho(im) = dlnkap_dlnrho
            dlnkappa_dlnt(im) = dlnkap_dlnt
! MHP 10/02 variable index error
            if (shell_diag%sesum(im).gt.0.0d0) then
!            IF(SESUM(I).GT.0.0D0)THEN
!               ETOT = SESUM(I)
!               EGNEUT = SEG(6,I)+SEG(7,I)
               total_energy_sum = shell_diag%sesum(im)
               neutrino_and_grav_sum = shell_diag%seg(6,im)+shell_diag%seg(7,im)
               neutrino_loss_fraction(im) = (total_energy_sum - &
                    neutrino_and_grav_sum)/total_energy_sum
            else
               neutrino_loss_fraction(im) = 0.0d0
            end if
            dlnepsilon_dlnrho(im) = zone_dlnepsilon_dlnrho
            dlnepsilon_dlnt(im) = zone_dlnepsilon_dlnt
         end if
! DBG PULSE
! MHP 8/25 unconditional: previously gated on pulsation_output_active
! (the legacy path-length-triggered OPAL pulsation writer's flag), but
! the GYRE-format pulsation writer (io/write_gyre_pulse.f90) is
! triggered independently (by model-number interval, checked in
! wrtout.f90) and needs these populated for every converged model, not
! just ones flagged by the older mechanism. All source locals here
! (dlnrho_dlnp, dlnkap_dlnrho/dlnt, zone_dlnepsilon_dlnrho/dlnt, etc.)
! are already computed unconditionally above, so this is just always
! copying them into the pulse1/mixing-length output arrays -- no
! change to any existing output (.short/.track/.store/.pmod/.penv/
! .patm) values, since none of those read from these arrays.
         pulse_diag%pulse_dlnrho_dlnp(im) = dlnrho_dlnp
         pulse_diag%pulse_dlneps_dlnrho(im) = zone_dlnepsilon_dlnrho
         pulse_diag%pulse_dlneps_dlnt(im) = zone_dlnepsilon_dlnt
         pulse_diag%pulse_dlnkap_dlnrho(im) = dlnkap_dlnrho
         pulse_diag%pulse_dlnkap_dlnt(im) = dlnkap_dlnt
         pulse_diag%pulse_specific_heat(im) = specific_heat_cp
         pulse_diag%pulse_mean_molecular_weight(im) = specific_gas_constant
         pulse_diag%pulse_electron_mean_molecular_weight(im) = &
              electron_mean_weight_inverse
         pulse_diag%pulse_dlnrho_dlnt(im) = dlnrho_dlnt
         valfmlt(im) = alfmlt
         vphmlt(im) = phmlt
         vcmxmlt(im) = cmxmlt
 30   continue

      return
end subroutine coefft

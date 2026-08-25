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
! equation of state (via eos_lib's eos_get), opacity (via kap_lib's
! kap_get), and
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
! Three dummy arguments below (star%pulse%qt, star%pulse%qp, star%pulse%qtl -- see common/pulse2/) are
! simultaneously common-block storage: COMMON/PULSE2/ was first
! declared (unused, as placeholders) in an earlier-converted file
! (wrtmod.f90) with generic member names kept close to the original;
! this file is the first to actually USE those three slots, as the
! live temperature/pressure Henyey-equation scratch terms (they also
! double as the eq_t_val/eq_p_val/dqt_dl arguments passed to reduce),
! so per the project's COMMON-block-reuse rule they keep the
! wrtmod.f90 names (star%pulse%qt/star%pulse%qp/star%pulse%qtl) here too, despite those names no
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
     kinetic_energy_rot_old, envelope_zone_index, log_teff, ierr)

      use net_lib
      use star_info_lib, only: star, i_eps_grav, i_eps_neu, i_grad_actual, i_grad_ad, i_grad_rad, json
      use luout_lib
      use phys_const_lib
      use eos_lib
      use kap_lib
      use burn_lib
      implicit none

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
! 7/91 MHP
! ZERO OUT SOLAR NEUTRINO FLUXES.
! FLUXTOT = TOTAL FLUX OF EACH OF THE NEUTRINOS PRODUCED IN THE SUN
      integer, intent(out) :: ierr

      ierr = 0

      if (star%ctrl%lsnu) then
         do j = 1,10
            star%neutrino_flux_total(j) = 0.0d0
         end do
      end if
! MHP 10/02 QFPR,QFTR NOT USED - OMIT
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
       luminosity_terms(j) = 0.0d0
      end do
      do im = 1,num_points
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
         call eos_get(zone_log_temperature, temperature, &
              zone_log_pressure, pressure, zone_log10_density, density, &
              hydrogen_fraction, metal_fraction, beta, beta_inverse, &
              beta14, ion_fraction, specific_gas_constant, &
              ion_mean_weight_inverse, electron_mean_weight_inverse, &
              electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
              specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
              dlnrho_dlnp_dt, adiabatic_gradient_dt, &
              adiabatic_gradient_dp, specific_heat_cp_dt, &
              specific_heat_cp_dp, want_derivatives, in_atmosphere, &
              saha_state, composition_at_zone=composition(:,im))
! DBG 12/95 GET OPACITY
         call kap_get(zone_log10_density, zone_log_temperature, &
              hydrogen_fraction, metal_fraction, opacity, log10_opacity, &
              dlnkap_dlnrho, dlnkap_dlnt, ion_fraction)
         star%iovim = im
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
              temperature_rotation_factor, log_teff, ierr)
         if (ierr /= 0) return
       log10_density(im) = zone_log10_density
! COMPUTE DERIVATIVES
!       IF(LROT) THEN
!  CALCULATE D(LOG FP)/D(LOG R) AND D(LOG FT)/D(LOG R)
!     *                 (CLN*(HR(IM+1) - HR(IM-1)))
!              QFTR = (DLOG(FT(IM+1)) - DLOG(FT(IM-1)))/
!     *                 (CLN*(HR(IM+1) - HR(IM-1)))
!             ELSE
!              QFPR = (DLOG(FP(M)) - DLOG(FP(M-1)))/
!     *                 (CLN*(HR(M) - HR(M-1)))
!              QFTR = (DLOG(FT(M)) - DLOG(FT(M-1)))/
!     *                 (CLN*(HR(M) - HR(M-1)))
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
       star%pulse%qp =-dexp(ln10*(cgl + zone_log_mass + zone_log_mass - &
            zone_log_pressure - qtemp - zone_log_radius ))*rotation_p_factor(im)
!       QPR = -QP - QP - QP - QP*(1.0D0 - QFPR)
       dqp_dr = -star%pulse%qp - star%pulse%qp - star%pulse%qp - star%pulse%qp
       dqp_dp = -star%pulse%qp
       convective_flag(im) = is_convective
       star%pulse%qt = actual_gradient*star%pulse%qp
       dqt_dr = -star%pulse%qt - star%pulse%qt - star%pulse%qt - star%pulse%qt
!       QTR = -QT - QT - QT - QT*(1.0D0 - QFTR)
       if (.not.is_convective) then
! TEMPERATURE GRADIENT IS RADIATIVE
          star%pulse%qtl = clni*star%pulse%qt/zone_luminosity_lsun
          dqt_dp = star%pulse%qt*dlnkap_dlnrho*dlnrho_dlnp
          dqt_dt = star%pulse%qt*(-4.0d0 + dlnkap_dlnt + dlnkap_dlnrho*dlnrho_dlnt)
       else
! TEMPERATURE GRADIENT IS CONVECTIVE
          star%pulse%qtl = 0.0d0
          dqt_dp = star%pulse%qt*(-1.0d0 + dgrad_dp_component)
          dqt_dt = star%pulse%qt*dgrad_dt_component
          dqt_dr = dqt_dr + star%pulse%qt*dgrad_dr_component
       end if
       eq_l_val = 0.0d0
       dql_dt = 0.0d0
       dql_dp = 0.0d0
       if (zone_log_temperature.gt.star%ctrl%tcut(1)) then
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
            energy_gen_component(6) = star%neutrino_loss_rate
            alpha_capture_energy_zone = star%alpha_capture_energy
! 7/91 MHP
! CONVERT NEUTRINO FLUX RATES (UNITS 10**10 ERGS PER GM)
! TO UNITS OF 10**10 ERGS BY MULTIPLYING BY THE MASS.
            if (star%ctrl%lsnu) then
               do j = 1,10
                  star%neutrino_flux_total(j) = star%neutrino_flux_total(j) + &
                       star%neutrino_flux(j)*shell_mass(im)
               end do
            end if
            do j = 1,6
               luminosity_terms(j) = luminosity_terms(j) + &
                    (shell_mass(im)/star%solar_luminosity_cgs)* &
                    energy_gen_component(j)
               zone_energy_luminosity = zone_energy_luminosity + &
                    (shell_mass(im)/star%solar_luminosity_cgs)* &
                    energy_gen_component(j)
            end do
! JVS 10/11 Calculate the He3+He3 and sum of He3+He3 and He3+He4 luminosity
            star%he3_he3_luminosity_zone(im) = (shell_mass(im)/ &
                 star%solar_luminosity_cgs)*star%he3_he3_energy_rate
            star%he3_burning_luminosity_zone(im) = (shell_mass(im)/ &
                 star%solar_luminosity_cgs)*star%he3_burning_energy_rate
! JVS end
            luminosity_terms(8)=luminosity_terms(8)+(shell_mass(im)/ &
                 star%solar_luminosity_cgs)*alpha_capture_energy_zone
            zone_energy_luminosity = zone_energy_luminosity + &
                 (shell_mass(im)/star%solar_luminosity_cgs)* &
                 alpha_capture_energy_zone
            eq_l_val = energy_gen_rate
            dql_dt = dql_dt + zone_dlnepsilon_dlnt + &
                 zone_dlnepsilon_dlnrho*dlnrho_dlnt
            dql_dp = dql_dp + zone_dlnepsilon_dlnrho*dlnrho_dlnp
         end if
         if (compute_entropy_term) then
! SET UP ENTROPY TERMS
            zone_dt = delta_time_inv
            if (star%job%use_mass_accretion.and.star%ctrl%mass_accretion_rate.gt.0.0d0) then
               if (im.ge.envelope_zone_index) then
                  zone_log_temperature_delta = log_temperature_delta(im)+ &
                       star%rot%delta_log_temperature
                  zone_log_pressure_delta = log_pressure_delta(im)+ &
                       star%rot%delta_log_pressure
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
                 star%solar_luminosity_cgs)*egrav
            eq_l_val = eq_l_val + egrav
            dql_dp = dql_dp + zone_dt*(entropy_term*(1.0d0-dlnrho_dlnp+ &
                 dlnrho_dlnp_dt)-entropy_term1 - entropy_term3* &
                 adiabatic_gradient_dp)
            dql_dt = dql_dt + zone_dt*(entropy_term*(-dlnrho_dlnt+ &
                 dlnrho_dlnt_dt) + entropy_term2 - entropy_term3* &
                 adiabatic_gradient_dt)
! 7/92 INCLUDE CHANGE IN ROTATIONAL KINETIC ENERGY IN ENERGY EQUATION.
            if (star%job%rotation_active) then
               star%rot%rotational_energy_term(im) = zone_dt*(kinetic_energy_rot(im)- &
                    kinetic_energy_rot_old(im))/shell_mass(im)
               eq_l_val = eq_l_val - star%rot%rotational_energy_term(im)
            end if
! ADD CHANGE IN ENTROPY FROM ACCRETED MATERIAL
         end if
         cccql = star%ln_solar_luminosity*mass_weight_ln(im)
         eq_l_val = cccql*eq_l_val
         dql_dp = cccql*dql_dp
         dql_dt = cccql*dql_dt
         if (im.gt.1) then
! REDUCE MATRIX FOR PAIR OF POINTS (IM-1,IM)
            im1 = im
            call reduce(im1,elim_coeff,elim_rhs,luminosity_lsun,max_residual, &
                 log_pressure,log_radius,log_mass,log_temperature, &
                 eq_p_val0,star%pulse%qp,dqp_dr0,dqp_dr,dqp_dp0,dqp_dp,eq_t_val0,star%pulse%qt, &
                 dqt_dr0,dqt_dr,dqt_dl0,star%pulse%qtl,dqt_dp0,dqt_dp,dqt_dt0,dqt_dt, &
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
         eq_p_val0 = star%pulse%qp
         dqp_dr0 = dqp_dr
         dqp_dp0 = dqp_dp
         eq_t_val0 = star%pulse%qt
         dqt_dr0 = dqt_dr
         dqt_dl0 = star%pulse%qtl
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
         star%svel(im) = convective_velocity
!         ENDIF
!  STORE VARIABLES FOR OUTPUT IN SCRIB2 IF MODEL IS TO BE PRINTED OUT
! DBG PULSE STORE VARIABLES FOR PULSATION OUPUT
! DBG 3/91 CHANGED TO ALWAYS EXECUTE THIS STUFF
!         LONG = MOD(MODEL,NPRT2).EQ.0 .OR. LROT
! MHP 10/02 LSHORT NOT USED, OMIT
!         LSHORT = .NOT.LONG .AND. MOD(MODEL,NPRT1).EQ.0
!  ZERO OUT NUCLEAR ENERGY TERMS IF T < NUCLEAR CUTOFF.
         if (log_temperature(im).lt.star%ctrl%tcut(1)) then
            star%sesum(im) = 0.0d0
            star%seg(i_eps_grav,im) = gravitational_luminosity(im)
            do j = 1,6
               star%seg(j,im) = 0.0d0
           end do
         else
!         ELSE IF(LONG) THEN
!  LONG OUTPUT NEEDED
            star%sesum(im) = energy_gen_component(1)+energy_gen_component(2)+ &
                 energy_gen_component(3)+energy_gen_component(4)+ &
                 energy_gen_component(5)
            star%seg(i_eps_neu,im) = energy_gen_component(6)
            star%seg(i_eps_grav,im) = gravitational_luminosity(im)
            if (star%sesum(im).gt.1.0d-22) then
               energy_sum_inverse = 1.0d0/star%sesum(im)
            else
               energy_sum_inverse = 0.0d0
            end if
            do j = 1,5
               star%seg(j,im) = energy_gen_component(j)*energy_sum_inverse
              end do
!  SHORT OUTPUT ONLY
         end if
         star%sbeta(im) = beta
         star%seta(im) = electron_degeneracy_parameter
         star%locons(im) = conductive_opacity_flag
         star%so(im) = opacity
         star%del_grad(i_grad_rad,im) = radiative_gradient
         star%del_grad(i_grad_actual,im) = actual_gradient
         star%del_grad(i_grad_ad,im) = adiabatic_gradient
         do j = 1,3
            star%sfxion(j,im) = ion_fraction(j)
         end do
         star%svel(im) = convective_velocity
         star%scp(im) = specific_heat_cp
! MHP 02/12 COMMENTED CODE OUT, AS REPLICATED BELOW
!         IF(LSOUND) THEN
! MHP 7/96 CALCULATION OF GAMMA1 FROM GUENTHER 1995 P.C.

! JVS 01/11 always want gamma:
            chi_rho = 1.0d0/dlnrho_dlnp
            chi_t = -chi_rho*dlnrho_dlnt
            specific_heat_cv = specific_heat_cp - exp(ln10*(log_pressure(im)- &
                 log10_density(im)-log_temperature(im)))*chi_t**2/chi_rho
            star%adiabatic_index_gamma1(im) = chi_rho*specific_heat_cp/ &
                 specific_heat_cv
            star%pulse_dlnrho_dlnp(im) = dlnrho_dlnp
            star%pulse_dlnrho_dlnt(im) = dlnrho_dlnt
! JVS END


         if (star%job%rotation_active) then
            star%rot%dlnkappa_dlnrho(im) = dlnkap_dlnrho
            star%rot%dlnkappa_dlnt(im) = dlnkap_dlnt
! MHP 10/02 variable index error
            if (star%sesum(im).gt.0.0d0) then
               total_energy_sum = star%sesum(im)
               neutrino_and_grav_sum = star%seg(i_eps_neu,im)+star%seg(i_eps_grav,im)
               star%rot%neutrino_loss_fraction(im) = (total_energy_sum - &
                    neutrino_and_grav_sum)/total_energy_sum
            else
               star%rot%neutrino_loss_fraction(im) = 0.0d0
            end if
            star%rot%dlnepsilon_dlnrho(im) = zone_dlnepsilon_dlnrho
            star%rot%dlnepsilon_dlnt(im) = zone_dlnepsilon_dlnt
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
         star%pulse_dlnrho_dlnp(im) = dlnrho_dlnp
         star%pulse_dlneps_dlnrho(im) = zone_dlnepsilon_dlnrho
         star%pulse_dlneps_dlnt(im) = zone_dlnepsilon_dlnt
         star%pulse_dlnkap_dlnrho(im) = dlnkap_dlnrho
         star%pulse_dlnkap_dlnt(im) = dlnkap_dlnt
         star%pulse_specific_heat(im) = specific_heat_cp
         star%pulse_mean_molecular_weight(im) = specific_gas_constant
         star%pulse_electron_mean_molecular_weight(im) = &
              electron_mean_weight_inverse
         star%pulse_dlnrho_dlnt(im) = dlnrho_dlnt
         star%rot%valfmlt(im) = star%rot%alfmlt
         star%rot%vphmlt(im) = star%rot%phmlt
         star%rot%vcmxmlt(im) = star%rot%cmxmlt
      end do

      return
end subroutine coefft

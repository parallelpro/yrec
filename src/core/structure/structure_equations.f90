!----------------------------------------------------------------------
! structure_equations
!----------------------------------------------------------------------
! The structure equations of the Henyey solve. compute_henyey_coefficients
! evaluates the physics at every shell, forms the four equations' terms
! and their partials, and hands each pair of shells to henyey_eliminate
! (henyey_elimination.f90) as it goes.
module structure_equations
      implicit none
      private
      public :: compute_henyey_coefficients

contains

!----------------------------------------------------------------------
! compute_henyey_coefficients (formerly coefft)
!----------------------------------------------------------------------
! Linearizes the structure equations at every shell, from the centre
! out, and eliminates as it goes. At each shell im it calls the equation
! of state (eos_lib's eos_get), the opacity (kap_lib's kap_get) and the
! temperature gradients (compute_temperature_gradients); above the
! burning cutoff the nuclear energy generation (burn_step's
! compute_zone_energy); and, when the step has a duration, the
! gravitational (entropy-change) energy term. It forms the four terms q
! and their partials in a henyey_shell_terms record (cur; prev holds the
! shell below) and calls henyey_eliminate for the pair (im-1, im), or
! sets the central boundary conditions at im = 1. The results are
! star%elim_coeff and star%elim_rhs. It also stores per-zone physics for
! output: star%logRho, eps_total and eps_channels, gradr/gradT/grada,
! adiabatic_index_gamma1, the pulse_* arrays, the luminosity channels and
! the neutrino fluxes.
!
! Variables y = (log P, log T, log R, L/Lsun), coordinate log m (log10,
! cgs). Each term is q = dy/dlog m:
!   qp = dlog P/dlog m = -G m^2/(4 pi r^4 P) fp
!   qt = dlog T/dlog m = grad qp,   grad = dln T/dln P
!   qr = dlog R/dlog m = m/(4 pi r^3 rho)
!   ql = dL/dlog m     = (ln10 m/Lsun) eps,
!        eps = eps_nuc + eps_grav - d(e_rot)/dt   (erg/g/s)
! with fp = star%fp_rot the rotational correction to hydrostatic
! equilibrium. The partials q_dp, q_dt, q_dr are dq/dln P, dq/dln T,
! dq/dln R (natural log) at fixed values of the other variables, and
! qt_dl = (dqt/dL)/ln10. Since ql is linear in eps, ql_dp and ql_dt
! carry the absolute derivatives d eps/dln P and d eps/dln T (times
! ln10 m/Lsun). The full contract and the difference scheme are in
! henyey_elimination.f90.
!
! pt_scr%qp, pt_scr%qt and pt_scr%qtl are shared point-scratch slots
! (point_scratch_lib; the envelope integrands use them for P, T and
! log T). Here they hold qp, qt and qt_dl, which are copied into cur
! before the elimination.
subroutine compute_henyey_coefficients(delta_time, dlnrho_dlnt, dlnrho_dlnp, &
     saha_state, envelope_zone_index, ierr)
      use convection_lib, only: compute_temperature_gradients
      use henyey_elimination

      use net_lib, only: nuclear_energy_result
      use star_info_lib, only: star, i_eps_grav, i_eps_neu, &
           i_h1, i_he4, i_metals, i_he3, i_c12, i_c13, i_n14, i_n15, &
           i_o16, i_o17, i_o18, i_h2, i_li6, i_li7, i_be9, &
           i_lum_he_c, n_lum_channels, n_nu_fluxes
      use point_scratch_lib
      use physical_constants
      use eos_lib
      use kap_lib
      use burn_step, only: compute_zone_energy
      use math_lib
      implicit none

      double precision, intent(in) :: delta_time
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
! The eos and opacity results at the shell (as in shell_physics). The
! dummies dlnrho_dlnt and dlnrho_dlnp are this routine's outputs to the
! driver; they are copied into eos_res before each eos call and back out
! after it, so each call starts from the previous shell's values.
      type(eos_result) :: eos_res
      double precision :: kap_res(num_kap_results)
! Every mesh point is interior (not atmosphere), and the eos and
! gradient derivatives are always wanted.
      logical, parameter :: want_derivatives = .true., in_atmosphere = .false.
      double precision :: actual_gradient, radiative_gradient, &
           dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
           convective_velocity
      logical :: is_convective
      double precision :: qtemp
! The equation terms at the current shell (cur) and, carried from the
! previous pass of the loop, at the shell below (prev). qp, qt and qt_dl
! are formed in pt_scr and copied into cur just before the elimination.
      type(henyey_shell_terms) :: prev, cur
! The shell's energy record from compute_zone_energy, unpacked right
! after the call into energy_gen_rate, energy_gen_component and
! zone_deps_dlnrho/dlnt.
      type(nuclear_energy_result) :: burn
      double precision :: zone_deps_dlnrho, zone_deps_dlnt
! The change of the rotational kinetic energy per unit mass and time.
      double precision :: rotational_energy_term
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
! The locals of the disabled rotation block in store_zone_thermodynamics:
!      double precision :: total_energy_sum, neutrino_and_grav_sum
      integer, intent(out) :: ierr

      ierr = 0

      call start_energy_bookkeeping()
      do im = 1,star%nz
       call gather_shell_inputs()
       call evaluate_eos_opacity_gradients()
       if (ierr /= 0) return
       call form_radius_pressure_temperature_terms()
       call form_nuclear_energy_terms()
       call form_gravitational_energy_terms()
       call eliminate_shell_pair()
       call store_zone_energy_channels()
       call store_zone_thermodynamics()
       call store_pulse_derivatives()
      end do

      return

contains

! ---------------------------------------------------------------
! Zero the solar-neutrino fluxes (star%neutrino_flux_total(j), the total
! flux of neutrino species j made in the star) and the luminosity
! channels, and set up the inverse time step of the entropy
! term when the step has a duration.
subroutine start_energy_bookkeeping()
      if (star%ctrl%calc_neutrinos) then
         do j = 1,n_nu_fluxes
            star%neutrino_flux_total(j) = 0.0d0
         end do
      end if
      compute_entropy_term = delta_time.gt.0.0d0
      if (compute_entropy_term) then
       delta_time_inv = 1.0d0/delta_time
! One year in seconds and its inverse: the entropy term's floor on the
! time step (form_gravitational_energy_terms).
       one_year_sec = 3.1558d7
       one_year_sec_inv = 3.1688d-8
      end if
      do j = 1,n_lum_channels
       star%luminosity_breakdown(j) = 0.0d0
      end do
end subroutine start_energy_bookkeeping

! ---------------------------------------------------------------
! Copy shell im's mass, temperature, pressure, radius, luminosity,
! composition, density and rotation factors into the locals.
subroutine gather_shell_inputs()
       zone_energy_luminosity = 0.0d0
       zone_log_mass = star%log_mass(im)
       zone_log_temperature = star%logT(im)
       zone_log_pressure = star%logP(im)
       zone_log_radius = star%logR(im)
       zone_luminosity_lsun = star%luminosity_lsun(im)
       hydrogen_fraction = star%xa(i_h1,im)
       helium_fraction = star%xa(i_he4,im)
       metal_fraction = star%xa(i_metals,im)
       he3_fraction = star%xa(i_he3,im)
       c12_fraction = star%xa(i_c12,im)
       c13_fraction = star%xa(i_c13,im)
       n14_fraction = star%xa(i_n14,im)
       n15_fraction = star%xa(i_n15,im)
       o16_fraction = star%xa(i_o16,im)
       o17_fraction = star%xa(i_o17,im)
       o18_fraction = star%xa(i_o18,im)
! The extended-composition species are always defined. Of these copies
! only hydrogen_fraction and metal_fraction are read (by the eos and the
! opacity); compute_zone_energy takes star%xa(:,im) itself.
       deuterium_fraction = star%xa(i_h2,im)
       li6_fraction = star%xa(i_li6,im)
       li7_fraction = star%xa(i_li7,im)
       be9_fraction = star%xa(i_be9,im)
       shell_index = im
       zone_log10_density = star%logRho(im)
       pressure_rotation_factor = star%fp_rot(im)
        temperature_rotation_factor = star%ft_rot(im)
end subroutine gather_shell_inputs

! ---------------------------------------------------------------
! The equation of state, the opacity and the temperature gradients at
! shell im; a failure leaves ierr nonzero.
subroutine evaluate_eos_opacity_gradients()
         eos_res%log10_density = zone_log10_density
         eos_res%dlnrho_dlnt = dlnrho_dlnt
         eos_res%dlnrho_dlnp = dlnrho_dlnp
         call eos_get(zone_log_temperature, zone_log_pressure, &
              hydrogen_fraction, metal_fraction, eos_res, &
              want_derivatives, in_atmosphere, saha_state, &
              composition_at_zone=star%xa(:,im), ierr=ierr)
         if (ierr /= 0) return
         dlnrho_dlnt = eos_res%dlnrho_dlnt
         dlnrho_dlnp = eos_res%dlnrho_dlnp
! The density the eos returns at (P, T) is the one everything below uses
! (qr, the energy generation, the central boundary condition).
         zone_log10_density = eos_res%log10_density
! The opacity at that density.
         call kap_get(eos_res%log10_density, zone_log_temperature, &
              hydrogen_fraction, metal_fraction, kap_res, &
              eos_res%ion_fraction, ierr=ierr)
         if (ierr /= 0) return
         call compute_temperature_gradients(zone_log_temperature, zone_log_pressure, &
              eos_res, kap_res, zone_log_radius, zone_log_mass, &
              zone_luminosity_lsun, actual_gradient, radiative_gradient, &
              dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
              convective_velocity, want_derivatives, is_convective, &
              pressure_rotation_factor, temperature_rotation_factor, &
              star%log_Teff, &
              im >= star%iov1 .and. im <= star%iov2, ierr)
         if (ierr /= 0) return
end subroutine evaluate_eos_opacity_gradients

! ---------------------------------------------------------------
! Store the shell's density; form the radius, pressure and temperature
! equations' terms and partials (radiative or convective gradient).
subroutine form_radius_pressure_temperature_terms()
       star%logRho(im) = eos_res%log10_density
! The radius equation, qr = m/(4 pi r^3 rho) (qtemp = log(4 pi r^3)):
!   qr_dr = -3 qr, qr_dp = -qr dln rho/dln P, qr_dt = -qr dln rho/dln T.
       qtemp = c4pil + zone_log_radius + zone_log_radius + zone_log_radius
       cur%qr =+exp(ln10*(zone_log_mass - zone_log10_density - qtemp))
       cur%qr_dr = - cur%qr - cur%qr - cur%qr
       cur%qr_dp = -cur%qr*dlnrho_dlnp
       cur%qr_dt = -cur%qr*dlnrho_dlnt
! The pressure equation, qp = -G m^2/(4 pi r^4 P) fp:
!   qp_dr = -4 qp, qp_dp = -qp.
! The rotation factors fp and ft are held fixed. A legacy variant
! subtracted their radius derivatives (legacy QFPR, QFTR):
! qp_dr = -qp (4 - QFPR), qt_dr = -qt (4 - QFTR).
       pt_scr%qp =-exp(ln10*(cgl + zone_log_mass + zone_log_mass - &
            zone_log_pressure - qtemp - zone_log_radius ))*star%fp_rot(im)
       cur%qp_dr = -pt_scr%qp - pt_scr%qp - pt_scr%qp - pt_scr%qp
       cur%qp_dp = -pt_scr%qp
       star%convective_flag(im) = is_convective
! The temperature equation, qt = grad qp, with grad = dln T/dln P from
! compute_temperature_gradients; qt_dr = -4 qt comes from qp.
       pt_scr%qt = actual_gradient*pt_scr%qp
       cur%qt_dr = -pt_scr%qt - pt_scr%qt - pt_scr%qt - pt_scr%qt
       if (.not.is_convective) then
! Radiative: qt is proportional to kappa L m/(r^4 T^4) (times ft), so
!   qt_dl = qt/(L ln10),
!   qt_dp = qt (dln kappa/dln rho)(dln rho/dln P),
!   qt_dt = qt (-4 + dln kappa/dln T + (dln kappa/dln rho)(dln rho/dln T)).
          pt_scr%qtl = clni*pt_scr%qt/zone_luminosity_lsun
          cur%qt_dp = pt_scr%qt*kap_res(i_dlnkap_dlnrho)*dlnrho_dlnp
          cur%qt_dt = pt_scr%qt*(-4.0d0 + kap_res(i_dlnkap_dlnt) + kap_res(i_dlnkap_dlnrho)*dlnrho_dlnt)
       else
! Convective: with dgrad_dX_component = dln grad/dln X from the
! convection routine,
!   qt_dp = qt (-1 + dln grad/dln P)   (-1 from qp's 1/P),
!   qt_dt = qt dln grad/dln T,  qt_dr = qt (-4 + dln grad/dln r).
! The convective gradient's dependence on L is not linearized (qt_dl = 0).
          pt_scr%qtl = 0.0d0
          cur%qt_dp = pt_scr%qt*(-1.0d0 + dgrad_dp_component)
          cur%qt_dt = pt_scr%qt*dgrad_dt_component
          cur%qt_dr = cur%qt_dr + pt_scr%qt*dgrad_dr_component
       end if
end subroutine form_radius_pressure_temperature_terms

! ---------------------------------------------------------------
! The luminosity equation's nuclear term and partials above the burning
! cutoff (log T > nuclear_logT_cutoffs(1)), with the shell's share of the
! luminosity channels and neutrino fluxes. Before the ln10 m/Lsun scaling
! of eliminate_shell_pair,
!   ql = eps_nuc,
!   ql_dt = d eps/dln T + (d eps/dln rho)(dln rho/dln T),
!   ql_dp = (d eps/dln rho)(dln rho/dln P).
! Below the cutoff the burn locals (zone_deps_*, energy_gen_component)
! are zero. (The F77 coefft.f left them at the last burning shell's
! values, which reached the pulse output's d eps columns of every cooler
! shell, and eps_total of a shell exactly at the cutoff; fixed 2026-10.)
subroutine form_nuclear_energy_terms()
       cur%ql = 0.0d0
       cur%ql_dt = 0.0d0
       cur%ql_dp = 0.0d0
       zone_deps_dlnrho = 0.0d0
       zone_deps_dlnt = 0.0d0
       energy_gen_component = 0.0d0
       if (zone_log_temperature.gt.star%ctrl%nuclear_logT_cutoffs(1)) then
            call compute_zone_energy(zone_log10_density, &
                 zone_log_temperature, star%xa(:,im), shell_index, burn)
            zone_deps_dlnrho = burn%deps_dlnrho
            zone_deps_dlnt = burn%deps_dlnt
            energy_gen_rate = burn%total_energy_gen_rate
            energy_gen_component(1) = burn%pp_chain_energy_gen
            energy_gen_component(2) = burn%he3he4_be7_electron_energy_gen
            energy_gen_component(3) = burn%he3he4_be7_proton_energy_gen
            energy_gen_component(4) = burn%cno_cycle_energy_gen
            energy_gen_component(5) = burn%triple_alpha_energy_gen
            energy_gen_component(6) = burn%eps_neutrino_loss
            alpha_capture_energy_zone = burn%eps_alpha_capture
! Add the shell's solar-neutrino fluxes: burn%neutrino_flux is the flux
! at the Earth per gram of the shell (10^10 cm^-2 s^-1 g^-1,
! nuclear_energy.f90), so times dm it is the shell's flux.
            if (star%ctrl%calc_neutrinos) then
               do j = 1,n_nu_fluxes
                  star%neutrino_flux_total(j) = star%neutrino_flux_total(j) + &
                       burn%neutrino_flux(j)*star%dm(im)
               end do
            end if
! The shell's luminosity in each channel, dm eps_j/Lsun: pp1, pp2, pp3,
! CNO, triple alpha and the neutrino losses (zone_energy_luminosity sums
! them too, but nothing reads it).
            do j = 1,6
               star%luminosity_breakdown(j) = star%luminosity_breakdown(j) + &
                    (star%dm(im)/star%solar_luminosity_cgs)* &
                    energy_gen_component(j)
               zone_energy_luminosity = zone_energy_luminosity + &
                    (star%dm(im)/star%solar_luminosity_cgs)* &
                    energy_gen_component(j)
            end do
! The shell's He3+He3 luminosity, and its He3-burning luminosity
! (He3+He3 plus He3+He4), in Lsun.
            star%he3_he3_luminosity_zone(im) = (star%dm(im)/ &
                 star%solar_luminosity_cgs)*burn%eps_he3_he3
            star%he3_burning_luminosity_zone(im) = (star%dm(im)/ &
                 star%solar_luminosity_cgs)*burn%eps_he3_burning
! Alpha captures have their own channel, i_lum_he_c.
            star%luminosity_breakdown(i_lum_he_c)=star%luminosity_breakdown(i_lum_he_c)+(star%dm(im)/ &
                 star%solar_luminosity_cgs)*alpha_capture_energy_zone
            zone_energy_luminosity = zone_energy_luminosity + &
                 (star%dm(im)/star%solar_luminosity_cgs)* &
                 alpha_capture_energy_zone
            cur%ql = energy_gen_rate
            cur%ql_dt = cur%ql_dt + zone_deps_dlnt + &
                 zone_deps_dlnrho*dlnrho_dlnt
            cur%ql_dp = cur%ql_dp + zone_deps_dlnrho*dlnrho_dlnp
         end if
end subroutine form_nuclear_energy_terms

! ---------------------------------------------------------------
! Add the gravitational (entropy-change) energy term and its partials,
! and the change of the rotational kinetic energy, when the step has a
! duration. With chi = dln rho/dln T (negative) and grad_ad from the eos,
!   eps_grav = -T ds/dt = (P chi/rho) (Delta ln T/grad_ad - Delta ln P)/dt,
! where Delta is the change over the step (log_*_delta times ln10) and
! dt is the time step, or one year in zones with X > 0.01 when the step
! is shorter than that.
subroutine form_gravitational_energy_terms()
         if (compute_entropy_term) then
            zone_dt = delta_time_inv
! With accretion, the zones at or above envelope_zone_index add the
! accretion shifts of log T and log P (mass_change.f90) to their change
! over the step.
            if (star%job%use_mass_accretion.and.star%ctrl%mass_accretion_rate.gt.0.0d0) then
               if (im.ge.envelope_zone_index) then
                  zone_log_temperature_delta = star%log_temperature_delta(im)+ &
                       star%delta_log_temperature
                  zone_log_pressure_delta = star%log_pressure_delta(im)+ &
                       star%delta_log_pressure
               else
                  zone_log_temperature_delta = star%log_temperature_delta(im)
                  zone_log_pressure_delta = star%log_pressure_delta(im)
               end if
            else
               zone_log_temperature_delta = star%log_temperature_delta(im)
               zone_log_pressure_delta = star%log_pressure_delta(im)
            end if
            if (star%xa(i_h1,im).gt.0.01d0 .and. delta_time.lt.one_year_sec) &
                 zone_dt = one_year_sec_inv
! entropy_term1 = P chi/rho, entropy_term2 = entropy_term1/grad_ad,
! entropy_term = entropy_term2 Delta ln T - entropy_term1 Delta ln P,
! entropy_term3 = entropy_term2 Delta ln T, and eps_grav (egrav) =
! entropy_term/dt, stored in erg/g/s and added to channel 7 (i_lum_grav).
            entropy_term1 = eos_res%pressure*dlnrho_dlnt/eos_res%density
            entropy_term2 = entropy_term1/eos_res%adiabatic_gradient
            entropy_term = (entropy_term2*zone_log_temperature_delta - &
                 entropy_term1*zone_log_pressure_delta)*ln10
            entropy_term3 = entropy_term2*ln10*zone_log_temperature_delta
            egrav = zone_dt*entropy_term
            star%gravitational_luminosity(im) = egrav
            star%luminosity_breakdown(7) = star%luminosity_breakdown(7) + (star%dm(im)/ &
                 star%solar_luminosity_cgs)*egrav
! The partials differentiate each factor of eps_grav in turn: P/rho,
! chi, Delta ln P or Delta ln T, and 1/grad_ad (adiabatic_gradient_dp/_dt
! are dln grad_ad/dln P and dln grad_ad/dln T):
!   ql_dp += [entropy_term (1 - dln rho/dln P + dlnrho_dlnt_dp)
!             - entropy_term1 - entropy_term3 dln grad_ad/dln P]/dt
!   ql_dt += [entropy_term (-dln rho/dln T + dlnrho_dlnt_dt)
!             + entropy_term2 - entropy_term3 dln grad_ad/dln T]/dt
! The chi factor's terms use dlnrho_dlnt_dp/_dt, which the eos defines as
! d chi/dln P and d chi/dln T (eos_results.f90); the product rule wants
! them divided by chi. Only the Jacobian is affected.
            cur%ql = cur%ql + egrav
            cur%ql_dp = cur%ql_dp + zone_dt*(entropy_term*(1.0d0-dlnrho_dlnp+ &
                 eos_res%dlnrho_dlnt_dp)-entropy_term1 - entropy_term3* &
                 eos_res%adiabatic_gradient_dp)
            cur%ql_dt = cur%ql_dt + zone_dt*(entropy_term*(-dlnrho_dlnt+ &
                 eos_res%dlnrho_dlnt_dt) + entropy_term2 - entropy_term3* &
                 eos_res%adiabatic_gradient_dt)
! The change of the rotational kinetic energy is a sink:
!   ql -= (E_rot - E_rot_old)/(dm dt),
! with E_rot the shell's rotational kinetic energy. It is not
! linearized.
            if (star%job%rotation_active) then
               rotational_energy_term = zone_dt*(star%kinetic_energy_rot(im)- &
                    star%kinetic_energy_rot_old(im))/star%dm(im)
               cur%ql = cur%ql - rotational_energy_term
            end if
! There is no separate term for the entropy that accreted material
! brings in; accretion enters only through the shifted changes above.
         end if
end subroutine form_gravitational_energy_terms

! ---------------------------------------------------------------
! Scale the luminosity terms by cccql = ln10 m/Lsun, so ql = dL/dlog m
! with L in Lsun; eliminate the pair (im-1, im), or set the central
! boundary conditions at im = 1; carry cur into prev.
subroutine eliminate_shell_pair()
         cccql = star%ln10_over_lsun*star%m(im)
         cur%ql = cccql*cur%ql
         cur%ql_dp = cccql*cur%ql_dp
         cur%ql_dt = cccql*cur%ql_dt
! Copy qp, qt and qt_dl from the scratch slots into cur.
         cur%qp = pt_scr%qp
         cur%qt = pt_scr%qt
         cur%qt_dl = pt_scr%qtl
         if (im.gt.1) then
            im1 = im
            call henyey_eliminate(im1, prev, cur)
         else
! The central boundary conditions, rows 3 and 4 at shell 1 in the form
! of henyey_elimination.f90,
!   d y = elim_rhs - elim_coeff(:,1) d log P - elim_coeff(:,2) d log T.
! The innermost point is taken as a uniform sphere, m = (4 pi/3) r^3 rho:
!   log R = (log m - log(4 pi/3) - log rho)/3,
! so elim_rhs(3,1) is that minus log R and elim_coeff(3,:,1) is
! (dln rho/dln P, dln rho/dln T)/3. Its luminosity is L = eps m/Lsun =
! ql/ln10, so elim_rhs(4,1) = ql/ln10 - L and elim_coeff(4,:,1) =
! -(ql_dp, ql_dt).
            star%elim_coeff(3,1,1) = cc13*dlnrho_dlnp
            star%elim_coeff(3,2,1) = cc13*dlnrho_dlnt
            star%elim_rhs(3,1) = -cc13*(c4pi3l + zone_log10_density - &
                 zone_log_mass) - zone_log_radius
            star%elim_coeff(4,1,1) = -cur%ql_dp
            star%elim_coeff(4,2,1) = -cur%ql_dt
            star%elim_rhs(4,1) = clni*cur%ql - zone_luminosity_lsun
         end if
! Carry this shell's terms into the next pair's elimination.
         prev = cur
end subroutine eliminate_shell_pair

! ---------------------------------------------------------------
! Store the shell's energy generation for output, for every model:
!   eps_total = eps_pp1 + eps_pp2 + eps_pp3 + eps_cno + eps_3alpha (erg/g/s),
!   eps_channels(1:5) = the fraction of eps_total in each of those,
!   eps_channels(i_eps_neu) = eps_nu and eps_channels(i_eps_grav) =
!   eps_grav (erg/g/s).
! Below the burning cutoff only eps_grav is nonzero.
subroutine store_zone_energy_channels()
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
end subroutine store_zone_energy_channels

! ---------------------------------------------------------------
! Store the shell's eos, opacity and gradient values and Gamma1 for output.
subroutine store_zone_thermodynamics()
         star%beta(im) = eos_res%beta
         star%eta(im) = eos_res%electron_degeneracy_parameter
         star%opacity_zone(im) = kap_res(i_kap)
         star%gradr(im) = radiative_gradient
         star%gradT(im) = actual_gradient
         star%grada(im) = eos_res%adiabatic_gradient
         do j = 1,3
            star%fxion_zone(j,im) = eos_res%ion_fraction(j)
         end do
         star%conv_vel(im) = convective_velocity
         star%scp(im) = eos_res%specific_heat_cp
! Gamma1 = chi_rho cp/cv, with chi_rho = 1/(dln rho/dln P),
! chi_T = -chi_rho dln rho/dln T and cv = cp - (P/(rho T)) chi_T^2/chi_rho
! (Guenther 1995, private communication).
            chi_rho = 1.0d0/dlnrho_dlnp
            chi_t = -chi_rho*dlnrho_dlnt
            specific_heat_cv = eos_res%specific_heat_cp - exp(ln10*(star%logP(im)- &
                 star%logRho(im)-star%logT(im)))*chi_t**2/chi_rho
            star%adiabatic_index_gamma1(im) = chi_rho*eos_res%specific_heat_cp/ &
                 specific_heat_cv

! Disabled with rotation's advection-diffusion transport (2026-10-12),
! their only reader: the copies of dln kappa/dln T, d eps/dln T and
! (eps - eps_nu - eps_grav)/eps for rotation_stability_setup. They come
! back with that transport.
!         if (star%job%rotation_active) then
!            rot_scr%dlnkappa_dlnt(im) = kap_res(i_dlnkap_dlnt)
!            if (star%eps_total(im).gt.0.0d0) then
!               total_energy_sum = star%eps_total(im)
!               neutrino_and_grav_sum = star%eps_channels(i_eps_neu,im)+star%eps_channels(i_eps_grav,im)
!               rot_scr%neutrino_loss_fraction(im) = (total_energy_sum - &
!                    neutrino_and_grav_sum)/total_energy_sum
!            else
!               rot_scr%neutrino_loss_fraction(im) = 0.0d0
!            end if
!            rot_scr%deps_dlnt(im) = zone_deps_dlnt
!         end if
end subroutine store_zone_thermodynamics

! ---------------------------------------------------------------
! Store the shell's derivatives and mean molecular weights for the
! pulsation output, for every model: the GYRE-format writer
! (pulsation_output's write_gyre_pulse) runs every pulse_gyre_interval
! models, independently of the older OPAL pulsation writer's flag, and
! needs them (MHP 8/25). No other output file reads these arrays. Below
! the burning cutoff pulse_deps_dlnrho/dlnt get the last burning shell's
! values (see form_nuclear_energy_terms).
subroutine store_pulse_derivatives()
         star%pulse_dlnrho_dlnp(im) = dlnrho_dlnp
         star%pulse_deps_dlnrho(im) = zone_deps_dlnrho
         star%pulse_deps_dlnt(im) = zone_deps_dlnt
         star%pulse_dlnkap_dlnrho(im) = kap_res(i_dlnkap_dlnrho)
         star%pulse_dlnkap_dlnt(im) = kap_res(i_dlnkap_dlnt)
         star%pulse_specific_heat(im) = eos_res%specific_heat_cp
! eos_res%specific_gas_constant is R/mu, so mu = R/specific_gas_constant;
! the electron slot is 1/mu_e straight from the eos.
         if (eos_res%specific_gas_constant > 0.0d0) then
            star%pulse_mean_molecular_weight(im) = gas_constant/eos_res%specific_gas_constant
         else
            star%pulse_mean_molecular_weight(im) = 0.0d0
         end if
         star%pulse_electron_mean_weight_inverse(im) = eos_res%electron_mean_weight_inverse
         star%pulse_dlnrho_dlnt(im) = dlnrho_dlnt
! These MLT pulse-output arrays are always zero: the alfmlt/phmlt/cmxmlt
! scalars they once copied were only ever zeroed, and are gone. Whether
! to keep the arrays awaits an author decision (audit SUMMARY.md 1.1 #10).
         star%valfmlt(im) = 0.0d0
         star%vphmlt(im) = 0.0d0
         star%vcmxmlt(im) = 0.0d0
end subroutine store_pulse_derivatives

end subroutine compute_henyey_coefficients

end module structure_equations

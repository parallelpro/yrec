!----------------------------------------------------------------------
! shell_equations
!----------------------------------------------------------------------
! The structure equations at one shell: a star-blind kernel. Given the
! shell's state (shell_state) and the run's settings
! (structure_settings), evaluate_shell_equations calls the equation of
! state, the opacity, the temperature gradients and the nuclear energy
! generation, and returns the four equations' terms and their partials
! (henyey_shell_terms) and the shell's values for output
! (shell_diagnostics). It reads nothing from the star and keeps nothing
! between calls; compute_henyey_coefficients (structure_equations.f90)
! gathers each shell from the star, calls it, eliminates, and stores.
!
! Variables y = (log P, log T, log R, L/Lsun), coordinate log m (log10,
! cgs). Each term is q = dy/dlog m:
!   qp = dlog P/dlog m = -G m^2/(4 pi r^4 P) fp
!   qt = dlog T/dlog m = grad qp,   grad = dln T/dln P
!   qr = dlog R/dlog m = m/(4 pi r^3 rho)
!   ql = dL/dlog m     = (ln10 m/Lsun) eps,
!        eps = eps_nuc + eps_grav - d(e_rot)/dt + eps_extra   (erg/g/s)
! with fp the rotational correction to hydrostatic equilibrium and
! eps_extra the sum of the registered extra energy sources
! (energy_sources.f90; zero, and not evaluated, when there are none). The
! partials q_dp, q_dt, q_dr are dq/dln P, dq/dln T, dq/dln R (natural
! log) at fixed values of the other variables, and qt_dl = (dqt/dL)/ln10.
! Since ql is linear in eps, ql_dp and ql_dt carry the absolute
! derivatives d eps/dln P and d eps/dln T (times ln10 m/Lsun). The full
! contract and the difference scheme are in henyey_elimination.f90.
! core/test/test_structure checks every partial against finite
! differences of these terms.
module shell_equations
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: i_h1, i_metals, n_species_extended
      use physical_constants, only: log10_newton_g_default
      use eos_lib, only: eos_result, eos_get
      use kap_lib, only: kap_get, num_kap_results, i_kap, i_dlnkap_dlnrho, &
           i_dlnkap_dlnt
      use convection_lib, only: compute_temperature_gradients
      use net_lib, only: compute_nuclear_energy, nuclear_energy_result
      use energy_sources, only: energy_point, num_energy_sources, &
           evaluate_energy_sources
      implicit none
      private
      public :: henyey_shell_terms, shell_state, structure_settings, &
           shell_diagnostics, evaluate_shell_equations

! The structure-equation terms at one shell: for each of the four
! equations (pressure, temperature, radius, luminosity), the term q and
! its nonzero partials.
      type :: henyey_shell_terms
         real(dp) :: qp, qp_dr, qp_dp
         real(dp) :: qt, qt_dr, qt_dl, qt_dp, qt_dt
         real(dp) :: qr, qr_dr, qr_dp, qr_dt
         real(dp) :: ql, ql_dp, ql_dt
      end type henyey_shell_terms

! One shell's state: its mass coordinate (log m and m in g), the four
! structure variables, the density the eos starts from, the composition
! (composition_layout's slots), the rotation factors fp and ft, the
! changes of log P and log T over the step (with any accretion shift
! already added), the rotational kinetic energy now and at the start of
! the step (erg), the shell's mass dm (g), and two flags: inside the
! overshoot region, and inside the surface convection zone (where
! deuterium burning is capped by the turnover time).
      type :: shell_state
         real(dp) :: log_mass = 0.0d0, mass = 0.0d0, shell_mass = 0.0d0
         real(dp) :: log_temperature = 0.0d0, log_pressure = 0.0d0, &
              log_radius = 0.0d0, luminosity_lsun = 0.0d0
         real(dp) :: log_density = 0.0d0
         real(dp) :: xa(n_species_extended) = 0.0d0
         real(dp) :: fp = 1.0d0, ft = 1.0d0
         real(dp) :: log_temperature_delta = 0.0d0, &
              log_pressure_delta = 0.0d0
         real(dp) :: rotational_ke = 0.0d0, rotational_ke_start = 0.0d0
         logical :: in_overshoot = .false., in_surface_cz = .false.
      end type shell_state

! The settings the shells share: the time step dt (s; zero or negative
! means no time derivatives), the model's log Teff (the depth-varying
! spot contrast), the convective turnover time (s, the deuterium cap),
! ln10/Lsun, whether the star rotates, and log10 G (cgs).
      type :: structure_settings
         real(dp) :: dt = 0.0d0
         real(dp) :: log_teff = 0.0d0
         real(dp) :: turnover_timescale = 0.0d0
         real(dp) :: ln10_over_lsun = 0.0d0
         logical :: rotation_active = .false.
         real(dp) :: log10_gravitational_constant = log10_newton_g_default
      end type structure_settings

! The shell's values for output: the eos density, whether the shell is
! convective, the radiative and actual gradients and the convective
! velocity, the opacity results, Gamma1; the nuclear energy record (its
! chains' energies, the neutrino loss; burn%burning above net's burning
! cutoff, everything zero below) and its derivatives;
! and, when the step has a duration (has_eps_grav), eps_grav (erg/g/s);
! eps_extra, the extra energy sources' sum (erg/g/s; zero without any).
      type :: shell_diagnostics
         real(dp) :: log_density = 0.0d0
         logical :: convective = .false.
         real(dp) :: radiative_gradient = 0.0d0, actual_gradient = 0.0d0, &
              convective_velocity = 0.0d0
         real(dp) :: kap_res(num_kap_results) = 0.0d0
         real(dp) :: gamma1 = 0.0d0
         type(nuclear_energy_result) :: burn
         real(dp) :: deps_dlnrho = 0.0d0, deps_dlnt = 0.0d0
         logical :: has_eps_grav = .false.
         real(dp) :: eps_grav = 0.0d0
         real(dp) :: eps_extra = 0.0d0
      end type shell_diagnostics

contains

!----------------------------------------------------------------------
! evaluate_shell_equations
!----------------------------------------------------------------------
! The terms of the four structure equations and their partials at one
! shell (terms), and its values for output (diag), in named steps:
!   evaluate_eos_opacity_gradients   eos, opacity, temperature gradients
!   form_radius_pressure_temperature_terms
!   form_nuclear_energy_terms        the nuclear term (zero below net's cutoff)
!   form_gravitational_energy_terms  when the step has a duration
!   form_extra_energy_terms          the registered extra sources, if any
!   scale_luminosity_terms           ql to dL/dlog m in Lsun
!   compute_gamma1
! eos is the eos result: on entry the previous shell's, which the eos
! takes as its starting point (with the density set to sh%log_density);
! on return this shell's. saha_iterations is the eos's Saha iteration counter,
! carried from shell to shell. ierr is nonzero when the eos, the opacity,
! the gradients or an extra energy source fail.
subroutine evaluate_shell_equations(sh, cfg, eos, saha_iterations, terms, diag, ierr)
      use physical_constants, only: c4pil, clni, ln10
      use math_lib
      type(shell_state), intent(in) :: sh
      type(structure_settings), intent(in) :: cfg
      type(eos_result), intent(inout) :: eos
      integer, intent(inout) :: saha_iterations
      type(henyey_shell_terms), intent(out) :: terms
      type(shell_diagnostics), intent(out) :: diag
      integer, intent(out) :: ierr

! Every mesh point is interior (not atmosphere), and the eos and
! gradient derivatives are always wanted.
      logical, parameter :: want_derivatives = .true., in_atmosphere = .false.
! One year in seconds and the inverse the entropy term uses, as the F77
! COEFFT wrote them (TIME, TIMEI; the inverse is 3.16877e-8 to five
! digits, physical_constants' seconds_per_year is the year: replacing
! either moves results, a B item): the entropy term's floor on the time
! step while hydrogen is above eps_grav_floor_min_x
! (form_gravitational_energy_terms).
      real(dp), parameter :: one_year_sec = 3.1558d7, &
           one_year_sec_inv = 3.1688d-8, eps_grav_floor_min_x = 0.01d0
      real(dp) :: hydrogen_fraction, metal_fraction
      real(dp) :: dgrad_dt_component, dgrad_dp_component, &
           dgrad_dr_component
      real(dp) :: qtemp, qp, qt, qt_dl
      real(dp) :: zone_dt, entropy_term1, entropy_term2, &
           entropy_term, entropy_term3, rotational_energy_term, cccql
      real(dp) :: chi_rho, chi_t, specific_heat_cv

      ierr = 0
      hydrogen_fraction = sh%xa(i_h1)
      metal_fraction = sh%xa(i_metals)
      call evaluate_eos_opacity_gradients()
      if (ierr /= 0) return
      call form_radius_pressure_temperature_terms()
      call form_nuclear_energy_terms()
      call form_gravitational_energy_terms()
      call form_extra_energy_terms()
      if (ierr /= 0) return
      call scale_luminosity_terms()
      call compute_gamma1()

contains

! ---------------------------------------------------------------
! The equation of state, the opacity and the temperature gradients at
! the shell; a failure leaves ierr nonzero.
subroutine evaluate_eos_opacity_gradients()
         eos%log10_density = sh%log_density
         call eos_get(sh%log_temperature, sh%log_pressure, &
              hydrogen_fraction, metal_fraction, eos, &
              want_derivatives, in_atmosphere, saha_iterations, &
              composition_at_zone=sh%xa, ierr=ierr)
         if (ierr /= 0) return
! The density the eos returns at (P, T) is the one everything below uses
! (qr, the energy generation, the central boundary condition).
         diag%log_density = eos%log10_density
! The opacity at that density.
         call kap_get(eos%log10_density, sh%log_temperature, &
              hydrogen_fraction, metal_fraction, diag%kap_res, &
              eos%ion_fraction, ierr=ierr)
         if (ierr /= 0) return
         call compute_temperature_gradients(sh%log_temperature, sh%log_pressure, &
              eos, diag%kap_res, sh%log_radius, sh%log_mass, &
              sh%luminosity_lsun, diag%actual_gradient, diag%radiative_gradient, &
              dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
              diag%convective_velocity, want_derivatives, diag%convective, &
              sh%fp, sh%ft, cfg%log_teff, sh%in_overshoot, ierr)
         if (ierr /= 0) return
end subroutine evaluate_eos_opacity_gradients

! ---------------------------------------------------------------
! The radius, pressure and temperature equations' terms and partials
! (radiative or convective gradient).
subroutine form_radius_pressure_temperature_terms()
! The radius equation, qr = m/(4 pi r^3 rho) (qtemp = log(4 pi r^3)):
!   qr_dr = -3 qr, qr_dp = -qr dln rho/dln P, qr_dt = -qr dln rho/dln T.
       qtemp = c4pil + sh%log_radius + sh%log_radius + sh%log_radius
       terms%qr =+exp(ln10*(sh%log_mass - diag%log_density - qtemp))
       terms%qr_dr = - terms%qr - terms%qr - terms%qr
       terms%qr_dp = -terms%qr*eos%dlnrho_dlnp
       terms%qr_dt = -terms%qr*eos%dlnrho_dlnt
! The pressure equation, qp = -G m^2/(4 pi r^4 P) fp:
!   qp_dr = -4 qp, qp_dp = -qp.
! The rotation factors fp and ft are held fixed. A legacy variant
! subtracted their radius derivatives (legacy QFPR, QFTR):
! qp_dr = -qp (4 - QFPR), qt_dr = -qt (4 - QFTR).
       qp =-exp(ln10*(cfg%log10_gravitational_constant + sh%log_mass + sh%log_mass - &
            sh%log_pressure - qtemp - sh%log_radius ))*sh%fp
       terms%qp_dr = -qp - qp - qp - qp
       terms%qp_dp = -qp
! The temperature equation, qt = grad qp, with grad = dln T/dln P from
! compute_temperature_gradients; qt_dr = -4 qt comes from qp.
       qt = diag%actual_gradient*qp
       terms%qt_dr = -qt - qt - qt - qt
       if (.not.diag%convective) then
! Radiative: qt is proportional to kappa L m/(r^4 T^4) (times ft), so
!   qt_dl = qt/(L ln10),
!   qt_dp = qt (dln kappa/dln rho)(dln rho/dln P),
!   qt_dt = qt (-4 + dln kappa/dln T + (dln kappa/dln rho)(dln rho/dln T)).
          qt_dl = clni*qt/sh%luminosity_lsun
          terms%qt_dp = qt*diag%kap_res(i_dlnkap_dlnrho)*eos%dlnrho_dlnp
          terms%qt_dt = qt*(-4.0d0 + diag%kap_res(i_dlnkap_dlnt) + &
               diag%kap_res(i_dlnkap_dlnrho)*eos%dlnrho_dlnt)
       else
! Convective: with dgrad_dX_component = dln grad/dln X from the
! convection routine,
!   qt_dp = qt (-1 + dln grad/dln P)   (-1 from qp's 1/P),
!   qt_dt = qt dln grad/dln T,  qt_dr = qt (-4 + dln grad/dln r).
! The convective gradient's dependence on L is not linearized (qt_dl = 0).
          qt_dl = 0.0d0
          terms%qt_dp = qt*(-1.0d0 + dgrad_dp_component)
          terms%qt_dt = qt*dgrad_dt_component
          terms%qt_dr = terms%qt_dr + qt*dgrad_dr_component
       end if
end subroutine form_radius_pressure_temperature_terms

! ---------------------------------------------------------------
! The luminosity equation's nuclear term and partials. Before the
! ln10 m/Lsun scaling of scale_luminosity_terms,
!   ql = eps_nuc,
!   ql_dt = d eps/dln T + (d eps/dln rho)(dln rho/dln T),
!   ql_dp = (d eps/dln rho)(dln rho/dln P).
! At or below net's burning cutoff (nuclear_logT_cutoffs(1)) the record
! burn has burning = .false. and every value zero, and so are the terms.
! (The F77 coefft.f left them at the last burning shell's values,
! which reached the pulse output's d eps columns of every cooler shell,
! and eps_total of a shell exactly at the cutoff; fixed 2026-10.)
subroutine form_nuclear_energy_terms()
       terms%ql = 0.0d0
       terms%ql_dt = 0.0d0
       terms%ql_dp = 0.0d0
       diag%deps_dlnrho = 0.0d0
       diag%deps_dlnt = 0.0d0
       call compute_nuclear_energy(diag%log_density, sh%log_temperature, &
            sh%xa, sh%in_surface_cz, cfg%turnover_timescale, diag%burn)
       if (diag%burn%burning) then
            diag%deps_dlnrho = diag%burn%deps_dlnrho
            diag%deps_dlnt = diag%burn%deps_dlnt
            terms%ql = diag%burn%total_energy_gen_rate
            terms%ql_dt = terms%ql_dt + diag%deps_dlnt + &
                 diag%deps_dlnrho*eos%dlnrho_dlnt
            terms%ql_dp = terms%ql_dp + diag%deps_dlnrho*eos%dlnrho_dlnp
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
         diag%has_eps_grav = cfg%dt > 0.0d0
         if (diag%has_eps_grav) then
            zone_dt = 1.0d0/cfg%dt
            if (sh%xa(i_h1) > eps_grav_floor_min_x .and. cfg%dt < one_year_sec) &
                 zone_dt = one_year_sec_inv
! entropy_term1 = P chi/rho, entropy_term2 = entropy_term1/grad_ad,
! entropy_term = entropy_term2 Delta ln T - entropy_term1 Delta ln P,
! entropy_term3 = entropy_term2 Delta ln T, and eps_grav =
! entropy_term/dt (erg/g/s).
            entropy_term1 = eos%pressure*eos%dlnrho_dlnt/eos%density
            entropy_term2 = entropy_term1/eos%adiabatic_gradient
            entropy_term = (entropy_term2*sh%log_temperature_delta - &
                 entropy_term1*sh%log_pressure_delta)*ln10
            entropy_term3 = entropy_term2*ln10*sh%log_temperature_delta
            diag%eps_grav = zone_dt*entropy_term
! The partials differentiate each factor of eps_grav in turn: P/rho,
! chi, Delta ln P or Delta ln T, and 1/grad_ad (adiabatic_gradient_dp/_dt
! are dln grad_ad/dln P and dln grad_ad/dln T):
!   ql_dp += [entropy_term (1 - dln rho/dln P + dln chi/dln P)
!             - entropy_term1 - entropy_term3 dln grad_ad/dln P]/dt
!   ql_dt += [entropy_term (-dln rho/dln T + dln chi/dln T)
!             + entropy_term2 - entropy_term3 dln grad_ad/dln T]/dt
! The eos gives d chi/dln P and d chi/dln T (dlnrho_dlnt_dp/_dt,
! eos_results.f90), so dln chi = d chi/chi. (They were used without the
! 1/chi until 2026-10, as in coefft.f; only the Jacobian, so the Newton
! convergence, was affected.)
            terms%ql = terms%ql + diag%eps_grav
            terms%ql_dp = terms%ql_dp + zone_dt*(entropy_term*(1.0d0-eos%dlnrho_dlnp+ &
                 eos%dlnrho_dlnt_dp/eos%dlnrho_dlnt)-entropy_term1 - entropy_term3* &
                 eos%adiabatic_gradient_dp)
            terms%ql_dt = terms%ql_dt + zone_dt*(entropy_term*(-eos%dlnrho_dlnt+ &
                 eos%dlnrho_dlnt_dt/eos%dlnrho_dlnt) + entropy_term2 - entropy_term3* &
                 eos%adiabatic_gradient_dt)
! The change of the rotational kinetic energy is a sink:
!   ql -= (E_rot - E_rot_old)/(dm dt),
! with E_rot the shell's rotational kinetic energy. It is not
! linearized.
            if (cfg%rotation_active) then
               rotational_energy_term = zone_dt*(sh%rotational_ke- &
                    sh%rotational_ke_start)/sh%shell_mass
               terms%ql = terms%ql - rotational_energy_term
            end if
! There is no separate term for the entropy that accreted material
! brings in; accretion enters only through the shifted changes in
! sh%log_*_delta.
         end if
end subroutine form_gravitational_energy_terms

! ---------------------------------------------------------------
! Add the registered extra energy sources (energy_sources.f90) and their
! partials, as the nuclear term's are added: before the ln10 m/Lsun
! scaling,
!   ql += eps_extra,
!   ql_dt += d eps_extra/dln T + (d eps_extra/dln rho)(dln rho/dln T),
!   ql_dp += (d eps_extra/dln rho)(dln rho/dln P).
! With no source registered nothing is evaluated or added.
subroutine form_extra_energy_terms()
         type(energy_point) :: pt
         real(dp) :: eps_extra, deps_dlnrho, deps_dlnt
         diag%eps_extra = 0.0d0
         if (num_energy_sources() == 0) return
         pt%mass = sh%mass
         pt%shell_mass = sh%shell_mass
         pt%log_temperature = sh%log_temperature
         pt%log_pressure = sh%log_pressure
         pt%log_density = diag%log_density
         pt%log_radius = sh%log_radius
         pt%luminosity_lsun = sh%luminosity_lsun
         pt%xa = sh%xa
         pt%dt = cfg%dt
         call evaluate_energy_sources(pt, eps_extra, deps_dlnrho, deps_dlnt, ierr)
         if (ierr /= 0) return
         diag%eps_extra = eps_extra
         terms%ql = terms%ql + eps_extra
         terms%ql_dt = terms%ql_dt + deps_dlnt + deps_dlnrho*eos%dlnrho_dlnt
         terms%ql_dp = terms%ql_dp + deps_dlnrho*eos%dlnrho_dlnp
end subroutine form_extra_energy_terms

! ---------------------------------------------------------------
! Scale the luminosity terms by cccql = ln10 m/Lsun, so ql = dL/dlog m
! with L in Lsun, and set qp, qt and qt_dl.
subroutine scale_luminosity_terms()
         cccql = cfg%ln10_over_lsun*sh%mass
         terms%ql = cccql*terms%ql
         terms%ql_dp = cccql*terms%ql_dp
         terms%ql_dt = cccql*terms%ql_dt
         terms%qp = qp
         terms%qt = qt
         terms%qt_dl = qt_dl
end subroutine scale_luminosity_terms

! ---------------------------------------------------------------
! Gamma1 = chi_rho cp/cv, with chi_rho = 1/(dln rho/dln P),
! chi_T = -chi_rho dln rho/dln T and cv = cp - (P/(rho T)) chi_T^2/chi_rho
! (Guenther 1995, private communication).
subroutine compute_gamma1()
            chi_rho = 1.0d0/eos%dlnrho_dlnp
            chi_t = -chi_rho*eos%dlnrho_dlnt
            specific_heat_cv = eos%specific_heat_cp - exp(ln10*(sh%log_pressure- &
                 diag%log_density-sh%log_temperature))*chi_t**2/chi_rho
            diag%gamma1 = chi_rho*eos%specific_heat_cp/specific_heat_cv
end subroutine compute_gamma1

end subroutine evaluate_shell_equations

end module shell_equations

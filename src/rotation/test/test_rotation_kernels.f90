!----------------------------------------------------------------------
! test_rotation_kernels
!----------------------------------------------------------------------
! Rotation's kernels against the laws they keep. The model is synthetic
! (no deck, no tables): 60 zones of a solar-like star, r(k) = R (0.002 +
! 0.998 f), f = (k-1)/59, R = 6.9598e10 cm, with log rho and the shell
! masses as test_rotation builds them, a convective envelope above
! f = 0.71, and omega about ten times the solar rate. Each expected
! value is formed here from the model, not from the kernels. Checks:
!   inertia      zone_moments_of_inertia at omega = 0 (no distortion,
!                eta_2 = 0, r0 = r): I = (2/3) r^2 dm;
!   shape        rotation_shape_factors at omega = 1e-12 rad/s: fp = ft
!                = 1; from omega = 3.2e-6 to 2 omega, 1 - fp grows as
!                omega^2 (ratio 4) and 1 - ft as omega^4 (ratio 16: the
!                level surface's area and <g><1/g> change only at second
!                order in the distortion), each to 1e-3;
!   omega_from_j each radiative zone carries its own angular momentum,
!                I omega = j dm, and the convective envelope its total,
!                sum I omega = sum j dm, at one omega (a solid body);
!   AM diffusion am_diffusion_coeffs + solve_angular_momentum_transport
!                on 20 points (diffusion numbers dt D/(dchi I) of order
!                1, so that round-off stays round-off): without wind, sum I omega is unchanged
!                and the change sumdj is zero; with a wind, sum I omega
!                drops by exactly the surface term, -(1/2) (dOmega_exp +
!                dOmega_imp) I_N, and sumdj equals that drop; a solid
!                body stays one;
!   mixing       composition_diffusion_coeffs + solve_tridiagonal on 20
!                points (diffusion numbers of order 1): sum m X is unchanged, and a uniform X stays;
!   wind law     compute_solid_body_braking, Kawaler: dOmega ~ Omega^n
!                below saturation (ratio 2^n from Omega to 2 Omega),
!                ~ Omega Omega_sat^(n-1) above it (ratio 2);
!   spin-down    brake_convection_zone (Kawaler) over T in N = 2..16
!                calls, against the analytic laws: saturated, Omega_0
!                exp(-K Omega_sat^(n-1) T), with the implicit midpoint
!                rule's leading error (KT)^3/(12 N^2) to 1%; unsaturated
!                (n = 3), Omega_0 (1 + 2 K Omega_0^2 T)^(-1/2); both
!                errors fall as 1/N^2 (the order to 0.02); a step long
!                enough to be split into 4 sub-steps equals 4 calls of
!                a quarter step;
!   decay        am_diffusion_coeffs + solve_angular_momentum_transport
!                with I_i = h and D = 1 (dOmega/dt = d^2 Omega/dchi^2
!                on [0, 1], no flux at the ends): the cosine modes 1
!                and 3 are the implicit step's eigenvectors, with its
!                factor 1/(1 + 4 r sin^2(k pi/2m)), to round-off; the
!                grid's decay rate converges to pi^2 as h^2 (the order
!                to 0.01), and the amplitude at T = 1/pi^2 to exp(-1)
!                as 1/M in M steps (backward Euler; within 1% at M =
!                64, the order to 0.05);
!   velocities   instability_velocities on the radiative interfaces
!                2..40, omega rising outward (a shear):
!                Eddington-Sweet: U_ES ~ omega^2 (ratio 4 from omega to
!                2 omega, the quadrupole moment scaled with omega^2, the
!                omega^4 coefficient zero), and a solid body has no
!                shear term (U_ES does not depend on es_shear_coeff);
!                the mu gradient: a uniform mu leaves U_ES alone, and a
!                stabilizing one (mu rising inward) multiplies it by
!                1/2 to 1 (grad_mu enters twice; the author's decision,
!                docs/HOW_YREC_WORKS.md 4.6);
!                GSF (no inhibition): zero for a solid body, ~ omega^2
!                at a fixed shape of omega(r) (ratio 4); with viscous
!                inhibition, zero below 2 sqrt(nu/K) of the dynamical
!                limit;
!                shear: zero for a solid body and below the secular
!                threshold; above it, quadratic in the excess (ratio 4
!                from an excess e to 2 e); zero where a mu gradient's
!                limit exceeds the shear (mu falling linearly outward);
!                where the shear is dynamical, its correction lowers
!                U_GSF even when the trigger's centred gradient exceeds
!                the plain difference quotient (zones uneven in ln r);
!   regions      next_transport_region on a synthetic U: the regions
!                (2,4), (6,7) and (8,12), then none;
!   grids        on the region of zones 6..40, whose first zone is the
!                top of a convection zone (zones 1..6): equal_spaced_grid's points are
!                equally spaced in chi from the region's first chi to
!                its last; am_transport_grid's grid zone masses add up to
!                the region's mass, half-shells at its ends and the
!                convection zone below included (exactly: the sum
!                telescopes; the grid's J is a point quadrature, a few
!                percent off the model's here, which the map back's
!                fractional changes do not feel); equal_grid_to_model conserves J
!                in a pure redistribution (to round-off: it then always
!                rescales), matches a surface loss within its documented
!                1e-5, and gives the convection zone below the region the
!                region's first fractional change. (When the region's
!                first zone is radiative, the grid still counts a
!                convection zone below it but the map back leaves it
!                alone: docs/bug_sweep.md O13, open.)
!   region       diffuse_region_angular_momentum on the grids' region:
!                J conserved without a wind (to round-off) while j
!                moves; the wind's surface term, -(1/2) (dOmega_exp +
!                dOmega_imp) I_n, taken within the map back's 1e-5; a
!                solid body unchanged;
!   wind terms   outer_point_wind_terms: none below the surface; under
!                disk locking the point held and no loss; with the wind,
!                compute_solid_body_braking's terms at the point's I and
!                omega, which become the surface state, the implicit term
!                averaged with the surface region's earlier one only
!                (fresh when there is none); none with the wind torque
!                off;
!   Matt forms   the Matt law's solid-body and convection-zone forms
!                agree for one star over a short step (to 1e-5);
!   stability    rotation_stability_setup on a structure whose logs of
!                P, T and rho and whose other inputs are cubics in r, so
!                that its 4-point Lagrange interpolation to the interfaces
!                is exact: each coefficient against its formula
!                (stability_coefficients' header) evaluated here at the
!                interface: the GSF coefficient (32/3) H_P r^3/(G m
!                Delta), the Eddington-Sweet coefficients (8/3) r^3/(G m)
!                P/(g Delta c_P rho T) (L/m - eps) and its omega^4 part,
!                the shear coefficient, the dynamical-shear limit 2 g
!                sqrt((rho/P) (grad_ad - grad)), the mu-gradient
!                Richardson term, and the margins Delta and
!                del_grad_diff with the end interfaces' floors;
!   mechanisms   compute_mechanism_velocity on the stability check's
!                structure with a sheared omega: each slot's mechanism
!                writes its own slot and the one it declares it corrects
!                (transport_mechanism_corrects), and no other; a
!                corrected slot comes earlier; a slot out of range does
!                nothing;
!   checkj       check_angular_momentum: converged when the largest
!                |dj/j| since the previous iterate is within
!                convergence_tolerance, not converged above it; the
!                previous iterate becomes the new j; omega is the new j's
!                (I omega = j dm in the radiative zones); a zone with
!                j <= 0 stops the run (ierr /= 0);
!   checkc       check_composition: an interior abundance below zero by
!                less than 1e-5 of its surface value becomes its
!                neighbours' mean, with no cut; a larger one halves dt
!                and asks for a redo; the fourth cut stops the run; from
!                the second iteration, mu of a fully ionized gas: A/2
!                for H1, A/3 for He4, 12/7 for C12.
! Each check prints its relative error and ok (below its tolerance) or
! FAIL; the program stops with status 1 on any FAIL.
program test_rotation_kernels
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: json
      use rotation_configuration, only: rotation_settings
      use physical_constants, only: cgl, log10_radius_from_l_teff
      use omega_profile, only: omega_from_j
      use moments_of_inertia, only: zone_moments_of_inertia
      use rotating_shape, only: rotation_shape_factors
      use am_transport_coefficients, only: am_diffusion_coeffs
      use am_transport_coefficients, only: solve_angular_momentum_transport
      use rotational_mixing, only: composition_diffusion_coeffs
      use numerics_lib, only: solve_tridiagonal
      use magnetic_braking, only: compute_solid_body_braking, wind_conditions
      use magnetic_braking, only: brake_convection_zone, outer_point_wind_terms
      use region_transport, only: diffuse_region_angular_momentum
      use instability_velocities, only: compute_eddington_sweet_velocity, &
           inhibit_by_mu_gradient, compute_gsf_velocity, compute_shear_velocities
      use transport_regions, only: next_transport_region
      use transport_mechanisms, only: compute_mechanism_velocity
      use transport_mechanism_layout, only: n_transport_mechanisms, transport_mechanism_corrects
      use circulation, only: compute_quadrupole
      use rotation_scratch_lib, only: transport_grid_state, substep_physics_state
      use equal_grid, only: equal_spaced_grid
      use am_grid, only: am_transport_grid, equal_grid_to_model
      use stability_coefficients, only: rotation_stability_setup
      use rotation_scratch_lib, only: stability_state
      use transport_checks, only: check_angular_momentum, check_composition
      use composition_layout, only: i_h1, i_he4, i_metals, i_he3
      use star_info_lib, only: star_info
      implicit none

      integer, parameter :: n = 60
      real(dp), parameter :: rsun = 6.9598d10
      real(dp) :: log_rho(json), log_r(json), log_m(json), dm(json), r(json), &
           rho(json), f
      logical :: convective(json)
      integer :: i, failures

      failures = 0
! G and the settings the shape and the wind laws read (the controls'
! defaults where a run would take them).
      cgl = -7.17571d0
      rotation_settings%shape_outer_iterations = 5
      rotation_settings%shape_inner_iterations = 20
      rotation_settings%shape_tolerance = 1.0d-36
      rotation_settings%cz_rotation_law_exponent = 0.0d0
      rotation_settings%solar_radius_cgs = rsun
      rotation_settings%solar_mass_cgs = 1.989d33
      rotation_settings%log10_solar_luminosity = log10(3.828d33)
      do i = 1, n
         f = dble(i - 1)/dble(n - 1)
         r(i) = rsun*(0.002d0 + 0.998d0*f)
         log_r(i) = log10(r(i))
         log_rho(i) = 2.18d0 - 3.8d0*f**1.1d0 - 5.0d0*f**30
         rho(i) = 10.0d0**log_rho(i)
         convective(i) = f > 0.71d0
      end do
      dm(1) = 4.0d0/3.0d0*acos(-1.0d0)*r(1)**3*rho(1)
      log_m(1) = log10(dm(1))
      do i = 2, n
         dm(i) = 2.0d0*acos(-1.0d0)*(r(i)**2*rho(i) + r(i-1)**2*rho(i-1))*(r(i) - r(i-1))
         log_m(i) = log10(10.0d0**log_m(i-1) + dm(i))
      end do

      call inertia_check()
      call shape_check()
      call omega_from_j_check()
      call am_diffusion_check()
      call mixing_check()
      call wind_law_check()
      call spindown_check()
      call decay_checks()
      call velocity_checks()
      call dynamical_shear_check()
      call region_check()
      call grid_checks()
      call region_transport_checks()
      call wind_terms_checks()
      call matt_forms_check()
      call stability_check()
      call mechanism_slot_check()
      call checkj_check()
      call checkc_check()
      if (failures > 0) then
         write(*,'(i0,a)') failures, " check(s) failed"
         stop 1
      end if

contains

! one check: its relative error against its tolerance
subroutine report(label, rel_err, tolerance)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: rel_err, tolerance
      if (rel_err <= tolerance) then
         write(*,'(a,t56,1pe10.2,a)') label, rel_err, "  ok"
      else
         write(*,'(a,t56,1pe10.2,a)') label, rel_err, "  FAIL"
         failures = failures + 1
      end if
end subroutine report

subroutine inertia_check()
      real(dp) :: w(json), eta2(json), r0(json), inertia(json), di_dw(json), err
      w = 0.0d0
      eta2 = 0.0d0
      r0 = r
      call zone_moments_of_inertia(eta2, log_r, log_m, dm, 1, n, w, r0, inertia, di_dw)
      err = maxval(abs(inertia(1:n) - (2.0d0/3.0d0)*r(1:n)**2*dm(1:n))/ &
           ((2.0d0/3.0d0)*r(1:n)**2*dm(1:n)))
      call report("inertia: I = (2/3) r^2 dm at omega = 0", err, 1.0d-13)
end subroutine inertia_check

subroutine shape_check()
      real(dp) :: w(json), eta2(json), r0(json), fp(json), ft(json), gbar(json), &
           fp1(json), ft1(json), err, ratio_fp, ratio_ft
      integer :: ierr
      w = 1.0d-12
      call rotation_shape_factors(log_rho, log_r, log_m, n, w, eta2, fp, ft, gbar, r0, ierr)
      err = max(maxval(abs(fp(1:n) - 1.0d0)), maxval(abs(ft(1:n) - 1.0d0)))
      call report("shape: fp = ft = 1 as omega -> 0", err, 1.0d-12)
      w = 3.2d-6
      call rotation_shape_factors(log_rho, log_r, log_m, n, w, eta2, fp1, ft1, gbar, r0, ierr)
      w = 6.4d-6
      call rotation_shape_factors(log_rho, log_r, log_m, n, w, eta2, fp, ft, gbar, r0, ierr)
! zone 45, in the envelope, where the distortion is well above round-off
      ratio_fp = (1.0d0 - fp(45))/(1.0d0 - fp1(45))
      ratio_ft = (1.0d0 - ft(45))/(1.0d0 - ft1(45))
      call report("shape: 1 - fp ~ omega^2 (ratio 4)", abs(ratio_fp/4.0d0 - 1.0d0), 1.0d-3)
      call report("shape: 1 - ft ~ omega^4 (ratio 16)", abs(ratio_ft/16.0d0 - 1.0d0), 1.0d-3)
end subroutine shape_check

subroutine omega_from_j_check()
      real(dp) :: w(json), eta2(json), r0(json), inertia(json), di_dw(json), &
           jm(json), qiw(json), j_before(json), err, j_cz, i_omega_cz
      integer :: k, k_cz
      eta2 = 0.0d0
      r0 = r
      qiw = 0.0d0
      do k = 1, n
         w(k) = 2.9d-5*(1.0d0 + 0.5d0*(1.0d0 - dble(k - 1)/dble(n - 1))**2)
      end do
      call zone_moments_of_inertia(eta2, log_r, log_m, dm, 1, n, w, r0, inertia, di_dw)
      jm(1:n) = inertia(1:n)*w(1:n)/dm(1:n)
! the envelope spins faster than the solid body it becomes
      do k = 1, n
         if (convective(k)) jm(k) = 1.2d0*jm(k)
      end do
      j_before(1:n) = jm(1:n)*dm(1:n)
      call omega_from_j(log_rho, jm, log_r, log_m, dm, convective, n, eta2, inertia, &
           w, qiw, r0, 0.0d0)
      err = 0.0d0
      k_cz = n + 1
      do k = 1, n
         if (convective(k)) then
            k_cz = min(k_cz, k)
         else
            err = max(err, abs(inertia(k)*w(k) - j_before(k))/j_before(k))
         end if
      end do
      call report("omega_from_j: radiative zones keep I omega = j dm", err, 1.0d-12)
      j_cz = sum(j_before(k_cz:n))
      i_omega_cz = sum(inertia(k_cz:n)*w(k_cz:n))
      call report("omega_from_j: the envelope keeps its J", abs(i_omega_cz - j_cz)/j_cz, 1.0d-12)
      err = maxval(abs(w(k_cz:n) - w(n)))/w(n)
      call report("omega_from_j: the envelope rotates as a solid body", err, 1.0d-15)
end subroutine omega_from_j_check

subroutine am_diffusion_check()
      integer, parameter :: m = 20
      real(dp) :: d(json), ei(json), ej(json), ew(json), sub(json), diag(json), &
           sup(json), rhs(json), sol(json), dj(json), sumdj, wind_term, &
           j0, j1, dwe, dwi
      integer :: k, ierr
      do k = 1, m
         ei(k) = 1.0d50*(1.0d0 + 0.3d0*sin(1.7d0*k))**2
         ew(k) = 3.0d-6*(1.0d0 + 0.8d0*cos(0.9d0*k))
         ej(k) = ei(k)*ew(k)
         d(k) = 1.0d34*(1.0d0 + 0.5d0*sin(2.3d0*k))
      end do
      j0 = sum(ei(1:m)*ew(1:m))
! no wind
      call am_diffusion_coeffs(d, 0.05d0, 3.0d14, ei, ej, ew, m, 0.0d0, 0.0d0, &
           .false., sub, diag, sup, rhs, wind_term)
      call solve_angular_momentum_transport(m, ei, dj, sumdj, sub, diag, sup, rhs, &
           sol, wind_term, ierr)
      j1 = sum(ei(1:m)*sol(1:m))
      call report("AM diffusion: J conserved without wind", abs(j1 - j0)/j0, 1.0d-13)
      call report("AM diffusion: sumdj = 0 without wind", abs(sumdj)/j0, 1.0d-13)
! a wind: the surface loses -(1/2)(dOmega_exp + dOmega_imp) I_N
      dwe = 1.0d-3*ew(m)
      dwi = 2.0d-3*ew(m)
      call am_diffusion_coeffs(d, 0.05d0, 3.0d14, ei, ej, ew, m, dwe, dwi, &
           .false., sub, diag, sup, rhs, wind_term)
      call solve_angular_momentum_transport(m, ei, dj, sumdj, sub, diag, sup, rhs, &
           sol, wind_term, ierr)
      j1 = sum(ei(1:m)*sol(1:m))
      call report("AM diffusion: the wind takes its surface term", &
           abs((j1 - j0) - (-0.5d0*(dwe + dwi)*ei(m)))/abs(0.5d0*(dwe + dwi)*ei(m)), 1.0d-9)
      call report("AM diffusion: sumdj = the wind's loss", &
           abs(sumdj - (j1 - j0))/abs(j1 - j0), 1.0d-9)
! a solid body
      ew(1:m) = 3.0d-6
      ej(1:m) = ei(1:m)*ew(1:m)
      call am_diffusion_coeffs(d, 0.05d0, 3.0d14, ei, ej, ew, m, 0.0d0, 0.0d0, &
           .false., sub, diag, sup, rhs, wind_term)
      call solve_angular_momentum_transport(m, ei, dj, sumdj, sub, diag, sup, rhs, &
           sol, wind_term, ierr)
      call report("AM diffusion: a solid body stays one", &
           maxval(abs(sol(1:m) - 3.0d-6))/3.0d-6, 1.0d-14)
end subroutine am_diffusion_check

subroutine mixing_check()
      integer, parameter :: m = 20
      real(dp) :: d(json), x(json), em(json), sub(json), diag(json), sup(json), &
           rhs(json), x1(json), err
      integer :: k, ierr
      do k = 1, m
         em(k) = 1.0d31*(1.0d0 + 0.4d0*sin(1.3d0*k))
         x(k) = 0.7d0 - 0.3d0*exp(-(dble(k)/5.0d0)**2)
         d(k) = 1.0d15*(1.0d0 + 0.5d0*cos(2.1d0*k))
      end do
      call composition_diffusion_coeffs(d, 0.05d0, 3.0d14, x, em, m, sub, diag, sup, rhs)
      call solve_tridiagonal(sub(1:m), diag(1:m), sup(1:m), rhs(1:m), m, x1(1:m), ierr)
      err = abs(sum(em(1:m)*x1(1:m)) - sum(em(1:m)*x(1:m)))/sum(em(1:m)*x(1:m))
      call report("mixing: sum m X conserved", err, 1.0d-14)
      x(1:m) = 0.7d0
      call composition_diffusion_coeffs(d, 0.05d0, 3.0d14, x, em, m, sub, diag, sup, rhs)
      call solve_tridiagonal(sub(1:m), diag(1:m), sup(1:m), rhs(1:m), m, x1(1:m), ierr)
      call report("mixing: a uniform X stays", maxval(abs(x1(1:m) - 0.7d0))/0.7d0, 1.0d-14)
end subroutine mixing_check

subroutine wind_law_check()
      real(dp) :: d1, d2, dend, n_exp
      integer :: ierr
      type(wind_conditions) :: no_rossby
      rotation_settings%use_pmm_wind_law = .false.
      rotation_settings%wind_constant = 2.7d47
      rotation_settings%wind_law_omega_exponent = 3.0d0
      rotation_settings%wind_mass_loss_exponent = 0.5d0
      rotation_settings%wind_radius_exponent = 2.0d0
      rotation_settings%wind_mass_exponent = -0.5d0
      rotation_settings%wind_saturation_omega = 3.0d-5
      n_exp = rotation_settings%wind_law_omega_exponent
      dend = 0.0d0
      call compute_solid_body_braking(0.0d0, 1.0d15, 1.0d54, 1, 1.0d-6, 1.0d0, 3.762d0, &
           1.0d-6, d1, dend, no_rossby, ierr)
      call compute_solid_body_braking(0.0d0, 1.0d15, 1.0d54, 1, 2.0d-6, 1.0d0, 3.762d0, &
           2.0d-6, d2, dend, no_rossby, ierr)
      call report("wind law: dOmega ~ Omega^n below saturation", abs(d2/d1/2.0d0**n_exp - 1.0d0), 1.0d-13)
      call compute_solid_body_braking(0.0d0, 1.0d15, 1.0d54, 1, 1.0d-4, 1.0d0, 3.762d0, &
           1.0d-4, d1, dend, no_rossby, ierr)
      call compute_solid_body_braking(0.0d0, 1.0d15, 1.0d54, 1, 2.0d-4, 1.0d0, 3.762d0, &
           2.0d-4, d2, dend, no_rossby, ierr)
      call report("wind law: dOmega ~ Omega above saturation", abs(d2/d1/2.0d0 - 1.0d0), 1.0d-13)
end subroutine wind_law_check

! The wind's spin-down over T in N calls of brake_convection_zone
! (Kawaler law, the settings of wind_law_check) against the analytic
! solutions, with K = C 2^a (R/Rsun)^b M^c / I_cz (Mdot = 2e-14 Msun/yr,
! M = 1): saturated, dOmega/dt = -K Omega_sat^(n-1) Omega, so Omega =
! Omega_0 exp(-K Omega_sat^(n-1) T); unsaturated (n = 3), dOmega/dt =
! -K Omega^3, so Omega = Omega_0 (1 + 2 K Omega_0^2 T)^(-1/2). Each call
! solves the implicit midpoint rule (to 1e-6 between iterates), so the
! error falls as 1/N^2.
subroutine spindown_check()
      real(dp), parameter :: log_teff = 3.762d0
      real(dp) :: k_wind, radius, t_total, omega0, exact, err(4), w_split, w_steps
      integer :: level
      logical :: transport_saved
      transport_saved = rotation_settings%instability_transport_active
      rotation_settings%instability_transport_active = .false.
      radius = 10.0d0**log10_radius_from_l_teff(0.0d0, log_teff, &
           rotation_settings%log10_solar_luminosity)
      k_wind = rotation_settings%wind_constant*2.0d0**rotation_settings%wind_mass_loss_exponent &
           *(radius/rotation_settings%solar_radius_cgs)**rotation_settings%wind_radius_exponent &
           /1.0d54
! saturated: K Omega_sat^2 T = 1/2, from Omega_0 = 1e-4 (> Omega_sat
! throughout); N = 2, 4, 8, 16
      omega0 = 1.0d-4
      t_total = 0.5d0/(k_wind*rotation_settings%wind_saturation_omega**2)
      exact = omega0*exp(-0.5d0)
      do level = 1, 4
         err(level) = abs(spun_down(omega0, t_total, 2**level, log_teff) - exact)/exact
      end do
! the trapezoid's leading error for a linear decay: (K Omega_sat^2 T)^3
! /(12 N^2), at N = 16
      call report("spin-down: saturated, error = (KT)^3/(12 N^2)", &
           abs(err(4)/(0.5d0**3/(12.0d0*16.0d0**2)) - 1.0d0), 1.0d-2)
      call report("spin-down: saturated, error ~ 1/N^2", &
           abs(log(err(3)/err(4))/log(2.0d0) - 2.0d0), 2.0d-2)
! unsaturated: 2 K Omega_0^2 T = 1, from Omega_0 = 1e-6, to Omega_0/sqrt 2
      omega0 = 1.0d-6
      t_total = 0.5d0/(k_wind*omega0**2)
      exact = omega0/sqrt(2.0d0)
      do level = 1, 4
         err(level) = abs(spun_down(omega0, t_total, 2**level, log_teff) - exact)/exact
      end do
      call report("spin-down: unsaturated, Omega_0 (1+2K Omega_0^2 T)^-1/2", err(4), 1.0d-4)
      call report("spin-down: unsaturated, error ~ 1/N^2", &
           abs(log(err(3)/err(4))/log(2.0d0) - 2.0d0), 2.0d-2)
! a step whose first-guess loss is 3 Omega is split into 4 sub-steps:
! the same as 4 calls of a quarter step
      omega0 = 1.0d-4
      t_total = 3.0d0/(k_wind*rotation_settings%wind_saturation_omega**2)
      w_split = spun_down(omega0, t_total, 1, log_teff)
      w_steps = spun_down(omega0, t_total, 4, log_teff)
      call report("spin-down: a long step = its sub-steps", abs(w_split - w_steps)/w_steps, 1.0d-12)
      rotation_settings%instability_transport_active = transport_saved
end subroutine spindown_check

! Omega after n_calls calls of brake_convection_zone over t_total/n_calls
! each, on a one-zone convection zone (I_cz = 1e54, M_cz = 1e32), from
! Omega_0; each call's Omega is the zone's j M_cz/I_cz.
real(dp) function spun_down(omega0, t_total, n_calls, log_teff)
      real(dp), intent(in) :: omega0, t_total, log_teff
      integer, intent(in) :: n_calls
      real(dp), parameter :: i_cz = 1.0d54, m_cz = 1.0d32
      real(dp) :: jm(json)
      type(wind_conditions) :: no_rossby
      integer :: call_idx, ierr
      jm(1) = omega0*i_cz/m_cz
      do call_idx = 1, n_calls
         call brake_convection_zone(0.0d0, t_total/dble(n_calls), 0.0d0, m_cz, 1, 1, .true., &
              jm(1)*m_cz/i_cz, 1.0d0, log_teff, i_cz, jm, no_rossby, ierr)
      end do
      spun_down = jm(1)*m_cz/i_cz
end function spun_down

! Constant-coefficient angular momentum diffusion against its
! eigenmodes: shells of I_i = h on m points of spacing h = L/m, D =
! kappa, so that dOmega/dt = kappa d^2 Omega/dchi^2 on [0, L] with no
! flux at either end. On the grid the mode cos(k pi (i - 1/2)/m) is an
! eigenvector of the implicit step, with the factor 1/(1 + 4 r
! sin^2(k pi/(2 m))), r = dt kappa/h^2, exactly. In the continuum
! cos(pi chi/L) decays as exp(-kappa (pi/L)^2 t): the grid's rate,
! (1/a - 1)/dt for a step's factor a, converges to kappa (pi/L)^2 as
! h^2, and the amplitude after a time T in M steps (backward Euler) as
! 1/M.
subroutine decay_checks()
      real(dp), parameter :: pi = acos(-1.0d0)
      real(dp) :: a, h, r, expected, rate, err(4), lambda
      integer :: k_mode, level, m
! the grid's eigenvectors, modes 1 and 3 on 20 points, r about 1.3
      m = 20
      h = 1.0d0/dble(m)
      do k_mode = 1, 3, 2
         r = 1.3d0
         a = mode_amplitude(m, k_mode, r*h**2, 1)
         expected = 1.0d0/(1.0d0 + 4.0d0*r*sin(dble(k_mode)*pi/(2.0d0*m))**2)
         if (k_mode == 1) then
            call report("decay: mode 1 is the step's eigenvector", abs(a/expected - 1.0d0), 1.0d-13)
         else
            call report("decay: mode 3 is the step's eigenvector", abs(a/expected - 1.0d0), 1.0d-13)
         end if
      end do
! the rate against kappa (pi/L)^2 (kappa = L = 1), on m = 10..80 points,
! one step of dt = 1/lambda each
      lambda = pi**2
      do level = 1, 4
         m = 10*2**(level - 1)
         a = mode_amplitude(m, 1, 1.0d0/lambda, 1)
         rate = (1.0d0/a - 1.0d0)*lambda
         err(level) = abs(rate/lambda - 1.0d0)
      end do
      call report("decay: grid rate -> kappa (pi/L)^2 as h^2", &
           abs(log(err(3)/err(4))/log(2.0d0) - 2.0d0), 1.0d-2)
! the amplitude at T = 1/lambda on 400 points (rate error 5e-6), in
! M = 8..64 steps, against exp(-1)
      do level = 1, 4
         a = mode_amplitude(400, 1, 1.0d0/(lambda*dble(4*2**level)), 4*2**level)
         err(level) = abs(a/exp(-1.0d0) - 1.0d0)
      end do
      call report("decay: amplitude at T within 1% of exp(-1) (M = 64)", err(4), 1.0d-2)
      call report("decay: amplitude error ~ 1/M", &
           abs(log(err(3)/err(4))/log(2.0d0) - 1.0d0), 5.0d-2)
end subroutine decay_checks

! The amplitude of mode k_mode after n_steps steps of dt, from omega = 1
! + cos(k pi (i - 1/2)/m)/2 on m points of I_i = h = 1/m, D = 1: its
! projection on the mode over the starting one's.
real(dp) function mode_amplitude(m, k_mode, dt, n_steps)
      integer, intent(in) :: m, k_mode, n_steps
      real(dp), intent(in) :: dt
      real(dp), parameter :: pi = acos(-1.0d0)
      real(dp) :: d(json), ei(json), ej(json), ew(json), sub(json), diag(json), &
           sup(json), rhs(json), sol(json), dj(json), mode(json), sumdj, wind_term, h
      integer :: i, step, ierr
      h = 1.0d0/dble(m)
      do i = 1, m
         mode(i) = cos(dble(k_mode)*pi*(dble(i) - 0.5d0)/dble(m))
         ei(i) = h
         ew(i) = 1.0d0 + 0.5d0*mode(i)
         d(i) = 1.0d0
      end do
      do step = 1, n_steps
         ej(1:m) = ei(1:m)*ew(1:m)
         call am_diffusion_coeffs(d, h, dt, ei, ej, ew, m, 0.0d0, 0.0d0, &
              .false., sub, diag, sup, rhs, wind_term)
         call solve_angular_momentum_transport(m, ei, dj, sumdj, sub, diag, sup, rhs, &
              sol, wind_term, ierr)
         ew(1:m) = sol(1:m)
      end do
      mode_amplitude = sum((ew(1:m) - 1.0d0)*mode(1:m))/sum(0.5d0*mode(1:m)**2)
end function mode_amplitude

! ---------------------------------------------------------------
! The instability velocities on interfaces z1..z2 of the radiative
! interior (see the header).
subroutine velocity_checks()
      integer, parameter :: z1 = 2, z2 = 40
      real(dp), parameter :: g_cgs = 6.674d-8, ln10 = 2.302585092994046d0
      real(dp) :: w(json), wi(json), ri(json), grav(json), phi(json), c1(json), &
           c2(json), c3(json), ratio(json), u1(json), u2(json), mu(json), &
           mu_flat(json), mu_lin(json), mui_lin(json), lp(json), del(json), dw(json), qdyn(json), mui(json), &
           nu(json), kt(json), gf(json), kgsf(json), rmu(json), u0(json), q
      logical :: conv(json)
      integer :: k
      conv = .false.
      c2 = 0.0d0
      do k = 1, n
         w(k) = 3.0d-5*(1.0d0 + 0.5d0*r(k)/rsun)
         grav(k) = g_cgs*10.0d0**log_m(k)/r(k)**2
         phi(k) = 0.1d0*w(k)**2*r(k)**2
         c1(k) = 1.0d10*(1.0d0 + r(k)/rsun)
         c3(k) = 1.0d11
         lp(k) = 17.0d0 - 8.0d0*r(k)/rsun
         mu(k) = 0.62d0 + 0.25d0*exp(-(r(k)/(0.12d0*rsun))**2)
         mu_flat(k) = 0.62d0
         mu_lin(k) = 0.8d0 - 0.15d0*r(k)/rsun
         nu(k) = 10.0d0
         kt(k) = 1.0d7
         gf(k) = grav(k)
         kgsf(k) = 1.0d-3
         rmu(k) = 1.0d-4
         qdyn(k) = 1.0d-4
      end do
      do k = 2, n
         ri(k) = 0.5d0*(r(k) + r(k-1))
         wi(k) = 0.5d0*(w(k) + w(k-1))
         mui(k) = 0.5d0*(mu(k) + mu(k-1))
         mui_lin(k) = 0.5d0*(mu_lin(k) + mu_lin(k-1))
      end do
! Eddington-Sweet: omega -> 2 omega with the quadrupole moment ~ omega^2
      u1 = 0.0d0; u2 = 0.0d0
      call compute_eddington_sweet_velocity(z1, z2, conv, r, w, wi, ri, grav, phi, &
           c1, c2, c3, ratio, u1)
      call compute_eddington_sweet_velocity(z1, z2, conv, r, 2.0d0*w, 2.0d0*wi, ri, grav, &
           4.0d0*phi, c1, c2, c3, ratio, u2)
      call report("velocities: U_ES ~ omega^2 (ratio 4)", &
           maxval(abs(u2(z1:z2)/u1(z1:z2) - 4.0d0))/4.0d0, 1.0d-13)
! a solid body: no shear term
      u1 = 0.0d0; u2 = 0.0d0
      call compute_eddington_sweet_velocity(z1, z2, conv, r, 0.0d0*w + 3.0d-5, 0.0d0*wi + 3.0d-5, &
           ri, grav, phi, c1, c2, c3, ratio, u1)
      call compute_eddington_sweet_velocity(z1, z2, conv, r, 0.0d0*w + 3.0d-5, 0.0d0*wi + 3.0d-5, &
           ri, grav, phi, c1, c2, 0.0d0*c3, ratio, u2)
      call report("velocities: solid body, U_ES without a shear term", &
           maxval(abs(u2(z1:z2) - u1(z1:z2))/u1(z1:z2)), 0.0d0)
! the mu gradient: uniform mu, then a stabilizing one with
! del_grad_diff = 0.1 + grad_mu
      do k = 2, n
         del(k) = 0.1d0 + (log(mu(k)) - log(mu(k-1)))/(ln10*(lp(k) - lp(k-1)))
      end do
      u0 = 1.0d0; u1 = 1.0d0
      call inhibit_by_mu_gradient(z1, z2, mu_flat, lp, del, u1)
      call report("velocities: a uniform mu leaves U_ES alone", &
           maxval(abs(u1(z1:z2) - 1.0d0)), 0.0d0)
      u1 = 1.0d0
      call inhibit_by_mu_gradient(z1, z2, mu, lp, del, u1)
      q = max(0.0d0, 0.5d0 - minval(u1(z1:z2))) + max(0.0d0, maxval(u1(z1:z2)) - 1.0d0)
      if (minval(u1(z1:z2)) > 0.999d0) q = 1.0d0
      call report("velocities: a stabilizing mu cuts U_ES by 1/2 to 1", q, 0.0d0)
! GSF without inhibition: a solid body, then omega -> 2 omega
      do k = 2, n
         dw(k) = (w(k) - w(k-1))/(log(r(k)) - log(r(k-1)))
      end do
      u1 = 0.0d0
      call compute_gsf_velocity(z1, z2, 1, conv, r, 0.0d0*w + 3.0d-5, 0.0d0*dw, qdyn, mu, mui, &
           0.0d0*wi + 3.0d-5, nu, kt, gf, kgsf, u1)
      call report("velocities: U_GSF = 0 for a solid body", maxval(abs(u1(z1:z2))), 0.0d0)
      u1 = 0.0d0; u2 = 0.0d0
      call compute_gsf_velocity(z1, z2, 1, conv, r, w, dw, qdyn, mu, mui, wi, nu, kt, gf, kgsf, u1)
      call compute_gsf_velocity(z1, z2, 1, conv, r, 2.0d0*w, 2.0d0*dw, qdyn, mu, mui, 2.0d0*wi, &
           nu, kt, gf, kgsf, u2)
! (d ln omega/dr is a difference of logs of nearby omegas, which loses
! about three digits to cancellation, hence the looser tolerance)
      call report("velocities: U_GSF ~ omega^2 (ratio 4)", &
           maxval(abs(u2(z1:z2)/u1(z1:z2) - 4.0d0))/4.0d0, 1.0d-11)
! viscous inhibition (mode 2): nu/K so large that every |d omega/d ln r|
! is below 2 sqrt(nu/K) q_dyn
      u1 = 1.0d0
      call compute_gsf_velocity(z1, z2, 2, conv, r, w, dw, qdyn, mu, mui, wi, &
           0.0d0*nu + 1.0d7, kt, gf, kgsf, u1)
      call report("velocities: viscosity stabilizes the GSF", maxval(abs(u1(z1:z2))), 0.0d0)
! shear: a solid body, below the threshold, the quadratic onset, a mu
! gradient's limit. The threshold is sqrt(nu/K Re_c/8) q_dyn.
      q = sqrt(nu(1)/kt(1)*0.125d0*1.0d3)*qdyn(1)
      u1 = 0.0d0; u2 = 0.0d0
      call compute_shear_velocities(z1, z2, 1.0d3, log_r, r, w, 0.0d0*dw, qdyn, mu_flat, &
           0.0d0*mui + 0.62d0, ri, nu, kt, rmu, u1, u2)
      call report("velocities: no shear velocity for a solid body", maxval(abs(u1(z1:z2))), 0.0d0)
      u1 = 0.0d0
      call compute_shear_velocities(z1, z2, 1.0d3, log_r, r, w, 0.0d0*dw + 0.5d0*q, qdyn, &
           mu_flat, 0.0d0*mui + 0.62d0, ri, nu, kt, rmu, u1, u2)
      call report("velocities: no shear velocity below the threshold", maxval(abs(u1(z1:z2))), 0.0d0)
      u1 = 0.0d0; u0 = 0.0d0
      call compute_shear_velocities(z1, z2, 1.0d3, log_r, r, w, 0.0d0*dw + q*(1.0d0 + 1.0d-2), &
           qdyn, mu_flat, 0.0d0*mui + 0.62d0, ri, nu, kt, rmu, u1, u2)
      call compute_shear_velocities(z1, z2, 1.0d3, log_r, r, w, 0.0d0*dw + q*(1.0d0 + 2.0d-2), &
           qdyn, mu_flat, 0.0d0*mui + 0.62d0, ri, nu, kt, rmu, u0, u2)
      call report("velocities: shear onset quadratic (ratio 4)", &
           maxval(abs(u0(z1:z2)/u1(z1:z2) - 4.0d0))/4.0d0, 1.0d-10)
      u1 = 0.0d0
      call compute_shear_velocities(z1, z2, 1.0d3, log_r, r, w, 0.0d0*dw + q*(1.0d0 + 1.0d-2), &
           qdyn, mu_lin, mui_lin, ri, nu, kt, 0.0d0*rmu + 1.0d0, u1, u2)
      call report("velocities: a mu gradient stabilizes the shear", maxval(abs(u1(z1:z2))), 0.0d0)
end subroutine velocity_checks

! ---------------------------------------------------------------
! The dynamical shear's correction of U_GSF (compute_shear_velocities,
! docs/bug_sweep.md O11): on zones unevenly spaced in ln r, with each
! interface midway in r, the centred d omega/d ln r exceeds the plain
! difference quotient (Delta omega)/(Delta ln r); with the limit set
! between the two, the trigger fires while Delta omega is below limit
! Delta ln r. The corrected inner omega must still lie between the two
! zones' values, so the correction lowers U_GSF at every such interface
! (an inner omega beyond omega(i-1) would raise it).
subroutine dynamical_shear_check()
      integer, parameter :: m = 12
      real(dp) :: rr(m), lr(m), w(m), dw(m), lim(m), mu(m), mui(m), rint(m), visc(m), &
           kt(m), rich(m), uss(m), gsf(m), lnc, plain, centred
      integer :: k, bad
      do k = 1, m
         rr(k) = 1.0d10*(1.0d0 + 0.3d0*dble(k))**2
         lr(k) = log10(rr(k))
         w(k) = 1.0d-5*(1.0d0 + 0.2d0*dble(k))
      end do
      dw = 0.0d0; lim = 1.0d0; rint = 0.0d0
      do k = 2, m
         rint(k) = 0.5d0*(rr(k) + rr(k-1))
         lnc = log(rint(k))
         centred = 0.25d0*(w(k) - w(k-1))*(1.0d0/(log(rr(k)) - lnc) + 1.0d0/(lnc - log(rr(k-1))))
         plain = (w(k) - w(k-1))/(log(rr(k)) - log(rr(k-1)))
         dw(k) = centred
         lim(k) = 0.5d0*(abs(centred) + abs(plain))
      end do
      mu = 0.6d0; mui = 0.6d0; visc = 1.0d0; kt = 1.0d7; rich = 0.0d0
      uss = 0.0d0; gsf = 1.0d0
      call compute_shear_velocities(2, m, 1.0d3, lr, rr, w, dw, lim, mu, mui, rint, visc, &
           kt, rich, uss, gsf)
      bad = count(gsf(2:m) >= 1.0d0 .or. gsf(2:m) <= 0.0d0)
      call report("velocities: the dynamical correction lowers U_GSF", dble(bad), 0.0d0)
end subroutine dynamical_shear_check

! ---------------------------------------------------------------
! next_transport_region on U = (0, 0, 1, 1, 0, 0, 2, 0, 3, 3, 3, 3) at
! interfaces 2..12: regions (2,4), (6,7), (8,12), then none; and none
! for U = 0.
subroutine region_check()
      real(dp) :: u(12)
      integer :: expected(2,3), got(2,3), scan_start, zone_begin, zone_end, k, mismatches
      logical :: in_region, found
      u = [0.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 2.0d0, 0.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0]
      expected = reshape([2, 4, 6, 7, 8, 12], [2, 3])
      got = 0
      scan_start = 2; in_region = .false.; zone_begin = 0; zone_end = 0
      mismatches = 0
      do k = 1, 3
         call next_transport_region(u, 12, scan_start, in_region, zone_begin, zone_end, found)
         if (.not.found) mismatches = mismatches + 1
         got(:,k) = [zone_begin, zone_end]
      end do
      mismatches = mismatches + count(got /= expected)
      if (scan_start <= 12) mismatches = mismatches + 1
      call next_transport_region(u, 12, scan_start, in_region, zone_begin, zone_end, found)
      if (found) mismatches = mismatches + 1
      u = 0.0d0; scan_start = 2; in_region = .false.
      call next_transport_region(u, 12, scan_start, in_region, zone_begin, zone_end, found)
      if (found) mismatches = mismatches + 1
      call report("regions: (2,4), (6,7), (8,12), then none", dble(mismatches), 0.0d0)
end subroutine region_check

! ---------------------------------------------------------------
! The transport grid and the map back (see the header).
subroutine grid_checks()
      integer, parameter :: zb = 6, ze = 40
      type(transport_grid_state), save :: chi_grid
      type(substep_physics_state), save :: sub_physics
      real(dp) :: lum(json), lp(json), m(json), inertia(json), jm(json), w(json), &
           d_am(json), d_mix(json), dj(json), jm1(json), spacing, region_mass, &
           grid_j, s_loss, err, ratio_cz, ratio_first
      logical :: conv(json)
      integer :: k
      rotation_settings%chi_mass_scale = 8.0d-2
      rotation_settings%chi_luminosity_scale = 2.0d-2
      rotation_settings%chi_pressure_scale = 5.0d-2
      rotation_settings%solar_luminosity_cgs = 3.828d33
      do k = 1, n
         f = dble(k - 1)/dble(n - 1)
         lum(k) = 1.0d0 - exp(-(f/0.1d0)**2) + 1.0d-6
         lp(k) = 17.0d0 - 8.0d0*f
         m(k) = 10.0d0**log_m(k)
         inertia(k) = 2.0d0/3.0d0*r(k)**2*dm(k)
         w(k) = 3.0d-5*(1.0d0 + 0.5d0*f)
         jm(k) = w(k)*inertia(k)/dm(k)
         d_am(k) = 1.0d24*(1.0d0 + f)
         d_mix(k) = 1.0d23*(1.0d0 + f)
         conv(k) = k <= zb
         sub_physics%eps_plus_neu_grav(k) = 1.0d0
      end do
! the chi grid
      call equal_spaced_grid(chi_grid, lum, lp, log_m, zb, ze, n)
      spacing = chi_grid%eq_chi_spacing
      err = 0.0d0
      do k = 2, chi_grid%num_eq_points
         err = max(err, abs(chi_grid%eq_chi(k) - chi_grid%eq_chi(k-1) - spacing)/abs(spacing))
      end do
      err = max(err, abs(chi_grid%eq_chi(chi_grid%num_eq_points) - chi_grid%chi(ze - zb + 1))/ &
           abs(chi_grid%chi(ze - zb + 1) - chi_grid%chi(1)))
      call report("grids: chi points equally spaced, end to end", err, 1.0d-12)
! the momentum grid's zone masses and J
      call am_transport_grid(chi_grid, sub_physics, d_am, d_mix, log_rho, inertia, jm, lum, lp, &
           log_r, r, log_m, m, dm, log_m(n), zb, ze, conv, n, w)
      region_mass = 0.5d0*(m(ze) + m(ze+1)) - 0.5d0*(m(zb) + m(zb-1)) + sum(dm(1:zb-1))
      call report("grids: the grid masses add up to the region's", &
           abs(sum(chi_grid%eq_mass(1:chi_grid%num_eq_points)) - region_mass)/region_mass, 1.0d-13)
      grid_j = sum(chi_grid%eq_angular_momentum(1:chi_grid%num_eq_points))
! the map back: a redistribution, then a loss of 1e-3 of the grid's J
      do k = 1, chi_grid%num_eq_points
         dj(k) = 1.0d-3*chi_grid%eq_angular_momentum(k)*sin(6.0d0*dble(k)/dble(chi_grid%num_eq_points))
      end do
      dj(chi_grid%num_eq_points) = dj(chi_grid%num_eq_points) - sum(dj(1:chi_grid%num_eq_points))
      jm1 = jm
      call equal_grid_to_model(chi_grid, dj, chi_grid%eq_angular_momentum, dm, zb, ze, conv, n, 0.0d0, jm1)
      call report("grids: the map back conserves J in a redistribution", &
           abs(sum((jm1(1:n) - jm(1:n))*dm(1:n)))/sum(jm(1:n)*dm(1:n)), 1.0d-13)
      ratio_first = jm1(zb)/jm(zb)
      ratio_cz = maxval(abs(jm1(1:zb-1)/jm(1:zb-1) - ratio_first))
      call report("grids: the convection zone below takes the first change", ratio_cz, 1.0d-15)
      s_loss = -1.0d-3*grid_j
      do k = 1, chi_grid%num_eq_points
         dj(k) = s_loss*chi_grid%eq_angular_momentum(k)/grid_j
      end do
      jm1 = jm
      call equal_grid_to_model(chi_grid, dj, chi_grid%eq_angular_momentum, dm, zb, ze, conv, n, s_loss, jm1)
      call report("grids: the map back takes a loss to 1e-5", &
           abs(sum((jm1(1:n) - jm(1:n))*dm(1:n)) - s_loss)/abs(s_loss), 1.0d-5)
end subroutine grid_checks

! ---------------------------------------------------------------
! diffuse_region_angular_momentum on the region of grid_checks (zones
! 6..40 above a convection zone): without a wind the region's J is
! unchanged (to round-off: the map back rescales to the grid's total,
! zero) while j moves; with the wind's terms it loses (1/2) (dOmega_exp
! + dOmega_imp) I_n, within the map back's 1e-5; a solid body stays one.
! D is chosen for diffusion numbers dt D/(dchi I) of order 1 or below
! (j moves by up to 1e-3), so that round-off stays round-off.
subroutine region_transport_checks()
      integer, parameter :: zb = 6, ze = 40
      type(transport_grid_state), save :: chi_grid
      type(substep_physics_state), save :: sub_physics
      real(dp) :: lum(json), lp(json), m(json), inertia(json), jm(json), w(json), &
           d_am(json), d_mix(json), jm1(json), j0, err, loss, dwe, dwi
      logical :: conv(json)
      integer :: k, ierr
      do k = 1, n
         f = dble(k - 1)/dble(n - 1)
         lum(k) = 1.0d0 - exp(-(f/0.1d0)**2) + 1.0d-6
         lp(k) = 17.0d0 - 8.0d0*f
         m(k) = 10.0d0**log_m(k)
         inertia(k) = 2.0d0/3.0d0*r(k)**2*dm(k)
         w(k) = 3.0d-5*(1.0d0 + 0.5d0*f)
         jm(k) = w(k)*inertia(k)/dm(k)
         d_am(k) = 1.0d3*(1.0d0 + f)
         d_mix(k) = 1.0d23*(1.0d0 + f)
         conv(k) = k <= zb
         sub_physics%eps_plus_neu_grav(k) = 1.0d0
      end do
      j0 = sum(jm(1:n)*dm(1:n))
      call equal_spaced_grid(chi_grid, lum, lp, log_m, zb, ze, n)
      call am_transport_grid(chi_grid, sub_physics, d_am, d_mix, log_rho, inertia, jm, lum, lp, &
           log_r, r, log_m, m, dm, log_m(n), zb, ze, conv, n, w)
! no wind
      jm1 = jm
      call diffuse_region_angular_momentum(chi_grid, 1.0d13, 0.0d0, 0.0d0, .false., &
           dm, zb, ze, conv, n, jm1, ierr)
      err = abs(sum((jm1(1:n) - jm(1:n))*dm(1:n)))/j0
      if (ierr /= 0 .or. maxval(abs(jm1(zb:ze)/jm(zb:ze) - 1.0d0)) < 1.0d-4) err = 1.0d0
      call report("region: J conserved without a wind, j moves", err, 1.0d-13)
! the wind's terms: 1e-3 and 2e-3 of the outermost point's omega
      dwe = 1.0d-3*chi_grid%eq_omega(chi_grid%num_eq_points)
      dwi = 2.0d-3*chi_grid%eq_omega(chi_grid%num_eq_points)
      loss = -0.5d0*(dwe + dwi)*chi_grid%eq_moment_of_inertia(chi_grid%num_eq_points)
      jm1 = jm
      call diffuse_region_angular_momentum(chi_grid, 1.0d13, dwe, dwi, .false., &
           dm, zb, ze, conv, n, jm1, ierr)
      err = abs(sum((jm1(1:n) - jm(1:n))*dm(1:n)) - loss)/abs(loss)
      if (ierr /= 0) err = 1.0d0
      call report("region: the wind's surface term, within 1e-5", err, 1.0d-5)
! a solid body
      do k = 1, n
         w(k) = 3.0d-5
         jm(k) = w(k)*inertia(k)/dm(k)
      end do
      call am_transport_grid(chi_grid, sub_physics, d_am, d_mix, log_rho, inertia, jm, lum, lp, &
           log_r, r, log_m, m, dm, log_m(n), zb, ze, conv, n, w)
      jm1 = jm
      call diffuse_region_angular_momentum(chi_grid, 1.0d13, 0.0d0, 0.0d0, .false., &
           dm, zb, ze, conv, n, jm1, ierr)
      err = maxval(abs(jm1(1:n)/jm(1:n) - 1.0d0))
      if (ierr /= 0) err = 1.0d0
      call report("region: a solid body stays one", err, 1.0d-13)
end subroutine region_transport_checks

! ---------------------------------------------------------------
! outer_point_wind_terms (Kawaler law, the settings of wind_law_check):
! a region below the surface gets no terms (disk locking or not), and
! leaves the surface region's stored implicit term alone; at the surface
! under disk locking the point is held and loses nothing; with the wind,
! the terms are compute_solid_body_braking's for the point's inertia and
! omega, which become the surface state, the implicit term averaged with
! the stored one when an earlier iteration set it and fresh otherwise
! (docs/bug_sweep.md O27), and stored; with use_wind_torque off, none.
subroutine wind_terms_checks()
      real(dp), parameter :: inertia_n = 1.0d52, w_now = 2.0d-5, w_start = 2.1d-5
      real(dp) :: e, i_term, cz_i, w_s, e_ref, i_ref, err, stored
      logical :: fixed, stored_set
      integer :: ierr
      type(wind_conditions) :: no_rossby
      rotation_settings%use_wind_torque = .true.
! below the surface, even under disk locking
      cz_i = -1.0d0; w_s = -1.0d0; e = 9.0d0; i_term = 9.0d0
      stored = -3.0d0; stored_set = .true.
      call outer_point_wind_terms(.false., .true., 2, 0.0d0, 1.0d13, inertia_n, w_now, &
           1.0d0, 3.762d0, w_start, no_rossby, stored, stored_set, e, i_term, fixed, &
           cz_i, w_s, ierr)
      err = 0.0d0
      if (e /= 0.0d0 .or. i_term /= 0.0d0 .or. fixed .or. cz_i /= -1.0d0 .or. w_s /= -1.0d0 &
           .or. stored /= -3.0d0 .or. .not. stored_set) err = 1.0d0
      call report("wind terms: none below the surface", err, 0.0d0)
! disk locking
      e = 9.0d0; i_term = 9.0d0
      call outer_point_wind_terms(.true., .true., 2, 0.0d0, 1.0d13, inertia_n, w_now, &
           1.0d0, 3.762d0, w_start, no_rossby, stored, stored_set, e, i_term, fixed, &
           cz_i, w_s, ierr)
      err = 0.0d0
      if (e /= 0.0d0 .or. i_term /= 0.0d0 .or. .not. fixed .or. cz_i /= -1.0d0 &
           .or. stored /= -3.0d0) err = 1.0d0
      call report("wind terms: disk locking holds the point", err, 0.0d0)
! the wind, at the second iteration, with a stored implicit term
      i_ref = 1.0d-9
      call compute_solid_body_braking(0.0d0, 1.0d13, inertia_n, 2, w_now, 1.0d0, 3.762d0, &
           w_start, e_ref, i_ref, no_rossby, ierr)
      stored = 1.0d-9; stored_set = .true.
      call outer_point_wind_terms(.true., .false., 2, 0.0d0, 1.0d13, inertia_n, w_now, &
           1.0d0, 3.762d0, w_start, no_rossby, stored, stored_set, e, i_term, fixed, &
           cz_i, w_s, ierr)
      err = 0.0d0
      if (e /= e_ref .or. i_term /= i_ref .or. fixed .or. cz_i /= inertia_n .or. w_s /= w_now &
           .or. e_ref <= 0.0d0 .or. stored /= i_ref) err = 1.0d0
      call report("wind terms: the braking's, at the point's I and omega", err, 0.0d0)
! the wind at the second iteration with no stored term: fresh, as the
! first iteration's
      call compute_solid_body_braking(0.0d0, 1.0d13, inertia_n, 1, w_now, 1.0d0, 3.762d0, &
           w_start, e_ref, i_ref, no_rossby, ierr)
      stored = 7.0d0; stored_set = .false.
      call outer_point_wind_terms(.true., .false., 2, 0.0d0, 1.0d13, inertia_n, w_now, &
           1.0d0, 3.762d0, w_start, no_rossby, stored, stored_set, e, i_term, fixed, &
           cz_i, w_s, ierr)
      err = 0.0d0
      if (i_term /= i_ref .or. stored /= i_ref .or. .not. stored_set) err = 1.0d0
      call report("wind terms: a first surface term starts afresh", err, 0.0d0)
! no wind torque
      rotation_settings%use_wind_torque = .false.
      e = 9.0d0; i_term = 9.0d0; stored = -3.0d0
      call outer_point_wind_terms(.true., .false., 1, 0.0d0, 1.0d13, inertia_n, w_now, &
           1.0d0, 3.762d0, w_start, no_rossby, stored, stored_set, e, i_term, fixed, &
           cz_i, w_s, ierr)
      err = 0.0d0
      if (e /= 0.0d0 .or. i_term /= 0.0d0 .or. fixed .or. stored /= -3.0d0) err = 1.0d0
      call report("wind terms: none without the wind torque", err, 0.0d0)
end subroutine wind_terms_checks

! ---------------------------------------------------------------
! The Matt law's two forms (magnetic_braking) for one star: over a short
! step, the solid-body form's drop of omega at the step's start
! (compute_solid_body_braking) equals the convection zone's loss
! (brake_convection_zone, from its change of j) up to the implicit
! midpoint rule's first-order term. The solid-body form once lacked the
! 1/G of the solar centrifugal term (docs/bug_sweep.md O17), 0.17% here.
subroutine matt_forms_check()
! dt for a drop of about 2e-6 of omega (with these settings the torque
! is weak), where the midpoint and explicit losses differ by 1.5 times
! that (the torque goes as omega^3)
      real(dp), parameter :: i_cz = 1.0d54, m_cz = 1.0d32, omega = 1.0d-5, dt = 2.0d17
      real(dp) :: jm(json), d_start, d_end, d_cz
      type(wind_conditions) :: cond
      integer :: ierr
      rotation_settings%use_pmm_wind_law = .true.
      rotation_settings%scale_by_rossby_number = .false.
      rotation_settings%scale_by_b_field = .false.
      rotation_settings%wind_constant = 3.4d47
      rotation_settings%wind_law_omega_exponent = 3.0d0
      rotation_settings%wind_mass_exponent = -0.5d0
      rotation_settings%wind_radius_exponent = 3.1d0
      rotation_settings%wind_luminosity_exponent = 0.0d0
      rotation_settings%wind_pressure_exponent = 0.44d0
      rotation_settings%wind_turnover_exponent = 0.0d0
      rotation_settings%wind_centrifugal_exponent = 0.5d0
      rotation_settings%wind_centrifugal_constant = 0.0506d0
      rotation_settings%pmm_solar_omega = 2.6d-6
      rotation_settings%pmm_solar_pressure = 4.9752d0
      rotation_settings%pmm_solar_turnover_timescale = 9.6132d5
      rotation_settings%wind_saturation_omega = 3.0d-5
      cond = wind_conditions(turnover_timescale_start=1.0d6, turnover_timescale_end=1.0d6, &
           log_pphot_start=4.9752d0, log_pphot_end=4.9752d0, step_fraction=0.5d0)
      d_end = 0.0d0
      call compute_solid_body_braking(0.0d0, dt, i_cz, 1, omega, 1.0d0, 3.762d0, omega, &
           d_start, d_end, cond, ierr)
      jm(1) = omega*i_cz/m_cz
      call brake_convection_zone(0.0d0, dt, 0.0d0, m_cz, 1, 1, .true., omega, 1.0d0, 3.762d0, &
           i_cz, jm, cond, ierr)
      d_cz = omega - jm(1)*m_cz/i_cz
      call report("wind law: the Matt law's two forms agree", abs(d_cz/d_start - 1.0d0), 1.0d-5)
      rotation_settings%use_pmm_wind_law = .false.
end subroutine matt_forms_check

! ---------------------------------------------------------------
! rotation_stability_setup against its formulas (see the header).
subroutine stability_check()
      use physical_constants, only: cpi, radiation_constant_over_3
      integer, parameter :: zmin = 5, zmax = 50
      real(dp), parameter :: ln10 = 2.302585092994046d0, lsun = 3.828d33
      type(stability_state), save :: stab
      type(substep_physics_state), save :: phys
      real(dp) :: lrho(json), lt(json), lp(json), grav(json), lum(json), mass(json), &
           radius(json), shear_limit(json), x, xc, err(8), gconst, p, tt, rho, g, &
           mm, ll, eps, cp, grad, grada, dlnmu, floor_, delta, hp, thermal, b, rc, &
           expect, got
      integer :: k, i
      rotation_settings%solar_luminosity_cgs = lsun
      radiation_constant_over_3 = 7.5657d-15/3.0d0
      gconst = 10.0d0**cgl
      do k = 1, n
         x = r(k)/rsun
         lp(k) = poly(x, [17.2d0, -3.0d0, 1.5d0, -2.0d0])
         lt(k) = poly(x, [7.2d0, -1.2d0, 0.3d0, -0.2d0])
         lrho(k) = poly(x, [2.2d0, -4.0d0, 1.5d0, -0.5d0])
         grav(k) = poly(x, [3.0d4, 6.0d4, -3.0d4, 1.5d4])
         lum(k) = poly(x, [0.1d0, 0.8d0, 0.2d0, -0.1d0])
         mass(k) = poly(x, [1.0d32, 1.8d33, 4.0d32, -3.0d32])
         phys%gradr(k) = poly(x, [0.30d0, 0.05d0, 0.02d0, -0.01d0])
         phys%grada(k) = poly(x, [0.40d0, -0.02d0, 0.01d0, 0.005d0])
         phys%gradT(k) = min(phys%gradr(k), phys%grada(k))
         phys%qdt(k) = poly(x, [-1.0d0, -0.1d0, 0.0d0, 0.0d0])
         phys%eps_total(k) = poly(x, [5.0d0, -4.0d0, 1.0d0, 0.5d0])
         phys%cp(k) = poly(x, [3.0d8, 1.5d8, 3.0d7, -1.5d7])
         phys%mu(k) = 0.62d0 + 0.2d0*(1.0d0 - x)**3
      end do
      call rotation_stability_setup(stab, phys, lrho, grav, lum, lp, log_r, mass, lt, &
           zmin, zmax, n)
      radius = stab%radius
      shear_limit = stab%dynamical_shear_omega_limit
      err = 0.0d0
      do i = zmin, zmax
         rc = 0.5d0*(r(i) + r(i-1)); xc = rc/rsun
         p = 10.0d0**poly(xc, [17.2d0, -3.0d0, 1.5d0, -2.0d0])
         tt = 10.0d0**poly(xc, [7.2d0, -1.2d0, 0.3d0, -0.2d0])
         rho = 10.0d0**poly(xc, [2.2d0, -4.0d0, 1.5d0, -0.5d0])
         g = poly(xc, [3.0d4, 6.0d4, -3.0d4, 1.5d4])
         ll = lsun*poly(xc, [0.1d0, 0.8d0, 0.2d0, -0.1d0])
         mm = poly(xc, [1.0d32, 1.8d33, 4.0d32, -3.0d32])
         eps = poly(xc, [5.0d0, -4.0d0, 1.0d0, 0.5d0])
         cp = poly(xc, [3.0d8, 1.5d8, 3.0d7, -1.5d7])
         grada = poly(xc, [0.40d0, -0.02d0, 0.01d0, 0.005d0])
         grad = min(poly(xc, [0.30d0, 0.05d0, 0.02d0, -0.01d0]), grada)
         dlnmu = (log10(phys%mu(i)) - log10(phys%mu(i-1)))/(lp(i) - lp(i-1))
         floor_ = 1.0d-3
         if (i == zmin) floor_ = max(1.0d-3, 0.5d0*(phys%grada(zmin+1) - phys%gradT(zmin+1)) + dlnmu)
         if (i == zmax) floor_ = max(1.0d-3, 0.5d0*(phys%grada(zmax-1) - phys%gradT(zmax-1)) + dlnmu)
         delta = max(grada - grad, floor_)
         hp = p*rc**2/(rho*gconst*mm)
         thermal = p/(g*delta*cp*rho*tt)
         b = radiation_constant_over_3*tt**4/p
         call compare_into(err(1), stab%gsf_kippenhahn_coeff(i), 32.0d0/3.0d0*hp*rc**3/(gconst*mm*delta))
         call compare_into(err(2), stab%es_velocity_coeff1(i), 8.0d0/3.0d0*rc**3/(gconst*mm)*thermal*(ll/mm - eps))
         call compare_into(err(3), stab%es_velocity_coeff2(i), &
              -8.0d0/3.0d0*rc**3/(gconst*mm)*thermal*(ll/mm)/(2.0d0*cpi*gconst*rho))
         call compare_into(err(4), stab%es_shear_coeff(i), (ll/mm)*thermal/3.0d0/(cpi*gconst)*rc/rho)
         call compare_into(err(5), shear_limit(i), 2.0d0*g*sqrt(rho/p*max(grada - grad, 1.0d-3)))
         call compare_into(err(6), stab%mu_gradient_richardson_coeff(i), &
              rho/p*g**2/(ln10*(lp(i) - lp(i-1)))*(b - 1.0d0)/(1.0d0 + 3.0d0*b))
         call compare_into(err(7), stab%del_grad_diff_interface(i), max(grada - grad + dlnmu, floor_))
         call compare_into(err(8), stab%radius_interface(i), rc)
      end do
      call report("stability: GSF coefficient (32/3) H_P r^3/(G m Delta)", err(1), 1.0d-11)
      call report("stability: Eddington-Sweet coefficient", err(2), 1.0d-11)
      call report("stability: Eddington-Sweet omega^4 coefficient", err(3), 1.0d-11)
      call report("stability: Eddington-Sweet shear coefficient", err(4), 1.0d-11)
      call report("stability: dynamical-shear limit", err(5), 1.0d-11)
      call report("stability: mu-gradient Richardson term", err(6), 1.0d-11)
      call report("stability: del_grad_diff with the end floors", err(7), 1.0d-11)
      call report("stability: interfaces midway in r", err(8), 1.0d-14)
end subroutine stability_check

! ---------------------------------------------------------------
! check_angular_momentum's contract (see the header).
subroutine checkj_check()
      type(star_info), save :: st
      real(dp) :: j0(json), jprev(json), jnew(json), eta2(json), inertia(json), &
           w(json), qiw(json), r0(json), dt, err
      logical :: conv(json), converged, redo
      integer :: k, cuts, ierr
      rotation_settings%convergence_tolerance = 1.0d-5
      rotation_settings%max_diffusion_iters = 10
      rotation_settings%disk_locking_active = .false.
      conv = .false.
      st%disk_gate_age_gyr = 1.0d0
      do k = 1, n
         inertia(k) = 2.0d0/3.0d0*r(k)**2*dm(k)
         j0(k) = 3.0d-6*(1.0d0 + 0.3d0*dble(k)/dble(n))*inertia(k)/dm(k)
         st%old_omega(k) = 3.0d-6
      end do
      eta2 = 0.0d0; r0 = r
! a change of 1e-6 (within the tolerance)
      jprev = j0; jnew = j0*(1.0d0 + 1.0d-6); dt = 1.0d13; cuts = 0
      call check_angular_momentum(st, log_rho, jprev, j0, log_r, log_m, dm, 2, conv, n, dt, &
           eta2, inertia, jnew, cuts, converged, redo, w, qiw, r0, .false., ierr)
      err = 0.0d0
      if (.not. converged .or. redo .or. ierr /= 0) err = 1.0d0
      err = err + maxval(abs(jprev(1:n) - jnew(1:n))/jnew(1:n))
      do k = 1, n
         err = max(err, relerr(inertia(k)*w(k), jnew(k)*dm(k)))
      end do
      call report("checkj: within the tolerance: converged, omega from j", err, 1.0d-12)
! a change of 1e-3 (above it)
      jprev = j0; jnew = j0*(1.0d0 + 1.0d-3*sin(dble([(k, k = 1, json)])))
      call check_angular_momentum(st, log_rho, jprev, j0, log_r, log_m, dm, 2, conv, n, dt, &
           eta2, inertia, jnew, cuts, converged, redo, w, qiw, r0, .false., ierr)
      err = 0.0d0
      if (converged .or. redo .or. ierr /= 0) err = 1.0d0
      call report("checkj: above the tolerance: not converged", err, 0.0d0)
! a negative j stops the run
      jprev = j0; jnew = j0; jnew(20) = -jnew(20)
      call check_angular_momentum(st, log_rho, jprev, j0, log_r, log_m, dm, 2, conv, n, dt, &
           eta2, inertia, jnew, cuts, converged, redo, w, qiw, r0, .false., ierr)
      err = 0.0d0
      if (ierr == 0) err = 1.0d0
      call report("checkj: a negative j stops the run", err, 0.0d0)
end subroutine checkj_check

! ---------------------------------------------------------------
! compute_mechanism_velocity's slots (see the header): every element of
! mechanism_velocity starts at a sentinel no kernel writes, and after
! slot m's call only rows m and transport_mechanism_corrects(m) may
! differ from it; slots 0 and n + 1 write nothing.
subroutine mechanism_slot_check()
      integer, parameter :: zmin = 5, zmax = 50
      real(dp), parameter :: sentinel = -123.456d0
      type(stability_state), save :: stab
      type(substep_physics_state), save :: phys
      real(dp) :: lrho(json), lt(json), lp(json), grav(json), lum(json), mass(json), &
           w(json), v(n_transport_mechanisms, json), x
      logical :: conv(json)
      integer :: k, m, row, mismatches, written, ierr
      do k = 1, n
         x = r(k)/rsun
         lp(k) = poly(x, [17.2d0, -3.0d0, 1.5d0, -2.0d0])
         lt(k) = poly(x, [7.2d0, -1.2d0, 0.3d0, -0.2d0])
         lrho(k) = poly(x, [2.2d0, -4.0d0, 1.5d0, -0.5d0])
         grav(k) = poly(x, [3.0d4, 6.0d4, -3.0d4, 1.5d4])
         lum(k) = poly(x, [0.1d0, 0.8d0, 0.2d0, -0.1d0])
         mass(k) = poly(x, [1.0d32, 1.8d33, 4.0d32, -3.0d32])
         phys%gradr(k) = poly(x, [0.30d0, 0.05d0, 0.02d0, -0.01d0])
         phys%grada(k) = poly(x, [0.40d0, -0.02d0, 0.01d0, 0.005d0])
         phys%gradT(k) = min(phys%gradr(k), phys%grada(k))
         phys%qdt(k) = poly(x, [-1.0d0, -0.1d0, 0.0d0, 0.0d0])
         phys%eps_total(k) = poly(x, [5.0d0, -4.0d0, 1.0d0, 0.5d0])
         phys%cp(k) = poly(x, [3.0d8, 1.5d8, 3.0d7, -1.5d7])
! a uniform mu, so that no mu gradient stabilizes the shear: all three
! mechanisms act somewhere
         phys%mu(k) = 0.62d0
         w(k) = 3.0d-5*(1.0d0 + 2.0d0*x**2)
         conv(k) = .false.
      end do
      call rotation_stability_setup(stab, phys, lrho, grav, lum, lp, log_r, mass, lt, &
           zmin, zmax, n)
      call compute_quadrupole(stab, lrho, grav, w, n, ierr)
      do k = 2, n
         stab%omega_interface(k) = 0.5d0*(w(k) + w(k-1))
         stab%mean_molecular_weight_interface(k) = 0.5d0*(phys%mu(k) + phys%mu(k-1))
         stab%kinematic_viscosity_interface(k) = 30.0d0
         stab%thermal_diffusivity_interface(k) = 1.0d7
         stab%domega_dlnr(k) = (w(k) - w(k-1))/(log(r(k)) - log(r(k-1)))
      end do
! two interfaces dynamically unstable, where the shear corrects the GSF
! velocity
      stab%dynamical_shear_omega_limit(30:31) = 1.0d-3*stab%dynamical_shear_omega_limit(30:31)
      mismatches = 0
      written = 0
      do m = 0, n_transport_mechanisms + 1
         v = sentinel
         call compute_mechanism_velocity(m, zmin, zmax, stab, phys, w, log_r, lp, conv, v)
         do row = 1, n_transport_mechanisms
            if (m >= 1 .and. m <= n_transport_mechanisms) then
               if (row == m) then
                  if (any(v(row, zmin:zmax) /= sentinel)) written = written + 1
                  cycle
               end if
               if (row == transport_mechanism_corrects(m)) cycle
            end if
            if (any(v(row, 1:n) /= sentinel)) mismatches = mismatches + 1
         end do
      end do
! every mechanism wrote its own slot (here all three act somewhere)
      if (written /= n_transport_mechanisms) mismatches = mismatches + 1
      call report("mechanisms: each writes its slot and what it corrects", dble(mismatches), 0.0d0)
      mismatches = 0
      do m = 1, n_transport_mechanisms
         if (transport_mechanism_corrects(m) >= m .or. transport_mechanism_corrects(m) < 0) &
              mismatches = mismatches + 1
      end do
      call report("mechanisms: a corrected slot comes earlier", dble(mismatches), 0.0d0)
end subroutine mechanism_slot_check

! check_composition (see the header): a small negative interior
! abundance repaired, a larger one a cut, the fourth cut a stop, and mu
! for a fully ionized gas.
subroutine checkc_check()
      type(star_info), save :: st
      type(substep_physics_state), save :: phys
      real(dp), parameter :: a_h = 1.007825d0, a_he = 4.002603d0, a_he3 = 3.01603d0
      real(dp) :: comp(15, json), dt, err
      logical :: converged, redo
      integer :: k, cuts, ierr
      integer, parameter :: nz = 30
      rotation_settings%max_diffusion_iters = 10
      do k = 1, nz
         comp(:, k) = 0.0d0
         comp(i_h1, k) = 0.7d0; comp(i_he4, k) = 0.28d0
         comp(i_metals, k) = 0.02d0; comp(i_he3, k) = 1.0d-5
      end do
! He3 at -1e-11 in zone 10, below 1e-5 of its surface value in size:
! the mean of its neighbours, no cut
      comp(i_he3, 10) = -1.0d-11
      dt = 1.0d13; cuts = 0; converged = .true.
      call check_composition(st, phys, comp, 1, nz, dt, cuts, converged, redo, ierr)
      err = abs(comp(i_he3, 10) - 1.0d-5)/1.0d-5
      if (redo .or. cuts /= 0 .or. dt /= 1.0d13 .or. ierr /= 0) err = 1.0d0
      call report("checkc: a small negative X becomes its neighbours' mean", err, 1.0d-15)
! H at -1e-3 in zone 10: a cut, dt halved, a redo, not converged
      comp(i_h1, 10) = -1.0d-3
      call check_composition(st, phys, comp, 1, nz, dt, cuts, converged, redo, ierr)
      err = 0.0d0
      if (.not. redo .or. converged .or. cuts /= 1 .or. dt /= 0.5d13 .or. ierr /= 0) err = 1.0d0
      call report("checkc: a negative X halves the sub-step for a redo", err, 0.0d0)
! the fourth cut of the sub-step stops the run
      cuts = 3
      call check_composition(st, phys, comp, 1, nz, dt, cuts, converged, redo, ierr)
      err = 0.0d0
      if (ierr == 0) err = 1.0d0
      call report("checkc: the fourth cut stops the run", err, 0.0d0)
! mu at the second iteration: pure H1, He4 and metals (as C12) in zones
! 2, 3 and 4, the envelope mixture elsewhere; star%amuenv is that
! mixture's sum X_i/A_i. A fully ionized gas: mu = A/(1 + Z_nuc), so
! A(H)/2, A(He4)/3 and 12/7.
      comp(i_h1, 10) = 0.7d0
      st%envelope_hydrogen_fraction = 0.7d0; st%envelope_helium_fraction = 0.28d0
      st%envelope_metal_fraction = 0.02d0; st%envelope_he3_fraction = 1.0d-5
      st%amuenv = 0.7d0/a_h + 0.28d0/a_he + 0.02d0/12.0d0 + 1.0d-5/a_he3
      comp(1:4, 2) = 0.0d0; comp(i_h1, 2) = 1.0d0
      comp(1:4, 3) = 0.0d0; comp(i_he4, 3) = 1.0d0
      comp(1:4, 4) = 0.0d0; comp(i_metals, 4) = 1.0d0
      phys%mu = 0.0d0
      call check_composition(st, phys, comp, 2, nz, dt, cuts, converged, redo, ierr)
      err = max(relerr(phys%mu(2), a_h/2.0d0), relerr(phys%mu(3), a_he/3.0d0), &
           relerr(phys%mu(4), 12.0d0/7.0d0))
      if (redo .or. ierr /= 0) err = 1.0d0
      call report("checkc: mu = A/(1+Z) for ionized H, He4, C12", err, 1.0d-14)
end subroutine checkc_check

! a cubic c(1) + c(2) x + c(3) x^2 + c(4) x^3
pure real(dp) function poly(x, c)
      real(dp), intent(in) :: x, c(4)
      poly = c(1) + x*(c(2) + x*(c(3) + x*c(4)))
end function poly

! the relative difference of a from e
pure real(dp) function relerr(a, e)
      real(dp), intent(in) :: a, e
      relerr = abs(a - e)/max(abs(e), tiny(1.0d0))
end function relerr

! keep the largest relative difference in slot
subroutine compare_into(slot, a, e)
      real(dp), intent(inout) :: slot
      real(dp), intent(in) :: a, e
      slot = max(slot, relerr(a, e))
end subroutine compare_into

end program test_rotation_kernels

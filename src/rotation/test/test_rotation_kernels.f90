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
! Each check prints its relative error and ok (below its tolerance) or
! FAIL; the program stops with status 1 on any FAIL.
program test_rotation_kernels
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: json
      use rotation_configuration, only: rotation_settings
      use physical_constants, only: cgl
      use omega_profile, only: omega_from_j
      use moments_of_inertia, only: zone_moments_of_inertia
      use rotating_shape, only: rotation_shape_factors
      use am_transport_coefficients, only: am_diffusion_coeffs
      use am_transport_coefficients, only: solve_angular_momentum_transport
      use rotational_mixing, only: composition_diffusion_coeffs
      use numerics_lib, only: solve_tridiagonal
      use magnetic_braking, only: compute_solid_body_braking, wind_conditions
      use instability_velocities, only: compute_eddington_sweet_velocity, &
           inhibit_by_mu_gradient, compute_gsf_velocity, compute_shear_velocities
      use transport_regions, only: next_transport_region
      use rotation_scratch_lib, only: transport_grid_state, substep_physics_state
      use equal_grid, only: equal_spaced_grid
      use am_grid, only: am_transport_grid, equal_grid_to_model
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
      call velocity_checks()
      call region_check()
      call grid_checks()
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

end program test_rotation_kernels

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
!                ~ Omega Omega_sat^(n-1) above it (ratio 2).
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

end program test_rotation_kernels

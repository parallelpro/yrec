!----------------------------------------------------------------------
! test_rotation
!----------------------------------------------------------------------
! Standalone test of rotation's self-contained numerics, printed for
! byte-comparison against expected_test_rotation.out:
!   magnetic braking, the wind torque on the surface convection zone
!   (the controls and model quantities the laws read are set here);
!   the two quadratures of the shape factors (the trapezoid refinement
!   of the distortion integral, the Gauss-Legendre average over an
!   equipotential);
!   the tridiagonal solve of angular-momentum transport;
!   the rotating structure on a synthetic model: the moments of inertia
!   (and the imposed convection-zone law), omega from the specific
!   angular momentum, the shape factors fp and ft, and the viscosity
!   and thermal diffusivity (viscos).
! The braking entry points:
!   brake_convection_zone       the spin-down of the convection zone
!                               over a step: its j/m afterwards. The
!                               Matt (PMM) law with and without Rossby
!                               scaling, with B-field scaling,
!                               saturated, sub-stepped; the Kawaler law
!                               (use_pmm_wind_law false) with and
!                               without saturation scaled by tau_cz,
!                               inactive, and its error path (tau_cz
!                               undefined)
!   compute_solid_body_braking  the solid-body form: d(omega) at the
!                               start and end of the step, on a first
!                               and a later iteration: the Matt law
!                               with Rossby scaling, with B-field
!                               scaling too, with B-field scaling
!                               alone; the Kawaler law with and without
!                               saturation scaled by tau_cz
! The values are plausible solar-like ones, chosen to reach every
! branch.
program test_rotation
      use moments_of_inertia, only: zone_moments_of_inertia
      use rotating_shape, only: rotation_shape_factors
      use star_info_lib, only: star, json
      use physical_constants, only: csig, csigl, cgl
      use magnetic_braking, only: brake_convection_zone, compute_solid_body_braking
      use shell_quadrature, only: refine_distortion_integral, average_over_equipotential
      use secular_transport, only: solve_angular_momentum_transport
      implicit none
      double precision, parameter :: msun = 1.989d33
      double precision :: j(json), omega, dt, domega_start, domega_end
      integer :: ierr, k

! the constants setups computes at run time (its default G)
      csig = 5.67051d-5
      csigl = log10(csig)
      cgl = -7.17571d0
      star%solar_mass_cgs = msun
      star%solar_radius_cgs = 6.957d10
      star%log10_solar_luminosity = log10(3.828d33)
! a Matt-type law: tau_cz and P_phot at their solar values, exponents
! of order unity
      star%ctrl%constfactor = 3.4d47
      star%ctrl%wind_law_omega_exponent = 3.0d0
      star%ctrl%exm = -0.5d0
      star%ctrl%exr = 3.1d0
      star%ctrl%exl = 0.0d0
      star%ctrl%expr = 0.44d0
      star%ctrl%extau = 0.0d0
      star%ctrl%exmd = 0.5d0
      star%ctrl%excen = 0.5d0
      star%ctrl%c_2 = 0.0506d0
      star%ctrl%pmm_solar_omega = 2.6d-6
      star%ctrl%pmm_solar_pressure = 4.9752d0
      star%ctrl%pmm_solar_turnover_timescale = 9.6132d5
      star%job%wind_saturation_omega = 3.0d-5
      star%convective_turnover_timescale_old = 9.0d5
      star%convective_turnover_timescale = 1.1d6
      star%fracstep = 0.5d0
      star%pphot0 = 4.95d0
      star%pphot = 5.00d0

      write(*,'(a)') "# brake_convection_zone (j/m of the convection zone after the step)"
      star%ctrl%use_pmm_wind_law = .true.
      star%ctrl%scale_by_rossby_number = .true.
      star%ctrl%scale_by_b_field = .false.
      call cz_case("pmm rossby, unsaturated, 1e17 s", 2.6d-6, 1.0d17, .true.)
      call cz_case("pmm rossby, saturated, 1e15 s", 1.0d-4, 1.0d15, .true.)
      call cz_case("pmm rossby, sub-stepped, 1e18 s", 2.6d-6, 1.0d18, .true.)
      star%ctrl%scale_by_b_field = .true.
      call cz_case("pmm rossby + B field, 1e17 s", 2.6d-6, 1.0d17, .true.)
      call cz_case("pmm rossby + B field, saturated", 1.0d-4, 1.0d15, .true.)
! B-field scaling alone: the iteration scales omega by a tau_cz that
! is only set with Rossby scaling (zero), the first guess does not
      star%ctrl%scale_by_rossby_number = .false.
      call cz_case("pmm B field, no rossby, 1e17 s", 2.6d-6, 1.0d17, .true.)
      star%ctrl%scale_by_b_field = .false.
      call cz_case("pmm, no scaling, 1e17 s", 2.6d-6, 1.0d17, .true.)
      call cz_case("pmm, no scaling, wind flag off", 2.6d-6, 1.0d17, .false.)
      star%ctrl%use_pmm_wind_law = .false.
      star%job%instability_transport_active = .false.
      call cz_case("kawaler, 1e17 s", 2.6d-6, 1.0d17, .true.)
      call cz_case("kawaler, saturated, 1e15 s", 1.0d-4, 1.0d15, .true.)
      call cz_case("kawaler, sub-stepped, 1e18 s", 2.6d-6, 1.0d18, .true.)
      call cz_case("kawaler, wind flag off", 2.6d-6, 1.0d17, .false.)
      star%job%instability_transport_active = .true.
      star%job%wind_saturation_omega = 30.0d0
      call cz_case("kawaler, w_sat/tau_cz, 1e17 s", 2.6d-6, 1.0d17, .true.)
      star%convective_turnover_timescale = 0.5d0
      call cz_case("kawaler, tau_cz undefined", 2.6d-6, 1.0d17, .true.)
      star%convective_turnover_timescale = 1.1d6
      star%job%wind_saturation_omega = 3.0d-5

      write(*,'(a)') "# compute_solid_body_braking (d omega at the start and end of the step)"
      do k = 1, 5
         select case (k)
         case (1)
            star%ctrl%use_pmm_wind_law = .true.
            star%ctrl%scale_by_rossby_number = .true.
            star%ctrl%scale_by_b_field = .false.
         case (2)
            star%ctrl%scale_by_b_field = .true.
         case (3)
            star%ctrl%use_pmm_wind_law = .false.
            star%job%wind_saturation_omega = 30.0d0
         case (4)
            star%job%wind_saturation_omega = 3.0d-5
         case (5)
            star%ctrl%use_pmm_wind_law = .true.
            star%ctrl%scale_by_rossby_number = .false.
         end select
         omega = 2.6d-6
         dt = 1.0d16
         domega_end = 0.0d0
         call compute_solid_body_braking(0.0d0, dt, 1.0d54, 1, 0.9d0*omega, 1.0d0, &
              3.762d0, omega, domega_start, domega_end, ierr)
         write(*,'(a,i2,a,i2,2(1pe24.15))') "case", k, " iteration 1 ierr", ierr, &
              domega_start, domega_end
         call compute_solid_body_braking(0.0d0, dt, 1.0d54, 2, 0.95d0*omega, 1.0d0, &
              3.762d0, omega, domega_start, domega_end, ierr)
         write(*,'(a,i2,a,i2,2(1pe24.15))') "case", k, " iteration 2 ierr", ierr, &
              domega_start, domega_end
      end do
      call shape_quadrature_check()
      call angular_momentum_solve_check()
      call structure_check()
      write(*,'(a)') "test_rotation: done"

contains

! The convection zone: zones 3-6 of the model, 0.02 Msun with moment of
! inertia 1e54 g cm^2, rotating at omega; prints its j/m after the
! braking.
subroutine cz_case(label, omega_surface, timestep, wind_flag)
      character(len=*), intent(in) :: label
      double precision, intent(in) :: omega_surface, timestep
      logical, intent(in) :: wind_flag
      integer :: zone
      j = 0.0d0
      do zone = 1, 8
         j(zone) = 1.0d16*zone
      end do
      call brake_convection_zone(0.0d0, timestep, 0.98d0*msun, 1.0d0*msun, 3, 6, wind_flag, &
           omega_surface, 1.0d0, 3.762d0, 1.0d54, j, ierr)
      write(*,'(a,a,a,i2)') "cz ", label, "  ierr", ierr
      write(*,'(4(1pe24.15))') j(2:7)
end subroutine cz_case

! The trapezoid refinement of the distortion integral for one shell
! (n = 1..4) and the next (n restarts), then the Gauss-Legendre
! averages of g, 1/g and the area over the equipotential of shell 3.
subroutine shape_quadrature_check()
      double precision :: s, q, qp, r0(json), log_mass(json), g0g, ginvg, sphig
      integer :: n
      write(*,'(a)') "# shape quadratures (trapezoid refinement, equipotential average)"
      qp = 0.02d0
      do n = 1, 4
         call refine_distortion_integral(0.4d10, 0.6d10, s, n, 1.5d0, 2.0d0, 0.8d33, 0.6d33, 1.0d-11, &
              1.2d-11, 2.3d0, 2.25d0, q, qp)
         write(*,'(a,i2,2(1pe24.15))') "trapezoid n s q", n, s, q
      end do
      do n = 1, 2
         call refine_distortion_integral(0.6d10, 0.7d10, s, n, 1.2d0, 1.5d0, 0.9d33, 0.8d33, 0.9d-11, &
              1.0d-11, 2.35d0, 2.3d0, q, qp)
         write(*,'(a,i2,2(1pe24.15))') "trapezoid shell 2 n s q", n, s, q
      end do
      r0 = 0.0d0
      log_mass = 0.0d0
      r0(3) = 0.5d10
      log_mass(3) = 32.9d0
      call average_over_equipotential(g0g, ginvg, sphig, 0.125d0*4.0d0*3.1415926535898d0, &
           r0, log_mass, 1.0d33*s, q, 1.0d-11, 1.0d-3, 3)
      write(*,'(a,3(1pe24.15))') "equipotential g0g ginvg sphig", g0g, ginvg, sphig
end subroutine shape_quadrature_check

! The angular-momentum tridiagonal solve on six points: the new
! omega, the change in j weighted by the moment of inertia (dj, its
! surface element seeded with the wind loss) and their sum; then a
! singular system.
subroutine angular_momentum_solve_check()
      integer, parameter :: n = 6
      double precision :: sub(json), diag(json), sup(json), rhs(json), solution(json), &
           inertia(json), dj(json), sumdj
      integer :: i
      write(*,'(a)') "# angular-momentum tridiagonal solve"
      sub = 0.0d0
      diag = 0.0d0
      sup = 0.0d0
      rhs = 0.0d0
      do i = 1, n
         sub(i) = -1.0d0 + 0.1d0*i
         diag(i) = 4.0d0 + 0.5d0*i
         sup(i) = -1.5d0 + 0.05d0*i*i
         rhs(i) = 1.0d0 + i*(7.0d0 - i)
         inertia(i) = 0.2d0*i
      end do
      call solve_angular_momentum_transport(n, inertia, dj, sumdj, sub, diag, sup, rhs, solution, 0.125d0, ierr)
      write(*,'(a)') "omega"
      write(*,'(3(1pe24.15))') solution(1:n)
      write(*,'(a)') "dj"
      write(*,'(3(1pe24.15))') dj(1:n)
      write(*,'(a,1pe24.15,i3)') "sumdj ierr", sumdj, ierr
      diag(2) = sub(2)*sup(1)/diag(1)
      call solve_angular_momentum_transport(n, inertia, dj, sumdj, sub, diag, sup, rhs, solution, 0.125d0, ierr)
      write(*,'(a,i3)') "singular: ierr", ierr
end subroutine angular_momentum_solve_check

! The rotating structure on a synthetic 60-zone solar-like model with a
! convective envelope above 0.71 R, rotating differentially at about
! ten times the solar rate (so the distortion shows), in
! rezone's order: the moments of inertia, omega from the specific
! angular momentum (convection zones rotate as solid bodies), the shape
! factors fp and ft with eta2 and the mean radius, the moments again,
! and the imposed convection-zone law (walpcz /= 0: thin shells). Then
! the radiative and molecular viscosity and the thermal diffusivity.
subroutine structure_check()
      use star_info_lib, only: i_h1, i_he4, i_metals, i_he3
      integer, parameter :: n = 60
      integer, parameter :: shown(6) = [1, 12, 30, 43, 50, 60]
      double precision, parameter :: rsun = 6.9598d10
      double precision :: log_rho(json), log_r(json), log_m(json), log_t(json), &
           dm(json), w(json), eta2(json), r0(json), fp(json), ft(json), gbar(json), &
           inertia(json), di_dw(json), jm(json), qiw(json), f, r(json), rho(json)
      double precision, allocatable :: xa(:,:)
      logical :: convective(json)
      integer :: i, s
      write(*,'(a)') "# rotating structure"
      allocate(xa(15,json))
      xa = 0.0d0
      do i = 1, n
         f = dble(i - 1)/dble(n - 1)
         r(i) = rsun*(0.002d0 + 0.998d0*f)
         log_r(i) = log10(r(i))
         log_t(i) = 7.196d0 - 0.85d0*f**1.3d0 - 2.6d0*f**25
         log_rho(i) = 2.18d0 - 3.8d0*f**1.1d0 - 5.0d0*f**30
         rho(i) = 10.0d0**log_rho(i)
         convective(i) = f > 0.71d0
         w(i) = 2.9d-5*(1.0d0 + 0.5d0*(1.0d0 - f)**2)
         xa(i_h1,i) = 0.71d0 - 0.35d0*exp(-(f/0.12d0)**2)
         xa(i_metals,i) = 0.018d0
         xa(i_he3,i) = 1.0d-5
         xa(i_he4,i) = 1.0d0 - xa(i_h1,i) - xa(i_metals,i) - xa(i_he3,i)
      end do
      dm(1) = 4.0d0/3.0d0*acos(-1.0d0)*r(1)**3*rho(1)
      log_m(1) = log10(dm(1))
      do i = 2, n
         dm(i) = 2.0d0*acos(-1.0d0)*(r(i)**2*rho(i) + r(i-1)**2*rho(i-1))*(r(i) - r(i-1))
         log_m(i) = log10(10.0d0**log_m(i-1) + dm(i))
      end do
      eta2 = 0.0d0
      r0 = r
      qiw = 0.0d0
      star%ctrl%walpcz = 0.0d0
      call zone_moments_of_inertia(eta2, log_r, log_m, dm, 1, n, w, r0, inertia, di_dw)
      jm(1:n) = inertia(1:n)*w(1:n)/dm(1:n)
      call omega_from_j(log_rho, jm, log_r, log_m, dm, convective, n, eta2, inertia, w, qiw, r0)
      write(*,'(a)') "zone  omega from j/m"
      write(*,'(i4,1pe24.15)') (shown(s), w(shown(s)), s = 1, size(shown))
      call rotation_shape_factors(log_rho, log_r, log_m, n, w, eta2, fp, ft, gbar, r0, ierr)
      write(*,'(a,i2)') "shape factors: ierr", ierr
      write(*,'(a)') "zone  fp ft eta2 r0 gbar"
      do s = 1, size(shown)
         i = shown(s)
         write(*,'(i4,5(1pe24.15))') i, fp(i), ft(i), eta2(i), r0(i), gbar(i)
      end do
      call zone_moments_of_inertia(eta2, log_r, log_m, dm, 1, n, w, r0, inertia, di_dw)
      write(*,'(a)') "zone  I dI/domega"
      do s = 1, size(shown)
         write(*,'(i4,2(1pe24.15))') shown(s), inertia(shown(s)), di_dw(shown(s))
      end do
      star%ctrl%walpcz = -1.0d0
      call zone_moments_of_inertia(eta2, log_r, log_m, dm, 1, n, w, r0, inertia, di_dw)
      write(*,'(a)') "walpcz = -1: zone  I dI/domega"
      do s = 1, size(shown)
         write(*,'(i4,2(1pe24.15))') shown(s), inertia(shown(s)), di_dw(shown(s))
      end do
      star%ctrl%walpcz = 0.0d0
      do i = 1, n
         star%cp(i) = 3.0d8*(1.0d0 + dble(i)/dble(n))
         star%opacity_zone(i) = 1.0d0 + 100.0d0*(dble(i)/dble(n))**4
      end do
      call viscos(xa, log_rho, log_t, n)
      write(*,'(a)') "zone  viscosity thermal-diffusivity"
      do s = 1, size(shown)
         write(*,'(i4,2(1pe24.15))') shown(s), star%visc(shown(s)), star%thdif(shown(s))
      end do
end subroutine structure_check

end program test_rotation

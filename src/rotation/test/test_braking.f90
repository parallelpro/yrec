!----------------------------------------------------------------------
! test_braking
!----------------------------------------------------------------------
! Standalone test of magnetic braking, the wind torque on the surface
! convection zone. It sets the controls and the model quantities the
! braking laws read, then prints what each returns:
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
! branch. Results print for byte-comparison against
! expected_test_braking.out.
program test_braking
      use star_info_lib, only: star, json
      use physical_constants, only: csig, csigl, cgl
      use magnetic_braking, only: brake_convection_zone, compute_solid_body_braking
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

      write(*,'(a)') "# test_braking: brake_convection_zone (j/m of the convection zone after the step)"
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

      write(*,'(a)') "# test_braking: compute_solid_body_braking (d omega at the start and end of the step)"
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
      write(*,'(a)') "test_braking: done"

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

end program test_braking

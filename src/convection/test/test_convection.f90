!----------------------------------------------------------------------
! test_convection
!----------------------------------------------------------------------
! Standalone test of the convection physics, printed for
! byte-comparison against expected_test_convection.out:
!   compute_temperature_gradients  on hand-set EOS and opacity results:
!       a radiative point (and the adiabatic-overshoot gate), a
!       convective core above tgcut, mixing-length convection with and
!       without derivatives, starspots (constant and depth-dependent
!       contrast), and the error path (an unstable square root)
!   find_convection_zones  on a synthetic 40-zone model: no convection
!       zone, a core, an envelope, three zones, a fully convective
!       star; then overshoot (core, envelope, intermediate, the lovmax
!       limit, the rotating-model rule) with the pressure scale height
!       from the Yale EOS, including two zones that merge
!   homogenize_convection_zones  basic and extended composition, and
!       the renormalization floor
! Then the ok/FAIL checks, each confirmed to fail under a mutation
! (test_convection.mutations.toml; convection/README.md section 4 lists
! them): the gradient (del_rad's formula, the radiative and adiabatic
! cases, the mixing-length solution against one computed here, del
! between del_ad and del_rad, the efficient and inefficient limits'
! scalings, the derivatives against centred differences, with spots,
! the error contract), the zones (the scan, twelve zones and a
! thirteenth), overshoot (H_P, the depth, each kind of region's extent
! against alpha H_P, the rotating rule, the merge), homogenization (each
! species' mass kept, uniform zones, the floors) and semiconvection on a
! model with a convective core, with the OPAL95 opacities (an
! extension, the distance limit, the merge, a stable edge). The program
! exits 1 when a check fails.
program test_convection
      use star_info_lib, only: json, max_convective_zones, max_radiative_zones
      use eos_lib, only: eos_result, eos_init, eos_settings_t, eos_table_paths_t, eos_get
      use kap_lib, only: num_kap_results, i_kap, i_dlnkap_dlnrho, i_dlnkap_dlnt, &
           kap_init, kap_get, kap_settings_t, kap_table_paths_t
      use convection_lib, only: convection_settings_t, set_convection_settings, &
           compute_temperature_gradients, find_convection_zones, homogenize_convection_zones, &
           extend_by_semiconvection, pressure_scale_height, envelope_overshoot_depth, &
           convection_unstable_radicand, convection_too_many_zones
      use physical_constants, only: ln10, cdelrl, csig
      use run_messages, only: run_log_unit
      implicit none
! the settings the cases change, installed before each use
      type(convection_settings_t) :: cs
      integer :: ov_bottom, ov_top
      integer, parameter :: nz = 40
      double precision :: xa(15,json), log_rho(json), log_p(json), log_r(json), &
           log_m(json), log_t(json), dm(json)
      logical :: convective(json)
! the semiconvection model (build_semiconvective_model): its core is
! shells 1..sc_core, Z = sc_metals everywhere
      integer, parameter :: sc_core = 10
      double precision, parameter :: sc_metals = 0.016232d0
      double precision :: sc_xa(15,json), sc_log_rho(json), sc_log_p(json), &
           sc_log_r(json), sc_log_m(json), sc_log_t(json), sc_lum(json), sc_target(json)
      integer :: num_failed = 0
      integer :: rad_bounds(max_radiative_zones,2), mix_bounds(max_convective_zones,2), &
           mix_bounds_no_os(max_convective_zones,2)
      integer :: core_edge, env_edge, n_rad, n_mix, n_mix_no_os, ierr, i
      character(len=256) :: yrec_input

      call get_environment_variable("YREC_INPUT", yrec_input)
      if (len_trim(yrec_input) == 0) then
         write(*,'(a)') "test_convection: FAIL (YREC_INPUT not set)"
         stop 1
      end if
      call boot_yale_eos()
      run_log_unit = 6

      call gradient_cases()
      call build_model()
      call zone_cases()
      call homogenize_cases()
      call gradient_checks()
      call zone_checks()
      call overshoot_checks()
      call homogenize_checks()
      call semiconvection_checks()
      write(*,'(a)') "test_convection: done"
      if (num_failed > 0) then
         write(*,'(a,i0,a)') "test_convection: ", num_failed, " check(s) FAILED"
         stop 1
      end if

contains

! The Yale EOS alone (the Fermi table, the Saha cutoff at log T = 6),
! for the pressure scale heights of the overshoot cases.
subroutine boot_yale_eos()
      type(eos_settings_t) :: eos_cfg
      type(eos_table_paths_t) :: eos_paths
      run_log_unit = 20
      open(run_log_unit, file="test_convection.short", status="replace")
      eos_cfg%saha_log10t_cutoff = 6.0d0
      eos_paths%fermi = trim(yrec_input)//"/eos/yale/FERMI.TAB"
      call eos_init(eos_cfg, eos_paths, ierr)
      if (ierr /= 0) then
         write(*,'(a)') "test_convection: FAIL (eos_init error)"
         stop 1
      end if
      close(run_log_unit)
end subroutine boot_yale_eos

! One gradient evaluation at (log T, log P, log r, log m, L/Lsun) with
! the given EOS and opacity values; prints every output.
subroutine gradient(label, log_temperature, log_pressure, log_radius, log_mass, &
     luminosity, density, cp, grad_ad, dlnrho_dlnt, opacity, derivs, zone)
      character(len=*), intent(in) :: label
      double precision, intent(in) :: log_temperature, log_pressure, log_radius, &
           log_mass, luminosity, density, cp, grad_ad, dlnrho_dlnt, opacity
      logical, intent(in) :: derivs
      integer, intent(in) :: zone
      type(eos_result) :: eos_res
      double precision :: kap_res(num_kap_results), actual, radiative, dt, dp, dr, vconv
      logical :: is_convective
      eos_res%temperature = 10.0d0**log_temperature
      eos_res%pressure = 10.0d0**log_pressure
      eos_res%density = density
      eos_res%dlnrho_dlnt = dlnrho_dlnt
      eos_res%dlnrho_dlnp = 0.98d0
      eos_res%specific_heat_cp = cp
      eos_res%adiabatic_gradient = grad_ad
      eos_res%dlnrho_dlnt_dt = 0.05d0
      eos_res%dlnrho_dlnt_dp = -0.02d0
      eos_res%adiabatic_gradient_dt = -0.03d0
      eos_res%adiabatic_gradient_dp = 0.01d0
      eos_res%specific_heat_cp_dt = 0.4d0
      eos_res%specific_heat_cp_dp = -0.1d0
      kap_res = 0.0d0
      kap_res(i_kap) = opacity
      kap_res(i_dlnkap_dlnt) = -3.5d0
      kap_res(i_dlnkap_dlnrho) = 1.0d0
      actual = -1.0d0
      radiative = -1.0d0
      dt = -1.0d0
      dp = -1.0d0
      dr = -1.0d0
      vconv = -1.0d0
      call compute_temperature_gradients(log_temperature, log_pressure, eos_res, kap_res, &
           log_radius, log_mass, luminosity, actual, radiative, dt, dp, dr, vconv, &
           derivs, is_convective, 1.02d0, 0.99d0, 3.76d0, zone >= 10 .and. zone <= 14, ierr)
! only what the path sets: the gradients always; the velocity unless
! the square root failed; the derivatives for a convective point asked
! for them
      write(*,'(a,l2,i3)') label//": convective ierr", is_convective, ierr
      write(*,'(2(1pe24.15))') actual, radiative
      if (ierr == 0) write(*,'(a,1pe24.15)') "  velocity", vconv
      if (ierr == 0 .and. is_convective .and. derivs) &
           write(*,'(a,3(1pe24.15))') "  derivatives", dt, dp, dr
end subroutine gradient

subroutine gradient_cases()
      write(*,'(a)') "# compute_temperature_gradients"
! (the overshoot region below the envelope is zones 10-14)
      cs%mixing_length_alpha = 1.9d0
      cs%log10_solar_luminosity = log10(3.828d33)
      cs%log10_t_cut = 6.8d0
      cs%adiabatic_overshoot = .false.
      cs%spot_filling_factor = 0.0d0
      call set_convection_settings(cs)
! radiative: a small luminosity and opacity
      call gradient("radiative", 6.5d0, 16.0d0, 10.2d0, 33.2d0, 0.2d0, &
           5.0d0, 3.0d8, 0.4d0, -1.0d0, 1.0d0, .true., 12)
      cs%adiabatic_overshoot = .true.
      call set_convection_settings(cs)
      call gradient("radiative, adiabatic overshoot zone", 6.5d0, 16.0d0, 10.2d0, 33.2d0, 0.2d0, &
           5.0d0, 3.0d8, 0.4d0, -1.0d0, 1.0d0, .true., 12)
      call gradient("radiative, outside the overshoot zones", 6.5d0, 16.0d0, 10.2d0, 33.2d0, 0.2d0, &
           5.0d0, 3.0d8, 0.4d0, -1.0d0, 1.0d0, .true., 20)
      cs%adiabatic_overshoot = .false.
      call set_convection_settings(cs)
! convective core: log T above tgcut
      call gradient("core above tgcut", 7.1d0, 17.0d0, 9.8d0, 32.8d0, 2.5d0, &
           80.0d0, 2.0d8, 0.4d0, -1.0d0, 3.0d0, .true., 3)
! the surface convection zone: mixing-length theory
      call gradient("mixing length, derivatives", 4.8d0, 9.5d0, 10.82d0, 33.298d0, 1.0d0, &
           1.0d-5, 5.0d9, 0.25d0, -8.0d0, 5.0d3, .true., 38)
      call gradient("mixing length, no derivatives", 4.8d0, 9.5d0, 10.82d0, 33.298d0, 1.0d0, &
           1.0d-5, 5.0d9, 0.25d0, -8.0d0, 5.0d3, .false., 38)
! starspots
      cs%spot_filling_factor = 0.3d0
      cs%spot_temp_contrast = 0.8d0
      cs%spot_depth_varies = .false.
      call set_convection_settings(cs)
      call gradient("spots, constant contrast", 4.8d0, 9.5d0, 10.82d0, 33.298d0, 1.0d0, &
           1.0d-5, 5.0d9, 0.25d0, -8.0d0, 5.0d3, .true., 38)
      cs%spot_depth_varies = .true.
      call set_convection_settings(cs)
      call gradient("spots, contrast varying with depth", 4.8d0, 9.5d0, 10.82d0, 33.298d0, 1.0d0, &
           1.0d-5, 5.0d9, 0.25d0, -8.0d0, 5.0d3, .true., 38)
      cs%spot_filling_factor = 0.0d0
      call set_convection_settings(cs)
! a superadiabatic point whose density rises with T: no real root
      call gradient("unstable square root", 4.8d0, 9.5d0, 10.82d0, 33.298d0, 1.0d0, &
           1.0d-5, 5.0d9, 0.25d0, 0.5d0, 5.0d3, .true., 38)
end subroutine gradient_cases

! A 40-zone model from the centre (zone 1) out, solar-like: log T 7.2
! -> 3.8, log P 17.3 -> 5, X = 0.7, Z = 0.02.
subroutine build_model()
      double precision :: f
      xa = 0.0d0
      do i = 1, nz
         f = dble(i - 1)/dble(nz - 1)
         log_t(i) = 7.2d0 - 3.4d0*f**1.5d0
         log_p(i) = 17.3d0 - 12.3d0*f**1.2d0
         log_r(i) = 9.0d0 + 1.84d0*f**0.5d0
         log_m(i) = 31.0d0 + 2.3d0*f**0.3d0
         log_rho(i) = 2.2d0 - 9.0d0*f**1.3d0
         dm(i) = 1.0d31*(1.0d0 + f)
         xa(1,i) = 0.70d0 - 0.3d0*(1.0d0 - f)**4
         xa(3,i) = 0.02d0
         xa(2,i) = 1.0d0 - xa(1,i) - xa(3,i)
         xa(9,i) = 1.0d-5*(1.0d0 + f)
         xa(12,i) = 3.0d-6*f
      end do
end subroutine build_model

! convective(i) true for the listed ranges, then find_convection_zones
subroutine zones(label, ranges)
      character(len=*), intent(in) :: label
      integer, intent(in) :: ranges(:)
      integer :: k
      convective = .false.
      do k = 1, size(ranges), 2
         convective(ranges(k):ranges(k+1)) = .true.
      end do
      mix_bounds = 0
      mix_bounds_no_os = 0
      rad_bounds = 0
      ov_bottom = -7
      ov_top = -7
      call set_convection_settings(cs)
      call find_convection_zones(xa, log_rho, log_p, log_r, log_m, log_t, convective, nz, &
           rad_bounds, mix_bounds, mix_bounds_no_os, core_edge, env_edge, n_rad, n_mix, &
           n_mix_no_os, ov_bottom, ov_top, ierr)
      write(*,'(a,6i4)') label//": core_edge env_edge n_rad n_mix n_mix_no_os ierr", &
           core_edge, env_edge, n_rad, n_mix, n_mix_no_os, ierr
      write(*,'(a,24i4)') "  mixed", (mix_bounds(k,1), mix_bounds(k,2), k = 1, max(n_mix, 1))
      write(*,'(a,24i4)') "  mixed without overshoot", (mix_bounds_no_os(k,1), &
           mix_bounds_no_os(k,2), k = 1, max(n_mix_no_os, 1))
      write(*,'(a,26i4)') "  radiative", (rad_bounds(k,1), rad_bounds(k,2), k = 1, max(n_rad, 1))
      write(*,'(a,2i4)') "  iov1 iov2", ov_bottom, ov_top
end subroutine zones

subroutine zone_cases()
      write(*,'(a)') "# find_convection_zones"
      cs%core_overshoot = .false.
      cs%envelope_overshoot = .false.
      cs%intermediate_overshoot = .false.
      cs%rotation_active = .false.
      cs%instability_transport_active = .false.
      call zones("no convection", [1, 0])
      call zones("core", [1, 6])
      call zones("envelope", [31, 40])
      call zones("core, intermediate, envelope", [1, 5, 18, 21, 31, 40])
      call zones("fully convective", [1, 40])
      write(*,'(a)') "# find_convection_zones with overshoot (Yale EOS scale heights)"
      cs%core_overshoot = .true.
      cs%envelope_overshoot = .true.
      cs%intermediate_overshoot = .true.
      cs%overshoot_alpha_core = 0.6d0
      cs%overshoot_alpha_envelope = 0.6d0
      cs%overshoot_alpha_intermediate = 0.8d0
      cs%limit_core_overshoot = .false.
      cs%core_overshoot_limit = 0.1d0
      call zones("core, intermediate, envelope", [1, 5, 18, 21, 31, 40])
      call zones("fully convective", [1, 40])
      cs%limit_core_overshoot = .true.
      call zones("core, lovmax limit", [1, 5, 31, 40])
      cs%limit_core_overshoot = .false.
      cs%overshoot_alpha_core = 2.5d0
      call zones("core and intermediate merging", [1, 5, 9, 12, 31, 40])
      cs%overshoot_alpha_core = 0.6d0
      cs%overshoot_alpha_envelope = 1.0d-4
      cs%rotation_active = .true.
      cs%instability_transport_active = .true.
      call zones("tiny envelope overshoot, rotating", [1, 5, 31, 40])
      cs%rotation_active = .false.
      cs%instability_transport_active = .false.
end subroutine zone_cases

subroutine homogenize_cases()
      double precision :: comp(15,json)
      write(*,'(a)') "# homogenize_convection_zones"
      cs%envelope_metals_no_cno = 0.015d0
      call set_convection_settings(cs)
      do i = 1, nz
         comp(:,i) = xa(:,i)
         comp(13,i) = 1.0d-7*i
      end do
! the renormalization floors, in radiative zones
      comp(1,20) = -1.0d-9
      comp(9,21) = 0.0d0
      convective = .false.
      convective(1:6) = .true.
      convective(31:40) = .true.
      cs%use_extended_composition = .false.
      call set_convection_settings(cs)
      call homogenize_convection_zones(comp, dm, convective, nz)
      write(*,'(a)') "basic composition, every zone (X, Y, Z, x9, x13)"
      write(*,'(5(1pe24.15))') (comp(1,i), comp(2,i), comp(3,i), comp(9,i), comp(13,i), &
           i = 1, nz)
      do i = 1, nz
         comp(:,i) = xa(:,i)
         comp(13,i) = 1.0d-7*i
      end do
      cs%use_extended_composition = .true.
      call set_convection_settings(cs)
      call homogenize_convection_zones(comp, dm, convective, nz)
      write(*,'(a)') "extended composition: x13 in zones 1, 6, 31, 40"
      write(*,'(4(1pe24.15))') comp(13,1), comp(13,6), comp(13,31), comp(13,40)
      cs%use_extended_composition = .false.
end subroutine homogenize_cases

!======================================================================
! The ok/FAIL checks. Each prints its label, the measured error (or 0)
! and ok or FAIL; test_convection.py compares the verdicts, not the
! errors' digits (round-off differs between platforms). Every check has
! a mutation in convection/test/test_convection.mutations.toml.
!======================================================================

! One check: ok when error <= tolerance (a NaN fails).
subroutine check(label, error, tolerance)
      character(len=*), intent(in) :: label
      double precision, intent(in) :: error, tolerance
      character(len=72) :: padded
      padded = label
      if (error <= tolerance) then
         write(*,'(a,1x,1pe10.2,a)') padded, error, "  ok"
      else
         write(*,'(a,1x,1pe10.2,a)') padded, error, "  FAIL"
         num_failed = num_failed + 1
      end if
end subroutine check

! A check of a condition: its error is 0 when it holds, 1 when not.
subroutine check_that(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      call check(label, merge(0.0d0, 1.0d0, condition), 0.5d0)
end subroutine check_that

! The settings of the gradient checks: alpha 1.9, Lsun 3.828e33, the
! cut at log T 6.8, no spots, no adiabatic overshoot.
subroutine gradient_settings()
      cs = convection_settings_t()
      cs%mixing_length_alpha = 1.9d0
      cs%log10_solar_luminosity = log10(3.828d33)
      cs%log10_t_cut = 6.8d0
      call set_convection_settings(cs)
end subroutine gradient_settings

! A synthetic EOS and opacity about the surface convection zone's point
! of gradient_cases (T0 = 10^4.8 K, P0 = 10^9.5, rho0 = rho_scale),
! power laws with constant second derivatives, so that the derivatives
! compute_temperature_gradients takes from them are exact:
!   ln rho = ln rho0 + 0.98 p + (-8) t + 0.25 t^2 + 0.3 t p,
!   del_ad = 0.25 exp(-0.03 t + 0.01 p), c_P = 5e9 exp(0.4 t - 0.1 p),
!   kappa = kappa0 (rho/rho0) (T/T0)^-3.5,
! t = ln(T/T0), p = ln(P/P0).
subroutine synthetic_point(log_t, log_p, rho_scale, kappa0, eos_res, kap_res)
      double precision, intent(in) :: log_t, log_p, rho_scale, kappa0
      type(eos_result), intent(out) :: eos_res
      double precision, intent(out) :: kap_res(num_kap_results)
      double precision :: t, p
      t = ln10*(log_t - 4.8d0)
      p = ln10*(log_p - 9.5d0)
      eos_res%temperature = 10.0d0**log_t
      eos_res%pressure = 10.0d0**log_p
      eos_res%density = rho_scale*exp(0.98d0*p - 8.0d0*t + 0.25d0*t*t + 0.3d0*t*p)
      eos_res%dlnrho_dlnt = -8.0d0 + 0.5d0*t + 0.3d0*p
      eos_res%dlnrho_dlnp = 0.98d0 + 0.3d0*t
      eos_res%dlnrho_dlnt_dt = 0.5d0
      eos_res%dlnrho_dlnt_dp = 0.3d0
      eos_res%adiabatic_gradient = 0.25d0*exp(-0.03d0*t + 0.01d0*p)
      eos_res%adiabatic_gradient_dt = -0.03d0
      eos_res%adiabatic_gradient_dp = 0.01d0
      eos_res%specific_heat_cp = 5.0d9*exp(0.4d0*t - 0.1d0*p)
      eos_res%specific_heat_cp_dt = 0.4d0
      eos_res%specific_heat_cp_dp = -0.1d0
      kap_res = 0.0d0
      kap_res(i_kap) = kappa0*(eos_res%density/rho_scale)*exp(-3.5d0*t)
      kap_res(i_dlnkap_dlnrho) = 1.0d0
      kap_res(i_dlnkap_dlnt) = -3.5d0
end subroutine synthetic_point

! The gradient at a synthetic point (log r 10.82, log m 33.298, L = 1
! Lsun, f_P = 1.02, f_T = 0.99, log Teff 3.76, outside the overshoot
! region): del, del_rad, the three derivatives, v_c and ierr.
subroutine synthetic_gradient(log_t, log_p, log_r, rho_scale, kappa0, grad, grad_rad, &
     grad_ad, derivs, velocity, jerr)
      double precision, intent(in) :: log_t, log_p, log_r, rho_scale, kappa0
      double precision, intent(out) :: grad, grad_rad, grad_ad, derivs(3), velocity
      integer, intent(out) :: jerr
      type(eos_result) :: eos_res
      double precision :: kap_res(num_kap_results)
      logical :: is_convective
      call synthetic_point(log_t, log_p, rho_scale, kappa0, eos_res, kap_res)
      derivs = 0.0d0
      call compute_temperature_gradients(log_t, log_p, eos_res, kap_res, log_r, &
           33.298d0, 1.0d0, grad, grad_rad, derivs(1), derivs(2), derivs(3), velocity, &
           .true., is_convective, 1.02d0, 0.99d0, 3.76d0, .false., jerr)
      grad_ad = eos_res%adiabatic_gradient
end subroutine synthetic_gradient

! The mixing-length solution computed here from the equations of
! temperature_gradients' header, with the cubic solved by bisection: del
! and v_c at a synthetic point (spots as installed in cs).
subroutine independent_mlt(log_t, log_p, log_r, rho_scale, kappa0, grad, velocity)
      double precision, intent(in) :: log_t, log_p, log_r, rho_scale, kappa0
      double precision, intent(out) :: grad, velocity
      type(eos_result) :: eos_res
      double precision :: kap_res(num_kap_results), kappa, grad_rad, excess, g, h_p, phi, &
           pf, a1, a3, v, lo, hi, s, x, log_teff_u
      integer :: k
      call synthetic_point(log_t, log_p, rho_scale, kappa0, eos_res, kap_res)
      kappa = kap_res(i_kap)
      grad_rad = 10.0d0**cdelrl*kappa*1.0d0*10.0d0**cs%log10_solar_luminosity* &
           eos_res%pressure/(10.0d0**cs%log10_gravitational_constant*10.0d0**33.298d0* &
           eos_res%temperature**4)*0.99d0/1.02d0
      s = 1.0d0
      if (cs%spot_filling_factor /= 0.0d0) then
         x = cs%spot_temp_contrast
         if (cs%spot_depth_varies) then
            log_teff_u = 3.76d0 - 0.25d0*log10(cs%spot_filling_factor*x**4 + 1.0d0 - &
                 cs%spot_filling_factor)
            x = 1.0d0 - (1.0d0 - cs%spot_temp_contrast)*10.0d0**log_teff_u/eos_res%temperature
         end if
         s = cs%spot_filling_factor*x**4 + 1.0d0 - cs%spot_filling_factor
      end if
      excess = grad_rad/s - eos_res%adiabatic_gradient
      g = 10.0d0**cs%log10_gravitational_constant*10.0d0**33.298d0/(10.0d0**log_r)**2
      h_p = eos_res%pressure/(eos_res%density*g)
      phi = cs%mixing_length_alpha*kappa*eos_res%density*h_p
      pf = 1.0d0/(1.0d0 + phi*phi/3.0d0)
      a1 = 16.0d0*sqrt(2.0d0)*csig*eos_res%temperature**3*kappa*pf/ &
           (eos_res%specific_heat_cp*sqrt(g*(-eos_res%dlnrho_dlnt)*excess/h_p))
      a3 = 0.75d0*phi*phi*pf/a1
      lo = 0.0d0
      hi = 1.0d0
      do k = 1, 200
         v = 0.5d0*(lo + hi)
         if (a3*v**3 + v*v + a1*v - 1.0d0 > 0.0d0) then
            hi = v
         else
            lo = v
         end if
      end do
      v = 0.5d0*(lo + hi)
      grad = eos_res%adiabatic_gradient + excess*v*(v + a1)
      velocity = v*sqrt(g*(-eos_res%dlnrho_dlnt)*cs%mixing_length_alpha**2*h_p*excess/8.0d0)
end subroutine independent_mlt

! The errors of the three derivatives against centred differences of
! ln del in ln T, ln P and ln r (step h = 1e-5), each relative to
! max(|derivative|, 1e-3), at a synthetic point.
subroutine derivative_errors(rho_scale, kappa0, errors)
      double precision, intent(in) :: rho_scale, kappa0
      double precision, intent(out) :: errors(3)
      double precision, parameter :: h = 1.0d-5, log_t0 = 4.8d0, log_p0 = 9.5d0, &
           log_r0 = 10.82d0
      double precision :: grad, grad_rad, grad_ad, derivs(3), velocity, gp, gm, &
           dummy(3), numeric
      integer :: jerr, k
      double precision :: shift(3)
      call synthetic_gradient(log_t0, log_p0, log_r0, rho_scale, kappa0, grad, grad_rad, &
           grad_ad, derivs, velocity, jerr)
      do k = 1, 3
         shift = 0.0d0
         shift(k) = h/ln10
         call synthetic_gradient(log_t0 + shift(1), log_p0 + shift(2), log_r0 + shift(3), &
              rho_scale, kappa0, gp, grad_rad, grad_ad, dummy, velocity, jerr)
         call synthetic_gradient(log_t0 - shift(1), log_p0 - shift(2), log_r0 - shift(3), &
              rho_scale, kappa0, gm, grad_rad, grad_ad, dummy, velocity, jerr)
         numeric = (log(gp) - log(gm))/(2.0d0*h)
         errors(k) = abs(numeric - derivs(k))/max(abs(derivs(k)), 1.0d-3)
      end do
end subroutine derivative_errors

subroutine gradient_checks()
      double precision :: grad, grad_rad, grad_ad, derivs(3), velocity, expected, &
           expected_velocity, errors(3), worst, ratio(2), log_rho_point
      type(eos_result) :: eos_res
      double precision :: kap_res(num_kap_results), actual, radiative, dt, dp, dr, vconv
      logical :: is_convective, inside
      integer :: jerr, k
      write(*,'(a)') "# compute_temperature_gradients: checks"
      call gradient_settings()
! del_rad from its formula
      call synthetic_gradient(4.8d0, 9.5d0, 10.82d0, 1.0d-5, 5.0d3, grad, grad_rad, grad_ad, &
           derivs, velocity, jerr)
      expected = 10.0d0**cdelrl*5.0d3*10.0d0**cs%log10_solar_luminosity*10.0d0**9.5d0/ &
           (10.0d0**cs%log10_gravitational_constant*10.0d0**33.298d0*10.0d0**(4.0d0*4.8d0))* &
           0.99d0/1.02d0
      call check("del_rad = 10^cdelrl kappa L P/(G m T^4) f_T/f_P", &
           abs(grad_rad/expected - 1.0d0), 1.0d-13)
! a radiative shell, and the adiabatic overshoot region
      call synthetic_point(6.5d0, 16.0d0, 5.0d0, 1.0d0, eos_res, kap_res)
      kap_res(i_kap) = 1.0d-3
      cs%adiabatic_overshoot = .true.
      call set_convection_settings(cs)
      do k = 1, 2
         inside = (k == 1)
         call compute_temperature_gradients(6.5d0, 16.0d0, eos_res, kap_res, 10.2d0, 33.2d0, &
              0.2d0, actual, radiative, dt, dp, dr, vconv, .true., is_convective, 1.0d0, &
              1.0d0, 3.76d0, inside, jerr)
         if (inside) then
            call check_that("radiative shell in the overshoot region: del = del_ad", &
                 .not. is_convective .and. actual == eos_res%adiabatic_gradient .and. &
                 radiative < eos_res%adiabatic_gradient .and. vconv == 0.0d0)
         else
            call check_that("radiative shell outside it: del = del_rad, v_c = 0", &
                 .not. is_convective .and. actual == radiative .and. vconv == 0.0d0)
         end if
      end do
      call gradient_settings()
! above log10_t_cut
      call synthetic_point(7.1d0, 17.0d0, 80.0d0, 3.0d3, eos_res, kap_res)
      kap_res(i_kap) = 30.0d0
      call compute_temperature_gradients(7.1d0, 17.0d0, eos_res, kap_res, 9.8d0, 32.8d0, &
           2.5d0, actual, radiative, dt, dp, dr, vconv, .true., is_convective, 1.0d0, &
           1.0d0, 3.76d0, .false., jerr)
      call check_that("above log10_t_cut: del = del_ad, del_ad's derivatives", &
           is_convective .and. actual == eos_res%adiabatic_gradient .and. &
           dt == eos_res%adiabatic_gradient_dt .and. dp == eos_res%adiabatic_gradient_dp &
           .and. dr == 0.0d0 .and. jerr == 0)
! the mixing-length solution against one computed here
      worst = 0.0d0
      do k = 0, 3
         log_rho_point = -21.0d0 + 6.0d0*k
         call synthetic_gradient(4.8d0, 9.5d0, 10.82d0, 10.0d0**log_rho_point, 1.0d0, grad, &
              grad_rad, grad_ad, derivs, velocity, jerr)
         call independent_mlt(4.8d0, 9.5d0, 10.82d0, 10.0d0**log_rho_point, 1.0d0, expected, &
              expected_velocity)
         worst = max(worst, abs(grad/expected - 1.0d0), abs(velocity/expected_velocity - 1.0d0))
      end do
      call check("MLT: del and v_c as the equations give, cubic by bisection", worst, 1.0d-10)
! del between del_ad and del_rad, over 24 decades of density
      worst = 0.0d0
      do k = 0, 24
         log_rho_point = -21.0d0 + k
         call synthetic_gradient(4.8d0, 9.5d0, 10.82d0, 10.0d0**log_rho_point, 1.0d0, grad, &
              grad_rad, grad_ad, derivs, velocity, jerr)
         if (grad < grad_ad .or. grad > grad_rad .or. jerr /= 0) worst = 1.0d0
      end do
      call check("MLT: del_ad <= del <= del_rad from rho = 1e-21 to 1e3", worst, 0.0d0)
! the efficient limit: (del - del_ad)/Delta ~ rho^(-1/3)
      do k = 1, 2
         log_rho_point = 2.0d0 + k
         call synthetic_gradient(4.8d0, 9.5d0, 10.82d0, 10.0d0**log_rho_point, 1.0d0, grad, &
              grad_rad, grad_ad, derivs, velocity, jerr)
         ratio(k) = (grad - grad_ad)/(grad_rad - grad_ad)
      end do
      call check("MLT efficient limit: (del - del_ad)/Delta ~ rho^-1/3", &
           abs(log10(ratio(2)/ratio(1)) + 1.0d0/3.0d0), 1.0d-2)
! the inefficient limit: (del_rad - del)/Delta ~ rho^2
      do k = 1, 2
         log_rho_point = -22.0d0 + k
         call synthetic_gradient(4.8d0, 9.5d0, 10.82d0, 10.0d0**log_rho_point, 1.0d0, grad, &
              grad_rad, grad_ad, derivs, velocity, jerr)
         ratio(k) = (grad_rad - grad)/(grad_rad - grad_ad)
      end do
      call check("MLT inefficient limit: (del_rad - del)/Delta ~ rho^2", &
           abs(log10(ratio(2)/ratio(1)) - 2.0d0), 1.0d-2)
! the derivatives against centred differences
      call derivative_errors(1.0d-17, 1.0d0, errors)
      call check("derivatives vs centred differences: dln del/dln T", errors(1), 1.0d-6)
      call check("derivatives vs centred differences: dln del/dln P", errors(2), 1.0d-6)
      call check("derivatives vs centred differences: dln del/dln r", errors(3), 1.0d-6)
! with spots: the r derivative holds; the T and P ones miss the spot
! factor (M242), so their errors are printed, not checked
      cs%spot_filling_factor = 0.3d0
      cs%spot_temp_contrast = 0.8d0
      do k = 1, 2
         cs%spot_depth_varies = (k == 2)
         call set_convection_settings(cs)
         call synthetic_gradient(4.8d0, 9.5d0, 10.82d0, 1.0d-17, 1.0d0, grad, grad_rad, &
              grad_ad, derivs, velocity, jerr)
         call independent_mlt(4.8d0, 9.5d0, 10.82d0, 1.0d-17, 1.0d0, expected, &
              expected_velocity)
         call derivative_errors(1.0d-17, 1.0d0, errors)
         if (k == 1) then
            call check("spots, constant contrast: del and v_c as the equations give", &
                 max(abs(grad/expected - 1.0d0), abs(velocity/expected_velocity - 1.0d0)), &
                 1.0d-10)
            call check("spots, constant contrast: dln del/dln r vs centred difference", &
                 errors(3), 1.0d-6)
            write(*,'(a,2(1pe10.2))') "  known gap (M242), dln del/dln T and dln P:", &
                 errors(1), errors(2)
         else
            call check("spots, depth-varying contrast: del and v_c as the equations give", &
                 max(abs(grad/expected - 1.0d0), abs(velocity/expected_velocity - 1.0d0)), &
                 1.0d-10)
            call check("spots, depth-varying contrast: dln del/dln r vs centred difference", &
                 errors(3), 1.0d-6)
            write(*,'(a,2(1pe10.2))') "  known gap (M242), dln del/dln T and dln P:", &
                 errors(1), errors(2)
         end if
      end do
      call gradient_settings()
! the error contract: a density rising with T has no real root
      call synthetic_point(4.8d0, 9.5d0, 1.0d-5, 5.0d3, eos_res, kap_res)
      eos_res%dlnrho_dlnt = 0.5d0
      call compute_temperature_gradients(4.8d0, 9.5d0, eos_res, kap_res, 10.82d0, 33.298d0, &
           1.0d0, actual, radiative, dt, dp, dr, vconv, .true., is_convective, 1.0d0, &
           1.0d0, 3.76d0, .false., jerr)
      call check_that("a negative radicand returns convection_unstable_radicand", &
           jerr == convection_unstable_radicand)
end subroutine gradient_checks

! find_convection_zones on the 40-zone model with the flags of ranges,
! silently (the overshoot cases log their regions to the run log).
subroutine find_zones(ranges)
      integer, intent(in) :: ranges(:)
      integer :: k
      convective = .false.
      do k = 1, size(ranges), 2
         convective(ranges(k):ranges(k+1)) = .true.
      end do
      mix_bounds = 0
      mix_bounds_no_os = 0
      rad_bounds = 0
      ov_bottom = -7
      ov_top = -7
      call set_convection_settings(cs)
      call find_convection_zones(xa, log_rho, log_p, log_r, log_m, log_t, convective, nz, &
           rad_bounds, mix_bounds, mix_bounds_no_os, core_edge, env_edge, n_rad, n_mix, &
           n_mix_no_os, ov_bottom, ov_top, ierr)
end subroutine find_zones

! Every shell 1..nz lies in exactly one region, mixed or radiative.
logical function regions_tile_the_model()
      integer :: covered(nz), k
      covered = 0
      do k = 1, n_mix
         covered(mix_bounds(k,1):mix_bounds(k,2)) = covered(mix_bounds(k,1):mix_bounds(k,2)) + 1
      end do
      do k = 1, n_rad
         covered(rad_bounds(k,1):rad_bounds(k,2)) = covered(rad_bounds(k,1):rad_bounds(k,2)) + 1
      end do
      regions_tile_the_model = all(covered == 1)
end function regions_tile_the_model

subroutine zone_checks()
      integer :: ranges(26), k
      write(*,'(a)') "# find_convection_zones: checks"
      cs%core_overshoot = .false.
      cs%envelope_overshoot = .false.
      cs%intermediate_overshoot = .false.
      cs%rotation_active = .false.
      cs%instability_transport_active = .false.
! a single convective shell is not a zone
      call find_zones([1, 5, 15, 15, 18, 21, 31, 40])
      call check_that("zones: runs of more than one shell; a single shell is not one", &
           ierr == 0 .and. n_mix == 3 .and. all(mix_bounds(1:3,1) == [1, 18, 31]) .and. &
           all(mix_bounds(1:3,2) == [5, 21, 40]) .and. core_edge == 5 .and. env_edge == 31)
      call check_that("zones: the mixed and radiative regions tile the model", &
           regions_tile_the_model())
! twelve zones fill the twelve rows; a thirteenth is an error
      do k = 1, 12
         ranges(2*k-1) = 3*k - 2
         ranges(2*k) = 3*k - 1
      end do
      call find_zones(ranges(1:24))
      call check_that("zones: twelve zones fill the bounds' twelve rows", &
           ierr == 0 .and. n_mix == 12 .and. mix_bounds(12,1) == 34 .and. &
           regions_tile_the_model())
      ranges(25) = 37
      ranges(26) = 38
      call find_zones(ranges)
      call check_that("zones: a thirteenth zone returns convection_too_many_zones", &
           ierr == convection_too_many_zones)
end subroutine zone_checks

! H_P at shell k with the Yale EOS's density for its P, T and
! composition, P r^2/(G m rho), computed here.
double precision function eos_scale_height(k)
      integer, intent(in) :: k
      type(eos_result) :: eos_res
      integer :: saha_iterations, jerr
      saha_iterations = 0
      eos_res = eos_result()
      eos_res%log10_density = log_rho(k)
      call eos_get(log_t(k), log_p(k), xa(1,k), xa(3,k), eos_res, .false., .true., &
           saha_iterations, composition_at_zone=xa(:,k), ierr=jerr)
      eos_scale_height = 10.0d0**log_p(k)*(10.0d0**log_r(k))**2/(10.0d0**( &
           cs%log10_gravitational_constant + log_m(k))*10.0d0**eos_res%log10_density)
end function eos_scale_height

! An upper edge moved from old to new by distance d: new is the last
! shell within d above old (r(new) <= r(old) + d < r(new+1)), or the
! surface.
logical function moved_up_by(old, new, d)
      integer, intent(in) :: old, new
      double precision, intent(in) :: d
      double precision :: target
      target = 10.0d0**log_r(old) + d
      moved_up_by = 10.0d0**log_r(new) <= target*(1.0d0 + 1.0d-12)
      if (new < nz) moved_up_by = moved_up_by .and. 10.0d0**log_r(new+1) > target
end function moved_up_by

! A lower edge moved down by d: new is the deepest shell within d below
! old (r(new-1) < r(old) - d <= r(new)), or the centre.
logical function moved_down_by(old, new, d)
      integer, intent(in) :: old, new
      double precision, intent(in) :: d
      double precision :: target
      target = 10.0d0**log_r(old) - d
      moved_down_by = 10.0d0**log_r(new) >= target*(1.0d0 - 1.0d-12)
      if (new > 1) moved_down_by = moved_down_by .and. 10.0d0**log_r(new-1) < target
end function moved_down_by

subroutine overshoot_checks()
      double precision :: h_p, expected, d
      write(*,'(a)') "# overshoot: checks"
      cs%core_overshoot = .false.
      cs%envelope_overshoot = .false.
      cs%intermediate_overshoot = .false.
      cs%limit_core_overshoot = .false.
      cs%rotation_active = .false.
      cs%instability_transport_active = .false.
      cs%overshoot_alpha_envelope = 0.7d0
      call set_convection_settings(cs)
! H_P and the envelope overshoot depth
      h_p = pressure_scale_height(log_p(31), log_r(31), log_rho(31), log_m(31))
      expected = 10.0d0**log_p(31)*(10.0d0**log_r(31))**2/(10.0d0**( &
           cs%log10_gravitational_constant + log_m(31))*10.0d0**log_rho(31))
      call check("pressure_scale_height = P r^2/(G m rho)", abs(h_p/expected - 1.0d0), 1.0d-13)
      call check_that("envelope_overshoot_depth = 0 with envelope overshoot off", &
           envelope_overshoot_depth(log_p, log_r, log_rho, log_m, 31) == 0.0d0)
      cs%envelope_overshoot = .true.
      call set_convection_settings(cs)
      call check("envelope_overshoot_depth = alpha_env H_P at the base", abs( &
           envelope_overshoot_depth(log_p, log_r, log_rho, log_m, 31)/(0.7d0*h_p) - 1.0d0), &
           1.0d-14)
! the core: up by alpha_core H_P
      cs%envelope_overshoot = .false.
      cs%core_overshoot = .true.
      cs%overshoot_alpha_core = 2.5d0
      call find_zones([1, 5, 31, 40])
      d = 2.5d0*eos_scale_height(5)
      call check_that("overshoot, core: the last shell within alpha_core H_P above", &
           ierr == 0 .and. mix_bounds(1,2) > 5 .and. moved_up_by(5, mix_bounds(1,2), d) &
           .and. mix_bounds(2,1) == 31)
! the core with the Woo & Demarque limit, min(alpha H_P, beta r)
      cs%limit_core_overshoot = .true.
      cs%core_overshoot_limit = 0.5d0
      call find_zones([1, 5, 31, 40])
      d = min(2.5d0*eos_scale_height(5), 0.5d0*10.0d0**log_r(5))
      call check_that("overshoot, core: limited to beta r (Woo & Demarque)", &
           ierr == 0 .and. d < 2.5d0*eos_scale_height(5) .and. &
           moved_up_by(5, mix_bounds(1,2), d))
      cs%limit_core_overshoot = .false.
! the envelope: down by alpha_env H_P; the overshoot region below it
! (H_P at shell 31 is a tenth of a shell's width: alpha_env 30)
      cs%core_overshoot = .false.
      cs%envelope_overshoot = .true.
      cs%overshoot_alpha_envelope = 30.0d0
      call find_zones([1, 5, 31, 40])
      d = 30.0d0*eos_scale_height(31)
      call check_that("overshoot, envelope: the deepest shell within alpha_env H_P", &
           ierr == 0 .and. mix_bounds(2,1) < 31 .and. moved_down_by(31, mix_bounds(2,1), d) &
           .and. mix_bounds(1,2) == 5)
      call check_that("overshoot, envelope: the region below it, bottom..old edge", &
           ov_bottom == mix_bounds(2,1) .and. ov_top == 31)
! an intermediate zone: both ways by alpha_intermediate H_P
      cs%envelope_overshoot = .false.
      cs%intermediate_overshoot = .true.
      cs%overshoot_alpha_intermediate = 3.0d0
      call find_zones([1, 5, 18, 21, 31, 40])
      call check_that("overshoot, intermediate: both edges by alpha_int H_P", ierr == 0 &
           .and. mix_bounds(2,1) < 18 .and. mix_bounds(2,2) > 21 .and. &
           moved_down_by(18, mix_bounds(2,1), 3.0d0*eos_scale_height(18)) .and. &
           moved_up_by(21, mix_bounds(2,2), 3.0d0*eos_scale_height(21)))
! the rotating rule: at least one shell
      cs%intermediate_overshoot = .false.
      cs%envelope_overshoot = .true.
      cs%overshoot_alpha_envelope = 1.0d-4
      cs%rotation_active = .true.
      cs%instability_transport_active = .true.
      call find_zones([1, 5, 31, 40])
      call check_that("overshoot, rotating: at least one shell below the envelope", &
           ierr == 0 .and. mix_bounds(2,1) == 30)
      cs%rotation_active = .false.
      cs%instability_transport_active = .false.
! the merge after overshoot: the regions that meet become one
      cs%envelope_overshoot = .false.
      cs%core_overshoot = .true.
      cs%intermediate_overshoot = .true.
      cs%overshoot_alpha_core = 2.5d0
      cs%overshoot_alpha_intermediate = 0.8d0
      call find_zones([1, 5, 9, 12, 31, 40])
      call check_that("overshoot merge: met regions become one, first to last shell", &
           ierr == 0 .and. n_mix == 2 .and. n_mix_no_os == 3 .and. mix_bounds(1,1) == 1 &
           .and. mix_bounds(1,2) == 13 .and. mix_bounds(2,1) == 31 .and. &
           regions_tile_the_model())
      cs%core_overshoot = .false.
      cs%intermediate_overshoot = .false.
      cs%overshoot_alpha_core = 0.6d0
      call set_convection_settings(cs)
end subroutine overshoot_checks

subroutine homogenize_checks()
      double precision :: comp(15,json), saved(15,json), before(15,2), after(15,2), worst
      integer :: zone(2,2), j, s
      write(*,'(a)') "# homogenize_convection_zones: checks"
      cs%envelope_metals_no_cno = 0.015d0
      cs%use_extended_composition = .true.
      call set_convection_settings(cs)
      zone(:,1) = [1, 6]
      zone(:,2) = [31, 40]
      do i = 1, nz
         comp(:,i) = xa(:,i)
         comp(13,i) = 1.0d-7*i
      end do
      saved = comp
      convective = .false.
      convective(1:6) = .true.
      convective(31:40) = .true.
      do j = 1, 2
         do s = 1, 15
            before(s,j) = sum(comp(s,zone(1,j):zone(2,j))*dm(zone(1,j):zone(2,j)))
         end do
      end do
      call homogenize_convection_zones(comp, dm, convective, nz)
      worst = 0.0d0
      do j = 1, 2
         do s = 1, 15
            after(s,j) = sum(comp(s,zone(1,j):zone(2,j))*dm(zone(1,j):zone(2,j)))
            if (before(s,j) /= 0.0d0) worst = max(worst, abs(after(s,j)/before(s,j) - 1.0d0))
         end do
      end do
      call check("homogenize: each species' mass in each zone is kept", worst, 1.0d-14)
      worst = 0.0d0
      do j = 1, 2
         do i = zone(1,j), zone(2,j)
            worst = max(worst, maxval(abs(comp(:,i) - comp(:,zone(1,j)))))
         end do
      end do
      call check("homogenize: each zone is uniform", worst, 0.0d0)
      call check_that("homogenize: the radiative shells are untouched", &
           all(comp(:,7:30) == saved(:,7:30)))
! a uniform composition stays (to round-off)
      do i = 1, nz
         comp(:,i) = [0.7d0, 0.28d0, 0.02d0, 1.0d-5, 3.0d-3, 4.0d-5, 1.0d-3, 4.0d-6, &
              9.0d-3, 4.0d-6, 2.0d-5, 2.0d-5, 1.0d-9, 1.0d-9, 1.0d-10]
      end do
      call homogenize_convection_zones(comp, dm, convective, nz)
      worst = 0.0d0
      do i = 1, nz
         worst = max(worst, maxval(abs(comp(:,i)/[0.7d0, 0.28d0, 0.02d0, 1.0d-5, 3.0d-3, &
              4.0d-5, 1.0d-3, 4.0d-6, 9.0d-3, 4.0d-6, 2.0d-5, 2.0d-5, 1.0d-9, 1.0d-9, &
              1.0d-10] - 1.0d0)))
      end do
      call check("homogenize: a uniform composition stays (relative)", worst, 1.0d-15)
! the floors
      do i = 1, nz
         comp(:,i) = xa(:,i)
      end do
      comp(1,20) = -1.0d-9
      comp(3,22) = 0.9d0
      comp(9,21) = 0.0d0
      call homogenize_convection_zones(comp, dm, convective, nz)
      call check_that("homogenize's floors: X >= 0, Z <= 1 - X, O16 >= 0.99e-3 CNO", &
           comp(1,20) == 0.0d0 .and. comp(3,22) == 1.0d0 - comp(1,22) .and. &
           abs(comp(9,21)/(0.99d-3*(comp(3,21) - 0.015d0)) - 1.0d0) < 1.0d-15)
      cs%use_extended_composition = .false.
      call set_convection_settings(cs)
end subroutine homogenize_checks

! The opacity tables semiconvection needs: OPAL95 (GS98) only, at Z =
! 0.016232, as test_kap loads them.
subroutine boot_opacity()
      type(kap_settings_t) :: kap_cfg
      type(kap_table_paths_t) :: kap_paths
      kap_cfg%use_opal95_tables = .true.
      kap_cfg%opal95_single_table_z = sc_metals
      kap_cfg%molecular_opacity_logt_min = 4.0d0
      kap_cfg%molecular_opacity_logt_max = 4.1d0
      kap_paths%opal95 = trim(yrec_input)//"/opacity/atomic/legacy/opal95/GS98.OP17"
      call kap_init(kap_cfg, kap_paths, 0.70d0, sc_metals, ierr)
      if (ierr /= 0) then
         write(*,'(a)') "test_convection: FAIL (kap_init error)"
         stop 1
      end if
end subroutine boot_opacity

! del_rad and del_ad at shell k of the semiconvection model, with the
! composition of shell comp_shell: the EOS, the opacity and the
! gradient, as extend_by_semiconvection evaluates a shell.
subroutine sc_gradients(k, comp_shell, luminosity, grad_rad, grad_ad, log_density)
      integer, intent(in) :: k, comp_shell
      double precision, intent(in) :: luminosity
      double precision, intent(out) :: grad_rad, grad_ad, log_density
      type(eos_result) :: eos_res
      double precision :: kap_res(num_kap_results), actual, dt, dp, dr, vconv
      integer :: saha_iterations, jerr
      logical :: is_convective
      saha_iterations = 0
      eos_res = eos_result()
      eos_res%log10_density = sc_log_rho(k)
      call eos_get(sc_log_t(k), sc_log_p(k), sc_xa(1,comp_shell), sc_xa(3,comp_shell), &
           eos_res, .false., .true., saha_iterations, &
           composition_at_zone=sc_xa(:,comp_shell), ierr=jerr)
      call kap_get(eos_res%log10_density, sc_log_t(k), sc_xa(1,comp_shell), &
           sc_xa(3,comp_shell), kap_res, eos_res%ion_fraction, ierr=jerr)
      call compute_temperature_gradients(sc_log_t(k), sc_log_p(k), eos_res, kap_res, &
           sc_log_r(k), sc_log_m(k), luminosity, actual, grad_rad, dt, dp, dr, vconv, &
           .false., is_convective, 1.0d0, 1.0d0, 3.9d0, .false., jerr)
      grad_ad = eos_res%adiabatic_gradient
      log_density = eos_res%log10_density
end subroutine sc_gradients

! A 40-shell model of a convective core (shells 1..sc_core, X = 0.30)
! under a hydrogen gradient (X from 0.45 up to 0.70), with L chosen so
! that each outer shell's del_rad, with its own composition, is
! sc_target(k) times its del_ad (times `scale`): 1.30, 1.20 and 1.05 in
! the first three shells outside the core, then 0.90 falling by 0.1 a
! shell to 0.3. With the core's composition the first shell's del_rad
! is f ~ 0.9 times its own, so semiconvection adds two shells: the third
! has del_rad > del_ad but f del_rad < del_ad.
subroutine build_semiconvective_model(scale)
      double precision, intent(in) :: scale
      double precision, parameter :: first_targets(3) = [1.30d0, 1.20d0, 1.05d0]
      double precision :: f, grad_rad, grad_ad, log_density
      integer :: k
      sc_xa = 0.0d0
      do k = 1, nz
         f = dble(k - 1)/dble(nz - 1)
         sc_log_t(k) = 7.50d0 - 0.40d0*f
         sc_log_p(k) = 17.0d0 - 2.0d0*f
         sc_log_r(k) = 9.8d0 + 0.8d0*sqrt(f)
         sc_log_m(k) = 33.0d0 + 0.8d0*sqrt(f)
         if (k <= sc_core) then
            sc_xa(1,k) = 0.30d0
         else
            sc_xa(1,k) = 0.45d0 + 0.25d0*dble(k - sc_core - 1)/dble(nz - sc_core - 1)
         end if
         sc_xa(3,k) = sc_metals
         sc_xa(2,k) = 1.0d0 - sc_xa(1,k) - sc_xa(3,k)
         sc_log_rho(k) = sc_log_p(k) - sc_log_t(k) - 7.92d0 + &
              log10(4.0d0/(3.0d0 + 5.0d0*sc_xa(1,k) - sc_xa(3,k)))
         call sc_gradients(k, k, 1.0d0, grad_rad, grad_ad, log_density)
         sc_log_rho(k) = log_density
         if (k <= sc_core) then
            sc_target(k) = 1.5d0
         else if (k <= sc_core + 3) then
            sc_target(k) = scale*first_targets(k - sc_core)
         else
            sc_target(k) = scale*max(0.90d0 - 0.1d0*dble(k - sc_core - 4), 0.3d0)
         end if
         sc_lum(k) = sc_target(k)*grad_ad/grad_rad
      end do
end subroutine build_semiconvective_model

! extend_by_semiconvection on the model with the given zone list and
! timestep (s); the model's densities are not touched (a copy).
subroutine semiconvect(bounds, num, timestep)
      integer, intent(inout) :: bounds(:,:), num
      double precision, intent(in) :: timestep
      double precision :: log_density(json)
      log_density = sc_log_rho
      call extend_by_semiconvection(timestep, sc_xa, log_density, sc_lum, sc_log_p, &
           sc_log_r, sc_log_m, sc_log_t, nz, bounds, num, 3.9d0, ierr)
end subroutine semiconvect

subroutine semiconvection_checks()
      integer :: bounds(max_convective_zones,2), num, edge, k
      double precision :: grad_mix, grad_local, grad_rad, grad_ad, log_density, ratio
      logical :: holds
      write(*,'(a)') "# extend_by_semiconvection: checks (OPAL95 opacities)"
      call boot_opacity()
      cs = convection_settings_t()
      cs%mixing_length_alpha = 1.9d0
      cs%log10_solar_luminosity = log10(3.828d33)
      cs%solar_luminosity_cgs = 3.828d33
      cs%log10_t_cut = 6.8d0
      call set_convection_settings(cs)
      call build_semiconvective_model(1.0d0)
! the core's edge moves out while f del_rad > del_ad
      bounds = 0
      bounds(1,:) = [1, sc_core]
      num = 1
      call semiconvect(bounds, num, 1.0d14)
      edge = bounds(1,2)
      call check_that("semiconvection: the core's upper edge moves out", &
           ierr == 0 .and. num == 1 .and. edge > sc_core + 1)
      call sc_gradients(sc_core + 1, sc_core, sc_lum(sc_core + 1), grad_mix, grad_ad, log_density)
      call sc_gradients(sc_core + 1, sc_core + 1, sc_lum(sc_core + 1), grad_local, grad_ad, &
           log_density)
      ratio = grad_mix/grad_local
      holds = ratio < 1.0d0 .and. edge < nz
      do k = sc_core + 2, min(edge + 1, nz)
         call sc_gradients(k, k, sc_lum(k), grad_rad, grad_ad, log_density)
         if (k <= edge) holds = holds .and. ratio*grad_rad >= grad_ad
         if (k == edge + 1) holds = holds .and. ratio*grad_rad < grad_ad
      end do
      call check_that("semiconvection: f del_rad >= del_ad in the shells added, not next", &
           holds)
! a timestep too short for the instability: only the first shell
      bounds = 0
      bounds(1,:) = [1, sc_core]
      num = 1
      call semiconvect(bounds, num, 1.0d0)
      call check_that("semiconvection: a short timestep adds only the first shell", &
           ierr == 0 .and. bounds(1,2) == sc_core + 1)
! the merge: a region the extension overlaps joins; one that only
! shares the new edge shell stays apart
      bounds = 0
      bounds(1,:) = [1, sc_core]
      bounds(2,:) = [edge - 1, edge + 5]
      num = 2
      call semiconvect(bounds, num, 1.0d14)
      call check_that("semiconvection merge: an overlapped region joins the core", &
           ierr == 0 .and. num == 1 .and. bounds(1,1) == 1 .and. bounds(1,2) == edge + 5)
      bounds = 0
      bounds(1,:) = [1, sc_core]
      bounds(2,:) = [edge, edge + 5]
      num = 2
      call semiconvect(bounds, num, 1.0d14)
      call check_that("semiconvection merge: a region sharing the edge shell stays", &
           ierr == 0 .and. num == 2 .and. bounds(1,2) == edge .and. bounds(2,1) == edge)
! a stable region stays
      call build_semiconvective_model(0.7d0)
      bounds = 0
      bounds(1,:) = [1, sc_core]
      num = 1
      call semiconvect(bounds, num, 1.0d14)
      call check_that("semiconvection: stable with the core's composition, it stays", &
           ierr == 0 .and. num == 1 .and. bounds(1,2) == sc_core)
end subroutine semiconvection_checks

end program test_convection

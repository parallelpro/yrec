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
! Semiconvection needs the opacity tables; the configuration matrix's
! semiconvection deck pins it.
program test_convection
      use run_setup, only: set_up_constants_and_tables
      use star_info_lib, only: star_info, json, max_convective_zones, max_radiative_zones
      use eos_lib, only: eos_result
      use kap_lib, only: num_kap_results, i_kap, i_dlnkap_dlnrho, i_dlnkap_dlnt
      use convection_lib, only: convection_settings_t, set_convection_settings, &
           compute_temperature_gradients, find_convection_zones, homogenize_convection_zones
      use io_units, only: run_log_unit
      implicit none
! the test's star (star_info_lib declares no instance); static, as it
! is too large for the stack
      type(star_info), save :: star
! the settings the cases change, installed before each use
      type(convection_settings_t) :: cs
      integer :: ov_bottom, ov_top
      integer, parameter :: nz = 40
      double precision :: xa(15,json), log_rho(json), log_p(json), log_r(json), &
           log_m(json), log_t(json), dm(json)
      logical :: convective(json)
      integer :: rad_bounds(max_radiative_zones,2), mix_bounds(max_convective_zones,2), &
           mix_bounds_no_os(max_convective_zones,2)
      integer :: core_edge, env_edge, n_rad, n_mix, n_mix_no_os, ierr, setups_ierr, i
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
      write(*,'(a)') "test_convection: done"

contains

! setups with the Yale EOS alone (Fermi table) and every table-based
! physics switched off; it sets the constants the gradients use.
subroutine boot_yale_eos()
      star%ctrl%fermi_unit = 15
      run_log_unit = 20
      open(run_log_unit, file="test_convection.short", status="replace")
      star%ctrl%use_mhd_eos = .false.
      star%ctrl%use_scv_eos = .false.
      star%ctrl%use_opal95_eos = .false.
      star%ctrl%use_opal2001_eos = .false.
      star%ctrl%use_opal2006_eos = .false.
      star%job%use_diffusion_z = .false.
      star%ctrl%use_numerical_derivatives = .false.
      star%ctrl%use_debye_huckel_correction = .false.
      star%ctrl%saha_log10t_cutoff = 6.0d0
      star%ctrl%use_opal95_tables = .false.
      star%ctrl%use_opal92_tables = .false.
      star%ctrl%use_laol89_tables = .false.
      star%ctrl%use_alex06_tables = .false.
      star%ctrl%use_alex95_tables = .false.
      star%ctrl%use_kurucz90_tables = .false.
      star%use_two_z_tables = .false.
      star%ctrl%use_pure_z_table = .false.
      star%ctrl%use_conductive_opacity = .false.
      star%job%atm_choice = 0
      star%job%fermi_table_path = trim(yrec_input)//"/eos/yale/FERMI.TAB"
      call set_up_constants_and_tables(star, setups_ierr)
      if (setups_ierr /= 0) then
         write(*,'(a)') "test_convection: FAIL (setups error)"
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

end program test_convection

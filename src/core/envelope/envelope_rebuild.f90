!----------------------------------------------------------------------
! envelope_rebuild
!----------------------------------------------------------------------
! Moving the fitting point outward when the star gains mass: accretion
! thickens the envelope above a fixed fitting point, so
! rebuild_envelope integrates a fresh envelope and turns its points
! between the old and the new fitting point into interior shells.
! (The legacy getnewenv's comment said "when the total mass is being
! decreased"; evolve_step calls it after change_mass accretes.)
! The arithmetic of the steps is envelope_refit_kernel's, star-free;
! this driver reads the star, integrates, calls the kernels with the
! star's arrays and reports.
module envelope_rebuild
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star_info
      implicit none
      private
      public :: rebuild_envelope

contains

! ---------------------------------------------------------------
! rebuild_envelope (formerly getnewenv)
! ---------------------------------------------------------------
! Moves the fitting point outward, to target_envelope_mass (a log(m/M),
! clamped to max_log_fit_point_mass_fraction: an envelope of at least
! about 2.3e-12 of the star), after the star has accreted
! (change_mass's new_atmosphere_fit_needed). In steps:
! 1. integrate the atmosphere and envelope from the surface down to the
!    current fitting point (integrate_envelope_from_surface);
! 2. shift the envelope onto the last interior shell
!    (shift_envelope_onto_shell);
! 3. append the envelope points inside the new fitting mass as interior
!    shells, the last one interpolated onto it (append_envelope_points;
!    a point exactly at the fitting mass is interpolated);
! 4. the new shells' masses (refit_shell_masses);
! 5. with rotation, a first guess at the new shells' rotation
!    (guess_new_shell_rotation), which omega_from_j makes consistent
!    (omega, the shape and the moments of inertia of every zone from j),
!    and the star's angular momentum and rotational energy
!    (regridding's compute_angular_momentum_totals).
! Updates star%log_fit_point_mass_fraction and nz, and asks for a new
! envelope triangle; called from evolve_step. ierr: the integration
! failed, the zones would overflow json, or the mass interpolation was
! degenerate.
! The rotation rates are star%omega itself: the new shells' guess starts
! from the rate of the old outermost shell, and omega_from_j's result is
! what the caller sees.
!   star reads: convective_flag, disk_gate_age_gyr, log_mass,
!     log_total_mass, logP, logR, logRho, logT, nz,
!     physics%cz_rotation_law_exponent, physics%rotation_active,
!     physics%use_extended_composition
!   star writes: dm, env_points, eta_squared, i_rot, j_rot,
!     kinetic_energy_rot, log_fit_point_mass_fraction, m, mean_radius,
!     omega, recompute_envelope_triangle, total_angular_momentum,
!     total_rotational_ke
subroutine rebuild_envelope(s, target_envelope_mass, ierr)
      use omega_profile, only: omega_from_j
      use envelope_integration_kernel, only: max_log_fit_point_mass_fraction
      use envelope_refit, only: append_envelope_points
      use envelope_refit_kernel, only: shift_envelope_onto_shell, refit_shell_masses, &
           guess_new_shell_rotation
      use regridding, only: compute_angular_momentum_totals
      use star_info_lib, only: json, n_species_basic, n_species_extended
      type(star_info), intent(inout) :: s
      real(dp), intent(inout) :: target_envelope_mass
      integer, intent(out) :: ierr
! di_domega of omega_from_j, which this routine does not use
      real(dp) :: qiw(json)
      integer :: last_copied_species, old_num_zones
      real(dp) :: old_fit_mass_fraction

      ierr = 0
      if (s%physics%use_extended_composition) then
         last_copied_species = n_species_extended
      else
         last_copied_species = n_species_basic
      end if
      target_envelope_mass = min(target_envelope_mass,max_log_fit_point_mass_fraction)
      old_num_zones = s%nz
      old_fit_mass_fraction = s%log_fit_point_mass_fraction

      call integrate_envelope_from_surface(s, ierr)
      if (ierr /= 0) return
      s%log_fit_point_mass_fraction = target_envelope_mass
      if (s%nz+s%env_points%num_env_points >= json) then
         call report_rebuild_error(' ERROR: envelope rebuild: the new envelope shells would' &
              //' overflow the json-point arrays')
         ierr = 1
         return
      end if
      call shift_envelope_onto_shell(s%env_points, s%logP(s%nz), s%logRho(s%nz), &
           s%logT(s%nz), s%logR(s%nz), s%log_mass(s%nz), s%log_total_mass)
      call append_envelope_points(s, .false., last_copied_species, ierr)
      if (ierr /= 0) then
         call report_rebuild_error(' ERROR: envelope rebuild: a degenerate interval in mass' &
              //' while appending the envelope shells')
         return
      end if
      call refit_shell_masses(old_num_zones, s%log_mass(1:s%nz), s%log_total_mass, &
           s%m(1:s%nz), s%dm(1:s%nz))
      if (s%physics%rotation_active) then
         call guess_new_shell_rotation(old_num_zones, s%physics%cz_rotation_law_exponent, &
              s%logR(1:s%nz), s%omega(1:s%nz), s%i_rot(1:s%nz), s%j_rot(1:s%nz))
         call omega_from_j(s%logRho,s%j_rot,s%logR,s%log_mass, &
              s%dm,s%convective_flag,s%nz,s%eta_squared, &
              s%i_rot,s%omega,qiw,s%mean_radius,s%disk_gate_age_gyr)
         call compute_angular_momentum_totals(s%j_rot(1:s%nz), s%omega(1:s%nz), &
              s%dm(1:s%nz), s%kinetic_energy_rot(1:s%nz), s%total_angular_momentum, &
              s%total_rotational_ke)
      end if
      call report_new_interior_points(s, old_num_zones)
      s%recompute_envelope_triangle = .true.
      call report_new_envelope_mass(old_fit_mass_fraction, s%log_fit_point_mass_fraction)
end subroutine rebuild_envelope

! ---------------------------------------------------------------
! One integration of the atmosphere and the envelope (atm_get) from the
! star's surface down to the current fitting point, into
! star%env_points: at the model's L and Teff, with R from
! L = 4 pi R^2 sigma Teff^4, g = G M / R^2, the X, Z and composition of
! the outermost shell, and no rotational distortion (f_P = f_T = 1). A
! spotted run with a convective surface starts from the Teff of the
! unspotted surface (ambient_log10_teff, legacy ATEFFL). The envelope
! takes a fixed step in log P, chi_grid_scale(ichi_dp_env_max), for
! this call only; the first step that would cross the pressure of the
! outermost shell stops at it. No vertex is fitted. ierr is atm_get's.
!   star reads: convective_flag, log_L_lsun, log_Teff, log_total_mass,
!     logP, nz, physics%chi_grid_scale, physics%spot_filling_factor,
!     physics%spot_temp_contrast, xa
!   star writes: none
subroutine integrate_envelope_from_surface(s, ierr)
      use envelope_integration, only: atm_get, envint_step_config, fixed_envint_step, &
           surface_log10_radius_cm, ambient_log10_teff
      use star_info_lib, only: i_h1, i_metals
      use control_slots, only: ichi_dp_env_max
      use physical_constants, only: cgl, ln10
      use math_lib
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      real(dp), parameter :: no_distortion = 1.0D0
      type(envint_step_config) :: env_steps
      integer :: saha_iterations
      real(dp) :: luminosity_linear, log_radius_surface, log_gravity_surface
      real(dp) :: hydrogen_fraction, metal_fraction, pressure_limit, log_teff

      env_steps = fixed_envint_step(s%physics%chi_grid_scale(ichi_dp_env_max))
      saha_iterations = 0
      luminosity_linear = exp(ln10*s%log_L_lsun)
      log_radius_surface = surface_log10_radius_cm(s, s%log_L_lsun, s%log_Teff)
      log_gravity_surface = cgl + s%log_total_mass - log_radius_surface - log_radius_surface
      hydrogen_fraction = s%xa(i_h1,s%nz)
      metal_fraction = s%xa(i_metals,s%nz)
      pressure_limit = s%logP(s%nz)
      if (s%convective_flag(s%nz).and.s%physics%spot_filling_factor /= 0.0.and. &
           s%physics%spot_temp_contrast /= 1.0) then
         log_teff = ambient_log10_teff(s, s%log_Teff)
      else
         log_teff = s%log_Teff
      end if
      call atm_get(s, luminosity_linear, no_distortion, no_distortion, &
           log_gravity_surface, log_radius_surface, log_teff, &
           hydrogen_fraction, metal_fraction, s%xa(:,s%nz), pressure_limit, &
           saha_iterations, ierr, env_steps=env_steps)
end subroutine integrate_envelope_from_surface

! ---------------------------------------------------------------
! An error's line, to the terminal and the run log.
subroutine report_rebuild_error(text)
      use core_messages, only: core_report, message_length, unwritten_line
      character(len=*), intent(in) :: text
      character(len=message_length) :: lines(1)
      lines = unwritten_line
      lines(1) = text
      call core_report(lines)
end subroutine report_rebuild_error

! ---------------------------------------------------------------
! The fitting mass before and after, to the terminal and the run log.
subroutine report_new_envelope_mass(old_fit_mass_fraction, new_fit_mass_fraction)
      use core_messages, only: core_report, message_length, unwritten_line
      real(dp), intent(in) :: old_fit_mass_fraction, new_fit_mass_fraction
      character(len=message_length) :: lines(2)
      lines = unwritten_line
      write(lines,597) old_fit_mass_fraction, new_fit_mass_fraction
 597  format(5X,'***** NEW ENVELOPE MASS CALCULATED *****'/8X, &
     'OLD SENV ',1PE22.13,'  NEW SENV',E22.13)
      call core_report(lines)
end subroutine report_new_envelope_mass

! ---------------------------------------------------------------
! The new interior shells, from the model's last zone before the
! rebuild (old_num_zones) to its last zone now, to the terminal and the
! run log: log rho, L/Lsun, log P, log r, log(m/M), log T and the
! convective flag. (The heading's "LOG L" and "LOG M" are the F77's.)
!   star reads: convective_flag, log_mass, log_total_mass, logP, logR,
!     logRho, logT, luminosity_lsun, nz
!   star writes: none
subroutine report_new_interior_points(s, old_num_zones)
      use core_messages, only: core_report, message_length, unwritten_line
      type(star_info), intent(in) :: s
      integer, intent(in) :: old_num_zones
      character(len=message_length) :: heading(2)
      character(len=message_length), allocatable :: rows(:)
      integer :: zone_index, last_zone
      last_zone = s%nz
      heading = unwritten_line
      write(heading,910)
 910  format(1X,'NEW INTERIOR POINTS FROM CHANGE IN ENVELOPE MASS'/ &
     ' J,LOG RHO, LOG L, LOG P, LOG R, LOG M, LOG T, CONV T/F')
      call core_report(heading)
      allocate(rows(max(last_zone-old_num_zones+1, 1)))
      rows = unwritten_line
      write(rows,911)(zone_index,s%logRho(zone_index),s%luminosity_lsun(zone_index), &
           s%logP(zone_index),s%logR(zone_index), &
           s%log_mass(zone_index)-s%log_total_mass, &
           s%logT(zone_index),s%convective_flag(zone_index), &
           zone_index = old_num_zones,last_zone)
 911  format(I5,1P6E16.8,L2)
      call core_report(rows)
end subroutine report_new_interior_points

end module envelope_rebuild

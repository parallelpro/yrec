!----------------------------------------------------------------------
! envelope_rebuild
!----------------------------------------------------------------------
! Moving the fitting point outward when the star gains mass: accretion
! thickens the envelope above a fixed fitting point, so
! rebuild_envelope integrates a fresh envelope and turns its points
! between the old and the new fitting point into interior shells.
! (The legacy getnewenv's comment said "when the total mass is being
! decreased"; evolve_step calls it after change_mass accretes.)
module envelope_rebuild
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star_info
      implicit none
      private
      public :: rebuild_envelope

contains

! ---------------------------------------------------------------
! rebuild_envelope (formerly getnewenv)
! --------------------------------------------------------------- Moves
! the fitting point outward, to target_envelope_mass (a log(m/M),
! clamped to max_log_fit_point_mass_fraction), after the star has
! accreted (change_mass's new_atmosphere_fit_needed). It integrates the
! atmosphere and envelope from the surface (atm_get at the current L and
! Teff) down to the current fitting point, shifts the envelope onto the
! last interior shell, and appends the envelope points inside the new
! fitting mass as interior shells (append_envelope_points). Then come
! the new shells' masses and, with rotation, a first guess at their
! rotation, which omega_from_j makes consistent, and the star's angular
! momentum and rotational energy. Updates
! star%log_fit_point_mass_fraction and num_zones; called from
! evolve_step. ierr: the integration failed, the zones would overflow
! json, or the mass interpolation was degenerate.
!   star reads: convective_flag, disk_gate_age_gyr, log_L_lsun,
!     log_mass, log_Teff, log_total_mass, logP, logR, logRho, logT,
!     luminosity_lsun, physics%chi_grid_scale,
!     physics%cz_rotation_law_exponent, physics%rotation_active,
!     physics%spot_filling_factor, physics%spot_temp_contrast,
!     physics%use_extended_composition, xa
!   star writes: dm, env_points, eta_squared, i_rot, j_rot,
!     kinetic_energy_rot, log_fit_point_mass_fraction, m, mean_radius,
!     nz, omega, recompute_envelope_triangle, total_angular_momentum,
!     total_rotational_ke
subroutine rebuild_envelope(s, target_envelope_mass, ierr)
      use omega_profile, only: omega_from_j
      use envelope_integration, only: atm_get, envint_step_config, fixed_envint_step, &
           surface_log10_radius_cm, ambient_log10_teff
      use envelope_integration_kernel, only: max_log_fit_point_mass_fraction
      use envelope_refit, only: append_envelope_points
      use star_info_lib, only: i_h1, i_metals, json, n_species_basic, &
           n_species_extended
      use control_slots, only: ichi_dp_env_max
      use physical_constants, only: cc23, cgl, ln10
      use math_lib
      use core_messages, only: core_stdout, message_length, unwritten_line
      implicit none
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      integer :: jerr_atm
      real(dp), intent(inout) :: target_envelope_mass
! di_domega of omega_from_j, which this routine does not use
      real(dp) :: qiw(json)

      integer :: species_end_index
      integer :: old_num_zones
      real(dp) :: envelope_mass_before
      type(envint_step_config) :: env_steps
      logical :: surface_bc_flag, print_flag
      integer :: saha_iterations
      real(dp) :: luminosity_linear
      real(dp) :: log_radius_surface, log_gravity_surface
      real(dp) :: hydrogen_fraction, metal_fraction
      real(dp) :: fp_surface, ft_surface
      real(dp) :: pressure_limit
      real(dp) :: spot_adjusted_log_teff
      real(dp) :: pressure_offset, density_offset, temperature_offset, &
           radius_offset
      real(dp) :: mass_test, mass_test2
      integer :: num_new_env_points
      integer :: zone_index
      real(dp) :: omega_ref
      real(dp) :: sum_angular_momentum, sum_rotational_ke
      real(dp) :: angular_momentum_shell
      real(dp) :: mass_at_base
! The rotation rates are star%omega itself: the new shells' guess starts
! from the rate of the old outermost shell, and omega_from_j's result
! is what the caller sees. (The legacy getnewenv read a local OMEGA
! that nothing had set, so new shells got omega = j = 0, and the
! rebuilt omega never reached the caller; bugsweep 2026.)

      ierr = 0
      call prepare_envelope_rebuild()
      call set_envelope_surface()
      call integrate_new_envelope()
      if (ierr /= 0) return
! Return an error if the new points could overflow the arrays (legacy
! stop 9999).
      if(s%nz+s%env_points%num_env_points >= json) then
         call write_rebuild_note(' rebuild_envelope: grid overflow (json too small)')
         ierr = 1
         return
      end if
      call shift_envelope_points()
! Append the envelope points inside the new fitting mass as shells,
! interpolating the last one onto it; a point exactly at the fitting
! mass is interpolated (accept_point_at_fit_mass = .false.).
      call append_envelope_points(s, .false., species_end_index, ierr)
! A degenerate mass interval returns an error (legacy stop 9998).
      if (ierr /= 0) then
         call write_rebuild_note(' rebuild_envelope: degenerate mass interpolation')
         return
      end if
      call compute_new_shell_masses()
! With rotation: guess omega and j of the new shells, find omega, the
! shape and the moments of inertia of every zone from j (omega_from_j),
! and sum the angular momentum and rotational energy.
      if(s%physics%rotation_active)then
         call guess_new_shell_rotation()
         call omega_from_j(s%logRho,s%j_rot,s%logR,s%log_mass, &
              s%dm,s%convective_flag,s%nz,s%eta_squared, &
              s%i_rot,s%omega,qiw,s%mean_radius,s%disk_gate_age_gyr)
         call sum_rotation_totals()
      end if
      call print_new_interior_points()
      s%recompute_envelope_triangle = .true.
      call report_new_envelope_mass()
      return

contains

! ---------------------------------------------------------------
! One line to standard output.
subroutine write_rebuild_note(text)
      character(len=*), intent(in) :: text
      character(len=message_length) :: lines(1)
      lines = unwritten_line
      lines(1) = text
      call core_stdout(lines)
end subroutine write_rebuild_note

! ---------------------------------------------------------------
! The fitting mass before and after, to standard output.
subroutine report_new_envelope_mass()
      character(len=message_length) :: lines(2)
      lines = unwritten_line
      write(lines,597)envelope_mass_before,s%log_fit_point_mass_fraction
 597  format(5X,'***** NEW ENVELOPE MASS CALCULATED *****'/8X, &
     'OLD SENV ',1PE22.13,'  NEW SENV',E22.13)
      call core_stdout(lines)
end subroutine report_new_envelope_mass

! ---------------------------------------------------------------
! The species count, the target envelope mass (clamped), the current
! number of zones and envelope mass, and the integration settings.
subroutine prepare_envelope_rebuild()
      if(s%physics%use_extended_composition)then
         species_end_index = n_species_extended
      else
         species_end_index = n_species_basic
      end if
! Keep log_fit_point_mass_fraction at or below
! max_log_fit_point_mass_fraction: an envelope of at least about 2.3e-12
! of the star.
      target_envelope_mass = min(target_envelope_mass,max_log_fit_point_mass_fraction)
! The zone count and the fitting mass before the rebuild.
      old_num_zones = s%nz
      envelope_mass_before = s%log_fit_point_mass_fraction
! The envelope integration takes a fixed step in log P,
! chi_grid_scale(ichi_dp_env_max), handed to atm_get for this call
! only. No boundary vertex is saved, and saha_iterations starts at 0.
      env_steps = fixed_envint_step(s%physics%chi_grid_scale(ichi_dp_env_max))
      surface_bc_flag = .false.
      print_flag = .true.
      saha_iterations = 0
end subroutine prepare_envelope_rebuild

! ---------------------------------------------------------------
! The surface of the envelope integration: L, R, g and the composition,
! no rotational distortion, the pressure limit and the eos composition.
subroutine set_envelope_surface()
! L in Lsun.
      luminosity_linear = exp(ln10*s%log_L_lsun)
! R from L = 4 pi R^2 sigma Teff^4.
      log_radius_surface = surface_log10_radius_cm(s, s%log_L_lsun, s%log_Teff)
! g = G M / R^2.
      log_gravity_surface = cgl + s%log_total_mass - log_radius_surface - log_radius_surface
! X and Z of the outermost shell.
      hydrogen_fraction = s%xa(i_h1,s%nz)
      metal_fraction = s%xa(i_metals,s%nz)
! No rotational distortion: f_P = f_T = 1.
      fp_surface = 1.0D0
      ft_surface = 1.0D0
! The first step that would cross the pressure of the outermost shell
! stops at it.
      pressure_limit = s%logP(s%nz)
end subroutine set_envelope_surface

! ---------------------------------------------------------------
! Integrate the atmosphere and envelope down to the current fitting
! point, then adopt the target envelope mass.
subroutine integrate_new_envelope()
! Integrate down to the current fitting point from the surface L and
! Teff. A spotted run with a convective surface starts from the Teff of
! the unspotted surface (ambient_log10_teff, legacy ATEFFL).
      if(s%convective_flag(s%nz).and.s%physics%spot_filling_factor /= 0.0.and. &
           s%physics%spot_temp_contrast /= 1.0)then
         spot_adjusted_log_teff = ambient_log10_teff(s, s%log_Teff)
      else
         spot_adjusted_log_teff = s%log_Teff
      end if
      call atm_get(s, luminosity_linear,fp_surface,ft_surface,log_gravity_surface, &
           s%log_total_mass,print_flag,surface_bc_flag,pressure_limit, &
           log_radius_surface, &
           spot_adjusted_log_teff,hydrogen_fraction,metal_fraction, &
           s%xa(:,s%nz),saha_iterations,ierr=jerr_atm,env_steps=env_steps)
! An integration failure goes back to evolve_step in ierr.
      if (jerr_atm /= 0) then
         ierr = jerr_atm
         return
      end if
      s%log_fit_point_mass_fraction = target_envelope_mass
end subroutine integrate_new_envelope

! ---------------------------------------------------------------
! Match the envelope points to the last interior point and keep, shifted
! down by one, those above its mass.
subroutine shift_envelope_points()
! Point 1 of the envelope lies at the mass of the last interior shell.
! Its differences from that shell in log P, log rho, log T and log R
! become offsets that make the envelope continuous with the interior.
! Every later point whose log m exceeds the last kept one by more than
! 1e-10 moves down one slot, with the offsets added. Log P, log rho,
! log R, log m, log T, X, Z and the convective flag move (the F77
! getnewenv.f left the flag unshifted, one point off; fixed 2026-10);
! the other arrays are not read for the appended points.
      pressure_offset = s%logP(s%nz) - s%env_points%env_log10_pressure(1)
      density_offset = s%logRho(s%nz) - s%env_points%env_log10_density(1)
      temperature_offset = s%logT(s%nz) - s%env_points%env_log10_temperature(1)
      radius_offset = s%logR(s%nz) - s%env_points%env_log10_radius(1)
      mass_test = s%log_mass(s%nz)
      num_new_env_points = 0
      do zone_index = 1,s%env_points%num_env_points - 1
         mass_test2 = s%log_total_mass+s%env_points%env_log10_mass(zone_index+1)
         if(mass_test2-mass_test > 1.0D-10)then
            s%env_points%env_log10_density(zone_index) = s%env_points%env_log10_density(zone_index+1)+density_offset
            s%env_points%env_log10_pressure(zone_index) = s%env_points%env_log10_pressure(zone_index+1)+pressure_offset
            s%env_points%env_log10_radius(zone_index) = s%env_points%env_log10_radius(zone_index+1)+radius_offset
            s%env_points%env_log10_mass(zone_index) = s%env_points%env_log10_mass(zone_index+1)
            s%env_points%env_log10_temperature(zone_index) = s%env_points%env_log10_temperature(zone_index+1)+temperature_offset
            s%env_points%env_hydrogen_fraction(zone_index) = s%env_points%env_hydrogen_fraction(zone_index+1)
            s%env_points%env_metal_fraction(zone_index) = s%env_points%env_metal_fraction(zone_index+1)
            s%env_points%env_convective_flag(zone_index) = s%env_points%env_convective_flag(zone_index+1)
            mass_test = mass_test2
            num_new_env_points = num_new_env_points + 1
         end if
      end do
      s%env_points%num_env_points = num_new_env_points
end subroutine shift_envelope_points

! ---------------------------------------------------------------
! The (unlogged) masses of the new shells and their mass contents.
subroutine compute_new_shell_masses()
! enclosed_mass (legacy HS1) = 10^log m, and the mass of each zone
! (shell_mass, legacy HS2) = (m(k+1) - m(k-1))/2; the outermost zone
! holds M - (m(nz) + m(nz-1))/2.
      do zone_index = old_num_zones,s%nz
         s%m(zone_index) = exp10(s%log_mass(zone_index))
      end do
      do zone_index = old_num_zones,s%nz-1
         s%dm(zone_index) = 0.5D0*(s%m(zone_index+1)-s%m(zone_index-1))
      end do
      mass_at_base = 0.5D0*(s%m(s%nz)+s%m(s%nz-1))
      s%dm(s%nz) = exp10(s%log_total_mass) - mass_at_base
end subroutine compute_new_shell_masses

! ---------------------------------------------------------------
! The first guess at the rotation of the new shells, from the outermost
! old shell by the convection-zone law walpcz: constant j, solid body,
! or omega ~ r^walpcz. Each shell's specific moment of inertia is
! the thin-shell i = (2/3) r^2, and j = omega i.
subroutine guess_new_shell_rotation()
         if(s%physics%cz_rotation_law_exponent <= -2.0D0)then
! walpcz <= -2: constant specific angular momentum j.
            do zone_index = old_num_zones+1,s%nz
               s%j_rot(zone_index) = s%j_rot(old_num_zones)
               s%i_rot(zone_index) = cc23*exp10((2.0D0*s%logR(zone_index)))
               s%omega(zone_index) = s%j_rot(zone_index)/s%i_rot(zone_index)
            end do
         else if(s%physics%cz_rotation_law_exponent >= 0.0D0)then
! walpcz >= 0: solid-body rotation, omega constant.
            do zone_index = old_num_zones+1,s%nz
               s%omega(zone_index) = s%omega(old_num_zones)
               s%i_rot(zone_index) = cc23*exp10((2.0D0*s%logR(zone_index)))
               s%j_rot(zone_index) = s%omega(old_num_zones)*s%i_rot(zone_index)
            end do
         else
! Otherwise omega ~ r^walpcz, the law enforce_rotation_profile
! (omega_profile) imposes: omega = omega_old (r/r_old)^walpcz.
            omega_ref = s%omega(old_num_zones)/exp10((s%logR(old_num_zones)*s%physics%cz_rotation_law_exponent))
            do zone_index = old_num_zones+1,s%nz
               s%omega(zone_index) = omega_ref*exp10((s%logR(zone_index)*s%physics%cz_rotation_law_exponent))
               s%i_rot(zone_index) = cc23*exp10((2.0D0*s%logR(zone_index)))
               s%j_rot(zone_index) = s%omega(zone_index)*s%i_rot(zone_index)
            end do
         end if
end subroutine guess_new_shell_rotation

! ---------------------------------------------------------------
! The angular momentum and rotational kinetic energy of each shell and
! of the star.
subroutine sum_rotation_totals()
! J_k = j_k dm_k and E_k = omega_k J_k / 2 for each shell, and their
! sums.
       sum_angular_momentum = 0.0D0
       sum_rotational_ke = 0.0D0
       do zone_index = 1,s%nz
          angular_momentum_shell = s%j_rot(zone_index)*s%dm(zone_index)
          s%kinetic_energy_rot(zone_index) = 0.5D0*s%omega(zone_index)*angular_momentum_shell
          sum_angular_momentum = sum_angular_momentum + angular_momentum_shell
          sum_rotational_ke = sum_rotational_ke + s%kinetic_energy_rot(zone_index)
       end do
       s%total_angular_momentum = sum_angular_momentum
       s%total_rotational_ke = sum_rotational_ke
end subroutine sum_rotation_totals

! ---------------------------------------------------------------
! Print the new interior shells.
subroutine print_new_interior_points()
      character(len=message_length) :: heading(2)
      character(len=message_length), allocatable :: rows(:)
      heading = unwritten_line
      write(heading,910)
 910  format(1X,'NEW INTERIOR POINTS FROM CHANGE IN ENVELOPE MASS'/ &
     ' J,LOG RHO, LOG L, LOG P, LOG R, LOG M, LOG T, CONV T/F')
      call core_stdout(heading)
      allocate(rows(max(s%nz-old_num_zones+1, 1)))
      rows = unwritten_line
      write(rows,911)(zone_index,s%logRho(zone_index),s%luminosity_lsun(zone_index), &
           s%logP(zone_index),s%logR(zone_index), &
           s%log_mass(zone_index)-s%log_total_mass, &
           s%logT(zone_index),s%convective_flag(zone_index), &
           zone_index = old_num_zones,s%nz)
 911  format(I5,1P6E16.8,L2)
      call core_stdout(rows)
end subroutine print_new_interior_points
end subroutine rebuild_envelope

end module envelope_rebuild

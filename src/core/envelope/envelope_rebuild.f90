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
      use star_info_lib, only: star_info
      implicit none
      private
      public :: rebuild_envelope

contains

! ---------------------------------------------------------------
! rebuild_envelope (formerly getnewenv)
! ---------------------------------------------------------------
! Moves the fitting point outward, to target_envelope_mass (a
! log(m/M), clamped to senv_thin_envelope), after the star has accreted
! (change_mass's new_atmosphere_fit_needed). It integrates the atmosphere and envelope from the surface
! (atm_get at the current L and Teff) down to the current fitting
! point, shifts the envelope onto the last interior shell, and appends
! the envelope points inside the new fitting mass as interior shells
! (append_envelope_points). Then come the new shells' masses and, with
! rotation, a first guess at their rotation, which omega_from_j makes
! consistent, and the star's angular momentum and rotational energy.
! Updates star%senv and num_zones; called from evolve_step. ierr: the
! integration failed, the zones would overflow json, or the mass
! interpolation was degenerate.
subroutine rebuild_envelope(s, target_envelope_mass, composition, log_density, &
     log_luminosity, log_pressure, log_radius, log_mass, enclosed_mass, &
     shell_mass, &
     log_total_mass, log_temperature, convective_flag, eta_squared, &
     moment_of_inertia, specific_angular_momentum, qiw, mean_radius, &
     rotational_kinetic_energy, log_luminosity_lsun, total_angular_momentum, &
     total_rotational_ke, log_teff, num_zones, new_points_added_flag, ierr)
      use omega_profile, only: omega_from_j
      use eos_lib, only: eos_set_debye_huckel_composition
      use envelope_integration, only: atm_get, envint_step_config, fixed_envint_step, &
           surface_log10_radius_cm, ambient_log10_teff
      use envelope_integration_kernel, only: senv_thin_envelope
      use envelope_refit, only: append_envelope_points
      use star_info_lib, only: i_h1, i_metals, json, n_species_basic, &
           n_species_extended
      use controls_buffer, only: ichi_dp_env_max
      use physical_constants
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      integer :: jerr_atm
      double precision, intent(inout) :: target_envelope_mass
      double precision, intent(inout) :: composition(n_species_extended,json), &
           log_density(json), log_luminosity(json), log_pressure(json), &
           log_radius(json), log_mass(json), enclosed_mass(json)
      double precision, intent(inout) :: shell_mass(json)
      double precision, intent(in) :: log_total_mass
      double precision, intent(inout) :: log_temperature(json)
      logical, intent(inout) :: convective_flag(json)
      double precision, intent(inout) :: eta_squared(json), &
           moment_of_inertia(json), specific_angular_momentum(json), &
           qiw(json), mean_radius(json)
      double precision, intent(inout) :: rotational_kinetic_energy(json)
      double precision, intent(in) :: log_luminosity_lsun
      double precision, intent(inout) :: total_angular_momentum, &
           total_rotational_ke
      double precision, intent(in) :: log_teff
      integer, intent(inout) :: num_zones
      logical, intent(out) :: new_points_added_flag

      integer :: species_end_index
      integer :: old_num_zones
      double precision :: envelope_mass_before
      type(envint_step_config) :: env_steps
      logical :: surface_bc_flag, print_flag
      integer :: ksaha
      double precision :: luminosity_linear
      double precision :: log_radius_surface, log_gravity_surface
      double precision :: hydrogen_fraction, metal_fraction
      double precision :: fp_surface, ft_surface
      double precision :: pressure_limit
      double precision :: spot_adjusted_log_teff
      double precision :: pressure_offset, density_offset, temperature_offset, &
           radius_offset
      double precision :: mass_test, mass_test2
      integer :: num_new_env_points
      integer :: zone_index
      double precision :: omega_ref
      double precision :: sum_angular_momentum, sum_rotational_ke
      double precision :: angular_momentum_shell
      double precision :: mass_at_base
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
      if(num_zones+s%env_points%num_env_points.ge.json) then
         write(*,*) 'rebuild_envelope: grid overflow (json too small)'
         ierr = 1
         return
      end if
      call shift_envelope_points()
! Append the envelope points inside the new fitting mass as shells,
! interpolating the last one onto it; a point exactly at the fitting
! mass is interpolated (accept_point_at_fit_mass = .false.).
      call append_envelope_points(s, .false., species_end_index, composition, &
           log_density, log_luminosity, log_pressure, log_radius, log_mass, &
           log_temperature, convective_flag, num_zones, ierr)
! A degenerate mass interval returns an error (legacy stop 9998).
      if (ierr /= 0) then
         write(*,*) 'rebuild_envelope: degenerate mass interpolation'
         return
      end if
      call compute_new_shell_masses()
! With rotation: guess omega and j of the new shells, find omega, the
! shape and the moments of inertia of every zone from j (omega_from_j),
! and sum the angular momentum and rotational energy.
      if(s%job%rotation_active)then
         call guess_new_shell_rotation()
         call omega_from_j(log_density,specific_angular_momentum,log_radius,log_mass, &
              shell_mass,convective_flag,num_zones,eta_squared, &
              moment_of_inertia,s%omega,qiw,mean_radius)
         call sum_rotation_totals()
      endif
      call print_new_interior_points()
      new_points_added_flag = .true.
       write(*,597)envelope_mass_before,s%senv
 597  format(5X,'***** NEW ENVELOPE MASS CALCULATED *****'/8X, &
     'OLD SENV ',1PE22.13,'  NEW SENV',E22.13)
      return

contains

! ---------------------------------------------------------------
! The species count, the target envelope mass (clamped), the current
! number of zones and envelope mass, and the integration settings.
subroutine prepare_envelope_rebuild()
      if(s%job%use_extended_composition)then
         species_end_index = n_species_extended
      else
         species_end_index = n_species_basic
      endif
! Keep senv at or below senv_thin_envelope: an envelope of at least
! about 2.3e-12 of the star.
      target_envelope_mass = min(target_envelope_mass,senv_thin_envelope)
! The zone count and the fitting mass before the rebuild.
      old_num_zones = num_zones
      envelope_mass_before = s%senv
! The envelope integration takes a fixed step in log P,
! chi_grid_scale(ichi_dp_env_max), handed to atm_get for this call
! only. No boundary vertex is saved, and saha_state starts at 0.
      env_steps = fixed_envint_step(s%ctrl%chi_grid_scale(ichi_dp_env_max))
      surface_bc_flag = .false.
      print_flag = .true.
      ksaha = 0
end subroutine prepare_envelope_rebuild

! ---------------------------------------------------------------
! The surface of the envelope integration: L, R, g and the composition,
! no rotational distortion, the pressure limit and the eos composition.
subroutine set_envelope_surface()
! L in Lsun.
      luminosity_linear = exp(ln10*log_luminosity_lsun)
! R from L = 4 pi R^2 sigma Teff^4.
      log_radius_surface = surface_log10_radius_cm(s, log_luminosity_lsun, log_teff)
! g = G M / R^2.
      log_gravity_surface = cgl + s%stotal - log_radius_surface - log_radius_surface
! X and Z of the outermost shell.
      hydrogen_fraction = composition(i_h1,num_zones)
      metal_fraction = composition(i_metals,num_zones)
! No rotational distortion: f_P = f_T = 1.
      fp_surface = 1.0D0
      ft_surface = 1.0D0
! The first step that would cross the pressure of the outermost shell
! stops at it.
      pressure_limit = log_pressure(num_zones)
! The composition of the outermost shell for the EOS.
      call eos_set_debye_huckel_composition(composition(:,num_zones))
end subroutine set_envelope_surface

! ---------------------------------------------------------------
! Integrate the atmosphere and envelope down to the current fitting
! point, then adopt the target envelope mass.
subroutine integrate_new_envelope()
! Integrate down to the current fitting point from the surface L and
! Teff. A spotted run with a convective surface starts from the Teff of
! the unspotted surface (ambient_log10_teff, legacy ATEFFL).
      if(convective_flag(num_zones).and.s%ctrl%spot_filling_factor.ne.0.0.and. &
           s%ctrl%spot_temp_contrast.ne.1.0)then
         spot_adjusted_log_teff = ambient_log10_teff(s, log_teff)
      else
         spot_adjusted_log_teff = log_teff
      endif
      call atm_get(s, luminosity_linear,fp_surface,ft_surface,log_gravity_surface, &
           log_total_mass,print_flag,surface_bc_flag,pressure_limit, &
           log_radius_surface, &
           spot_adjusted_log_teff,hydrogen_fraction,metal_fraction, &
           ksaha,ierr=jerr_atm,env_steps=env_steps)
! An integration failure goes back to evolve_step in ierr.
      if (jerr_atm /= 0) then
         ierr = jerr_atm
         return
      end if
      s%senv = target_envelope_mass
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
      pressure_offset = log_pressure(num_zones) - s%env_points%env_log10_pressure(1)
      density_offset = log_density(num_zones) - s%env_points%env_log10_density(1)
      temperature_offset = log_temperature(num_zones) - s%env_points%env_log10_temperature(1)
      radius_offset = log_radius(num_zones) - s%env_points%env_log10_radius(1)
      mass_test = log_mass(num_zones)
      num_new_env_points = 0
      do zone_index = 1,s%env_points%num_env_points - 1
         mass_test2 = s%stotal+s%env_points%env_log10_mass(zone_index+1)
         if(mass_test2-mass_test.gt.1.0D-10)then
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
         endif
      end do
      s%env_points%num_env_points = num_new_env_points
end subroutine shift_envelope_points

! ---------------------------------------------------------------
! The (unlogged) masses of the new shells and their mass contents.
subroutine compute_new_shell_masses()
! enclosed_mass (legacy HS1) = 10^log m, and the mass of each zone
! (shell_mass, legacy HS2) = (m(k+1) - m(k-1))/2; the outermost zone
! holds M - (m(nz) + m(nz-1))/2.
      do zone_index = old_num_zones,num_zones
         enclosed_mass(zone_index) = exp10(log_mass(zone_index))
      end do
      do zone_index = old_num_zones,num_zones-1
         shell_mass(zone_index) = 0.5D0*(enclosed_mass(zone_index+1)-enclosed_mass(zone_index-1))
      end do
      mass_at_base = 0.5D0*(enclosed_mass(num_zones)+enclosed_mass(num_zones-1))
      shell_mass(num_zones) = exp10(log_total_mass) - mass_at_base
end subroutine compute_new_shell_masses

! ---------------------------------------------------------------
! The first guess at the rotation of the new shells, from the outermost
! old shell by the convection-zone law walpcz: constant j, solid body,
! or omega ~ r^walpcz. Each shell's specific moment of inertia is
! the thin-shell i = (2/3) r^2, and j = omega i.
subroutine guess_new_shell_rotation()
         if(s%ctrl%walpcz.le.-2.0D0)then
! walpcz <= -2: constant specific angular momentum j.
            do zone_index = old_num_zones+1,num_zones
               specific_angular_momentum(zone_index) = specific_angular_momentum(old_num_zones)
               moment_of_inertia(zone_index) = cc23*exp10((2.0D0*log_radius(zone_index)))
               s%omega(zone_index) = specific_angular_momentum(zone_index)/moment_of_inertia(zone_index)
            end do
         else if(s%ctrl%walpcz.ge.0.0D0)then
! walpcz >= 0: solid-body rotation, omega constant.
            do zone_index = old_num_zones+1,num_zones
               s%omega(zone_index) = s%omega(old_num_zones)
               moment_of_inertia(zone_index) = cc23*exp10((2.0D0*log_radius(zone_index)))
               specific_angular_momentum(zone_index) = s%omega(old_num_zones)*moment_of_inertia(zone_index)
            end do
         else
! Otherwise omega ~ r^walpcz, the law enforce_rotation_profile
! (omega_profile) imposes: omega = omega_old (r/r_old)^walpcz.
            omega_ref = s%omega(old_num_zones)/exp10((log_radius(old_num_zones)*s%ctrl%walpcz))
            do zone_index = old_num_zones+1,num_zones
               s%omega(zone_index) = omega_ref*exp10((log_radius(zone_index)*s%ctrl%walpcz))
               moment_of_inertia(zone_index) = cc23*exp10((2.0D0*log_radius(zone_index)))
               specific_angular_momentum(zone_index) = s%omega(zone_index)*moment_of_inertia(zone_index)
            end do
         endif
end subroutine guess_new_shell_rotation

! ---------------------------------------------------------------
! The angular momentum and rotational kinetic energy of each shell and
! of the star.
subroutine sum_rotation_totals()
! J_k = j_k dm_k and E_k = omega_k J_k / 2 for each shell, and their
! sums.
       sum_angular_momentum = 0.0D0
       sum_rotational_ke = 0.0D0
       do zone_index = 1,num_zones
          angular_momentum_shell = specific_angular_momentum(zone_index)*shell_mass(zone_index)
          rotational_kinetic_energy(zone_index) = 0.5D0*s%omega(zone_index)*angular_momentum_shell
          sum_angular_momentum = sum_angular_momentum + angular_momentum_shell
          sum_rotational_ke = sum_rotational_ke + rotational_kinetic_energy(zone_index)
       end do
       total_angular_momentum = sum_angular_momentum
       total_rotational_ke = sum_rotational_ke
end subroutine sum_rotation_totals

! ---------------------------------------------------------------
! Print the new interior shells.
subroutine print_new_interior_points()
      write(*,910)
 910  format(1X,'NEW INTERIOR POINTS FROM CHANGE IN ENVELOPE MASS'/ &
     ' J,LOG RHO, LOG L, LOG P, LOG R, LOG M, LOG T, CONV T/F')
      write(*,911)(zone_index,log_density(zone_index),log_luminosity(zone_index), &
           log_pressure(zone_index),log_radius(zone_index), &
           log_mass(zone_index)-s%stotal, &
           log_temperature(zone_index),convective_flag(zone_index), &
           zone_index = old_num_zones,num_zones)
 911  format(I5,1P6E16.8,L2)
end subroutine print_new_interior_points
end subroutine rebuild_envelope

end module envelope_rebuild

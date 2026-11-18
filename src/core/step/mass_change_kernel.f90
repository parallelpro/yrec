!----------------------------------------------------------------------
! mass_change_kernel
!----------------------------------------------------------------------
! The arithmetic of the mass-change step as star-blind kernels: the
! longest step a rate allows and accretion's end at the disk-locking
! age, the masses of the convection zone after the step's mass change
! dM, the angular momentum a thermal wind carries off, and the accreted
! composition mixed into the zone. Each takes the arrays it works on as
! arguments and changes only its intent(inout) ones; nothing is read
! from the star or kept between calls. apply_mass_change
! (mass_change.f90) gathers the arrays from the star, calls these in
! order and stores the results; core/test/test_mass_change checks that
! they conserve the mass, the angular momentum and every species, and
! the timestep's limits against their formulas.
!
! Shells 1..nz run from the center to the fitting point; the surface
! convection zone is jenv..nz (its base jenv < nz). m is the mass inside
! each zone center (g), dm the shell masses (g), the last reaching the
! total mass M; j_rot is the specific angular momentum (cm^2/s).
module mass_change_kernel
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: n_species_basic
      use physical_constants, only: seconds_per_year
      use math_lib
      implicit none
      private
      public :: wind_angular_momentum, convection_zone_mass, &
           change_shell_masses, wind_specific_angular_momentum, &
           remove_wind_angular_momentum, mix_accreted_composition, &
           mass_change_timestep_limit, clip_step_to_disk_locking_age, &
           disk_has_dispersed

! A thermal wind's angular momentum over one step (g cm^2/s): lost, dJ =
! j_s dM, negative; zone_before, the convection zone's J_cz = sum of
! j_rot dm before the masses change; zone_after = J_cz + dJ, what the
! zone holds after.
      type :: wind_angular_momentum
         real(dp) :: lost = 0.0d0, zone_before = 0.0d0, zone_after = 0.0d0
      end type wind_angular_momentum

contains

!----------------------------------------------------------------------
! convection_zone_mass
!----------------------------------------------------------------------
! The mass of the surface convection zone, M_cz = sum of its shell
! masses dm (g), summed from its base up: zone_dm = dm(jenv:nz). The
! last shell reaches the total mass, so M_cz includes the mass outside
! the fitting point; with jenv = 1 it is the star's mass. The step's one
! M_cz: its timestep limit, the accretion's entropy shift and energy
! means, the compression and the mixing all take it.
pure function convection_zone_mass(zone_dm) result(cz_mass)
      real(dp), intent(in) :: zone_dm(:)
      real(dp) :: cz_mass
      integer :: zone_idx

      cz_mass = 0.0d0
      do zone_idx = 1, size(zone_dm)
         cz_mass = cz_mass + zone_dm(zone_idx)
      end do
end function convection_zone_mass

!----------------------------------------------------------------------
! mass_change_timestep_limit
!----------------------------------------------------------------------
! The longest step (s) the mass change allows at the rate rate_cgs =
! |Mdot| (g/s, positive), for loss and accretion alike: the mass changes
! by at most the fraction total_fraction of the star (star_mass_msun
! solar masses of solar_mass_cgs g), the fraction cz_fraction of the
! surface convection zone (cz_mass = M_cz, g) and half the mass outside
! the fitting point (outer_mass = M - m(nz), g):
!   dt_max = min(f_tot M, f_cz M_cz, (M - m(nz))/2)/|Mdot|.
pure function mass_change_timestep_limit(rate_cgs, total_fraction, &
     cz_fraction, star_mass_msun, solar_mass_cgs, cz_mass, outer_mass) &
     result(dt_max)
      real(dp), intent(in) :: rate_cgs, total_fraction, cz_fraction, &
           star_mass_msun, solar_mass_cgs, cz_mass, outer_mass
      real(dp) :: dt_max

      dt_max = min(total_fraction*star_mass_msun*solar_mass_cgs/rate_cgs, &
           cz_fraction*cz_mass/rate_cgs)
      dt_max = min(dt_max, 0.5d0*outer_mass/rate_cgs)
end function mass_change_timestep_limit

!----------------------------------------------------------------------
! disk_has_dispersed
!----------------------------------------------------------------------
! Whether a model of the age age_gyr is at or past the disk-locking age
! (both Gyr): there is no disk left to accrete from, so the caller
! switches the mass change off without a step (change_mass). The clip
! below is for a younger model only: it would hand an older one a step
! of negative length (docs/bug_sweep.md M319).
pure logical function disk_has_dispersed(age_gyr, disk_locking_age_gyr)
      real(dp), intent(in) :: age_gyr, disk_locking_age_gyr
      disk_has_dispersed = age_gyr >= disk_locking_age_gyr
end function disk_has_dispersed

!----------------------------------------------------------------------
! clip_step_to_disk_locking_age
!----------------------------------------------------------------------
! Accretion ends at the disk-locking age: a step dt (s) from the age
! age_gyr that would pass disk_locking_age_gyr (both Gyr) is clipped to
! end there, and disk_exhausted is set; a shorter step is left alone.
! For a model younger than the disk-locking age (disk_has_dispersed is
! false), so that the clipped step is positive.
pure subroutine clip_step_to_disk_locking_age(age_gyr, &
     disk_locking_age_gyr, dt, disk_exhausted)
      real(dp), intent(in) :: age_gyr, disk_locking_age_gyr
      real(dp), intent(inout) :: dt
      logical, intent(out) :: disk_exhausted

      disk_exhausted = .false.
      if (age_gyr + 1.0d-9*dt/seconds_per_year > disk_locking_age_gyr) then
         dt = (disk_locking_age_gyr - age_gyr)*1.0d9*seconds_per_year
         disk_exhausted = .true.
      end if
end subroutine clip_step_to_disk_locking_age

!----------------------------------------------------------------------
! change_shell_masses
!----------------------------------------------------------------------
! The masses of shells 1..nz after the step's mass change delta_mass =
! dM (g, negative for loss); total_mass is the new total mass M + dM
! (g). m and log_mass (log10 m) are the zone-center masses, dm the
! shell masses; only the convection zone jenv..nz changes.
!   - With mass loss, the zone-center masses above the base (and their
!     logs) compress about m(jenv):
!       m(k) = m(jenv) + f (m(k) - m(jenv))  for k > jenv,
!     with f = (M_cz + dM)/M_cz and M_cz = convection_zone_mass of the
!     shell masses dm(jenv..nz) on entry; the last shell, reaching the
!     new total, takes the rest of dM. Accreted mass is not spread over
!     the zones; it lands in dm(nz).
!   - The zone's shell masses are half the span between the neighbouring
!     zone centers, the last shell reaching the total mass:
!       dm(k) = (m(k+1) - m(k-1))/2,  dm(nz) = M + dM - (m(nz) + m(nz-1))/2.
!     A zone reaching the center (jenv = 1) starts from it: dm(1) =
!     (m(1) + m(2))/2, i.e. m_0 = -m_1 as in rezoning's
!     set_new_shell_masses.
pure subroutine change_shell_masses(jenv, delta_mass, total_mass, m, &
     log_mass, dm)
      integer, intent(in) :: jenv
      real(dp), intent(in) :: delta_mass, total_mass
      real(dp), intent(inout) :: m(:), log_mass(:), dm(:)
      integer :: nz, zone_idx
      real(dp) :: cz_mass, mass_scale_factor

      nz = size(m)
      if (delta_mass < 0.0d0) then
         cz_mass = convection_zone_mass(dm(jenv:nz))
         mass_scale_factor = (cz_mass + delta_mass)/cz_mass
         do zone_idx = jenv+1, nz
            m(zone_idx) = m(jenv) + mass_scale_factor*(m(zone_idx) - m(jenv))
            log_mass(zone_idx) = log10(m(zone_idx))
         end do
      end if
      do zone_idx = max(jenv,2), nz-1
         dm(zone_idx) = 0.5d0*(m(zone_idx+1) - m(zone_idx-1))
      end do
      if (jenv == 1) dm(1) = 0.5d0*(m(2) + m(1))
      dm(nz) = total_mass - 0.5d0*(m(nz) + m(nz-1))
end subroutine change_shell_masses

!----------------------------------------------------------------------
! wind_specific_angular_momentum
!----------------------------------------------------------------------
! The specific angular momentum j_s (cm^2/s) of the matter a thermal
! wind carries off: that of the surface, radius total_radius = R (cm),
! under the convection zone's rotation law walpcz, from the outer zone
! nz (omega_outer = omega(nz), rad/s; log_radius_outer = log10 r(nz),
! cm; j_rot_outer = j_rot(nz)). The moment of inertia per unit mass of
! a thin shell at the surface is (2/3) R^2.
!   walpcz >= 0, a solid-body zone:      j_s = (2/3) R^2 omega(nz);
!   walpcz <= -2, constant specific angular momentum: j_s = j_rot(nz);
!   otherwise a power law, omega ~ r^walpcz (the rotation module's law,
!   as in omega_profile), extrapolated from the outer zone to R:
!     j_s = (2/3) R^2 omega_s,  omega_s = omega(nz) (R/r(nz))^walpcz.
pure function wind_specific_angular_momentum(walpcz, total_radius, &
     omega_outer, log_radius_outer, j_rot_outer) result(j_surface)
      real(dp), intent(in) :: walpcz, total_radius, omega_outer, &
           log_radius_outer, j_rot_outer
      real(dp) :: j_surface
      real(dp) :: surface_moment_of_inertia_per_mass, surface_omega

      surface_moment_of_inertia_per_mass = 2.0d0/3.0d0*total_radius**2
      if (walpcz >= 0.0d0) then
         j_surface = omega_outer*surface_moment_of_inertia_per_mass
      else if (walpcz <= -2.0d0) then
         j_surface = j_rot_outer
      else
         surface_omega = omega_outer*pow(total_radius, walpcz)/ &
              exp10(log_radius_outer*walpcz)
         j_surface = surface_omega*surface_moment_of_inertia_per_mass
      end if
end function wind_specific_angular_momentum

!----------------------------------------------------------------------
! remove_wind_angular_momentum
!----------------------------------------------------------------------
! A thermal wind: the lost matter, delta_mass = dM < 0 (g), carries away
! only its own specific angular momentum j_surface
! (wind_specific_angular_momentum), dJ = j_s dM, and the convection
! zone keeps the rest. dm_before and dm_after are the zone's shell
! masses (jenv..nz) before and after change_shell_masses; j_rot, the
! zone's specific angular momentum, is scaled in place, and am returns
! dJ, J_cz and J_cz + dJ.
! J_cz = sum of j_rot dm_before is measured before the masses change.
! The shrunk shells at unchanged j_rot hold sum j_rot dm_after, already
! short of J_cz by about <j>|dM|; every j_rot is scaled by
! (J_cz + dJ)/sum j_rot dm_after, so the zone loses dJ once.
pure subroutine remove_wind_angular_momentum(j_surface, delta_mass, &
     dm_before, dm_after, j_rot, am)
      real(dp), intent(in) :: j_surface, delta_mass
      real(dp), intent(in) :: dm_before(:), dm_after(:)
      real(dp), intent(inout) :: j_rot(:)
      type(wind_angular_momentum), intent(out) :: am
      integer :: zone_idx
      real(dp) :: am_ratio

      am%lost = j_surface*delta_mass
      am%zone_before = zone_angular_momentum(j_rot, dm_before)
      am%zone_after = (am%zone_before + am%lost)
      am_ratio = am%zone_after/zone_angular_momentum(j_rot, dm_after)
      do zone_idx = 1, size(j_rot)
         j_rot(zone_idx) = j_rot(zone_idx)*am_ratio
      end do
end subroutine remove_wind_angular_momentum

! ---------------------------------------------------------------
! The angular momentum of the zone, sum of j_rot dm from its base up
! (g cm^2/s).
pure function zone_angular_momentum(j_rot, dm) result(am)
      real(dp), intent(in) :: j_rot(:), dm(:)
      real(dp) :: am
      integer :: zone_idx

      am = j_rot(1)*dm(1)
      do zone_idx = 2, size(j_rot)
         am = am + j_rot(zone_idx)*dm(zone_idx)
      end do
end function zone_angular_momentum

!----------------------------------------------------------------------
! mix_accreted_composition
!----------------------------------------------------------------------
! The accreted matter, delta_mass = dM > 0 (g) of composition
! x_accreted, mixed into the convection zone. zone_dm and xa are the
! zone's shell masses before the change (g) and its composition
! (species, shells jenv..nz). For each of the basic species
! (composition_layout's slots 1..n_species_basic), in every shell of the
! zone,
!   X = (sum of X dm + X_acc dM)/(M_cz + dM),  M_cz = sum of dm
! (convection_zone_mass), so the zone's mass of each species grows by
! X_acc dM exactly, whether or not the zone was uniform.
! The light elements of the extended set are not mixed here: burn_step
! takes the accreted gas's share of them as it burns them
! (surface_accretion_now, burn_li_be).
pure subroutine mix_accreted_composition(zone_dm, delta_mass, &
     x_accreted, xa)
      real(dp), intent(in) :: zone_dm(:), delta_mass
      real(dp), intent(in) :: x_accreted(:)
      real(dp), intent(inout) :: xa(:,:)
      integer :: species_idx, zone_idx
      real(dp) :: cz_mass, species_mass, mixed_abundance

      cz_mass = convection_zone_mass(zone_dm)
      do species_idx = 1, n_species_basic
         species_mass = 0.0d0
         do zone_idx = 1, size(zone_dm)
            species_mass = species_mass + xa(species_idx,zone_idx)*zone_dm(zone_idx)
         end do
         mixed_abundance = (species_mass + &
              x_accreted(species_idx)*delta_mass)/(delta_mass + cz_mass)
         do zone_idx = 1, size(xa,2)
            xa(species_idx,zone_idx) = mixed_abundance
         end do
      end do
end subroutine mix_accreted_composition

! ---------------------------------------------------------------
! Sources
!   The arithmetic of the F77's MDOT (wind/mdot, MHP 10/01), as kernels
!   since 2026-10. Two of its defects, fixed in 2026: for a zone reaching the
!   centre it read m(0), outside the array; and every j_rot was scaled by
!   (J_cz + dJ)/J_cz before the masses shrank, which removed <j>|dM| on top
!   of dJ. The timestep's limit and the disk-locking clip came out of
!   mass_change in 2026-11, with the gate before the clip
!   (docs/bug_sweep.md M319).
end module mass_change_kernel

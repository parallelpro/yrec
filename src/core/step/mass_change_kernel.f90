!----------------------------------------------------------------------
! mass_change_kernel
!----------------------------------------------------------------------
! The arithmetic of the mass-change step as star-blind kernels: the
! masses of the convection zone after the step's mass change dM, the
! angular momentum a thermal wind carries off, and the accreted
! composition mixed into the zone. Each takes the arrays it works on as
! arguments and changes only its intent(inout) ones; nothing is read
! from the star or kept between calls. apply_mass_change
! (mass_change.f90) gathers the arrays from the star, calls these in
! order and stores the results; core/test/test_mass_change checks that
! they conserve the mass, the angular momentum and every species.
!
! Shells 1..nz run from the center to the fitting point; the surface
! convection zone is jenv..nz (its base jenv < nz). m is the mass inside
! each zone center (g), dm the shell masses (g), the last reaching the
! total mass M; j_rot is the specific angular momentum (cm^2/s).
module mass_change_kernel
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: n_species_basic
      use math_lib
      implicit none
      private
      public :: wind_angular_momentum, convection_zone_mass, &
           change_shell_masses, wind_specific_angular_momentum, &
           remove_wind_angular_momentum, mix_accreted_composition

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
! The mass of the surface convection zone below the fitting point,
! M_cz = m(nz) - m(jenv) (g), or m(nz) when the zone reaches the center
! (jenv = 1). m: the zone-center masses of shells 1..nz.
pure function convection_zone_mass(jenv, m) result(cz_mass)
      integer, intent(in) :: jenv
      real(dp), intent(in) :: m(:)
      real(dp) :: cz_mass
      integer :: nz

      nz = size(m)
      if (jenv > 1) then
         cz_mass = m(nz) - m(jenv)
      else
         cz_mass = m(nz)
      end if
end function convection_zone_mass

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
!     with f = (M_cz + dM)/M_cz and M_cz = convection_zone_mass. Accreted
!     mass is not spread over the zones; it lands in dm(nz).
!   - The zone's shell masses are half the span between the neighbouring
!     zone centers, the last shell reaching the total mass:
!       dm(k) = (m(k+1) - m(k-1))/2,  dm(nz) = M + dM - (m(nz) + m(nz-1))/2.
!     A zone reaching the center (jenv = 1) starts from it: dm(1) =
!     (m(1) + m(2))/2, i.e. m_0 = -m_1 as in rezoning's
!     set_new_shell_masses. (It used to read m(0), outside the array.)
pure subroutine change_shell_masses(jenv, delta_mass, total_mass, m, &
     log_mass, dm)
      integer, intent(in) :: jenv
      real(dp), intent(in) :: delta_mass, total_mass
      real(dp), intent(inout) :: m(:), log_mass(:), dm(:)
      integer :: nz, zone_idx
      real(dp) :: cz_mass, mass_scale_factor

      nz = size(m)
      if (delta_mass < 0.0d0) then
         cz_mass = convection_zone_mass(jenv, m)
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
! (J_cz + dJ)/sum j_rot dm_after, so the zone loses dJ once. (It used
! to be scaled by (J_cz + dJ)/J_cz before the masses shrank, which
! removed <j>|dM| on top of dJ.)
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
! x_accreted, mixed into the convection zone of mass cz_mass = M_cz (g,
! the sum of its shell masses before the change). xa is the zone's
! composition (species, shells jenv..nz); the zone is taken to be mixed
! already, so its base column stands for it. For each of the basic
! species (composition_layout's slots 1..n_species_basic), in every
! shell of the zone,
!   X = (X(jenv) M_cz + X_acc dM)/(M_cz + dM).
! The light elements of the extended set are not mixed here: burn_step
! takes the accreted gas's share of them as it burns them
! (surface_accretion_now, burn_li_be).
pure subroutine mix_accreted_composition(cz_mass, delta_mass, &
     x_accreted, xa)
      real(dp), intent(in) :: cz_mass, delta_mass
      real(dp), intent(in) :: x_accreted(:)
      real(dp), intent(inout) :: xa(:,:)
      integer :: species_idx, zone_idx
      real(dp) :: mixed_abundance

      do species_idx = 1, n_species_basic
         mixed_abundance = (xa(species_idx,1)*cz_mass + &
              x_accreted(species_idx)*delta_mass)/(delta_mass + cz_mass)
         do zone_idx = 1, size(xa,2)
            xa(species_idx,zone_idx) = mixed_abundance
         end do
      end do
end subroutine mix_accreted_composition

end module mass_change_kernel

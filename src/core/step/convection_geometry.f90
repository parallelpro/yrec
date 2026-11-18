!----------------------------------------------------------------------
! convection_geometry
!----------------------------------------------------------------------
! Where the convection zones end, in one place. The zone flags
! (convection_lib's find_convection_zones) say which shells are
! convective; the routines here turn that into positions:
!   cz_base_fraction           where del_rad - del_ad changes sign between
!                              the last convective zone and the one below
!   first_zone_above           the lowest zone lying above a radius
!   locate_surface_mixing_base the base of the surface mixed region
!                              (convection zone plus envelope overshoot)
!                              at the start and end of the step, for the
!                              Li/Be burn
!   mixed_base_mass, mixed_base_fraction
!                              that base as a mass coordinate, and back
!                              to its zone's fraction
! Users: observables (the surface convection zone's mass, radius, T, rho,
! P and kappa), turnover_timescale (the radius of its base), burn_step
! (burn_li_be). Each caller keeps its own choice of gradients and its own
! switch for overshoot. The overshoot depth alpha_ov H_P and H_P itself
! are convection_lib's (envelope_overshoot_depth, pressure_scale_height).
module convection_geometry
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star_info
      implicit none
      private
      public :: cz_base_fraction, first_zone_above, locate_surface_mixing_base, &
           mixed_base_mass, mixed_base_fraction

contains

!----------------------------------------------------------------------
! cz_base_fraction: the fraction f of the way from the zone below the last
! convective zone up to that zone at which the gap del_rad - del_ad (or
! del_ad - del_rad: the sign cancels) vanishes, by linear interpolation
! between its value in the zone below (gap_below) and in the last
! convective zone (gap):
!   f = gap_below / (gap_below - gap).
real(dp) function cz_base_fraction(gap_below, gap)
      real(dp), intent(in) :: gap_below, gap
      cz_base_fraction = gap_below/(gap_below - gap)
end function cz_base_fraction

!----------------------------------------------------------------------
! first_zone_above: the lowest zone lying above radius_limit (cm).
! Searching down from start_zone, the first zone k with r(k) <
! radius_limit gives k + 1; the result is 1 if every zone down to the
! centre lies above the limit. log_radius is log10 r.
integer function first_zone_above(log_radius, start_zone, radius_limit)
      use physical_constants, only: ln10
      use math_lib
      real(dp), intent(in) :: log_radius(:), radius_limit
      integer, intent(in) :: start_zone
      integer :: k
      do k = start_zone, 1, -1
         if (exp(ln10*log_radius(k)) < radius_limit) then
            first_zone_above = k + 1
            return
         end if
      end do
      first_zone_above = 1
end function first_zone_above

!----------------------------------------------------------------------
! locate_surface_mixing_base (formerly locate_cz_base in burn_step)
! The base of the surface mixed region (the convection zone, plus the
! envelope overshoot when it is on) at the end of the step
! (cz_base_zone) and at the start (cz_base_zone_old), for burn_li_be.
! radius is log10 r of the end-of-step model; env_cz_zone and
! env_cz_zone_old are the convection zone's base zone at the end and at
! the start of the step. cz_base_frac places the true base within its
! zone k = cz_base_zone: its offset from the midpoint of r(k-1) and r(k),
! in units of r(k) - r(k-1), clipped to [-0.5, 0.5]; negative means the
! true base lies higher, positive lower. It is intent(inout) because it
! is assigned only when the base is an interior zone
! (1 < env_cz_zone < num_zones); the callers do not read it otherwise.
! cz_base_zone is kept in star%envelope_cz_base_zone_last_call, the starting
! zone of the next call within the step.
!   star reads: envelope_overshoot_depth_end,
!     envelope_overshoot_depth_start, logR_start,
!     physics%envelope_overshoot_active
!   star writes: cz_base_radius_last_call,
!     envelope_cz_base_zone_last_call
subroutine locate_surface_mixing_base(s, radius, grada, gradr, env_cz_zone, &
     env_cz_zone_old, num_zones, cz_base_zone, cz_base_zone_old, cz_base_frac)
      use physical_constants, only: ln10
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s

      real(dp), intent(in) :: radius(:)
! grada, gradr: the adiabatic and radiative gradients the base is located
! with (a rotation sub-step's interpolated model's, or the model's own)
      real(dp), intent(in) :: grada(:), gradr(:)
      integer, intent(in) :: env_cz_zone, env_cz_zone_old, num_zones
      integer, intent(out) :: cz_base_zone, cz_base_zone_old
      real(dp), intent(inout) :: cz_base_frac

      real(dp) :: del_diff, del_diff_below
      real(dp) :: search_radius, shell_radius, delta_radius, &
           cz_base_radius

! del_ad - del_rad at the last convective zone and at the one below it.
      if(env_cz_zone > 1.and.env_cz_zone < num_zones)then
         del_diff = grada(env_cz_zone) - gradr(env_cz_zone)
         del_diff_below = grada(env_cz_zone-1) - gradr(env_cz_zone-1)
! The true base (del_ad = del_rad) by linear interpolation, a fraction f
! (cz_base_fraction) of the way up from zone env_cz_zone-1; as an offset
! from the zone midpoint, cz_base_frac = 0.5 - f, clipped to
! [-0.5, 0.5]. Negative means the true base lies higher, positive lower.
         cz_base_frac = max(-0.5d0,0.5d0-cz_base_fraction(del_diff_below, &
              del_diff))
         cz_base_frac = min(0.5d0,cz_base_frac)
         if(.not.s%physics%envelope_overshoot_active)then
            cz_base_zone = env_cz_zone
            cz_base_zone_old = env_cz_zone_old
         else
! The starting depth. cz_base_radius_last_call = 0 (set by rezone_or_snapshot
! at every step) marks the first call of the step: the starting base
! radius is then the midpoint of zones env_cz_zone_old-1 and
! env_cz_zone_old in the start-of-step model, and the mixed region
! reaches the overshoot depth alpha_ov H_P (envelope_overshoot_depth_start)
! below it. Later calls in the step (the rotation sub-steps) start from
! the previous call's end-of-step zone.
            if(s%cz_base_radius_last_call == 0.0d0)then
               if (env_cz_zone_old > 1) then
                  s%cz_base_radius_last_call = 0.5d0*(exp(ln10*s%logR_start(env_cz_zone_old)) &
                           +exp(ln10*s%logR_start(env_cz_zone_old-1)))
                  search_radius = s%cz_base_radius_last_call - s%envelope_overshoot_depth_start
                  cz_base_zone_old = first_zone_above(s%logR_start, env_cz_zone_old-1, &
                       search_radius)
               else
! Fully convective at the start of the step: there is no zone below
! the base, which is zone 1 (the code used to read logR_start(0) here,
! out of bounds, and arrive at the same zone).
                  cz_base_zone_old = 1
               end if
            else
               cz_base_zone_old = s%envelope_cz_base_zone_last_call
            end if
! The ending depth, from the interpolated true base,
!   r_base = (r(k) + r(k-1))/2 - cz_base_frac (r(k) - r(k-1)),
! k = env_cz_zone; the mixed region reaches alpha_ov H_P
! (envelope_overshoot_depth_end) below it.
            delta_radius = exp(ln10*radius(env_cz_zone))-exp(ln10*radius(env_cz_zone-1))
            cz_base_radius = 0.5d0*(exp(ln10*radius(env_cz_zone))+exp(ln10*radius(env_cz_zone-1))) &
                    -cz_base_frac*delta_radius
            s%cz_base_radius_last_call = cz_base_radius
            search_radius = cz_base_radius - s%envelope_overshoot_depth_end
            cz_base_zone = first_zone_above(radius, env_cz_zone-1, search_radius)
! The overshoot edge's offset within its zone k = cz_base_zone, in the
! same convention: cz_base_frac = 0.5 - (r_edge - r(k-1))/(r(k) - r(k-1)),
! clipped to [-0.5, 0.5], both radii the current model's (the legacy
! code took r(k) from the start-of-step model until 2026-11, bug_sweep
! O10).
            if (cz_base_zone > 1) then
               shell_radius = exp(ln10*radius(cz_base_zone-1))
               delta_radius = exp(ln10*radius(cz_base_zone))-shell_radius
               cz_base_frac = 0.5d0-((search_radius-shell_radius)/delta_radius)
               cz_base_frac = max(-0.5d0,cz_base_frac)
               cz_base_frac = min(0.5d0,cz_base_frac)
            end if
         end if
      else
         cz_base_zone = env_cz_zone
         cz_base_zone_old = env_cz_zone_old
      end if
      s%envelope_cz_base_zone_last_call = cz_base_zone
      return
end subroutine locate_surface_mixing_base

!----------------------------------------------------------------------
! mixed_base_mass: the mass coordinate (g) of the base of a mixed region
! whose base zone is k = base_zone and whose edge sits base_frac of the
! way, in units of m(k) - m(k-1), below the zone's lower boundary (the
! convention of locate_surface_mixing_base and average_li_be_over_cz:
! the region holds shells k..num_zones and base_frac (m(k) - m(k-1))
! more):
!   m_base = sum of dm(1..k-1) - base_frac (m(k) - m(k-1)),
! with no fraction when k is 1 or num_zones. mass_coordinate is m, the
! mass inside each zone centre, and shell_mass dm.
pure function mixed_base_mass(mass_coordinate, shell_mass, base_zone, &
     base_frac, num_zones) result(base_mass)
      real(dp), intent(in) :: mass_coordinate(:), shell_mass(:), base_frac
      integer, intent(in) :: base_zone, num_zones
      real(dp) :: base_mass
      integer :: k

      base_mass = 0.0d0
      do k = 1, base_zone-1
         base_mass = base_mass + shell_mass(k)
      end do
      if (base_zone > 1 .and. base_zone < num_zones) base_mass = base_mass - &
           base_frac*(mass_coordinate(base_zone) - mass_coordinate(base_zone-1))
end function mixed_base_mass

!----------------------------------------------------------------------
! mixed_base_fraction: the edge fraction (mixed_base_mass's convention)
! that puts the base of a region whose base zone is k = base_zone at the
! mass coordinate base_mass (g):
!   f = (sum of dm(1..k-1) - base_mass)/(m(k) - m(k-1)),
! clipped to [-0.5, 0.5] as locate_surface_mixing_base clips its own;
! 0 when k is 1 or num_zones, or when base_mass is not positive (no base
! recorded).
pure function mixed_base_fraction(mass_coordinate, shell_mass, base_zone, &
     base_mass, num_zones) result(base_frac)
      real(dp), intent(in) :: mass_coordinate(:), shell_mass(:), base_mass
      integer, intent(in) :: base_zone, num_zones
      real(dp) :: base_frac

      base_frac = 0.0d0
      if (base_mass <= 0.0d0 .or. base_zone <= 1 .or. base_zone >= num_zones) return
      base_frac = (mixed_base_mass(mass_coordinate, shell_mass, base_zone, 0.0d0, &
           num_zones) - base_mass)/(mass_coordinate(base_zone) - mass_coordinate(base_zone-1))
      base_frac = min(0.5d0, max(-0.5d0, base_frac))
end function mixed_base_fraction

end module convection_geometry

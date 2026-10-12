!----------------------------------------------------------------------
! convection_geometry
!----------------------------------------------------------------------
! Where the convection zones end, in one place (2026-10). The zone flags
! (convection_lib's find_convection_zones) say which shells are convective; the
! routines here turn that into positions:
!   cz_base_fraction           where del_rad - del_ad changes sign between
!                              the last convective zone and the one below
!   pressure_scale_height      H_P = P r^2 / (rho G m), cm
!   first_zone_above           the lowest zone lying above a radius
!   locate_surface_mixing_base the base of the surface mixed region
!                              (convection zone plus envelope overshoot)
!                              at the start and end of the step, for the
!                              Li/Be burn
! Users: observables (the surface CZ's mass, radius, T, rho, P, kappa),
! turnover_timescale (the CZ-base radius), evolve_step and
! evolve_angular_momentum (the envelope overshoot depth, alpha * H_P),
! burn_step (burn_li_be, burn_li_be_midpoint). Each caller keeps its own
! choice of gradients and its own switch for overshoot.
module convection_geometry
      implicit none
      private
      public :: cz_base_fraction, pressure_scale_height, first_zone_above, &
           locate_surface_mixing_base

contains

!----------------------------------------------------------------------
! cz_base_fraction: gap = del_rad - del_ad (or del_ad - del_rad: the sign
! cancels) at the last convective zone (gap) and at the zone below it
! (gap_below). Returns the fraction of the way from the zone below to the
! last convective zone at which gap = 0, by linear interpolation.
double precision function cz_base_fraction(gap_below, gap)
      double precision, intent(in) :: gap_below, gap
      cz_base_fraction = gap_below/(gap_below - gap)
end function cz_base_fraction

!----------------------------------------------------------------------
! pressure_scale_height: H_P = P r^2 / (rho G m) in cm, from log10 P,
! r, rho and m.
double precision function pressure_scale_height(log_pressure, log_radius, &
     log_density, log_mass)
      use physical_constants, only: clndp, cgl
      use math_lib
      double precision, intent(in) :: log_pressure, log_radius, &
           log_density, log_mass
      pressure_scale_height = exp(clndp*(log_pressure+2.0D0*log_radius &
           -log_density-cgl-log_mass))
end function pressure_scale_height

!----------------------------------------------------------------------
! first_zone_above: searching down from start_zone, the first zone k with
! r(k) < radius_limit (cm) gives k + 1, the lowest zone lying above the
! limit; 1 if every zone down to the centre lies above it.
integer function first_zone_above(log_radius, start_zone, radius_limit)
      use star_info_lib, only: json
      use physical_constants, only: ln10
      use math_lib
      double precision, intent(in) :: log_radius(json), radius_limit
      integer, intent(in) :: start_zone
      integer :: k
      do k = start_zone, 1, -1
         if (exp(ln10*log_radius(k)).lt.radius_limit) then
            first_zone_above = k + 1
            return
         end if
      end do
      first_zone_above = 1
end function first_zone_above

!----------------------------------------------------------------------
! locate_surface_mixing_base (formerly locate_cz_base in burn_step; shared
! by burn_li_be and burn_li_be_midpoint, 2026 readability R3). Locates the true base of the
! surface convection zone at the end of the step (cz_base_zone, with
! cz_base_frac the offset of the true base from the zone midpoint in
! units of the zone width: negative = higher, positive = lower) and
! the base at the start of the step (cz_base_zone_old), including
! envelope overshoot when active, and records cz_base_zone in
! star%envelope_cz_base_zone_prev.
! cz_base_frac is intent(inout): it is assigned only when the CZ base
! is an interior zone (1 < env_cz_zone < num_zones), exactly as in the
! original inline block, and is never read by the callers otherwise.
!   star reads: grada, gradr, job%envelope_overshoot_active,
!     job%instability_transport_active, job%rotation_active, logR_start,
!     pressure_scale_height_end, pressure_scale_height_start
!   star writes: cz_base_radius_prev, envelope_cz_base_zone_prev
subroutine locate_surface_mixing_base(radius, env_cz_zone, env_cz_zone_old, num_zones, &
     cz_base_zone, cz_base_zone_old, cz_base_frac)
      use rotation_scratch_lib
      use star_info_lib, only: star, json
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: radius(json)
      integer, intent(in) :: env_cz_zone, env_cz_zone_old, num_zones
      integer, intent(out) :: cz_base_zone, cz_base_zone_old
      double precision, intent(inout) :: cz_base_frac

      double precision :: del_diff, del_diff_below
      double precision :: search_radius, shell_radius, delta_radius, &
           cz_base_radius

      if(env_cz_zone.gt.1.and.env_cz_zone.lt.num_zones)then
         if(star%job%rotation_active.and.star%job%instability_transport_active)then
            del_diff = mid_scr%grada(env_cz_zone) - &
                 mid_scr%gradr(env_cz_zone)
            del_diff_below = mid_scr%grada(env_cz_zone-1) - &
                 mid_scr%gradr(env_cz_zone-1)
         else
! EVALUATE DEL(AD) - DEL(RAD) AT THE LAST CONVECTIVE POINT AND THE ONE
! BELOW IT.
            del_diff = star%grada(env_cz_zone)-star%gradr(env_cz_zone)
            del_diff_below = star%grada(env_cz_zone-1)-star%gradr(env_cz_zone-1)
         endif
! USE LINEAR INTERPOLATION TO FIND THE DISTANCE OF THE TRUE LOCATION
! OF THE BASE FROM THE ZONE MIDPOINT. IF FX IS NEGATIVE,THEN THE TRUE
! BASE IS HIGHER; IF IT IS POSITIVE, THE TRUE BASE IS LOWER.
         cz_base_frac = max(-0.5d0,0.5d0-cz_base_fraction(del_diff_below, &
              del_diff))
         cz_base_frac = min(0.5d0,cz_base_frac)
         if(.not.star%job%envelope_overshoot_active)then
            cz_base_zone = env_cz_zone
            cz_base_zone_old = env_cz_zone_old
         else
! STARTING CZ DEPTH (cz_base_radius_prev = 0 means no previous step
! has stored one: locate it from the start-of-step model).
            if(star%cz_base_radius_prev.eq.0.0d0)then
               star%cz_base_radius_prev = 0.5d0*(exp(ln10*star%logR_start(env_cz_zone_old)) &
                        +exp(ln10*star%logR_start(env_cz_zone_old-1)))
               search_radius = star%cz_base_radius_prev - star%pressure_scale_height_start
               cz_base_zone_old = first_zone_above(star%logR_start, env_cz_zone_old-1, &
                    search_radius)
            else
               cz_base_zone_old = star%envelope_cz_base_zone_prev
            endif
! ENDING CZ DEPTH : DETERMINE OVERSHOOT FROM TRUE CZ BASE.
            delta_radius = exp(ln10*radius(env_cz_zone))-exp(ln10*radius(env_cz_zone-1))
            cz_base_radius = 0.5d0*(exp(ln10*radius(env_cz_zone))+exp(ln10*radius(env_cz_zone-1))) &
                    -cz_base_frac*delta_radius
            star%cz_base_radius_prev = cz_base_radius
            search_radius = cz_base_radius - star%pressure_scale_height_end
            cz_base_zone = first_zone_above(radius, env_cz_zone-1, search_radius)
! the true overshoot edge's offset within its zone. (The zone width uses
! the start-of-step radius above the edge and the current one below it,
! as it always has: see ROADMAP 13.)
            if (cz_base_zone.gt.1) then
               shell_radius = exp(ln10*radius(cz_base_zone-1))
               delta_radius = exp(ln10*star%logR_start(cz_base_zone))-shell_radius
               cz_base_frac = 0.5d0-((search_radius-shell_radius)/delta_radius)
               cz_base_frac = max(-0.5d0,cz_base_frac)
               cz_base_frac = min(0.5d0,cz_base_frac)
            endif
         endif
      else
         cz_base_zone = env_cz_zone
         cz_base_zone_old = env_cz_zone_old
      endif
      star%envelope_cz_base_zone_prev = cz_base_zone
      return
end subroutine locate_surface_mixing_base

end module convection_geometry

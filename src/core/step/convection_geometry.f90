!----------------------------------------------------------------------
! convection_geometry
!----------------------------------------------------------------------
! Where the convection zones end, in one place. The zone flags
! (convection_lib's find_convection_zones) say which shells are
! convective; the routines here turn that into positions:
!   cz_base_fraction           where del_rad - del_ad changes sign between
!                              the last convective zone and the one below
!   pressure_scale_height      H_P = P r^2 / (rho G m), cm
!   first_zone_above           the lowest zone lying above a radius
!   locate_surface_mixing_base the base of the surface mixed region
!                              (convection zone plus envelope overshoot)
!                              at the start and end of the step, for the
!                              Li/Be burn
! Users: observables (the surface convection zone's mass, radius, T, rho,
! P and kappa), turnover_timescale (the radius of its base), evolve_step
! and evolve_angular_momentum (the envelope overshoot depth, alpha_ov
! H_P), burn_step (burn_li_be). Each caller keeps its own choice of
! gradients and its own switch for overshoot.
module convection_geometry
      use star_info_lib, only: star_info
      implicit none
      private
      public :: cz_base_fraction, pressure_scale_height, first_zone_above, &
           locate_surface_mixing_base

contains

!----------------------------------------------------------------------
! cz_base_fraction: the fraction f of the way from the zone below the last
! convective zone up to that zone at which the gap del_rad - del_ad (or
! del_ad - del_rad: the sign cancels) vanishes, by linear interpolation
! between its value in the zone below (gap_below) and in the last
! convective zone (gap):
!   f = gap_below / (gap_below - gap).
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
! first_zone_above: the lowest zone lying above radius_limit (cm).
! Searching down from start_zone, the first zone k with r(k) <
! radius_limit gives k + 1; the result is 1 if every zone down to the
! centre lies above the limit. log_radius is log10 r.
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
! cz_base_zone is kept in star%envelope_cz_base_zone_prev, the starting
! zone of the next call within the step.
!   star reads: envelope_overshoot_depth_end,
!     envelope_overshoot_depth_start, job%envelope_overshoot_active,
!     logR_start
!   star writes: cz_base_radius_prev, envelope_cz_base_zone_prev
subroutine locate_surface_mixing_base(s, radius, grada, gradr, env_cz_zone, &
     env_cz_zone_old, num_zones, cz_base_zone, cz_base_zone_old, cz_base_frac)
      use star_info_lib, only: json
      use physical_constants
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s

      double precision, intent(in) :: radius(json)
! grada, gradr: the adiabatic and radiative gradients the base is located
! with (a rotation sub-step's interpolated model's, or the model's own)
      double precision, intent(in) :: grada(json), gradr(json)
      integer, intent(in) :: env_cz_zone, env_cz_zone_old, num_zones
      integer, intent(out) :: cz_base_zone, cz_base_zone_old
      double precision, intent(inout) :: cz_base_frac

      double precision :: del_diff, del_diff_below
      double precision :: search_radius, shell_radius, delta_radius, &
           cz_base_radius

! del_ad - del_rad at the last convective zone and at the one below it.
      if(env_cz_zone.gt.1.and.env_cz_zone.lt.num_zones)then
         del_diff = grada(env_cz_zone) - gradr(env_cz_zone)
         del_diff_below = grada(env_cz_zone-1) - gradr(env_cz_zone-1)
! The true base (del_ad = del_rad) by linear interpolation, a fraction f
! (cz_base_fraction) of the way up from zone env_cz_zone-1; as an offset
! from the zone midpoint, cz_base_frac = 0.5 - f, clipped to
! [-0.5, 0.5]. Negative means the true base lies higher, positive lower.
         cz_base_frac = max(-0.5d0,0.5d0-cz_base_fraction(del_diff_below, &
              del_diff))
         cz_base_frac = min(0.5d0,cz_base_frac)
         if(.not.s%job%envelope_overshoot_active)then
            cz_base_zone = env_cz_zone
            cz_base_zone_old = env_cz_zone_old
         else
! The starting depth. cz_base_radius_prev = 0 (set by rezone_or_snapshot
! at every step) marks the first call of the step: the starting base
! radius is then the midpoint of zones env_cz_zone_old-1 and
! env_cz_zone_old in the start-of-step model, and the mixed region
! reaches the overshoot depth alpha_ov H_P (pressure_scale_height_start)
! below it. Later calls in the step (the rotation sub-steps) start from
! the previous call's end-of-step zone.
            if(s%cz_base_radius_prev.eq.0.0d0)then
               if (env_cz_zone_old.gt.1) then
                  s%cz_base_radius_prev = 0.5d0*(exp(ln10*s%logR_start(env_cz_zone_old)) &
                           +exp(ln10*s%logR_start(env_cz_zone_old-1)))
                  search_radius = s%cz_base_radius_prev - s%envelope_overshoot_depth_start
                  cz_base_zone_old = first_zone_above(s%logR_start, env_cz_zone_old-1, &
                       search_radius)
               else
! Fully convective at the start of the step: there is no zone below
! the base, which is zone 1 (the code used to read logR_start(0) here,
! out of bounds, and arrive at the same zone).
                  cz_base_zone_old = 1
               end if
            else
               cz_base_zone_old = s%envelope_cz_base_zone_prev
            endif
! The ending depth, from the interpolated true base,
!   r_base = (r(k) + r(k-1))/2 - cz_base_frac (r(k) - r(k-1)),
! k = env_cz_zone; the mixed region reaches alpha_ov H_P
! (pressure_scale_height_end) below it.
            delta_radius = exp(ln10*radius(env_cz_zone))-exp(ln10*radius(env_cz_zone-1))
            cz_base_radius = 0.5d0*(exp(ln10*radius(env_cz_zone))+exp(ln10*radius(env_cz_zone-1))) &
                    -cz_base_frac*delta_radius
            s%cz_base_radius_prev = cz_base_radius
            search_radius = cz_base_radius - s%envelope_overshoot_depth_end
            cz_base_zone = first_zone_above(radius, env_cz_zone-1, search_radius)
! The overshoot edge's offset within its zone k = cz_base_zone, in the
! same convention: cz_base_frac = 0.5 - (r_edge - r(k-1))/(r(k) - r(k-1)),
! clipped to [-0.5, 0.5]. The zone width takes r(k) from the start-of-step
! model (logR_start) and r(k-1) from the current one, as the legacy code
! did (docs/ROADMAP.md section 13).
            if (cz_base_zone.gt.1) then
               shell_radius = exp(ln10*radius(cz_base_zone-1))
               delta_radius = exp(ln10*s%logR_start(cz_base_zone))-shell_radius
               cz_base_frac = 0.5d0-((search_radius-shell_radius)/delta_radius)
               cz_base_frac = max(-0.5d0,cz_base_frac)
               cz_base_frac = min(0.5d0,cz_base_frac)
            endif
         endif
      else
         cz_base_zone = env_cz_zone
         cz_base_zone_old = env_cz_zone_old
      endif
      s%envelope_cz_base_zone_prev = cz_base_zone
      return
end subroutine locate_surface_mixing_base

end module convection_geometry

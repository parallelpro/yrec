!----------------------------------------------------------------------
! momi
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original momi.f; only variable names, source form, and comment style
! were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Finds the moment of inertia (moment_of_inertia) and dI/d(omega)
! (di_domega) at each zone from zone_start to zone_end. momi assumes
! that SHAPE has already been called. The moment of inertia is
! calculated by the method used in Wai-Yuen Law's thesis (Yale, 1980),
! p.61: a power-series correction (evaluated by shell_inertia_integral.f90) to the
! spherical moment of inertia, in the rotation parameter
! omega**2*r**3/(GM), unless we are inside a convective zone forced to
! solid-body rotation (walpcz.ne.0), in which case the exact spherical
! formula is used directly.
subroutine zone_moments_of_inertia(eta_squared, log_radius, log_mass, shell_mass, zone_start, &
     zone_end, omega, mean_radius, moment_of_inertia, di_domega)
      use star_info_lib, only: star
      use star_info_lib, only: star, json
      use luout_lib
      use phys_const_lib
      implicit none

      double precision, intent(in) :: eta_squared(json), log_radius(json), &
           log_mass(json), shell_mass(json)
      integer, intent(in) :: zone_start, zone_end
      double precision, intent(in) :: omega(json), mean_radius(json)
      double precision, intent(out) :: moment_of_inertia(json), di_domega(json)
      integer :: zone_idx
      double precision :: prev_log_mean_radius, prev_log_true_radius, &
           rotation_param_const, mean_zone_radius, dlnr0_dlnr, &
           spherical_moment_of_inertia, mean_radius_cubed, &
           true_radius_cubed, r0_geom_factor, rotation_param, &
           eta_squared_zone, moment_of_inertia_per_mass, di_domega_per_mass

!  FIND THE MOMENT OF INERTIA (HI) AND DI/D(OMEGA) (QIW).
!  MOMI ASSUMES THAT SHAPE HAS ALREADY BEEN CALLED.
!  MOMENT OF INERTIA IS CALCULATED BY THE METHOD USED IN WAI-YUEN LAW'S
!  THESIS(YALE,1980) P.61.
! prev_log_mean_radius/prev_log_true_radius (originally R0P/RPHIP) are
! computed below at each zone but never subsequently read; preserved
! as dead code from the original.
      if (zone_start.eq.1) then
         prev_log_mean_radius = 0.0d0
         prev_log_true_radius = 0.0d0
      else
         prev_log_mean_radius = dlog(mean_radius(zone_start - 1))
         prev_log_true_radius = ln10*log_radius(zone_start - 1)
      end if
      rotation_param_const = cc13*5.0d0/dexp(ln10*cgl)
      if (star%ctrl%walpcz.ne.0.0d0) then
         do zone_idx = zone_start,zone_end
            moment_of_inertia(zone_idx) = cc23*shell_mass(zone_idx)* &
                 exp(ln10*2.0d0*log_radius(zone_idx))
            di_domega(zone_idx) = 0.0d0
         end do
         if (zone_start.eq.1) then
             mean_zone_radius = 0.5d0*(log_radius(1)+log_radius(2))
             moment_of_inertia(1) = 0.4d0*shell_mass(1)* &
                  exp(ln10*2.0d0*mean_zone_radius)
         end if
         return
      end if
      do zone_idx = zone_start,zone_end
!  QR0R = D LN R0/ D LN R
!
         dlnr0_dlnr = 1.0d0
! spherical_moment_of_inertia (originally H0) is computed here but
! never subsequently read; preserved as dead code from the original.
         spherical_moment_of_inertia = cc23*shell_mass(zone_idx)* &
              dexp(ln10*2.0d0*log_radius(zone_idx))
         mean_radius_cubed = mean_radius(zone_idx)**3
         true_radius_cubed = dexp(ln10*3.0d0*log_radius(zone_idx))
         r0_geom_factor = (mean_radius_cubed/true_radius_cubed)* &
              mean_radius(zone_idx)**2
         rotation_param = rotation_param_const*omega(zone_idx)**2* &
              mean_radius_cubed/(dexp(ln10*log_mass(zone_idx))* &
              (2.0d0+eta_squared(zone_idx)))
         eta_squared_zone = eta_squared(zone_idx)
!  EVALUATE INTEGRAL AT THE ZONE CENTER
         call shell_inertia_integral(rotation_param, eta_squared_zone, dlnr0_dlnr, &
              r0_geom_factor, moment_of_inertia_per_mass, di_domega_per_mass)
         moment_of_inertia(zone_idx) = moment_of_inertia_per_mass* &
              shell_mass(zone_idx)
         di_domega(zone_idx) = di_domega_per_mass*shell_mass(zone_idx)
         prev_log_mean_radius = dlog(mean_radius(zone_idx))
         prev_log_true_radius = ln10*log_radius(zone_idx)
      end do

      return
end subroutine zone_moments_of_inertia

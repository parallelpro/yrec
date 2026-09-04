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
! omega**2*r**3/(GM). When walpcz.ne.0 (imposed convection-zone
! rotation law) the thin spherical shell formula (2/3)*m*r**2 is used
! for every zone instead, with di_domega = 0, and the innermost zone
! (if zone_start = 1) gets the solid-sphere factor 2/5 at the mean of
! log_radius(1) and log_radius(2).
subroutine zone_moments_of_inertia(eta_squared, log_radius, log_mass, shell_mass, zone_start, &
     zone_end, omega, mean_radius, moment_of_inertia, di_domega)
      use star_info_lib, only: star, json
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: eta_squared(json), log_radius(json), &
           log_mass(json), shell_mass(json)
      integer, intent(in) :: zone_start, zone_end
      double precision, intent(in) :: omega(json), mean_radius(json)
      double precision, intent(out) :: moment_of_inertia(json), di_domega(json)
      integer :: zone_idx
      double precision :: rotation_param_const, mean_zone_radius, dlnr0_dlnr, &
           mean_radius_cubed, true_radius_cubed, r0_geom_factor, &
           rotation_param, eta_squared_zone, moment_of_inertia_per_mass, &
           di_domega_per_mass

      rotation_param_const = cc13*5.0d0/exp(ln10*cgl)
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
         mean_radius_cubed = mean_radius(zone_idx)**3
         true_radius_cubed = exp(ln10*3.0d0*log_radius(zone_idx))
         r0_geom_factor = (mean_radius_cubed/true_radius_cubed)* &
              mean_radius(zone_idx)**2
         rotation_param = rotation_param_const*omega(zone_idx)**2* &
              mean_radius_cubed/(exp(ln10*log_mass(zone_idx))* &
              (2.0d0+eta_squared(zone_idx)))
         eta_squared_zone = eta_squared(zone_idx)
!  EVALUATE INTEGRAL AT THE ZONE CENTER
         call shell_inertia_integral(rotation_param, eta_squared_zone, dlnr0_dlnr, &
              r0_geom_factor, moment_of_inertia_per_mass, di_domega_per_mass)
         moment_of_inertia(zone_idx) = moment_of_inertia_per_mass* &
              shell_mass(zone_idx)
         di_domega(zone_idx) = di_domega_per_mass*shell_mass(zone_idx)
      end do

      return
end subroutine zone_moments_of_inertia

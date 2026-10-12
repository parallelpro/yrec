!----------------------------------------------------------------------
! moments_of_inertia
!----------------------------------------------------------------------
! Each zone's moment of inertia and dI/d(omega), from the spherical value
! and a power-series correction in omega**2 r**3/(GM) (Law 1980).
module moments_of_inertia
      implicit none
      private
      public :: zone_moments_of_inertia

contains

! ---------------------------------------------------------------
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
      use physical_constants
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

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! intmom
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original intmom.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Returns the moment of inertia per unit mass (moment_of_inertia_per_mass,
! originally HIM) and dI/d(omega) per unit mass (di_domega_per_mass,
! originally QIWM) at a single zone center, given the local rotation
! parameter (rotation_param, originally A -- essentially
! omega**2*r0**3/(GM*(2+eta**2)), a measure of the centrifugal-to-
! gravitational force ratio used to expand the moment of inertia in
! powers of the rotational distortion) and the oblateness parameter
! eta_squared (originally ETA2X). Called once per zone by zone_moments_of_inertia.f90.
subroutine shell_inertia_integral(rotation_param, eta_squared, dlnr0_dlnr, r0_geom_factor, &
     moment_of_inertia_per_mass, di_domega_per_mass)
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: rotation_param, eta_squared, &
           dlnr0_dlnr, r0_geom_factor
      double precision, intent(out) :: moment_of_inertia_per_mass, &
           di_domega_per_mass


! series_coeff: coefficients of the power series (in rotation_param) used
! to correct the spherical moment of inertia for rotational distortion.
      double precision :: series_coeff(5)
      data series_coeff/0.6D0,0.72857142857143D0,0.50884115884115D0, &
           0.49728653699242D0,0.42852079651326D0/
      integer :: term_idx
      double precision :: leading_factor, series_sum, series_sum_domega

! RETURN MOMENT OF INERTIA PER UNIT MASS AT THE SPECIFIED POINT(HI)
! AND DI/D OMEGA PER UNIT MASS(QIW)
      leading_factor = cc23*dlnr0_dlnr*r0_geom_factor
      series_sum = 0.0d0
      series_sum_domega = 0.0d0
      do term_idx = 1,5
         series_sum = series_sum + series_coeff(term_idx)* &
              (term_idx*eta_squared + 5.0d0)*pow(rotation_param, term_idx)
         series_sum_domega = series_sum_domega + &
              term_idx*series_coeff(term_idx)*(term_idx*eta_squared + &
              5.0d0)*pow(rotation_param, term_idx)
      end do
      moment_of_inertia_per_mass = leading_factor*(1.0d0 + 1.5d-1*series_sum)
      di_domega_per_mass = leading_factor*(1.5d-1*series_sum_domega)

      return
end subroutine shell_inertia_integral

end module moments_of_inertia

!----------------------------------------------------------------------
! moments_of_inertia
!----------------------------------------------------------------------
! Each zone's moment of inertia I and its omega-derivative term
! di_domega: the spherical value times a power series in the distortion
! amplitude A = 5 omega^2 r0^3 / (3 G m (2 + eta_2)) (Law 1980).
module moments_of_inertia
      use rotation_configuration, only: rotation_settings
      implicit none
      private
      public :: zone_moments_of_inertia

contains

! ---------------------------------------------------------------
! zone_moments_of_inertia (legacy momi): I (moment_of_inertia) and
! di_domega of zones zone_start..zone_end, given r0 (mean_radius) and
! eta_2 (eta_squared) from shape. Method of Law (1980, PhD thesis, Yale,
! p. 61), per unit mass (shell_inertia_integral):
!   I/dm = (2/3) (r0^5/r^3) (dln r0/dln r)
!          [1 + 0.15 sum_{n=1..5} c_n (n eta_2 + 5) A^n],
!   di_domega/dm = (2/3) (r0^5/r^3) (dln r0/dln r)
!          0.15 sum_{n=1..5} n c_n (n eta_2 + 5) A^n,
! with r = 10**log_radius, A = 5 omega^2 r0^3 / (3 G m (2 + eta_2)) and
! dln r0/dln r taken as 1. Since A ~ omega^2, di_domega is
! (omega/2) dI/domega at fixed r0, not dI/domega.
! If rotation_settings%cz_rotation_law_exponent is not 0 (an imposed
! convection-zone law), every zone is instead a thin spherical shell,
! I = (2/3) dm r^2 with di_domega = 0, and zone 1 (when zone_start = 1)
! a sphere, I = (2/5) dm r^2 at r = sqrt(r(1) r(2)).
subroutine zone_moments_of_inertia(eta_squared, log_radius, log_mass, shell_mass, zone_start, &
     zone_end, omega, mean_radius, moment_of_inertia, di_domega)
      use star_info_lib, only: json
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
! An imposed convection-zone law: spherical shells, and a sphere for
! zone 1 (mean_zone_radius is log10 r there).
      if (rotation_settings%cz_rotation_law_exponent.ne.0.0d0) then
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
! dlnr0_dlnr = d ln r0/d ln r, taken as 1; r0_geom_factor = r0^5/r^3
! and rotation_param = A.
         dlnr0_dlnr = 1.0d0
         mean_radius_cubed = mean_radius(zone_idx)**3
         true_radius_cubed = exp(ln10*3.0d0*log_radius(zone_idx))
         r0_geom_factor = (mean_radius_cubed/true_radius_cubed)* &
              mean_radius(zone_idx)**2
         rotation_param = rotation_param_const*omega(zone_idx)**2* &
              mean_radius_cubed/(exp(ln10*log_mass(zone_idx))* &
              (2.0d0+eta_squared(zone_idx)))
         eta_squared_zone = eta_squared(zone_idx)
! The series at the zone, per unit mass.
         call shell_inertia_integral(rotation_param, eta_squared_zone, dlnr0_dlnr, &
              r0_geom_factor, moment_of_inertia_per_mass, di_domega_per_mass)
         moment_of_inertia(zone_idx) = moment_of_inertia_per_mass* &
              shell_mass(zone_idx)
         di_domega(zone_idx) = di_domega_per_mass*shell_mass(zone_idx)
      end do

      return
end subroutine zone_moments_of_inertia

! ---------------------------------------------------------------
! shell_inertia_integral: I and di_domega per unit mass at one zone
! (legacy HIM and QIWM), given A (rotation_param), eta_2 (eta_squared),
! dln r0/dln r (dlnr0_dlnr) and r0^5/r^3 (r0_geom_factor):
!   I/dm         = (2/3) (dln r0/dln r) (r0^5/r^3) (1 + 0.15 S0),
!   di_domega/dm = (2/3) (dln r0/dln r) (r0^5/r^3) 0.15 S1,
!   S0 = sum_{n=1..5} c_n (n eta_2 + 5) A^n,
!   S1 = sum_{n=1..5} n c_n (n eta_2 + 5) A^n.
subroutine shell_inertia_integral(rotation_param, eta_squared, dlnr0_dlnr, r0_geom_factor, &
     moment_of_inertia_per_mass, di_domega_per_mass)
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: rotation_param, eta_squared, &
           dlnr0_dlnr, r0_geom_factor
      double precision, intent(out) :: moment_of_inertia_per_mass, &
           di_domega_per_mass


! series_coeff: the c_n of Law's (1980) series in A that corrects the
! spherical moment of inertia for the rotational distortion.
      double precision :: series_coeff(5)
      data series_coeff/0.6D0,0.72857142857143D0,0.50884115884115D0, &
           0.49728653699242D0,0.42852079651326D0/
      integer :: term_idx
      double precision :: leading_factor, series_sum, series_sum_domega

! leading_factor = (2/3)(dln r0/dln r)(r0^5/r^3); series_sum = S0 and
! series_sum_domega = S1.
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

!----------------------------------------------------------------------
! rotating_shape
!----------------------------------------------------------------------
! The shape of a rotating star in the Kippenhahn & Thomas (1970) /
! Endal & Sofia (1976) formalism. Each level surface is
! r = r0 (1 - A P2(cos theta)); shape finds r0 and Radau's eta_2 zone by
! zone, and rotation_shape_factors averages the effective gravity over
! each surface to give the structure-equation factors fp and ft and the
! mean gravity <g>.
module rotating_shape
      use rotation_configuration, only: rotation_settings
      implicit none
      private
      public :: rotation_shape_factors, shape

contains

! ---------------------------------------------------------------
! rotation_shape_factors: the rotational corrections fp
! (pressure_rotation_factor) and ft (temperature_rotation_factor) to the
! hydrostatic and radiative-transfer equations (henyey_coefficients),
! for zones 1..num_points.
!   1. shape gives each zone's r0 and eta_2 (eta2).
!   2. The distortion integral
!        aint(r0) = int_0^r0 (rho/m) omega^2 r0'^7 (5+eta_2)/(2+eta_2) dr0'
!      is summed shell by shell. Across a shell it is Simpson's rule:
!      Richardson extrapolation (polint) of the one- and two-panel
!      trapezoid sums of refine_distortion_integral.
!   3. average_over_equipotential integrates g, 1/g and 1 over the
!      shell's level surface.
!   4. With r_psi = 10**log_radius (the radius of the sphere of the same
!      volume), S the area of the level surface and <x> = (1/S) int x dS,
!        fp = 4 pi r_psi^4 / (G m S <1/g>),
!        ft = (4 pi r_psi^2 / S)^2 / (<g> <1/g>),
!      and mean_gravity = <g> (cm/s^2).
! Also returns eta2 and r0 (cm); ierr comes from polint.
subroutine rotation_shape_factors(log_density, log_radius, log_mass, num_points, omega, &
     eta2, pressure_rotation_factor, temperature_rotation_factor, &
     mean_gravity, r0, ierr)
      use rotation_scratch_lib

      use composition_layout, only: json
      use physical_constants
      use numerics_lib, only: polint
      use shell_quadrature, only: refine_distortion_integral, average_over_equipotential
      use math_lib
      implicit none

      double precision, intent(in) :: log_density(json), log_radius(json), &
           log_mass(json)
      integer, intent(in) :: num_points
      double precision, intent(in) :: omega(json)
      double precision, intent(inout) :: eta2(json), r0(json)
      double precision, intent(out) :: pressure_rotation_factor(json), &
           temperature_rotation_factor(json), mean_gravity(json)
      integer, intent(out) :: ierr

      double precision :: extrap_step(20), xa(20), ya(20), aint0(10)
      double precision :: previous_shell_mass
      double precision :: b_coefficient, newton_g, eps
      integer :: jmax, k, km
      double precision :: prev_eta22, prev_omega_sq, prev_density, prev_r0, &
           prev_aint, prev_q
      integer :: i, j, n1, j1
      double precision :: density, g_times_mass, centrifugal_factor, &
           r0_cubed, distortion_a, shell_mass, eta22, omega_sq
      double precision :: aint1, dint, aint
      double precision :: g0, ginv0, sphi, ginv, f0, q

      ierr = 0
! r0 and eta_2 of every zone.
      call shape(log_density, log_radius, log_mass, 1, num_points, omega, &
           eta2, r0)
! The surface integrals run over colatitude 0 to b_coefficient = pi/2.
! The distortion integral takes jmax = 2 trapezoid refinements and a
! k = 2 point extrapolation, with relative tolerance eps. Because
! jmax = k, the result is always the j = 2 extrapolation (Simpson's
! rule); the eps test cannot change it.
      b_coefficient = 0.125d0*c4pi
      newton_g = exp(ln10*cgl)
      eps = 1.0d-6
      jmax = 2
      k = 2
      km = k-1
! The values at the inner edge of each shell (the previous zone), which
! refine_distortion_integral interpolates across the shell: m, eta_2 + 2,
! omega^2, rho, r0, aint and the integrand q. At the center: m = 0,
! eta_2 = 0, omega and rho of zone 1, r0 = 0, aint = 0 and q = 0.
      previous_shell_mass = 0.0d0
      prev_eta22 = 2.0d0
      prev_omega_sq = omega(1)**2
      prev_density = exp(ln10*log_density(1))
      prev_r0 = 0.0d0
      prev_aint = 0.0d0
      prev_q = 0.0d0
! fp and ft zone by zone, from eta_2 and r0.
      do i = 1,num_points
         density = exp(ln10*log_density(i))
         g_times_mass = newton_g*exp(ln10*log_mass(i))
! The distortion amplitude A = 5 omega^2 r0^3 / (3 G m (2 + eta_2))
! (distortion_a); centrifugal_factor is A / r0^3.
         centrifugal_factor = 5.0d0*cc13*omega(i)**2/ &
              (g_times_mass*(2.0d0+eta2(i)))
         r0_cubed = r0(i)**3
         distortion_a = centrifugal_factor*r0_cubed
         shell_mass = exp(ln10*log_mass(i))
         eta22 = eta2(i)+2.0d0
         omega_sq = omega(i)**2
         extrap_step(1) = 1.0d0
         do j = 1,jmax
! Trapezoid estimate j of the integral from r0(i-1) to r0(i); from
! j = k on, extrapolate the last k estimates to zero step (extrap_step
! is the squared step, quartered at each refinement).
            call refine_distortion_integral(prev_r0, r0(i), aint0(j), j, density, &
                 prev_density, shell_mass, previous_shell_mass, omega_sq, &
                 prev_omega_sq, eta22, prev_eta22, q, prev_q)
            if (j.ge.k) then
               n1 = j - km
               do j1 = 1,k
                  xa(j1) = extrap_step(n1)
                  ya(j1) = aint0(n1)
                  n1 = n1 + 1
               end do
               call polint(xa, ya, k, 0.0d0, aint1, dint, ierr)
               if (ierr /= 0) return
               if (dabs(dint).lt.eps*dabs(aint1)) exit
            end if
            aint0(j+1) = aint0(j)
            extrap_step(j+1) = 0.25d0*extrap_step(j)
         end do
! aint from the center to r0(i).
         aint = prev_aint + aint1
! Gauss-Legendre integrals over the level surface: g0 = int g dS,
! ginv0 = int dS/g and sphi = S.
         call average_over_equipotential(g0, ginv0, sphi, b_coefficient, r0, log_mass, &
              aint, q, omega_sq, distortion_a, i)
         mean_gravity(i) = g0/sphi
! f0 = log10(4 pi r_psi^4 / (S <1/g>)); fp and ft as in the header.
         ginv = log10(ginv0)
         f0 = c4pil + 4.0d0*log_radius(i) - ginv
         pressure_rotation_factor(i) = exp(ln10*(f0 - cgl - log_mass(i)))
         temperature_rotation_factor(i) = exp(ln10*(f0 + c4pil))/g0
! This zone becomes the inner edge of the next shell.
         prev_aint = aint
         prev_q = q
         previous_shell_mass = shell_mass
         prev_eta22 = eta22
         prev_omega_sq = omega_sq
         prev_density = density
         prev_r0 = r0(i)
      end do

      return
end subroutine rotation_shape_factors

! ---------------------------------------------------------------
! shape: the distortion of the level surfaces of zones
! zone_start..zone_end (Endal & Sofia). A level surface is
!   r = r0 (1 - A P2(cos theta)),  A = 5 omega^2 r0^3 / (3 G m (2 + eta_2)),
! with P2 the second Legendre polynomial. eta_2 (eta2; Radau's eta_2,
! not a square) solves Radau's equation
!   r d(eta_2)/dr + 6 (rho/rho_bar)(eta_2 + 1) + eta_2 (eta_2 - 1) = 6,
!   eta_2(0) = 0,  rho_bar = 3 m / (4 pi r^3),
! integrated outward by fourth-order Runge-Kutta (radau_rk4_step). r0
! makes the level surface enclose the volume of the sphere of radius
! r_psi = 10**log_radius,
!   r_psi^3 = r0^3 (1 + (3/5) A^2 - (2/35) A^3),
! solved by Newton's method in r0^3. As eta_2 depends on r0, the two
! are iterated in turn, at most rotation_settings%shape_outer_iterations
! times, each r0 solve taking at most shape_inner_iterations Newton
! steps; shape_tolerance is the relative tolerance on r0^3, and its cube
! root that on r0 and eta_2. For zone_start > 1 the integration starts
! from zone zone_start - 1, using its eta2 and r0 from an earlier call.
subroutine shape(log_density, log_radius, log_mass, zone_start, zone_end, &
     omega, eta2, r0)
      use composition_layout, only: json
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: log_density(json), log_radius(json), &
           log_mass(json)
      integer, intent(in) :: zone_start, zone_end
      double precision, intent(in) :: omega(json)
      double precision, intent(inout) :: eta2(json)
      double precision, intent(out) :: r0(json)
! --- locals ---
      double precision :: cg, c1, c2, c3, c4
      integer :: j_begin, i, j, k
      double precision :: rho_bar, r_phi, density, gm, r_phi_cubed, &
           r0_cubed, fact, delta_r0_cubed
      double precision :: rho_bar_prev, r_phi_prev, rho_prev
      double precision :: dr, rho_avg, rho_bar_avg
      double precision :: r0_estimate, acc_tol, err, eta2_temp
      double precision :: a_param

! G; c1 = 3/5 and c2 = 2/35 are the coefficients of the volume relation
! r_psi^3 = r0^3 (1 + c1 A^2 - c2 A^3), and c3 = 3 c1, c4 = 4 c2 those of
! its derivative with respect to r0^3, 1 + c3 A^2 - c4 A^3.
      cg = exp(ln10*cgl)
      c1 = 0.6d0
      c2 = 2.0d0/3.5d1
      c3 = 3.0d0*c1
      c4 = 4.0d0*c2
      j_begin = zone_start
      if (zone_start.eq.1) then
! Zone 1, from the center: eta_2 = 6 (1 - rho/rho_bar) there. Below,
! the suffix _prev marks a value at the previous zone.
         rho_bar = exp(ln10*(log_mass(1) - c4pi3l - 3.0d0*log_radius(1)))
         r_phi = exp(ln10*log_radius(1))
         density = exp(ln10*log_density(1))
         eta2(1) = 6.d0*(1.0d0 - density/rho_bar)
! Newton iteration for r0^3 at fixed r_psi and eta_2, as in the zone
! loop below. This one stops on a change <= shape_tolerance, the zone
! loop's on < shape_tolerance (as in the original).
         gm = cg*exp(ln10*log_mass(1))
         r_phi_cubed = r_phi**3
         r0_cubed = r_phi_cubed
         fact = 5.0d0*cc13*omega(1)**2/(gm*(2.0d0+eta2(1)))
         do j = 1,rotation_settings%shape_inner_iterations
            a_param = fact*r0_cubed
            delta_r0_cubed = (r_phi_cubed-r0_cubed*(1.0d0 + c1*a_param**2 - c2*a_param**3))/ &
            (1.0d0 + c3*a_param**2 - c4*a_param**3)
            r0_cubed = r0_cubed + delta_r0_cubed
            if(dabs(delta_r0_cubed/r0_cubed).le.rotation_settings%shape_tolerance)exit
         end do
         r0(1) = pow(r0_cubed, cc13)
         if (zone_end.eq.1) return
         rho_bar_prev = rho_bar
         r_phi_prev = r_phi
         rho_prev = density
         j_begin = 2
      else
! Start from the zone below zone_start.
         rho_bar_prev = exp(ln10*(log_mass(zone_start-1) - c4pi3l - &
              3.0d0*log_radius(zone_start-1)))
         r_phi_prev = exp(ln10*log_radius(zone_start-1))
         rho_prev = exp(ln10*log_density(zone_start-1))
      end if
! eta_2 and r0 of the remaining zones.
      do i = j_begin,zone_end
         rho_bar = exp(ln10*(log_mass(i) - c4pi3l - 3.0d0*log_radius(i)))
         density = exp(ln10*log_density(i))
         r_phi = exp(ln10*log_radius(i))
         gm = cg*exp(ln10*log_mass(i))
         dr = r_phi - r_phi_prev
! First guess for eta_2: integrate Radau's equation in r_psi from the
! previous zone, i.e. take r0 = r_psi. rho and rho_bar at the midpoint
! of the step are the means of their values at its two ends.
         rho_avg = 0.5d0*(density+rho_prev)
         rho_bar_avg = 0.5d0*(rho_bar+rho_bar_prev)
         eta2(i) = radau_rk4_step(dr, eta2(i-1), rho_prev, rho_bar_prev, &
              r_phi_prev, rho_avg, rho_bar_avg, density, rho_bar, r_phi)
         r_phi_cubed = r_phi**3
         r0_cubed = r_phi_cubed
         r0_estimate = r_phi
         acc_tol = pow(rotation_settings%shape_tolerance, cc13)
! Alternate between r0 at fixed eta_2 and eta_2 at the new r0.
         do k = 1,rotation_settings%shape_outer_iterations
            fact = 5.0d0*cc13*omega(i)**2/(gm*(2.0d0+eta2(i)))
! Newton iteration on x = r0^3 at fixed r_psi and eta_2:
!   x <- x + (r_psi^3 - x (1 + c1 A^2 - c2 A^3)) / (1 + c3 A^2 - c4 A^3),
! with A = fact x. (This one stops on < shape_tolerance, the zone-1
! loop above on <=, as in the original.)
            do j = 1,rotation_settings%shape_inner_iterations
               a_param = fact*r0_cubed
               delta_r0_cubed = (r_phi_cubed-r0_cubed*(1.0d0 + c1*a_param**2 - c2*a_param**3))/ &
               (1.0d0 + c3*a_param**2 - c4*a_param**3)
               r0_cubed = r0_cubed + delta_r0_cubed
               if(dabs(delta_r0_cubed/r0_cubed).lt.rotation_settings%shape_tolerance)exit
            end do
            r0(i) = pow(r0_cubed, cc13)
            err = r0(i) - r0_estimate
! Done if r0 differs from the radius eta_2 was found at by no more than
! acc_tol = shape_tolerance^(1/3), relative.
            if(dabs(err)/r0(i).le.acc_tol) exit
! Otherwise integrate Radau's equation again, now in r0 from r0(i-1) to
! r0(i), and stop once eta_2 changes by no more than acc_tol.
            eta2_temp = eta2(i)
            dr = r0(i) - r0(i-1)
            eta2(i) = radau_rk4_step(dr, eta2(i-1), rho_prev, rho_bar_prev, &
                 r0(i-1), rho_avg, rho_bar_avg, density, rho_bar, r0(i))
            err = eta2_temp - eta2(i)
            if(dabs(err).le.acc_tol) exit
            r0_estimate = r0(i)
         end do
         rho_prev = density
         rho_bar_prev = rho_bar
         r_phi_prev = r_phi
      end do
      return

contains

! One fourth-order Runge-Kutta step of Radau's equation,
!   d(eta)/dr = (6 - 6 (rho/rho_bar)(eta + 1) - eta (eta - 1)) / r,
! from (r_prev, eta_prev) to r = r_prev + dr, with rho and rho_bar given
! at the start, the midpoint (rho_avg, rho_bar_avg) and the end of the
! step.
      pure function radau_rk4_step(dr, eta_prev, rho_prev, rho_bar_prev, &
           r_prev, rho_avg, rho_bar_avg, rho, rho_bar, r) result(eta_new)
      double precision, intent(in) :: dr, eta_prev, rho_prev, rho_bar_prev, &
           r_prev, rho_avg, rho_bar_avg, rho, rho_bar, r
      double precision :: eta_new
      double precision :: deta1, deta2, deta3, deta4, eta_avg, r_avg
! The increment from the slope at the start of the step.
      deta1 = dr*(6.0d0 - 6.0d0*rho_prev*(eta_prev + 1.0d0)/rho_bar_prev &
      - eta_prev*(eta_prev - 1.0d0))/r_prev
      eta_avg = eta_prev + 0.5d0*deta1
      r_avg = r_prev + 0.5d0*dr
! The increment from the slope at the midpoint, at the eta of the first
! half step.
      deta2 = dr*(6.0d0 - 6.0d0*rho_avg*(eta_avg + 1.0d0)/rho_bar_avg &
      - eta_avg*(eta_avg - 1.0d0))/r_avg
      eta_avg = eta_prev + 0.5d0*deta2
! The midpoint again, at the eta of the refined half step.
      deta3 = dr*(6.0d0 - 6.0d0*rho_avg*(eta_avg + 1.0d0)/rho_bar_avg &
      - eta_avg*(eta_avg - 1.0d0))/r_avg
      eta_avg = eta_prev + deta3
! The end of the step, at the eta of a full step with the second
! midpoint slope.
      deta4 = dr*(6.0d0 - 6.0d0*rho*(eta_avg + 1.0d0)/rho_bar &
      - eta_avg*(eta_avg - 1.0d0))/r
! eta_new = eta_prev + (deta1 + 2 deta2 + 2 deta3 + deta4) / 6.
      eta_new = eta_prev+cc13*(0.5d0*deta1+deta2+deta3+0.5d0*deta4)
      end function radau_rk4_step

end subroutine shape

end module rotating_shape

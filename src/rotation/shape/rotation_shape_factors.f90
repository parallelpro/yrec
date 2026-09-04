!----------------------------------------------------------------------
! fpft
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original fpft.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Computes the rotational distortion correction factors pressure_
! rotation_factor/temperature_rotation_factor (FP/FT) applied to the
! hydrostatic/thermal structure equations in henyey_coefficients.f90 to account for
! centrifugal flattening (the Kippenhahn-Thomas / Endal-Sofia "shape"
! formalism): first calls shape (shape.f90) to get the run of the
! reference equipotential radius r0 and the distortion parameter eta2,
! then for each shell integrates the distorted potential (via
! trapzd/polint Richardson extrapolation) and averages the resulting
! local effective gravity over the shell (via qgauss with
! equipotential_integrand) to get mean_gravity, from which FP/FT and
! the diagnostic potential terms rot_scr%phisp/phirot/phidis are built.
subroutine rotation_shape_factors(log_density, log_radius, log_mass, num_points, omega, &
     eta2, pressure_rotation_factor, temperature_rotation_factor, &
     mean_gravity, r0, ierr)
      use rotation_scratch_lib

      use star_info_lib, only: json
      use phys_const_lib
      use numerics_lib
      use math_lib
      implicit none
! the shape integrand passed to numerics' qgauss
      external equipotential_integrand

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
           prev_aint, aintt, prev_q
      integer :: i, j, n1, j1
      double precision :: density, g_times_mass, centrifugal_factor, &
           r0_cubed, distortion_a, shell_mass, eta22, omega_sq
      double precision :: aint1, dint, aint
      double precision :: g0, ginv0, sphi, ginv, f0, q, rphi, rphi3

      ierr = 0
! FIND THE RUN OF R0 AND ETA2 FOR THE MODEL
      call shape(log_density, log_radius, log_mass, 1, num_points, omega, &
           eta2, r0)
      b_coefficient = 0.125d0*c4pi
      newton_g = exp(ln10*cgl)
      eps = 1.0d-6
      jmax = 2
      k = 2
      km = k-1
! THE VALUES OF THE PREVIOUS SHELL MPHI,ETA2+2,OMEGA**2,RHO, AND R0 ARE
! NEEDED FOR THE EVALUATION OF THE DISTORTED POTENTIAL.
! THE BELOW VALUES ARE THE CENTRAL VALUES OF THE ABOVE IN ORDER.
      previous_shell_mass = 0.0d0
      prev_eta22 = 2.0d0
      prev_omega_sq = omega(1)**2
      prev_density = exp(ln10*log_density(1))
      prev_r0 = 0.0d0
      prev_aint = 0.0d0
      aintt = 0.0d0
      prev_q = 0.0d0
! NOW CALCULATE FP AND FT USING THE ETA2 AND R0 VALUES
      do i = 1,num_points
         density = exp(ln10*log_density(i))
         g_times_mass = newton_g*exp(ln10*log_mass(i))
         centrifugal_factor = 5.0d0*cc13*omega(i)**2/ &
              (g_times_mass*(2.0d0+eta2(i)))
         r0_cubed = r0(i)**3
         distortion_a = centrifugal_factor*r0_cubed
         shell_mass = exp(ln10*log_mass(i))
         eta22 = eta2(i)+2.0d0
         omega_sq = omega(i)**2
         extrap_step(1) = 1.0d0
         do j = 1,jmax
! EVALUATE THE INTEGRAL AINT FROM 0 TO R0 USING THE TRAPEZOIDAL RULE
            call trapzd(prev_r0, r0(i), aint0(j), j, density, &
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
         aint = prev_aint + aint1
! FIND <G> AND <G-1> ACROSS THE SHELL BY GAUSSIAN QUADRATURE
! (func is passed as the integrand -- 2026, phase four step 2)
         call qgauss(equipotential_integrand, g0, ginv0, sphi, b_coefficient, r0, log_mass, &
              aint, q, omega_sq, distortion_a, i)
         mean_gravity(i) = g0/sphi
         ginv = log10(ginv0)
         f0 = c4pil + 4.0d0*log_radius(i) - ginv
         pressure_rotation_factor(i) = exp(ln10*(f0 - cgl - log_mass(i)))
         temperature_rotation_factor(i) = exp(ln10*(f0 + c4pil))/g0
! OUTPUT DATA
         rphi = exp(ln10*log_radius(i))
         rphi3 = rphi**3
         rot_scr%phisp(i) = g_times_mass/rphi
         rot_scr%phirot(i) =  omega_sq*rphi**2
         rot_scr%phidis(i) = c4pi*cc13*aint/rphi3
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

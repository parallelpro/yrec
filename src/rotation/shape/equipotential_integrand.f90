!----------------------------------------------------------------------
! func
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original func.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Evaluates, at colatitude x on the distorted equipotential surface of
! shell i (parameterized by r = r0(i)*(1 - a*P2(cos x)), the Roche/
! centrifugal-distortion approximation used throughout fpft/qgauss),
! the local effective gravity magnitude g = |grad(phi)| (local_gravity)
! and the surface-area-element factor s (area_element) needed by the
! <g> and <g^-1> integrals over the shell computed in qgauss. r0/
! log_mass/aint/q/w2/a/i are passed through unchanged from the caller
! (qgauss, which itself receives them unchanged from fpft).
subroutine equipotential_integrand(colatitude, local_gravity, area_element, r0, log_mass, &
     aint, q, w2, a, i)
      use star_info_lib, only: json

      use phys_const_lib
      implicit none



      double precision, intent(in) :: colatitude
      double precision, intent(out) :: local_gravity, area_element
      double precision, intent(in) :: r0(json), log_mass(json)
      double precision, intent(in) :: aint, q, w2, a
      integer, intent(in) :: i
      double precision :: cs, ss, p2, qr0r, qr0t, r, r2, const, qphir, qphith

      cs = dcos(colatitude)*dsin(colatitude)
      ss = dsin(colatitude)**2
! P2 = .5*(3COS(X)**2 - 1) = .5(3*(1-SIN(X)**2) - 1) = 1 - 1.5*SIN(X)**2
      p2 = 1.0d0 - 1.5d0*ss
! DERIVATIVES OF R0 WITH RESPECT TO R AND THETA
      qr0r = 1.0d0/(1.0d0 - 4.0d0*a*p2)
      qr0t = 3.0d0*cs*a*qr0r/(1.0d0 - a*p2)
! USE THE RELATION FOR R ON AN EQUIPOTENTIAL
      r = r0(i)*(1.0d0 - a*p2)
      r2 = r**2
      const = c4pi*cc13/(r2*r2)
! CALCULATE THE DERIVATIVES OF PHI WITH RESPECT TO R AND THETA
! D(PHI)/DR = GM/R2 + 12PI*G*P2*AINT/5R**4 -4PI*G*P2*Q*(DR0/DR)/5R**3
! -OMEGA**2*R*SIN(THETA)**2 - ASSUME DR0/DR = 1
      qphir = dexp(ln10*(cgl+log_mass(i)))/r2 - const*p2*(3.0d0*aint &
           - r*q*qr0r) - w2*r*ss
! D(PHI)/D(THETA) = 4PI*G*3DCOS(THETA)DSIN(THETA)*AINT/5R**4 -
! 4PI*G*P2*Q*DR0/D(THETA)/5R**4 - OMEGA**2*R*DCOS(THETA)SIN(THETA)
      qphith = const*(p2*q*qr0t-3.0d0*aint*cs) - w2*r*cs
      local_gravity = dsqrt(qphir**2 + qphith**2)
      area_element = c4pi*r*dsin(colatitude)*r0(i)* &
           dsqrt((1.0d0-a*p2)**2 + (3.0d0*cs*a)**2)

      return
end subroutine equipotential_integrand

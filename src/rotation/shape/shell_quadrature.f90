!----------------------------------------------------------------------
! shell_quadrature
!----------------------------------------------------------------------
! The two quadratures of the rotational shape factors
! (rotation_shape_factors), formerly in numerics_lib:
!   refine_distortion_integral (formerly trapzd): one refinement of the
!       trapezoid rule for the distortion integral
!       aint = int_0^r0 (rho/m) r0'**7 omega**2 (3+eta2+2)/(eta2+2) dr0'
!       across a shell, with rho, m, omega**2 and eta2+2 interpolated
!       linearly between the shell's two edges;
!   average_over_equipotential (formerly qgauss): the 10-point
!       Gauss-Legendre integrals of g*s, s/g and s over the colatitude
!       on a shell's equipotential, g and s from equipotential_integrand
!       (formerly func).
module shell_quadrature
      use star_info_lib, only: json
      use physical_constants
      use math_lib
      implicit none
      private
      public :: refine_distortion_integral, average_over_equipotential

contains

!----------------------------------------------------------------------
! refine_distortion_integral (formerly trapzd)
!----------------------------------------------------------------------
! Numerical Recipes' trapzd for the distortion integral over [b1,b2]:
! n = 1 evaluates the integrand q at r0 = b2 and returns the
! two-point estimate s; each later n (called in sequence 2, 3, ...)
! adds the 2**(n-2) new midpoints to s, with rho/sm/w2/eta22 linearly
! interpolated between their values at b1 (rhop, smp, w2p, eta22p) and
! at b2. (The number of midpoints was a SAVEd counter.)
subroutine refine_distortion_integral(b1, b2, s, n, rho, rhop, sm, smp, w2, w2p, eta22, &
     eta22p, q, qp)
      double precision, intent(in) :: b1, b2
      double precision, intent(inout) :: s
      integer, intent(in) :: n
      double precision, intent(in) :: rho, rhop, sm, smp, w2, w2p, &
           eta22, eta22p
      double precision, intent(inout) :: q
      double precision, intent(in) :: qp

! the number of midpoints added at refinement n
      integer :: it

      double precision :: r0, r03, tnm, dr, del, y, sum, drho, dm, &
           deta2, dw2, r03t, rhot, smt, w2t, eta22t, q0
      integer :: j

      r0 = b2
      r03 = r0**3
      if (n.eq.1) then
!  aint = int(0=>r0) (rho/m)*r0'**7*omega**2*(5+eta2)/(2+eta2) dr0'
!  q is the integrand (ro'**7,etc.) evaluated at r0(i)
!  aint and its derivatives w/r/to r and theta are needed to find <g>
       q = (rho*w2*r03*(3.0d0+eta22)/(sm*eta22))*r03*r0
!        q(i) = exp(cln*(hd(i)-hs(i)))*omega(i)**2*r0(i)**6
!    *   *(5.0d0+eta2(i))/(2.0d0+eta2(i))
       s = 0.5d0*(b2-b1)*(qp+q)
      else
       it = ishft(1, n-2)
       tnm = dfloat(it)
       dr = b2 - b1
       del = dr/tnm
       y = b1 + 0.5d0*del
       sum = 0.0d0
       drho = (rho - rhop)/dr
       dm = (sm - smp)/(b2**2 - b1**2)
       deta2 = (eta22 - eta22p)/dr
       dw2 = (w2 - w2p)/dr
       do j = 1, it
          r03t = y**3
! interpolate rho,m,omega,eta2+2 between shell i and shell i-1
! (Fixed in the 2026 bug sweep:) the linear interpolants must be evaluated
! at the offset (y - b1) of the current abscissa, as smt already is;
! the inherited `*del` used the constant sub-interval width instead.
          rhot = rhop+drho*(y - b1)
          smt = smp+dm*(y**2 - b1**2)
          w2t = w2p + dw2*(y - b1)
          eta22t = eta22p + deta2*(y - b1)
! calculate q between shells
          q0 = (rhot*w2t*r03t*(3.0d0+eta22t)/(smt*eta22t))*r03t*y
! q0 = rho*w2*r07t*(3.0d0+eta22)/(sm*eta22)
          sum = sum + q0
          y = y+del
       end do
       s = 0.5d0*(s+del*sum)
      end if

      return
end subroutine refine_distortion_integral

!----------------------------------------------------------------------
! average_over_equipotential (formerly qgauss)
!----------------------------------------------------------------------
! 10-point (5-point-symmetric) Gauss-Legendre quadrature over [0,b]:
! g0g = integral of g*s, ginvg = integral of s/g and sphig = integral
! of s, with g and s from equipotential_integrand at the two abscissas
! symmetric about the midpoint. r0, hs, aint, q, w2, a and i go to the
! integrand unchanged.
subroutine average_over_equipotential(g0g, ginvg, sphig, b, r0, hs, aint, q, w2, a, i)
      double precision, intent(out) :: g0g, ginvg, sphig
      double precision, intent(in) :: b
      double precision, intent(in) :: r0(json), hs(json)
      double precision, intent(in) :: aint, q, w2, a
      integer, intent(in) :: i

! the positive nodes on [-1,1] and their weights
      double precision, parameter :: x(5) = [.14887433898163d0, .43339539412925d0, &
           .67940956829902d0, .86506336668899d0, .97390652851717d0]
      double precision, parameter :: w(5) = [.29552422471475d0, .26926671931d0, &
           .21908636251598d0, .14945134915058d0, .06667134430869d0]

      double precision :: xm, xr, dx, g, s, g2, s2
      integer :: j

      xm = 0.5d0*b
      xr = xm
      g0g = 0.0d0
      ginvg = 0.0d0
      sphig = 0.0d0
      do j = 1, 5
       dx = xr*x(j)
       call equipotential_integrand(xm+dx, g, s, r0, hs, aint, q, w2, a, i)
       call equipotential_integrand(xm-dx, g2, s2, r0, hs, aint, q, w2, a, i)
       g0g = g0g+w(j)*(g*s+g2*s2)
       ginvg = ginvg+w(j)*(s/g+s2/g2)
       sphig = sphig+w(j)*(s+s2)
      end do
      g0g = g0g*xr
      ginvg = ginvg*xr
      sphig = sphig*xr
      return
end subroutine average_over_equipotential

!----------------------------------------------------------------------
! equipotential_integrand (formerly func)
!----------------------------------------------------------------------
! At colatitude x on the distorted equipotential of shell i,
! r = r0(i)*(1 - a*P2(cos x)) (the Roche/centrifugal-distortion
! approximation of the shape formalism): the effective gravity
! g = |grad(phi)| (local_gravity) and the area-element factor s
! (area_element).
subroutine equipotential_integrand(colatitude, local_gravity, area_element, r0, log_mass, &
     aint, q, w2, a, i)
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
      qphir = exp(ln10*(cgl+log_mass(i)))/r2 - const*p2*(3.0d0*aint &
           - r*q*qr0r) - w2*r*ss
! D(PHI)/D(THETA) = 4PI*G*3DCOS(THETA)DSIN(THETA)*AINT/5R**4 -
! 4PI*G*P2*Q*DR0/D(THETA)/5R**4 - OMEGA**2*R*DCOS(THETA)SIN(THETA)
      qphith = const*(p2*q*qr0t-3.0d0*aint*cs) - w2*r*cs
      local_gravity = dsqrt(qphir**2 + qphith**2)
      area_element = c4pi*r*dsin(colatitude)*r0(i)* &
           dsqrt((1.0d0-a*p2)**2 + (3.0d0*cs*a)**2)

      return
end subroutine equipotential_integrand

end module shell_quadrature

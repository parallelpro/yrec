!----------------------------------------------------------------------
! shell_quadrature
!----------------------------------------------------------------------
! The two quadratures of the rotational shape factors
! (rotation_shape_factors):
!   refine_distortion_integral (legacy trapzd): one refinement of the
!       trapezoid rule for a shell's part of the distortion integral
!       aint = int_0^r0 (rho/m) omega^2 r0'^7 (5 + eta_2)/(2 + eta_2) dr0',
!       with rho, omega^2 and eta_2 + 2 linear in r0' across the shell
!       and m linear in r0'^2;
!   average_over_equipotential (legacy qgauss): the 10-point
!       Gauss-Legendre integrals of g, 1/g and 1 over a shell's level
!       surface, with the effective gravity g and the area element from
!       equipotential_integrand (legacy func).
module shell_quadrature
      use star_info_lib, only: json
      use physical_constants
      use math_lib
      implicit none
      private
      public :: refine_distortion_integral, average_over_equipotential

contains

!----------------------------------------------------------------------
! refine_distortion_integral (legacy trapzd)
!----------------------------------------------------------------------
! Numerical Recipes' trapzd for the distortion integral over the shell
! [b1, b2] = [r0(i-1), r0(i)], called with n = 1, 2, ... in turn.
! n = 1 evaluates the integrand q at b2 and returns the one-panel sum
! s = (b2 - b1)(qp + q)/2, qp being the integrand at b1 (from the
! previous shell). Each later n adds the 2^(n-2) new midpoints to s.
! Between the edges, rho, omega^2 (w2) and eta_2 + 2 (eta22) are linear
! in r0 and the enclosed mass sm is linear in r0^2, from their values at
! b1 (rhop, w2p, eta22p, smp) and at b2.
subroutine refine_distortion_integral(b1, b2, s, n, rho, rhop, sm, smp, w2, w2p, eta22, &
     eta22p, q, qp)
      double precision, intent(in) :: b1, b2
      double precision, intent(inout) :: s
      integer, intent(in) :: n
      double precision, intent(in) :: rho, rhop, sm, smp, w2, w2p, &
           eta22, eta22p
      double precision, intent(inout) :: q
      double precision, intent(in) :: qp

! The number of midpoints added at refinement n.
      integer :: it

      double precision :: r0, r03, tnm, dr, del, y, sum, drho, dm, &
           deta2, dw2, r03t, rhot, smt, w2t, eta22t, q0
      integer :: j

      r0 = b2
      r03 = r0**3
      if (n.eq.1) then
! The integrand at r0 = b2,
!   q = rho omega^2 r0^7 (5 + eta_2) / (m (2 + eta_2)),
! with eta22 = eta_2 + 2. q is also d(aint)/dr0 at b2, which
! equipotential_integrand needs, and the next shell's qp.
       q = (rho*w2*r03*(3.0d0+eta22)/(sm*eta22))*r03*r0
       s = 0.5d0*(b2-b1)*(qp+q)
      else
! it = 2^(n-2) new midpoints, del apart, the first at b1 + del/2; the
! interpolation slopes are per unit r0 (per unit r0^2 for the mass).
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
! rho, m, omega^2 and eta_2 + 2 at the midpoint y. (Fixed in 2026: the
! legacy code used the sub-interval width del in place of y - b1.)
          rhot = rhop+drho*(y - b1)
          smt = smp+dm*(y**2 - b1**2)
          w2t = w2p + dw2*(y - b1)
          eta22t = eta22p + deta2*(y - b1)
! The integrand at y.
          q0 = (rhot*w2t*r03t*(3.0d0+eta22t)/(smt*eta22t))*r03t*y
          sum = sum + q0
          y = y+del
       end do
! The refined sum s_n = (s_(n-1) + del * sum of the new values) / 2.
       s = 0.5d0*(s+del*sum)
      end if

      return
end subroutine refine_distortion_integral

!----------------------------------------------------------------------
! average_over_equipotential (legacy qgauss)
!----------------------------------------------------------------------
! 10-point Gauss-Legendre quadrature over colatitude [0, b] on the level
! surface of shell i, evaluating the integrand at the pairs of nodes
! symmetric about b/2:
!   g0g = int g s dtheta,  ginvg = int s/g dtheta,  sphig = int s dtheta,
! with the effective gravity g and area element s from
! equipotential_integrand. With b = pi/2 (s covers both hemispheres)
! these are int g dS, int dS/g and the area S. r0, hs (log10 m), aint,
! q, w2 (omega^2), a (the distortion amplitude A) and i pass through to
! the integrand.
subroutine average_over_equipotential(g0g, ginvg, sphig, b, r0, hs, aint, q, w2, a, i)
      double precision, intent(out) :: g0g, ginvg, sphig
      double precision, intent(in) :: b
      double precision, intent(in) :: r0(json), hs(json)
      double precision, intent(in) :: aint, q, w2, a
      integer, intent(in) :: i

! The positive nodes on [-1,1] and their weights.
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
! equipotential_integrand (legacy func)
!----------------------------------------------------------------------
! At colatitude theta on the level surface of shell i,
! r = r0(i) (1 - a P2(cos theta)), the effective gravity g = |grad Psi|
! (local_gravity) and the area element s = dS/dtheta (area_element), for
!   Psi = -G m/r + (4 pi/3) aint(r0) P2/r^3 - omega^2 r^2 sin^2(theta)/2:
! the enclosed mass, the quadrupole of the distorted interior and the
! centrifugal potential (the shells outside r0 do not enter). aint
! depends on r0 with d(aint)/dr0 = q, and r0 on (r, theta) through the
! shape, with a proportional to r0^3.
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
! P2 = (3 cos^2(theta) - 1)/2 = 1 - 1.5 sin^2(theta).
      p2 = 1.0d0 - 1.5d0*ss
! The derivatives of r0: qr0r = dr0/dr = 1/(1 - 4 a P2) at fixed theta,
! and qr0t = -(1/r) dr0/dtheta = 3 a cos(theta) sin(theta) qr0r/(1 - a P2)
! at fixed r.
      qr0r = 1.0d0/(1.0d0 - 4.0d0*a*p2)
      qr0t = 3.0d0*cs*a*qr0r/(1.0d0 - a*p2)
! r on the level surface, and const = (4 pi/3)/r^4.
      r = r0(i)*(1.0d0 - a*p2)
      r2 = r**2
      const = c4pi*cc13/(r2*r2)
! dPsi/dr = G m/r^2 - (4 pi/3) P2 (3 aint - r q qr0r)/r^4
!           - omega^2 r sin^2(theta).
      qphir = exp(ln10*(cgl+log_mass(i)))/r2 - const*p2*(3.0d0*aint &
           - r*q*qr0r) - w2*r*ss
! (1/r) dPsi/dtheta = (4 pi/3)(P2 q qr0t - 3 aint cos(theta) sin(theta))/r^4
!                     - omega^2 r cos(theta) sin(theta).
! Unlike the radial r q qr0r term, the P2 q qr0t term carries no factor
! r (and (1/r) d/dtheta of aint(r0) would give -r P2 q qr0t).
      qphith = const*(p2*q*qr0t-3.0d0*aint*cs) - w2*r*cs
! g = |grad Psi|, and s = 4 pi r sin(theta) r0 sqrt((1 - a P2)^2
! + (3 a cos(theta) sin(theta))^2): the area per unit colatitude of the
! band at theta and its mirror below the equator, with
! dr/dtheta = 3 a r0 cos(theta) sin(theta) at fixed a.
      local_gravity = dsqrt(qphir**2 + qphith**2)
      area_element = c4pi*r*dsin(colatitude)*r0(i)* &
           dsqrt((1.0d0-a*p2)**2 + (3.0d0*cs*a)**2)

      return
end subroutine equipotential_integrand

end module shell_quadrature

!----------------------------------------------------------------------
! qgauss
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original qgauss.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! 10-point (5-point-symmetric) Gauss-Legendre quadrature over [0,b]
! of three related integrands evaluated by the external routine
! func: g0g = integral of g*s, ginvg = integral of s/g, and sphig =
! integral of s, where (g,s) and (g2,s2) are func's outputs at the
! two symmetric abscissas about the interval midpoint. r0, hs, aint,
! q, w2, a, and i are passed through unchanged to func on each call.
subroutine qgauss(g0g, ginvg, sphig, b, r0, hs, aint, q, w2, a, i)

      use const_lib
      implicit none
      integer, parameter :: json = 5000



      double precision, intent(out) :: g0g, ginvg, sphig
      double precision, intent(in) :: b
      double precision, intent(in) :: r0(json), hs(json)
      double precision, intent(in) :: aint, q, w2, a
      integer, intent(in) :: i

      double precision :: x(5), w(5)
      data x/.14887433898163d0,.43339539412925d0,.67940956829902d0, &
           .86506336668899d0,.97390652851717d0/
      data w/.29552422471475d0,.26926671931d0,.21908636251598d0, &
           .14945134915058d0,.06667134430869d0/

      double precision :: xm, xr, dx, g, s, g2, s2
      integer :: j
      save

      xm = 0.5d0*b
      xr = xm
      g0g = 0.0d0
      ginvg = 0.0d0
      sphig = 0.0d0
      do j = 1, 5
       dx = xr*x(j)
       call func(xm+dx, g, s, r0, hs, aint, q, w2, a, i)
       call func(xm-dx, g2, s2, r0, hs, aint, q, w2, a, i)
       g0g = g0g+w(j)*(g*s+g2*s2)
       ginvg = ginvg+w(j)*(s/g+s2/g2)
       sphig = sphig+w(j)*(s+s2)
      end do
      g0g = g0g*xr
      ginvg = ginvg*xr
      sphig = sphig*xr
      return
end subroutine qgauss

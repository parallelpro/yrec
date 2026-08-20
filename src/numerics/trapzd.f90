!----------------------------------------------------------------------
! trapzd
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original trapzd.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Numerical Recipes-style trapzd: on the n=1 call, evaluates the
! integrand q at the single point r0=b2 and returns the crude
! 2-point trapezoidal estimate s over [b1,b2]. On subsequent calls
! (n>1) it refines s by adding the contribution of the newly
! inserted midpoints, with the interpolated intermediate rho/sm/w2/
! eta22 values linearly interpolated between the (rhop,smp,w2p,
! eta22p) state at b1 and the (rho,sm,w2,eta22) state at b2. it is
! the number of subintervals from the previous call, held across
! calls via SAVE (the Numerical Recipes call-count doubling scheme;
! it used to be passed in explicitly as a dummy argument i -- see
! the commented-out original signature below).
!       SUBROUTINE TRAPZD(B1,B2,S,N,RHO,RHOP,SM,SMP,W2,W2P,ETA22,
!      *ETA22P,Q,QP,I)  ! KC 2025-05-31
subroutine trapzd(b1, b2, s, n, rho, rhop, sm, smp, w2, w2p, eta22, &
     eta22p, q, qp)

      implicit none

! common/const1/: not used in this file; declared only to preserve
! layout. Naming matches dburn.f90/eqburn.f90.
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
      common/const1/ ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi

! common/const2/: not used in this file; declared only to preserve
! layout. Naming matches eqstat2.f90/meqos.f90.
      double precision :: gas_constant, radiation_constant_over_3, ca3l, &
           csig, csigl, cgl, cmkh, cmkhn
      common/const2/ gas_constant, radiation_constant_over_3, ca3l, csig, &
           csigl, cgl, cmkh, cmkhn

      double precision, intent(in) :: b1, b2
      double precision, intent(inout) :: s
      integer, intent(in) :: n
      double precision, intent(in) :: rho, rhop, sm, smp, w2, w2p, &
           eta22, eta22p
      double precision, intent(inout) :: q
      double precision, intent(in) :: qp

      integer :: it
      save

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
!        q(i) = dexp(cln*(hd(i)-hs(i)))*omega(i)**2*r0(i)**6
!    *   *(5.0d0+eta2(i))/(2.0d0+eta2(i))
       s = 0.5d0*(b2-b1)*(qp+q)
       it = 1
      else
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
          rhot = rhop+drho*del
          smt = smp+dm*(y**2 - b1**2)
          w2t = w2p + dw2*del
          eta22t = eta22p + deta2*del
! calculate q between shells
          q0 = (rhot*w2t*r03t*(3.0d0+eta22t)/(smt*eta22t))*r03t*y
! q0 = rho*w2*r07t*(3.0d0+eta22)/(sm*eta22)
          sum = sum + q0
          y = y+del
       end do
       s = 0.5d0*(s+del*sum)
       it = it*2
      end if

      return
end subroutine trapzd

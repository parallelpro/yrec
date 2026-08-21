!----------------------------------------------------------------------
! deutrate
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original deutrate.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Compute the rate of nonequilibrium deuterium burning (excluding the
! abundance factor) at shell i, storing it in common/deuter/ for later
! use by dburn/dburnm. If the shell lies within (or below) the surface
! convection zone, the rate is capped so that deuterium burning cannot
! proceed faster than the local convective overturn timescale.
subroutine deutrate(dl,tl,x,i,itlvl)
      use light_burn_lib
      use turnover_lib
      use const_lib
      implicit none
      integer, parameter :: json=5000

      double precision, intent(in) :: dl, tl, x
      integer, intent(in) :: i, itlvl





      double precision :: c21
      data c21/5.240358E-8/
      save

! T9P13 IS THE TEMPERATURE IN UNITS OF 10^9 DEGREES K TO THE PLUS 1/3
!  POWER.  MINUS IS DENOTED BY M.  HERE T9 IS THE TEMPERATURE IN UNITS
!  OF 10^9 K, CONVERTED FROM THE LOG_10 (T) AND RHO IS THE DENSITY IN
!  CGS UNITS.
      double precision :: rho, t9, t9p13, t9p23, t9m13, t9m23, t9m1
      double precision :: z, tfacdeut, rdeut, rdeutmax, rdeut2

      rho=exp(ln10*dl)
      t9 = exp(ln10*(tl - 9.0d0))
      t9p13 = t9**cc13
      t9p23 = t9p13**2
      t9m13=1.0d0/t9p13
      t9m23=t9m13**2
      t9m1=1.0d0/t9
! MHP 5/02 ADD DEUTERIUM BURNING TERM TO THE CODE
! IF DEUTERIUM IS ABOVE A MINIMUM THRESHOLD VALUE.
! RDEUT IS THE RATE (EXCLUDING FACTORS OF THE
! ABUNDANCES) AND QRTDEUT IS THE DERIVATIVE W/R/T T.
! NOTE THAT SCREENING IS EXCLUDED - REASONABLE GIVEN
! THE LOW TEMPERATURES INVOLVED.
      z = -3.72d0*t9m13
      tfacdeut = 1.0d0+0.112d0*t9p13+3.38d0*t9p23+2.65d0*t9
! FACTOR OF 3.0115D23 REFLECTS AVAGADROS NUMBER DIVIDED BY THE
! MASS OF THE DEUTERON IN AMU
      rdeut = rho*2.240d3*t9p23*exp(z)*tfacdeut*3.0115d23
! NOW LIMIT DEUTERIUM BURNING IN A SURFACE CZ TO BE ON A TIME SCALE
! NO SHORTER THAN THE CONVECTIVE OVERTURN TIMESCALE.
      if(i.ge.light_burn%jcz .and. turnover%convective_turnover_timescale.gt.1.0d0)then
         rdeutmax = 3.0115d23/turnover%convective_turnover_timescale
         rdeut2 = rdeut*x
         if(i.eq.light_burn%jcz)then
! JVS 0712 Commented out write command
!            WRITE(*,911)RDEUT2,RDEUTMAX,TAUCZ
!  911        FORMAT(1P3E15.8)
         endif
         if(rdeut2.gt.rdeutmax)then
! JVS 0712 Commented out write command
!            WRITE(*,*)RDEUT2,RDEUTMAX
            if(x.gt.1.0d-6)then
               rdeut = rdeutmax/x
            endif
         endif
      endif
      light_burn%deuterium_burning_rate(i) = x*rdeut*c21
      if(itlvl.eq.1)then
         light_burn%deuterium_burning_rate_start(i) = light_burn%deuterium_burning_rate(i)
      endif
      return
end subroutine deutrate

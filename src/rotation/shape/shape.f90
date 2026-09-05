!----------------------------------------------------------------------
! shape
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original shape.f; only variable names, source form, and comment
! style were updated.
!
! SHAPE FINDS THE DISTORTION OF A GIVEN SHELL FROM SPHERICAL SYMMETRY;
! IT CALCULATES THE CHARACTERISTIC RADIUS R0 AND ASSOCIATED PARAMETER
! ETA2.  IN POLAR CO-ORDINATES THE RADIUS OF A SHELL IS
! R = R0*(1 - A*P2), WHERE A IS A FUNCTION OF ROTATION RATE AND P2 IS
! THE SECOND LEGENDRE POLYNOMIAL.  FOR MORE INFORMATION SEE ENDAL AND
! SOFIA.
subroutine shape(log_density, log_radius, log_mass, zone_start, zone_end, &
     omega, eta2, r0)
      use star_info_lib, only: star, json
      use phys_const_lib
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

      cg = exp(ln10*cgl)
      c1 = 0.6d0
      c2 = 2.0d0/3.5d1
      c3 = 3.0d0*c1
      c4 = 4.0d0*c2
      j_begin = zone_start
      if (zone_start.eq.1) then
!  CALCULATE FIRST POINT ETA2 USING CENTRAL B.C.
!  NOTE THAT THE SUFFIX 'P' DENOTES A QUANTITY FROM THE PREVIOUS SHELL.
!  E.G. RHOP IS THE DENSITY OF THE LAST SHELL AND RHO IS THE DENSITY OF
!  THE CURRENT ONE.
         rho_bar = exp(ln10*(log_mass(1) - c4pi3l - 3.0d0*log_radius(1)))
         r_phi = exp(ln10*log_radius(1))
         density = exp(ln10*log_density(1))
         eta2(1) = 6.d0*(1.0d0 - density/rho_bar)
! ITERATE FOR R0 GIVEN RPHI AND ETA2
! (Same Newton step as the per-zone loop below, but this one exits on
! .le. acfpft and that one on .lt. acfpft -- kept as in the original.)
         gm = cg*exp(ln10*log_mass(1))
         r_phi_cubed = r_phi**3
         r0_cubed = r_phi_cubed
         fact = 5.0d0*cc13*omega(1)**2/(gm*(2.0d0+eta2(1)))
         do j = 1,star%ctrl%itfp2
            a_param = fact*r0_cubed
            delta_r0_cubed = (r_phi_cubed-r0_cubed*(1.0d0 + c1*a_param**2 - c2*a_param**3))/ &
            (1.0d0 + c3*a_param**2 - c4*a_param**3)
            r0_cubed = r0_cubed + delta_r0_cubed
            if(dabs(delta_r0_cubed/r0_cubed).le.star%job%acfpft)exit
         end do
         r0(1) = pow(r0_cubed, cc13)
         if (zone_end.eq.1) return
         rho_bar_prev = rho_bar
         r_phi_prev = r_phi
         rho_prev = density
         j_begin = 2
      else
         rho_bar_prev = exp(ln10*(log_mass(zone_start-1) - c4pi3l - &
              3.0d0*log_radius(zone_start-1)))
         r_phi_prev = exp(ln10*log_radius(zone_start-1))
         rho_prev = exp(ln10*log_density(zone_start-1))
      end if
!  CALCULATE ETA2 AND R0 FOR REMAINING POINTS
      do i = j_begin,zone_end
         rho_bar = exp(ln10*(log_mass(i) - c4pi3l - 3.0d0*log_radius(i)))
         density = exp(ln10*log_density(i))
         r_phi = exp(ln10*log_radius(i))
         gm = cg*exp(ln10*log_mass(i))
         dr = r_phi - r_phi_prev
! FIND ETA2 USING 4-POINT RUNGE-KUTTE TECHNIQUE
! D(ETA2)/D(R0) IS COMPUTED USING RADAU'S EQUATION:
! R*(DETA2/DR)+6*RHO*(ETA2+1)/RHOBAR+ETA2*(ETA2-1) = 6 ,AND ETA2(0) = 0
! FOR A FIRST GUESS,R0 = RPHI IS ASSUMED.
! RHOA AND RHOBA ARE AVERAGES OF THE RHO,RHOBAR OF OLD AND NEW SHELLS
         rho_avg = 0.5d0*(density+rho_prev)
         rho_bar_avg = 0.5d0*(rho_bar+rho_bar_prev)
         eta2(i) = radau_rk4_step(dr, eta2(i-1), rho_prev, rho_bar_prev, &
              r_phi_prev, rho_avg, rho_bar_avg, density, rho_bar, r_phi)
         r_phi_cubed = r_phi**3
         r0_cubed = r_phi_cubed
         r0_estimate = r_phi
         acc_tol = pow(star%job%acfpft, cc13)
! ITERATE BETWEEN SOLUTION FOR ETA2 AND SOLUTION FOR R0 ITFP1 TIMES.
         do k = 1,star%ctrl%itfp1
            fact = 5.0d0*cc13*omega(i)**2/(gm*(2.0d0+eta2(i)))
! NOW ITERATE FOR R0 GIVEN RPHI AND ETA2, USING THE RELATION
! RPHI**3 = R0**3(1.0 + 3/5A**2 - 2/35A**3)
! WHERE A = OMEGA**2*R0**3*5/3GM(2+ETA2))
! (Same Newton step as the center loop above, but exits on .lt. acfpft
! where that one uses .le. -- kept as in the original.)
            do j = 1,star%ctrl%itfp2
               a_param = fact*r0_cubed
               delta_r0_cubed = (r_phi_cubed-r0_cubed*(1.0d0 + c1*a_param**2 - c2*a_param**3))/ &
               (1.0d0 + c3*a_param**2 - c4*a_param**3)
               r0_cubed = r0_cubed + delta_r0_cubed
               if(dabs(delta_r0_cubed/r0_cubed).lt.star%job%acfpft)exit
            end do
            r0(i) = pow(r0_cubed, cc13)
            err = r0(i) - r0_estimate
!  ETA2 IS A FUNCTION OF R0, AND R0=RPHI WAS USED TO CALCULATE ETA2
!  CORRECT ETA2 HERE IF DR/R0 > CUBE ROOT OF ACFPFT
            if(dabs(err)/r0(i).le.acc_tol) exit
! FIND ETA2 USING 4-POINT RUNGE-KUTTE TECHNIQUE AGAIN,BUT FINDING
! ETA2 AT R0(I) RATHER THAN ASSUMING R0 = RPHI.
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

! One 4th-order Runge-Kutta step of Radau's equation
!   R*(DETA2/DR) + 6*RHO*(ETA2+1)/RHOBAR + ETA2*(ETA2-1) = 6
! from (r_prev, eta_prev) to r = r_prev + dr, with rho/rho_bar given
! at the start, midpoint (rho_avg/rho_bar_avg) and end of the step.
! Body is token-identical to the two former inline copies.
      pure function radau_rk4_step(dr, eta_prev, rho_prev, rho_bar_prev, &
           r_prev, rho_avg, rho_bar_avg, rho, rho_bar, r) result(eta_new)
      double precision, intent(in) :: dr, eta_prev, rho_prev, rho_bar_prev, &
           r_prev, rho_avg, rho_bar_avg, rho, rho_bar, r
      double precision :: eta_new
      double precision :: deta1, deta2, deta3, deta4, eta_avg, r_avg
!    FIRST EVALUATE DETA/DR AT THE BEGINNING OF THE INTERVAL
      deta1 = dr*(6.0d0 - 6.0d0*rho_prev*(eta_prev + 1.0d0)/rho_bar_prev &
      - eta_prev*(eta_prev - 1.0d0))/r_prev
      eta_avg = eta_prev + 0.5d0*deta1
      r_avg = r_prev + 0.5d0*dr
!    USING THE ESTIMATED ETA2 AT THE MIDPOINT,FIND DETA/DR AT THE
!    MIDPOINT.
      deta2 = dr*(6.0d0 - 6.0d0*rho_avg*(eta_avg + 1.0d0)/rho_bar_avg &
      - eta_avg*(eta_avg - 1.0d0))/r_avg
      eta_avg = eta_prev + 0.5d0*deta2
!    USING THIS REFINED DERIVATIVE TO ESTIMATE ETA2 AT MIDPOINT,
!    FIND DETA/DR AT THE MIDPOINT AGAIN.
      deta3 = dr*(6.0d0 - 6.0d0*rho_avg*(eta_avg + 1.0d0)/rho_bar_avg &
      - eta_avg*(eta_avg - 1.0d0))/r_avg
      eta_avg = eta_prev + deta3
!    USING DETA/DR AT THE MIDPOINT TO ESTIMATE ETA2 AT THE END OF
!    THE INTERVAL,GET DETA/DR AT THE END OF THE INTERVAL.
      deta4 = dr*(6.0d0 - 6.0d0*rho*(eta_avg + 1.0d0)/rho_bar &
      - eta_avg*(eta_avg - 1.0d0))/r
!    PERFORM 4TH ORDER RUNGE-KUTTE INTEGRATION USING ABOVE 4 DERIVS.
      eta_new = eta_prev+cc13*(0.5d0*deta1+deta2+deta3+0.5d0*deta4)
      end function radau_rk4_step

end subroutine shape

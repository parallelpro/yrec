!----------------------------------------------------------------------
! instability_velocities
!----------------------------------------------------------------------
! The velocities of the rotational instabilities at the interfaces of a
! transport region (interface i lies between zones i-1 and i): the
! Eddington-Sweet circulation and its inhibition by the mu gradient, the
! GSF instability, and the secular and dynamical shear instabilities.
! Star-blind kernels (tools/check_boundaries.py STAR_BLIND_FILES): every
! input is an argument, from circulation_velocities, which interpolates
! the structure to the interfaces, and rotation_stability_setup, which
! forms the structural coefficients. test_rotation_kernels checks them
! against the laws they keep.
module instability_velocities
      use rotation_configuration, only: gradient_margin_floor
      implicit none
      private
      public :: compute_eddington_sweet_velocity, inhibit_by_mu_gradient, &
           compute_gsf_velocity, compute_shear_velocities

! The GSF velocity's mu-gradient inhibition is ignored where its limit
! is below gsf_mu_negligible_fraction of |d omega/d ln r| or below
! negligible_limit; the shear's mu gradient counts as absent below a
! fractional Delta mu of uniform_mu_fraction, and its Richardson term
! is floored at negligible_limit. (gradient_margin_floor, the margins'
! floor, is rotation_configuration's.)
      double precision, parameter :: gsf_mu_negligible_fraction = 1.0d-2
      double precision, parameter :: negligible_limit = 1.0d-20
      double precision, parameter :: uniform_mu_fraction = 1.0d-10
contains

! ---------------------------------------------------------------
! The Eddington-Sweet circulation velocity U_ES (cm/s) at the interfaces
! zone_min..zone_max, Zahn (1991) estimate, with the coefficients
! es_velocity_coeff1, es_velocity_coeff2 and es_shear_coeff from
! rotation_stability_setup. (The older KM1974 meridional circulation
! estimate is no longer offered; MHP 8/03.) Interfaces with convective
! zones on both sides are skipped. radius (cm), omega (rad/s),
! local_gravity and quadrupole_moment (compute_quadrupole) are at the
! zone centres, omega_interface and radius_interface at the interfaces.
! es_circulation_velocity and circulation_correction_ratio are set at
! the interfaces it treats; a negative ratio is reported.
subroutine compute_eddington_sweet_velocity(zone_min, zone_max, am_transport_convective_flag, &
     radius, omega, omega_interface, radius_interface, local_gravity, quadrupole_moment, &
     es_velocity_coeff1, es_velocity_coeff2, es_shear_coeff, &
     circulation_correction_ratio, es_circulation_velocity)
      use physical_constants, only: cc13, cc23
      implicit none
      integer, intent(in) :: zone_min, zone_max
      logical, intent(in) :: am_transport_convective_flag(:)
      double precision, intent(in) :: radius(:), omega(:), omega_interface(:), &
           radius_interface(:), local_gravity(:), quadrupole_moment(:), &
           es_velocity_coeff1(:), es_velocity_coeff2(:), es_shear_coeff(:)
      double precision, intent(inout) :: circulation_correction_ratio(:), &
           es_circulation_velocity(:)
      integer :: i
      double precision :: dr, gmid, phistd, phis2, phi2, qw
      do i = zone_min,zone_max
         if (am_transport_convective_flag(i).and. &
              am_transport_convective_flag(i-1)) cycle
! The uncorrected estimate, with dg/g = omega^2 r^3/(G m):
!   U_ES = omega^2 (es_velocity_coeff1 + omega^2 es_velocity_coeff2).
         es_circulation_velocity(i) = omega_interface(i)**2* &
              (es_velocity_coeff1(i)+omega_interface(i)**2*es_velocity_coeff2(i))
! The quadrupole correction, with dg/g as in Zahn (1992). The variation
! of gravity on a level surface is d/dr[(omega^2 r^2/3 - Phi)/g], with
! g the local mean gravity (~G m/r^2) and Phi the quadrupole moment
! (compute_quadrupole; it already contains omega^2). The classical
! expression neglects Phi and dm/dr:
!   phistd = (4/3) omega^2 r/g.
! U_ES is multiplied by the ratio of the full expression to it,
!   circulation_correction_ratio = (phis2 + phi2)/phistd,
!   phis2 = (omega^2/3) d(r^2/g)/dr,  phi2 = -d(Phi/g)/dr,
! and qw = omega d omega/dr. A negative ratio is reported.
         dr = radius(i) - radius(i-1)
         gmid = 0.5d0*(local_gravity(i)+local_gravity(i-1))
         phistd = 2.0d0*cc23*omega_interface(i)**2*radius_interface(i)/gmid
         phis2 = cc13*omega_interface(i)**2*(radius(i)**2/local_gravity(i)- &
                 radius(i-1)**2/local_gravity(i-1))/dr
         phi2 = (quadrupole_moment(i-1)/local_gravity(i-1)- &
              quadrupole_moment(i)/local_gravity(i))/dr
         qw = omega_interface(i)*(omega(i)-omega(i-1))/dr
         circulation_correction_ratio(i) = (phis2+phi2)/phistd
         if (circulation_correction_ratio(i).lt.0.0d0) then
            write(*,303) i,circulation_correction_ratio(i),phistd,phis2,phi2, &
            quadrupole_moment(i),quadrupole_moment(i-1),local_gravity(i), &
            local_gravity(i-1),radius(i), &
            radius(i-1)
 303        format(i5,' RAT ',1pe10.3,' POT SPH,CYL,QUA', &
                 3e10.3/' QUA,G,R ',6e12.3)
         end if
! U_ES = |circulation_correction_ratio U_ES - es_shear_coeff qw|.
         es_circulation_velocity(i) = &
              abs(circulation_correction_ratio(i)*es_circulation_velocity(i)- &
              es_shear_coeff(i)*qw)
! (commented out: the advection-diffusion branch)
! The advection-diffusion transport (Maeder & Zahn 1998; commented out,
! see secular_transport). theta_new is the d theta/dt term (MHP 6/02),
!   theta_new = es_relaxation_factor (theta_mean qw - theta_prev)/timestep.
! With ratio = circulation_correction_ratio, the advective part is
!   es_advective_velocity = (1/5) (ratio U_ES + theta_new)
! and the diffusive part, a diffusion coefficient (cm^2/s) for
! d omega/dr,
!   es_diffusive_velocity = (1/5) [(es_shear_coeff + difad_shear_coeff1)
!                                  omega^2 + difad_shear_coeff2 qw],
! both with the factor 1/5 (MHP 5/02). The total used for the transport
! regions is
!   U_ES = ratio U_ES + theta_new + (es_shear_coeff + difad_shear_coeff1) qw
!          + difad_shear_coeff2 (d omega/dr)^2,
! where the es_shear_coeff term has the opposite sign to the branch
! above. second_deriv_geom_factor_eqgrid is second_deriv_geom_factor
! times omega.
! Reviving it means passing these rot_scr arrays and the timestep in as
! arguments.
!            rot_scr%theta_new(i) = rot_scr%es_relaxation_factor(i)*(rot_scr%theta_mean(i)*qw- &
!                 rot_scr%theta_prev(i))/timestep
!            rot_scr%es_advective_velocity(i) = 0.2d0* &
!                 (rot_scr%circulation_correction_ratio(i)*es_circulation_velocity(i)+ &
!                 rot_scr%theta_new(i))
!            q1 = rot_scr%difad_shear_coeff1(i)*rot_scr%omega_interface(i)**2
!            q2 = rot_scr%difad_shear_coeff2(i)*qw
!            q0 = rot_scr%es_shear_coeff(i)*rot_scr%omega_interface(i)**2
!            rot_scr%es_diffusive_velocity(i) = 0.2d0*(q0+q1+q2)
!            es_circulation_velocity(i) = &
!                 rot_scr%circulation_correction_ratio(i)*es_circulation_velocity(i)+ &
!                 rot_scr%theta_new(i)+ &
!                 (rot_scr%es_shear_coeff(i)+rot_scr%difad_shear_coeff1(i))*qw+ &
!                 rot_scr%difad_shear_coeff2(i)*((omega(i)-omega(i-1))/dr)**2
!            rot_scr%second_deriv_geom_factor_eqgrid(i) = &
!                 rot_scr%second_deriv_geom_factor(i)*rot_scr%omega_interface(i)
      end do
end subroutine compute_eddington_sweet_velocity

! ---------------------------------------------------------------
! The inhibition of the circulation by the mu gradient (Maeder & Zahn
! 1998): Delta = grad_ad - grad in the denominator of U_ES becomes
! Delta + grad_mu, with grad_mu = d ln mu/d ln P = qmu/qp, so the
! velocities are multiplied by
!   qqq = max(del_grad_diff_interface, 1e-3)
!         / max(del_grad_diff_interface + grad_mu, 1e-3).
! del_grad_diff_interface already contains grad_mu
! (rotation_stability_setup), so grad_mu enters this factor twice and
! the mu gradient can cut U_ES by at most a half. This is the F77
! behaviour (setupv.f, vcirc.f), kept by the author's decision
! (2026-10-12; docs/HOW_YREC_WORKS.md section 4.6).
! mu and log_pressure (log10 P) are at the zone centres;
! es_circulation_velocity is scaled in place at the interfaces
! zone_min..zone_max.
pure subroutine inhibit_by_mu_gradient(zone_min, zone_max, mu, log_pressure, &
     del_grad_diff_interface, es_circulation_velocity)
      use physical_constants, only: ln10
      use math_lib
      implicit none
      integer, intent(in) :: zone_min, zone_max
      double precision, intent(in) :: mu(:), log_pressure(:), &
           del_grad_diff_interface(:)
      double precision, intent(inout) :: es_circulation_velocity(:)
      integer :: i
      double precision :: qmu, qp, ddel, ddtest, ddtest2, qqq
       do i = zone_min,zone_max
          qmu = log(mu(i))-log(mu(i-1))
          qp = ln10*(log_pressure(i)-log_pressure(i-1))
          ddel=del_grad_diff_interface(i)+ qmu/qp
          ddtest = max(del_grad_diff_interface(i),gradient_margin_floor)
          ddtest2 = max(ddel,gradient_margin_floor)
          qqq = ddtest/ddtest2
          es_circulation_velocity(i) = es_circulation_velocity(i)*qqq
! (commented out: the factor on the advection-diffusion transport's
! velocities and its coefficients that contain 1/Delta, set once per
! sub-step and so multiplied by qqq again at each iteration)
! Reviving it means passing these arrays in as arguments.
!          if (rotation_settings%use_diffusion_advection_transport) then
!             rot_scr%es_advective_velocity(i) = rot_scr%es_advective_velocity(i)*qqq
!             rot_scr%es_diffusive_velocity(i) = rot_scr%es_diffusive_velocity(i)*qqq
!          end if
!          rot_scr%es_relaxation_factor(i) = rot_scr%es_relaxation_factor(i)*qqq
!          rot_scr%velocity_coeff0(i) = rot_scr%velocity_coeff0(i)*qqq
       end do
end subroutine inhibit_by_mu_gradient

! ---------------------------------------------------------------
! The GSF instability velocity U_GSF (cm/s) at the interfaces
! zone_min..zone_max, Kippenhahn (1980) estimate, with the coefficient
! gsf_kippenhahn_coeff from rotation_stability_setup.
! gsf_inhibition_mode (legacy IGSF) sets what may inhibit it:
!   0: viscosity and mu gradients;
!   1: neither;
!   2: viscosity only.
! The related ABCD instability would not be inhibited by viscosity or mu
! gradients, so including it is equivalent to mode 1 (the original
! Kippenhahn 1980 estimate). The inhibition factors come from later work
! by Spruit.
! Interfaces with convective zones on both sides are skipped. radius
! (cm), omega (rad/s) and mu are at the zone centres; domega_dlnr is
! d omega/d ln r and dynamical_shear_omega_limit its dynamical-shear limit; the
! mean mu, omega, kinematic viscosity nu, thermal diffusivity K and
! gravity factor are at the interfaces. gsf_circulation_velocity is set
! at the interfaces it treats.
pure subroutine compute_gsf_velocity(zone_min, zone_max, gsf_inhibition_mode, &
     am_transport_convective_flag, radius, omega, domega_dlnr, dynamical_shear_omega_limit, &
     mu, mean_molecular_weight_interface, omega_interface, &
     kinematic_viscosity_interface, thermal_diffusivity_interface, &
     gravity_factor_interface, gsf_kippenhahn_coeff, gsf_circulation_velocity)
      use math_lib
      implicit none
      integer, intent(in) :: zone_min, zone_max, gsf_inhibition_mode
      logical, intent(in) :: am_transport_convective_flag(:)
      double precision, intent(in) :: radius(:), omega(:), domega_dlnr(:), &
           dynamical_shear_omega_limit(:), mu(:), mean_molecular_weight_interface(:), &
           omega_interface(:), kinematic_viscosity_interface(:), &
           thermal_diffusivity_interface(:), gravity_factor_interface(:), &
           gsf_kippenhahn_coeff(:)
      double precision, intent(inout) :: gsf_circulation_velocity(:)
      integer :: i
      double precision :: dr, qwrmx, fxx, qwrmx2, rmid, dlnwdr, dlnjmdr, fx
      do i = zone_min,zone_max
         if (am_transport_convective_flag(i).and. &
              am_transport_convective_flag(i-1)) cycle
! Viscous inhibition (modes 0 and 2; MHP 8/93, mu gradients neglected):
! stable while |d omega/d ln r| < q_nu = 2 sqrt(nu/K) dynamical_shear_omega_limit,
! with nu the kinematic viscosity and K the thermal diffusivity. Above
! it the velocity is multiplied by fxx = sqrt((|d omega/d ln r| - q_nu)
! /q_nu); without viscous inhibition fxx = 1.
         if (gsf_inhibition_mode.eq.2 .or. gsf_inhibition_mode.eq.0) then
            qwrmx = 2.0d0*sqrt(kinematic_viscosity_interface(i)/ &
                 thermal_diffusivity_interface(i))*dynamical_shear_omega_limit(i)
            if (abs(domega_dlnr(i)).lt.qwrmx) then
               gsf_circulation_velocity(i) = 0.0d0
               cycle
            else
              fxx = sqrt((abs(domega_dlnr(i))-qwrmx)/qwrmx)
            end if
         else
            fxx = 1.0d0
         end if
         dr = radius(i)-radius(i-1)
! mu-gradient inhibition (mode 0): stable while |d omega/d ln r| < q_mu
! = 2 sqrt(g |d mu/dr|/mu). Above it the factor is
! sqrt((|d omega/d ln r| - q_mu)/q_mu), or 1 when q_mu is below 1% of
! |d omega/d ln r| or below 1e-20 (which avoids dividing by zero;
! MHP 5/02), and fxx becomes the smaller of the two factors.
         if (gsf_inhibition_mode.eq.0) then
            qwrmx=2.0d0*sqrt(gravity_factor_interface(i)* &
                 abs(mu(i)-mu(i-1)) &
                 /dr/mean_molecular_weight_interface(i))
            if (abs(domega_dlnr(i)).lt.qwrmx) then
               gsf_circulation_velocity(i) = 0.0d0
               cycle
            else
               qwrmx2 = gsf_mu_negligible_fraction*domega_dlnr(i)
               if (abs(qwrmx).lt.abs(qwrmx2)) then
                  qwrmx2 = 1.0d0
               else
                  if (abs(qwrmx).lt.negligible_limit) then
                     qwrmx2 = 1.0d0
                  else
                     qwrmx2 = sqrt((abs(domega_dlnr(i))-qwrmx)/qwrmx)
                  end if
               end if
               fxx = min(fxx,qwrmx2)
            end if
         end if
         rmid = 0.5d0*(radius(i)+radius(i-1))
         dlnwdr = abs(log(omega(i))-log(omega(i-1)))/dr
! The GSF instability is triggered by d omega/dz /= 0 (rotation not on
! cylinders, true in general in these models) or by d ln j/dr < 0, with
! j = omega r^2 the specific angular momentum (legacy J/M), which
! sometimes occurs. Where j decreases outward both criteria are used
! and the larger wins; otherwise only the d omega/dz one:
!   U_GSF = gsf_kippenhahn_coeff K omega^2 max(2 |d ln j/dr|, dlnwdr/4)/r
! or
!   U_GSF = gsf_kippenhahn_coeff K omega^2 (dlnwdr/4)/r,
! with dlnwdr = |d ln omega/dr|, then multiplied by fxx. The 1/r is
! there because every circulation velocity is later multiplied by r (the
! length scale).
         if (omega(i)*radius(i)**2 .lt. omega(i-1)*radius(i-1)**2) then
            dlnjmdr = abs(2.0d0/rmid+(log(omega(i))- &
                      log(omega(i-1)))/dr)
            fx = max(2.0d0*dlnjmdr,0.25d0*dlnwdr)
            gsf_circulation_velocity(i) = gsf_kippenhahn_coeff(i)*fx* &
                 thermal_diffusivity_interface(i)*omega_interface(i)**2/rmid
         else
            gsf_circulation_velocity(i)=0.25d0*gsf_kippenhahn_coeff(i)* &
                 thermal_diffusivity_interface(i)*dlnwdr* &
                 omega_interface(i)**2/rmid
         end if
         gsf_circulation_velocity(i) = abs(fxx*gsf_circulation_velocity(i))
      end do
end subroutine compute_gsf_velocity

! ---------------------------------------------------------------
! The secular (diffusive) and dynamical shear instabilities (Endal &
! Sofia 1978, ApJ 220, 279, paper II), giving secular_shear_velocity
! U_SS (cm/s) at the interfaces zone_min..zone_max. Where the shear is
! dynamical it also corrects gsf_circulation_velocity
! (compute_gsf_velocity's) for the shear the instability leaves behind,
! and marks the interface (dynamical, with the inner omega at the limit
! in dynamical_inner_omega) for the caller's report. radius (cm), omega
! (rad/s) and mu are at the zone centres; critical_reynolds is the
! critical Reynolds number of the secular shear; the other inputs as in
! compute_gsf_velocity, with mu_gradient_richardson_coeff the
! mu-gradient Richardson term (rotation_stability_setup).
pure subroutine compute_shear_velocities(zone_min, zone_max, critical_reynolds, &
     radius, omega, domega_dlnr, dynamical_shear_omega_limit, &
     mu, mean_molecular_weight_interface, radius_interface, &
     kinematic_viscosity_interface, thermal_diffusivity_interface, &
     mu_gradient_richardson_coeff, secular_shear_velocity, gsf_circulation_velocity, &
     dynamical, dynamical_inner_omega)
      use math_lib
      implicit none
      integer, intent(in) :: zone_min, zone_max
      double precision, intent(in) :: critical_reynolds
      double precision, intent(in) :: radius(:), omega(:), &
           domega_dlnr(:), dynamical_shear_omega_limit(:), mu(:), &
           mean_molecular_weight_interface(:), radius_interface(:), &
           kinematic_viscosity_interface(:), thermal_diffusivity_interface(:), &
           mu_gradient_richardson_coeff(:)
      double precision, intent(inout) :: secular_shear_velocity(:), &
           gsf_circulation_velocity(:)
! the interfaces where the shear is dynamical, and the inner omega at
! the limit there (wmin), for the caller's report
      logical, intent(inout) :: dynamical(:)
      double precision, intent(inout) :: dynamical_inner_omega(:)
      integer :: i
      double precision :: qwr, dr, rmid, wmin, qwrmx, qwrmx2
      do i = zone_min,zone_max
         dynamical(i) = .false.
! Dynamical shear operates where |d omega/d ln r| > dynamical_shear_omega_limit
! (Ri < 1/4). The shear velocity is then the large dynamical value,
!   U_SS = (8/45) K (|d omega/d ln r|/dynamical_shear_omega_limit)^2/r.
         if (abs(domega_dlnr(i)).gt.dynamical_shear_omega_limit(i)) then
            qwr = abs(domega_dlnr(i))
            secular_shear_velocity(i)=8.0d0/4.5d1* &
                 thermal_diffusivity_interface(i)* &
                 (qwr/dynamical_shear_omega_limit(i))**2/radius_interface(i)
! Correct the GSF velocity for the shear the dynamical instability
! leaves behind: wmin is the omega of zone i-1 at which the trigger's
! own (centred) d omega/d ln r equals the limit,
!   wmin = omega(i) - (omega(i) - omega(i-1)) limit/|d omega/d ln r|,
! which lies strictly between omega(i) and omega(i-1) whenever the
! trigger fires (docs/bug_sweep.md O11), and U_GSF is scaled by its
! gradient term (Kippenhahn 1980, as compute_gsf_velocity forms it for
! each profile) evaluated with wmin over the same term evaluated with
! omega(i-1). Each such interface is marked.
            dr = radius(i)-radius(i-1)
            rmid = 0.5d0*(radius(i)+radius(i-1))
            wmin = omega(i) - (omega(i) - omega(i-1))* &
                 dynamical_shear_omega_limit(i)/qwr
            gsf_circulation_velocity(i) = gsf_circulation_velocity(i)* &
                 gsf_gradient_term(omega(i), wmin, radius(i), radius(i-1), dr, rmid)/ &
                 gsf_gradient_term(omega(i), omega(i-1), radius(i), radius(i-1), dr, rmid)
            dynamical(i) = .true.
            dynamical_inner_omega(i) = wmin
            cycle
         end if
! The largest |d omega/d ln r| stable against secular shear.
! dynamical_shear_omega_limit (legacy QWRMAX) is the dynamical-shear limit. With
! thermal diffusion the secular shear Richardson number is
! (Pr Re_c/8) Ri, with Pr = nu/K the Prandtl number and Re_c the critical
! Reynolds number (critical_reynolds), and the interface is unstable
! where it falls below 1/4, i.e. where |d omega/d ln r| exceeds
!   qwrmx = sqrt(Pr Re_c/8) dynamical_shear_omega_limit
! (the square root of Pr, not Pr; MHP 3/92). The velocity estimate is
! from Zahn (1991):
!   U_SS = (8/45) K ((|d omega/d ln r| - qwrmx)/dynamical_shear_omega_limit)^2/r.
         qwrmx = sqrt(kinematic_viscosity_interface(i)/ &
              thermal_diffusivity_interface(i)*1.25d-1*critical_reynolds)* &
              dynamical_shear_omega_limit(i)
         if (abs(domega_dlnr(i)).gt.qwrmx) then
! Unstable. Without a mu gradient (|Delta mu/mu| < 1e-10), U_SS as above.
            if (abs((mu(i)-mu(i-1))/mean_molecular_weight_interface(i)) &
                 .lt.uniform_mu_fraction) then
               qwrmx2 = 0.0d0
               qwr = abs(domega_dlnr(i)) - qwrmx
               secular_shear_velocity(i)=8.0d0/4.5d1* &
                    thermal_diffusivity_interface(i)* &
                    (qwr/dynamical_shear_omega_limit(i))**2/radius_interface(i)
            else
! With a mu gradient, its Richardson number
!   Ri_mu = (rho/P) (-d ln T/d ln mu) (d ln mu/d ln P) g^2
!           /(d omega/d ln r)^2,
! with -d ln T/d ln mu = (b - 1)/(1 + 3b), b = a T^4/(3P), from
! (1 - b) P = (R/mu) rho T. It reaches the dynamical criterion's 1/4
! (no Prandtl-number factor) at
!   qwrmx2 = 2 sqrt(mu_gradient_richardson_coeff |Delta mu/mu|).
               qwrmx2 = 2.0d0*sqrt(max(negligible_limit,mu_gradient_richardson_coeff(i)* &
                        abs((mu(i)-mu(i-1))/ &
                        mean_molecular_weight_interface(i))))
               if (abs(domega_dlnr(i)).gt.qwrmx2) then
! Unstable against both: U_SS from the excess over the larger limit (the
! mu limit where the mu gradient slows the transport of j). Stable
! against the mu gradient: U_SS stays zero.
                  qwrmx = max(qwrmx2,qwrmx)
                  qwr = abs(domega_dlnr(i)) - qwrmx
            secular_shear_velocity(i)=8.0d0/4.5d1* &
                 thermal_diffusivity_interface(i)* &
                 (qwr/dynamical_shear_omega_limit(i))**2/radius_interface(i)
               end if
            end if
         end if
      end do
end subroutine compute_shear_velocities

! ---------------------------------------------------------------
! The GSF velocity's gradient term (Kippenhahn 1980) across an interface
! with omega_top above and omega_bottom below, radii r_top and r_bottom,
! dr = r_top - r_bottom and rmid their mean: where the specific angular
! momentum omega r^2 decreases outward, max(2 |d ln j/dr|, |d ln
! omega/dr|/4); otherwise |d ln omega/dr|/4 (compute_gsf_velocity's two
! forms).
pure double precision function gsf_gradient_term(omega_top, omega_bottom, r_top, &
     r_bottom, dr, rmid)
      use math_lib
      double precision, intent(in) :: omega_top, omega_bottom, r_top, r_bottom, dr, rmid
      double precision :: dlnwdr, dlnjmdr
      dlnwdr = abs(log(omega_top)-log(omega_bottom))/dr
      if (omega_top*r_top**2 .lt. omega_bottom*r_bottom**2) then
         dlnjmdr = abs(2.0d0/rmid+(log(omega_top)-log(omega_bottom))/dr)
         gsf_gradient_term = max(2.0d0*dlnjmdr, 0.25d0*dlnwdr)
      else
         gsf_gradient_term = 0.25d0*dlnwdr
      end if
end function gsf_gradient_term

end module instability_velocities

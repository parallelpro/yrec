!----------------------------------------------------------------------
! circulation
!----------------------------------------------------------------------
! The velocities that drive rotational mixing: the circulation and
! instability velocities (circulation_velocities), the diffusion
! coefficients formed from them (diffusion_velocity_scales), the Zahn
! (1992) coupling factor (zahn_coupling_factor) and the quadrupole
! moment of the distorted star (compute_quadrupole).
module circulation
      use rotation_configuration, only: rotation_settings
      implicit none
      private
      public :: circulation_velocities, diffusion_velocity_scales, zahn_coupling_factor, compute_quadrupole

contains

! ---------------------------------------------------------------
! The circulation and instability velocities at each interface
! zone_min..zone_max of the transport region (interface i lies between
! zones i-1 and i). Each stability criterion is a structural coefficient
! (rotation_stability_setup, stability_coefficients) times a gradient
! across the interface, compared with a critical value; where the
! interface is unstable the mechanism gets a velocity:
!  - star%es_circulation_velocity, the Eddington-Sweet circulation (Zahn
!    1991, 1992), corrected for the quadrupole moment and reduced by the
!    mu gradient (Maeder & Zahn 1998);
!  - star%gsf_circulation_velocity, the GSF instability (Kippenhahn
!    1980), optionally inhibited by viscosity and mu gradients;
!  - star%secular_shear_velocity, the secular shear instability, or the
!    dynamical one where Ri < 1/4 (Endal & Sofia 1978).
! The Eddington-Sweet and GSF velocities are skipped at interfaces with
! convective zones on both sides. From the second iteration on, each
! velocity is averaged with the previous iteration's. Their sum is
! total_circulation_velocity (legacy HV); any_transport_active is set if
! it is positive anywhere. With no_am_transport_in_core the velocities
! are zeroed inside the energy-generating core (zero_transport_in_core).
! dlnomega_dlnr is d omega/d ln r (not
! d ln omega/d ln r) and dlnomega_dlnr_max its dynamical-shear limit.
! Each velocity is later multiplied by r (diffusion_velocity_scales) to
! give a diffusion coefficient.
subroutine circulation_velocities(log_radius, radius, zone_min, zone_max, iteration, &
     am_transport_convective_flag, num_zones, omega, any_transport_active, &
     dlnomega_dlnr, dlnomega_dlnr_max, total_circulation_velocity, &
     total_luminosity, log_pressure)
      use rotation_scratch_lib

      use star_info_lib, only: star, json
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: log_radius(json), radius(json)
      integer, intent(in) :: zone_min, zone_max, iteration
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(in) :: omega(json)
      logical, intent(out) :: any_transport_active
      double precision, intent(in) :: dlnomega_dlnr(json), &
           dlnomega_dlnr_max(json)
      double precision, intent(out) :: total_circulation_velocity(json)
      double precision, intent(in) :: total_luminosity, log_pressure(json)
! --- locals ---
      integer :: i, i0, i1
      double precision :: dr, gmid, phistd, phis2, phi2, qw, &
           qmu, qp, ddel, ddtest, ddtest2, qqq, &
           qwrmx, fxx, qwrmx2, qwr, rmid, dlnwdr, dlnjmdr, fx, qwrmx_dyn, &
           wmin, dlnwdr0, dlnjmdr0, dell

      any_transport_active = .false.
      if (iteration.gt.1) call store_previous_velocities()
      do i = 1,num_zones
         total_circulation_velocity(i) = 0.0d0
         star%es_circulation_velocity(i) = 0.0d0
         star%gsf_circulation_velocity(i) = 0.0d0
         star%secular_shear_velocity(i) = 0.0d0
      end do
      call interpolate_to_interfaces()
      call eddington_sweet_velocity()
      call inhibit_by_mu_gradient()
      call gsf_velocity()
      call shear_velocities()
      call apply_velocity_scales()
      if (iteration.gt.1) call average_with_previous()
      call sum_velocities()
      if (rotation_settings%no_am_transport_in_core) call zero_transport_in_core()

contains

! ---------------------------------------------------------------
! Keep the previous iteration's velocities for average_with_previous.
subroutine store_previous_velocities()
      do i = zone_min,zone_max
         circ_scr%gsf_circulation_velocity_prev(i) = star%gsf_circulation_velocity(i)
         circ_scr%es_circulation_velocity_prev(i) = star%es_circulation_velocity(i)
         circ_scr%secular_shear_velocity_prev(i) = star%secular_shear_velocity(i)
      end do
! (commented out: the advection-diffusion transport's velocities)
!      if (rotation_settings%use_diffusion_advection_transport) then
!         do i = zone_min,zone_max
!            rot_scr%es_advective_velocity_prev(i) = rot_scr%es_advective_velocity(i)
!            rot_scr%es_diffusive_velocity_prev(i) = rot_scr%es_diffusive_velocity(i)
!         end do
!      end if
end subroutine store_previous_velocities

! ---------------------------------------------------------------
! The mean molecular weight (and what depends on it) and omega change
! during a diffusion sub-step, so their values at the interfaces, where
! stability is evaluated, are found again each iteration. The kinematic
! viscosity and thermal diffusivity use 4-point Lagrange interpolation
! in ln (log_lagrange4) from zones i-2..i+1.
subroutine interpolate_to_interfaces()
! The first interface (2) uses zones 1-4.
      if (zone_min.lt.3) then
         rot_scr%kinematic_viscosity_interface(2)=log_lagrange4(mid_scr%visc, 1, 2)
         rot_scr%thermal_diffusivity_interface(2)=log_lagrange4(mid_scr%thdif, 1, 2)
         i0=3
      else
         i0=zone_min
      end if
! The last interface (num_zones) uses the last four zones.
      if (zone_max.eq.num_zones) then
         rot_scr%kinematic_viscosity_interface(num_zones)= &
              log_lagrange4(mid_scr%visc, num_zones-3, num_zones)
         rot_scr%thermal_diffusivity_interface(num_zones)= &
              log_lagrange4(mid_scr%thdif, num_zones-3, num_zones)
         i1=num_zones-1
      else
         i1=zone_max
      end if
! The other interfaces.
      do i = i0,i1
         rot_scr%kinematic_viscosity_interface(i)=log_lagrange4(mid_scr%visc, i-2, i)
         rot_scr%thermal_diffusivity_interface(i)=log_lagrange4(mid_scr%thdif, i-2, i)
      end do
! mu and omega: the mean of the two neighbouring zones.
      do i = 2,num_zones
         rot_scr%mean_molecular_weight_interface(i) = 0.5d0*(mid_scr%mu(i)+mid_scr%mu(i-1))
         rot_scr%omega_interface(i) = 0.5d0*(omega(i)+omega(i-1))
      end do
end subroutine interpolate_to_interfaces

! ---------------------------------------------------------------
! The Eddington-Sweet circulation velocity, Zahn (1991) estimate, with
! the coefficients from rotation_stability_setup. (The older KM1974
! meridional circulation estimate is no longer offered; MHP 8/03.)
subroutine eddington_sweet_velocity()
      do i = zone_min,zone_max
         if (am_transport_convective_flag(i).and. &
              am_transport_convective_flag(i-1)) cycle
! The uncorrected estimate, with dg/g = omega^2 r^3/(G m):
!   U_ES = omega^2 (es_velocity_coeff1 + omega^2 es_velocity_coeff2).
         star%es_circulation_velocity(i) = rot_scr%omega_interface(i)**2* &
              (rot_scr%es_velocity_coeff1(i)+rot_scr%omega_interface(i)**2*rot_scr%es_velocity_coeff2(i))
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
         gmid = 0.5d0*(rot_scr%local_gravity(i)+rot_scr%local_gravity(i-1))
         phistd = 2.0d0*cc23*rot_scr%omega_interface(i)**2*rot_scr%radius_interface(i)/gmid
         phis2 = cc13*rot_scr%omega_interface(i)**2*(radius(i)**2/rot_scr%local_gravity(i)- &
                 radius(i-1)**2/rot_scr%local_gravity(i-1))/dr
         phi2 = (rot_scr%quadrupole_moment(i-1)/rot_scr%local_gravity(i-1)- &
              rot_scr%quadrupole_moment(i)/rot_scr%local_gravity(i))/dr
         qw = rot_scr%omega_interface(i)*(omega(i)-omega(i-1))/dr
         rot_scr%circulation_correction_ratio(i) = (phis2+phi2)/phistd
         if (rot_scr%circulation_correction_ratio(i).lt.0.0d0) then
            write(*,303) i,rot_scr%circulation_correction_ratio(i),phistd,phis2,phi2, &
            rot_scr%quadrupole_moment(i),rot_scr%quadrupole_moment(i-1),rot_scr%local_gravity(i), &
            rot_scr%local_gravity(i-1),radius(i), &
            radius(i-1)
 303        format(i5,' RAT ',1pe10.3,' POT SPH,CYL,QUA', &
                 3e10.3/' QUA,G,R ',6e12.3)
         end if
! U_ES = |circulation_correction_ratio U_ES - es_shear_coeff qw|.
         star%es_circulation_velocity(i) = &
              abs(rot_scr%circulation_correction_ratio(i)*star%es_circulation_velocity(i)- &
              rot_scr%es_shear_coeff(i)*qw)
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
!            rot_scr%theta_new(i) = rot_scr%es_relaxation_factor(i)*(rot_scr%theta_mean(i)*qw- &
!                 rot_scr%theta_prev(i))/timestep
!            rot_scr%es_advective_velocity(i) = 0.2d0* &
!                 (rot_scr%circulation_correction_ratio(i)*star%es_circulation_velocity(i)+ &
!                 rot_scr%theta_new(i))
!            q1 = rot_scr%difad_shear_coeff1(i)*rot_scr%omega_interface(i)**2
!            q2 = rot_scr%difad_shear_coeff2(i)*qw
!            q0 = rot_scr%es_shear_coeff(i)*rot_scr%omega_interface(i)**2
!            rot_scr%es_diffusive_velocity(i) = 0.2d0*(q0+q1+q2)
!            star%es_circulation_velocity(i) = &
!                 rot_scr%circulation_correction_ratio(i)*star%es_circulation_velocity(i)+ &
!                 rot_scr%theta_new(i)+ &
!                 (rot_scr%es_shear_coeff(i)+rot_scr%difad_shear_coeff1(i))*qw+ &
!                 rot_scr%difad_shear_coeff2(i)*((omega(i)-omega(i-1))/dr)**2
!            rot_scr%second_deriv_geom_factor_eqgrid(i) = &
!                 rot_scr%second_deriv_geom_factor(i)*rot_scr%omega_interface(i)
      end do
end subroutine eddington_sweet_velocity

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
subroutine inhibit_by_mu_gradient()
       do i = zone_min,zone_max
          qmu = log(mid_scr%mu(i))-log(mid_scr%mu(i-1))
          qp = ln10*(log_pressure(i)-log_pressure(i-1))
          ddel=rot_scr%del_grad_diff_interface(i)+ qmu/qp
          ddtest = max(rot_scr%del_grad_diff_interface(i),1.0d-3)
          ddtest2 = max(ddel,1.0d-3)
          qqq = ddtest/ddtest2
          star%es_circulation_velocity(i) = star%es_circulation_velocity(i)*qqq
! (commented out: the factor on the advection-diffusion transport's
! velocities and its coefficients that contain 1/Delta, set once per
! sub-step and so multiplied by qqq again at each iteration)
!          if (rotation_settings%use_diffusion_advection_transport) then
!             rot_scr%es_advective_velocity(i) = rot_scr%es_advective_velocity(i)*qqq
!             rot_scr%es_diffusive_velocity(i) = rot_scr%es_diffusive_velocity(i)*qqq
!          end if
!          rot_scr%es_relaxation_factor(i) = rot_scr%es_relaxation_factor(i)*qqq
!          rot_scr%velocity_coeff0(i) = rot_scr%velocity_coeff0(i)*qqq
       end do
end subroutine inhibit_by_mu_gradient

! ---------------------------------------------------------------
! The GSF instability velocity, Kippenhahn (1980) estimate, with the
! coefficient gsf_kippenhahn_coeff from rotation_stability_setup.
! gsf_inhibition_mode (legacy IGSF) sets what may inhibit it:
!   0: viscosity and mu gradients;
!   1: neither;
!   2: viscosity only.
! The related ABCD instability would not be inhibited by viscosity or mu
! gradients, so including it is equivalent to mode 1 (the original
! Kippenhahn 1980 estimate). The inhibition factors come from later work
! by Spruit.
subroutine gsf_velocity()
      do i = zone_min,zone_max
         if (am_transport_convective_flag(i).and. &
              am_transport_convective_flag(i-1)) cycle
! Viscous inhibition (modes 0 and 2; MHP 8/93, mu gradients neglected):
! stable while |d omega/d ln r| < q_nu = 2 sqrt(nu/K) dlnomega_dlnr_max,
! with nu the kinematic viscosity and K the thermal diffusivity. Above
! it the velocity is multiplied by fxx = sqrt((|d omega/d ln r| - q_nu)
! /q_nu); without viscous inhibition fxx = 1.
         if (rotation_settings%gsf_inhibition_mode.eq.2 .or. rotation_settings%gsf_inhibition_mode.eq.0) then
            qwrmx = 2.0d0*sqrt(rot_scr%kinematic_viscosity_interface(i)/ &
                 rot_scr%thermal_diffusivity_interface(i))*dlnomega_dlnr_max(i)
            if (abs(dlnomega_dlnr(i)).lt.qwrmx) then
               star%gsf_circulation_velocity(i) = 0.0d0
               cycle
            else
              fxx = sqrt((abs(dlnomega_dlnr(i))-qwrmx)/qwrmx)
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
         if (rotation_settings%gsf_inhibition_mode.eq.0) then
            qwrmx=2.0d0*sqrt(rot_scr%gravity_factor_interface(i)* &
                 abs(mid_scr%mu(i)-mid_scr%mu(i-1)) &
                 /dr/rot_scr%mean_molecular_weight_interface(i))
            if (abs(dlnomega_dlnr(i)).lt.qwrmx) then
               star%gsf_circulation_velocity(i) = 0.0d0
               cycle
            else
               qwrmx2 = 1.0d-2*dlnomega_dlnr(i)
               if (abs(qwrmx).lt.abs(qwrmx2)) then
                  qwrmx2 = 1.0d0
               else
                  if (abs(qwrmx).lt.1.0d-20) then
                     qwrmx2 = 1.0d0
                  else
                     qwrmx2 = sqrt((abs(dlnomega_dlnr(i))-qwrmx)/qwrmx)
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
            star%gsf_circulation_velocity(i) = rot_scr%gsf_kippenhahn_coeff(i)*fx* &
                 rot_scr%thermal_diffusivity_interface(i)*rot_scr%omega_interface(i)**2/rmid
         else
            star%gsf_circulation_velocity(i)=0.25d0*rot_scr%gsf_kippenhahn_coeff(i)* &
                 rot_scr%thermal_diffusivity_interface(i)*dlnwdr* &
                 rot_scr%omega_interface(i)**2/rmid
         end if
         star%gsf_circulation_velocity(i) = abs(fxx*star%gsf_circulation_velocity(i))
      end do
end subroutine gsf_velocity

! ---------------------------------------------------------------
! The secular (diffusive) and dynamical shear instabilities (Endal &
! Sofia 1978, ApJ 220, 279, paper II), giving star%secular_shear_velocity.
subroutine shear_velocities()
      do i = zone_min,zone_max
! Dynamical shear operates where |d omega/d ln r| > dlnomega_dlnr_max
! (Ri < 1/4). The shear velocity is then the large dynamical value,
!   U_SS = (8/45) K (|d omega/d ln r|/dlnomega_dlnr_max)^2/r.
         if (abs(dlnomega_dlnr(i)).gt.dlnomega_dlnr_max(i)) then
            qwr = abs(dlnomega_dlnr(i))
            star%secular_shear_velocity(i)=8.0d0/4.5d1* &
                 rot_scr%thermal_diffusivity_interface(i)* &
                 (qwr/dlnomega_dlnr_max(i))**2/rot_scr%radius_interface(i)
! Correct the GSF velocity for the shear the dynamical instability
! leaves behind: wmin is the omega of zone i-1 that puts
! |d omega/d ln r| exactly at the limit, and U_GSF is scaled by its
! gradient term (Kippenhahn 1980) evaluated with wmin over the same
! term evaluated with omega(i-1). Each such interface is reported.
            dr = radius(i)-radius(i-1)
            rmid = 0.5d0*(radius(i)+radius(i-1))
            if (omega(i).lt.omega(i-1)) then
               wmin = omega(i) + dlnomega_dlnr_max(i)*ln10* &
                    (log_radius(i)-log_radius(i-1))
            else
               wmin = omega(i) - dlnomega_dlnr_max(i)*ln10* &
                    (log_radius(i)-log_radius(i-1))
            end if
            dlnwdr0 = abs(log(omega(i))-log(omega(i-1)))/dr
            dlnjmdr0 = abs(2.0d0/rmid+(log(omega(i))- &
                         log(omega(i-1)))/dr)
            dlnjmdr = abs(2.0d0/rmid+(log(omega(i))-log(wmin))/dr)
            dlnwdr = abs(log(omega(i))-log(wmin))/dr
            if (omega(i)*radius(i)**2 .lt. wmin*radius(i-1)**2) then
               fx = max(2.0d0*dlnjmdr,0.25d0*dlnwdr)
               qwrmx_dyn = max(2.0d0*dlnjmdr0,0.25d0*dlnwdr0)
               star%gsf_circulation_velocity(i) = star%gsf_circulation_velocity(i)* &
                    fx/qwrmx_dyn
            else
               star%gsf_circulation_velocity(i)=star%gsf_circulation_velocity(i)* &
                    dlnwdr/dlnwdr0
            end if
            write(6,9911) i,omega(i),omega(i-1),wmin
 9911       format(1x,'DYNAMICAL SHEAR-SHELL',i5,1p,' WTOP',e11.3, &
                 ' WBOT',e11.3,' LIMIT',e11.3)
            cycle
         end if
! The largest |d omega/d ln r| stable against secular shear.
! dlnomega_dlnr_max (legacy QWRMAX) is the dynamical-shear limit. With
! thermal diffusion the secular shear Richardson number is
! (Pr Re_c/8) Ri, with Pr = nu/K the Prandtl number and Re_c the critical
! Reynolds number (critical_reynolds), and the interface is unstable
! where it falls below 1/4, i.e. where |d omega/d ln r| exceeds
!   qwrmx = sqrt(Pr Re_c/8) dlnomega_dlnr_max
! (the square root of Pr, not Pr; MHP 3/92). The velocity estimate is
! from Zahn (1991):
!   U_SS = (8/45) K ((|d omega/d ln r| - qwrmx)/dlnomega_dlnr_max)^2/r.
         qwrmx = sqrt(rot_scr%kinematic_viscosity_interface(i)/ &
              rot_scr%thermal_diffusivity_interface(i)*1.25d-1*rotation_settings%critical_reynolds)* &
              dlnomega_dlnr_max(i)
         if (abs(dlnomega_dlnr(i)).gt.qwrmx) then
! Unstable. Without a mu gradient (|Delta mu/mu| < 1e-10), U_SS as above.
            if (abs((mid_scr%mu(i)-mid_scr%mu(i-1))/rot_scr%mean_molecular_weight_interface(i)) &
                 .lt.1.0d-10) then
               qwrmx2 = 0.0d0
               qwr = abs(dlnomega_dlnr(i)) - qwrmx
               star%secular_shear_velocity(i)=8.0d0/4.5d1* &
                    rot_scr%thermal_diffusivity_interface(i)* &
                    (qwr/dlnomega_dlnr_max(i))**2/rot_scr%radius_interface(i)
            else
! With a mu gradient, its Richardson number
!   Ri_mu = (rho/P) (-d ln T/d ln mu) (d ln mu/d ln P) g^2
!           /(d omega/d ln r)^2,
! with -d ln T/d ln mu = (b - 1)/(1 + 3b), b = a T^4/(3P), from
! (1 - b) P = (R/mu) rho T. It reaches the dynamical criterion's 1/4
! (no Prandtl-number factor) at
!   qwrmx2 = 2 sqrt(mu_gradient_richardson_coeff |Delta mu/mu|).
               qwrmx2 = 2.0d0*sqrt(max(1.0d-20,rot_scr%mu_gradient_richardson_coeff(i)* &
                        abs((mid_scr%mu(i)-mid_scr%mu(i-1))/ &
                        rot_scr%mean_molecular_weight_interface(i))))
               if (abs(dlnomega_dlnr(i)).gt.qwrmx2) then
! Unstable against both: U_SS from the excess over the larger limit (the
! mu limit where the mu gradient slows the transport of j). Stable
! against the mu gradient: U_SS stays zero.
                  qwrmx = max(qwrmx2,qwrmx)
                  qwr = abs(dlnomega_dlnr(i)) - qwrmx
            star%secular_shear_velocity(i)=8.0d0/4.5d1* &
                 rot_scr%thermal_diffusivity_interface(i)* &
                 (qwr/dlnomega_dlnr_max(i))**2/rot_scr%radius_interface(i)
               end if
            end if
         end if
      end do
end subroutine shear_velocities

! ---------------------------------------------------------------
! Multiply the velocities by the user's efficiency factors
! (es_velocity_scale, gsf_velocity_scale, secular_shear_velocity_scale;
! legacy FES, FGSF, FSS). The mu-gradient inhibition is already applied.
subroutine apply_velocity_scales()
      do i = zone_min,zone_max
         star%es_circulation_velocity(i)=abs(rotation_settings%es_velocity_scale* &
              star%es_circulation_velocity(i))
         star%gsf_circulation_velocity(i) = rotation_settings%gsf_velocity_scale* &
              star%gsf_circulation_velocity(i)
         star%secular_shear_velocity(i)= rotation_settings%secular_shear_velocity_scale* &
              star%secular_shear_velocity(i)
      end do
! (commented out: the scales of the advection-diffusion transport's
! parts, and their GSF and shear diffusion coefficients r |U|)
! The same for the advection-diffusion parts (MHP 11/94).
!      if (rotation_settings%use_diffusion_advection_transport) then
!         do i = zone_min,zone_max
!               rot_scr%es_diffusive_velocity(i)=abs(rotation_settings%es_velocity_scale* &
!                    rot_scr%es_diffusive_velocity(i))
!               rot_scr%es_advective_velocity(i)=rotation_settings%es_velocity_scale* &
!                    rot_scr%es_advective_velocity(i)
!         end do
!      end if
!      if (rotation_settings%use_diffusion_advection_transport) then
! In the advection-diffusion transport (MHP 5/02) the GSF and shear
! instabilities act as diffusion. In the Endal & Sofia form a velocity
! U (e.g. U_GSF ~ K dlnwdr omega^2/r) enters d/dr(con U r d omega/dr),
! so its diffusion coefficient is r U. In the Zahn form the equation is
! d/dr(con U omega), and a velocity proportional to d omega/dr is in
! effect diffusion with coefficient U omega/(d omega/dr), which is how
! es_diffusive_velocity is formed. For the GSF and shear the code uses
! the Endal & Sofia coefficient: r |U| of each is added to
! es_diffusive_velocity and kept in shear_diffusion_coeff and
! gsf_diffusion_coeff.
!         do i = zone_min,zone_max
!            rot_scr%es_diffusive_velocity(i) = rot_scr%es_diffusive_velocity(i)+ &
!                 rot_scr%radius_interface(i)*(abs(star%gsf_circulation_velocity(i))+ &
!                 abs(star%secular_shear_velocity(i)))
!            rot_scr%shear_diffusion_coeff(i) = rot_scr%radius_interface(i)* &
!                 abs(star%secular_shear_velocity(i))
!            rot_scr%gsf_diffusion_coeff(i) = rot_scr%radius_interface(i)* &
!                 abs(star%gsf_circulation_velocity(i))
!         end do
!      end if
end subroutine apply_velocity_scales

! ---------------------------------------------------------------
! From the second iteration on, the mean of the new velocities and the
! previous iteration's (MHP 11/94).
subroutine average_with_previous()
      do i = zone_min,zone_max
         star%gsf_circulation_velocity(i) = 0.5d0*(star%gsf_circulation_velocity(i) &
              + circ_scr%gsf_circulation_velocity_prev(i))
         star%es_circulation_velocity(i) = 0.5d0*(star%es_circulation_velocity(i) &
              + circ_scr%es_circulation_velocity_prev(i))
         star%secular_shear_velocity(i) = 0.5d0*(star%secular_shear_velocity(i) &
              + circ_scr%secular_shear_velocity_prev(i))
      end do
! (commented out: the same for the advection-diffusion parts)
!      if (rotation_settings%use_diffusion_advection_transport) then
!         do i = zone_min,zone_max
!            rot_scr%es_advective_velocity(i) = 0.5d0*(rot_scr%es_advective_velocity(i) &
!                 + rot_scr%es_advective_velocity_prev(i))
!            rot_scr%es_diffusive_velocity(i) = 0.5d0*(rot_scr%es_diffusive_velocity(i) &
!                 + rot_scr%es_diffusive_velocity_prev(i))
!         end do
!      end if
end subroutine average_with_previous

! ---------------------------------------------------------------
! total_circulation_velocity = U_GSF + U_ES + U_SS at each interface,
! zeroed below 1e-20 cm/s; any_transport_active if any is positive.
subroutine sum_velocities()
      do i =zone_min,zone_max
         total_circulation_velocity(i) = star%gsf_circulation_velocity(i) + &
              star%es_circulation_velocity(i) + star%secular_shear_velocity(i)
         if (total_circulation_velocity(i).lt.1.0d-20) &
              total_circulation_velocity(i)=0.0d0
         if (total_circulation_velocity(i).gt.0.0d0) any_transport_active=.true.
      end do
end subroutine sum_velocities

! ---------------------------------------------------------------
! With no_am_transport_in_core (mixing without angular momentum
! transport; 9/93): zero every velocity from zone_min outward up to the
! first interface where L reaches 99% of the total luminosity, to avoid
! numerical problems in the hydrogen-burning shell.
subroutine zero_transport_in_core()
      do i = zone_min,zone_max
         dell = rot_scr%luminosity_interface(i)/total_luminosity
         if (dell.lt.9.9d-1) then
            total_circulation_velocity(i) = 0.0d0
            star%gsf_circulation_velocity(i) = 0.0d0
            star%es_circulation_velocity(i) = 0.0d0
            star%secular_shear_velocity(i) = 0.0d0
!            rot_scr%es_advective_velocity(i) = 0.0d0
!            rot_scr%es_diffusive_velocity(i) = 0.0d0
         else
            exit
         end if
      end do
end subroutine zero_transport_in_core


! 4-point Lagrange interpolation of a(k..k+3) to interface i in ln a,
! with the weights rot_scr%lagrange_interp_weights(1:4,i) set by
! rotation_stability_setup (stability_coefficients).
      function log_lagrange4(a, k, i)
      double precision, intent(in) :: a(*)
      integer, intent(in) :: k, i
      double precision :: log_lagrange4
      log_lagrange4 = exp(log(a(k))* &
           rot_scr%lagrange_interp_weights(1,i)+log(a(k+1))* &
           rot_scr%lagrange_interp_weights(2,i)+ &
           log(a(k+2))*rot_scr%lagrange_interp_weights(3,i)+log(a(k+3))* &
           rot_scr%lagrange_interp_weights(4,i))
      end function log_lagrange4

end subroutine circulation_velocities

! ---------------------------------------------------------------
! The diffusion coefficients (velocity times the length scale r) of the
! angular momentum and composition equations, at the interfaces
! 2..num_zones; element i is the interface between zones i-1 and i.
! Inputs: radius_mid_prev, the unlogged zone radii (radius in
! circulation_velocities; legacy HRU), and the velocities in star%.
! Outputs:
!   radius_mid(i) = (radius_mid_prev(i) + radius_mid_prev(i-1))/2,
!     the interface radius r (radius_mid(1) = 0);
!   am_diffusion_coeff (legacy COD)
!     = 4 pi difad_velocity_scale (U_ES + U_GSF + U_SS) r;
!   mixing_diffusion_coeff (legacy COD2)
!     = 4 pi difad_velocity_scale mixing_velocity_scale
!       (es_mixing_scale U_ES + gsf_mixing_scale U_GSF
!        + secular_shear_mixing_scale U_SS) r,
!     each mechanism with its own mixing efficiency (2/91), or, with
!     use_zahn_coupling_factor, am_diffusion_coeff times star%vfc.
! difad_velocity_scale and mixing_velocity_scale are legacy FW and FC.
! With use_constant_background_diffusion, 4 pi difad_velocity_scale
! constant_background_diffusion_coeff is added to am_diffusion_coeff
! only (MHP 9/14, 8/13). With no_am_transport_in_core,
! am_diffusion_coeff is zero everywhere and only the mixing remains.
subroutine diffusion_velocity_scales(radius_mid_prev, num_zones, radius_mid, &
     am_diffusion_coeff, mixing_diffusion_coeff)
      use star_info_lib, only: star, json
      use physical_constants
      implicit none
      double precision, intent(in) :: radius_mid_prev(json)
      integer, intent(in) :: num_zones
      double precision, intent(out) :: radius_mid(json)
      double precision, intent(out) :: am_diffusion_coeff(json), &
           mixing_diffusion_coeff(json)
      integer :: i
      double precision :: con1, con2

      radius_mid(1) = 0.0d0
      do i = 2,num_zones
         radius_mid(i) = 0.5d0*(radius_mid_prev(i)+radius_mid_prev(i-1))
      end do
      con1 = c4pi*rotation_settings%difad_velocity_scale
      if (.not.rotation_settings%use_constant_background_diffusion) then
         if (.not.rotation_settings%use_zahn_coupling_factor) then
            con2 = con1*rotation_settings%mixing_velocity_scale
            do i = 2,num_zones
               am_diffusion_coeff(i)=con1*(star%es_circulation_velocity(i)+ &
                    star%gsf_circulation_velocity(i)+star%secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i)=con2*(rotation_settings%es_mixing_scale* &
                    star%es_circulation_velocity(i)+rotation_settings%gsf_mixing_scale* &
                    star%gsf_circulation_velocity(i)+rotation_settings%secular_shear_mixing_scale* &
                    star%secular_shear_velocity(i))*radius_mid(i)
            end do
         else
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*(star%es_circulation_velocity(i)+ &
                    star%gsf_circulation_velocity(i)+star%secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i) = am_diffusion_coeff(i)*star%vfc(i)
            end do
         end if
      else
         if (.not.rotation_settings%use_zahn_coupling_factor) then
            con2 = con1*rotation_settings%mixing_velocity_scale
            do i = 2,num_zones
               am_diffusion_coeff(i)=con1*(rotation_settings%constant_background_diffusion_coeff+ &
                    (star%es_circulation_velocity(i)+star%gsf_circulation_velocity(i)+ &
                    star%secular_shear_velocity(i))*radius_mid(i))
               mixing_diffusion_coeff(i)=con2*(rotation_settings%es_mixing_scale* &
                    star%es_circulation_velocity(i)+rotation_settings%gsf_mixing_scale* &
                    star%gsf_circulation_velocity(i)+rotation_settings%secular_shear_mixing_scale* &
                    star%secular_shear_velocity(i))*radius_mid(i)
            end do
         else
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*(star%es_circulation_velocity(i)+ &
                    star%gsf_circulation_velocity(i)+star%secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i) = am_diffusion_coeff(i)*star%vfc(i)
! The background coefficient is added after the mixing coefficient is
! formed, so it does not mix (MHP 8/13).
               am_diffusion_coeff(i) = am_diffusion_coeff(i)+ &
                    con1*rotation_settings%constant_background_diffusion_coeff
            end do
         end if
      end if
! Suppress angular momentum transport but keep the mixing.
      if (rotation_settings%no_am_transport_in_core) then
         do i = 2,num_zones
            am_diffusion_coeff(i) = 0.0d0
         end do
      end if
      return
end subroutine diffusion_velocity_scales

! ---------------------------------------------------------------
! The coupling factor star%vfc between the angular momentum and the
! composition diffusion coefficients, with the prescription of Zahn
! (1992), at the interfaces zone_min..zone_max:
!   vfc = min(1, f_c U/|2V - alpha U|),
! where U is the total circulation velocity (diffusion_velocity),
!   V = (1/(6 rho r)) d(rho r^2 U)/dr,
!   alpha = (1/2) d ln(r^2 omega)/d ln r,
! and f_c = mixing_velocity_scale (legacy FC) plays the role of C/30.
! vfc = 0 where U <= 0, and everywhere if f_c <= 0 (no mixing).
subroutine zahn_coupling_factor(log_density, radius, diffusion_velocity, zone_min, &
     zone_max, angular_velocity)
      use rotation_scratch_lib

      use star_info_lib, only: star, json
      use physical_constants
      use math_lib
      implicit none

! log_density: log10 rho (g/cm^3) at the zone centres.
! radius: the unlogged radius (cm) at the zone centres.
! diffusion_velocity: the total velocity U at the interfaces (zone edges).
! zone_min, zone_max: the first and last interface.
! angular_velocity: omega at the zone centres.
      double precision, intent(in) :: log_density(json), radius(json), &
           diffusion_velocity(json)
      integer, intent(in) :: zone_min, zone_max
      double precision, intent(in) :: angular_velocity(json)

      double precision :: half_dlnj_dlnr(json), circ_velocity(json)
      integer :: zone_index
      double precision :: omega_mid
      double precision :: interp_weight
      double precision :: v_minus, v_plus
      double precision :: density_r2_v_minus, density_r2_v_plus
      double precision :: density_radius, d_density_r2, du_dr
      double precision :: denom_test

! f_c <= 0: no mixing.

      if (rotation_settings%mixing_velocity_scale.le.0.0D0) then
         do zone_index = zone_min, zone_max
            star%vfc(zone_index) = 0.0D0
         end do
         return
      endif
! alpha = (1/2) d ln(r^2 omega)/d ln r
!       = (1/2) d(omega r^2)/dr/(omega r)
! at each interface, with omega the mean of the two zones. (The F77 code
! lacked the 1/2; corrected in 2026.)
      do zone_index = zone_min, zone_max
         omega_mid = 0.5D0*(angular_velocity(zone_index)+ &
              angular_velocity(zone_index-1))
         half_dlnj_dlnr(zone_index) = 0.5D0*(angular_velocity(zone_index)* &
              radius(zone_index)**2-angular_velocity(zone_index-1)* &
              radius(zone_index-1)**2)/ &
              (radius(zone_index)-radius(zone_index-1))/ &
              (omega_mid*rot_scr%radius_interface(zone_index))
      end do
! 2V = (1/(3 rho r)) d(rho r^2 U)/dr at interface i, from rho r^2 U at
! the zone centres i-1 and i, with U at a zone centre interpolated
! linearly between the interfaces on either side; rho and r at the
! interface are rot_scr's.
      interp_weight = (radius(zone_min)-rot_scr%radius_interface(zone_min))/ &
           (rot_scr%radius_interface(zone_min+1)-rot_scr%radius_interface(zone_min))
      v_minus = diffusion_velocity(zone_min)+interp_weight* &
           (diffusion_velocity(zone_min+1)-diffusion_velocity(zone_min))
      density_r2_v_minus = exp(ln10*log_density(zone_min))* &
           radius(zone_min)**2*v_minus
      do zone_index = zone_min + 1, zone_max - 1
         interp_weight = (radius(zone_index)-rot_scr%radius_interface(zone_index))/ &
              (rot_scr%radius_interface(zone_index+1)-rot_scr%radius_interface(zone_index))
         v_plus = diffusion_velocity(zone_index)+interp_weight* &
              (diffusion_velocity(zone_index+1)-diffusion_velocity(zone_index))
         density_r2_v_plus = exp(ln10*log_density(zone_index))* &
              radius(zone_index)**2*v_plus
         density_radius = rot_scr%density_interface(zone_index)*rot_scr%radius_interface(zone_index)
         circ_velocity(zone_index) = cc13*(density_r2_v_plus-density_r2_v_minus)/ &
              density_radius/(radius(zone_index)-radius(zone_index-1))
         density_r2_v_minus = density_r2_v_plus
      end do
! At the first and last interfaces,
!   2V = (1/3) [U d(rho r^2)/dr/(rho r) + r dU/dr],
! with d(rho r^2)/dr across the interface and dU/dr the one-sided
! difference to the neighbouring interface.
      density_radius = rot_scr%density_interface(zone_min)*rot_scr%radius_interface(zone_min)
      d_density_r2 = (exp(ln10*log_density(zone_min))*radius(zone_min)**2 - &
           exp(ln10*log_density(zone_min-1))*radius(zone_min-1)**2)/ &
           (radius(zone_min) - radius(zone_min-1))
      du_dr = (diffusion_velocity(zone_min+1)-diffusion_velocity(zone_min))/ &
           (rot_scr%radius_interface(zone_min+1)-rot_scr%radius_interface(zone_min))
      circ_velocity(zone_min) = cc13*(diffusion_velocity(zone_min)*d_density_r2/ &
           density_radius + rot_scr%radius_interface(zone_min)*du_dr)
      density_radius = rot_scr%density_interface(zone_max)*rot_scr%radius_interface(zone_max)
      d_density_r2 = (exp(ln10*log_density(zone_max))*radius(zone_max)**2 - &
           exp(ln10*log_density(zone_max-1))*radius(zone_max-1)**2)/ &
           (radius(zone_max) - radius(zone_max-1))
      du_dr = (diffusion_velocity(zone_max)-diffusion_velocity(zone_max-1))/ &
           (rot_scr%radius_interface(zone_max)-rot_scr%radius_interface(zone_max-1))
      circ_velocity(zone_max) = cc13*(diffusion_velocity(zone_max)*d_density_r2/ &
           density_radius + rot_scr%radius_interface(zone_max)*du_dr)
! vfc, the same for all mechanisms and at most 1:
!   vfc = f_c/max(|2V - alpha U|/U, f_c).
      do zone_index = zone_min, zone_max
         if (diffusion_velocity(zone_index).le.0.0D0) then
            star%vfc(zone_index) = 0.0D0
         else
            denom_test = max(abs(circ_velocity(zone_index)- &
                 half_dlnj_dlnr(zone_index)*diffusion_velocity(zone_index))/ &
                 diffusion_velocity(zone_index),rotation_settings%mixing_velocity_scale)
            star%vfc(zone_index) = rotation_settings%mixing_velocity_scale/denom_test
         endif
      end do

      return
end subroutine zahn_coupling_factor

! ---------------------------------------------------------------
! The quadrupole moment Phi of the rotating star at the zone centres
! (rot_scr%quadrupole_moment, legacy QUAD), with the Sweet (1950)
! formulation as given in Zahn (1992). Phi solves
!   (1/r) d^2(r Phi)/dr^2 - 6 Phi/r^2 + 4 pi G rho (d rho/dP) Phi
!     = -(4/3) pi G (r^2/g) d(rho omega^2)/dr,
! with g the local mean gravity (gravity, ~G m/r^2). With dP/dr =
! -rho g this is
!   Phi'' + (2/r) Phi' - 6 Phi/r^2 - 4 pi G (d rho/dr/g) Phi
!     = -(4/3) pi G (r^2/g) d(rho omega^2)/dr,
! which the code differences on the zone radii into the tridiagonal
! system
!   sub_diag Phi(i-1) + diag Phi(i) + super_diag Phi(i+1) = rhs.
! Boundary conditions: Phi = 0 at the centre; at the surface no mass
! lies outside, so Phi varies as 1/r^3 (the exterior l = 2 solution).
! gravity is copied to rot_scr%local_gravity for circulation_velocities.
! ierr = 1 if a pivot of the elimination is zero.
subroutine compute_quadrupole(log_density, gravity, radius, angular_velocity, num_zones, ierr)
      use rotation_scratch_lib
      use star_info_lib, only: json
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: log_density(json), gravity(json), &
           radius(json), angular_velocity(json)
      integer, intent(in) :: num_zones
      integer, intent(out) :: ierr

! The tridiagonal system and the work arrays of its solution (Thomas
! algorithm).
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), solution(json), gamma_elim(json)

      double precision :: density_omega2(json)
      integer :: zone_index, matrix_row
      double precision :: four_pi_g
      double precision :: dr, weight_plus, weight_minus, inv_dr2, dr_inv_r
      double precision :: drho_dr
      double precision :: dr_plus, dr_minus, d_density_omega2
      double precision :: radius_plus, outer_boundary_ratio
      double precision :: pivot

      ierr = 0
      do zone_index = 1,num_zones
         density_omega2(zone_index) = exp(ln10*log_density(zone_index))* &
              angular_velocity(zone_index)**2
      end do
      four_pi_g = exp(ln10*(c4pil+cgl))
! The derivatives are centred, with dr the smaller of the left and right
! intervals dr_minus, dr_plus and weights w_plus = dr/dr_plus,
! w_minus = dr/dr_minus:
!   Phi''  = [(Phi(i+1) - Phi(i))/dr_plus - (Phi(i) - Phi(i-1))/dr_minus]/dr,
!   Phi'   = [(Phi(i+1) - Phi(i))/dr_plus + (Phi(i) - Phi(i-1))/dr_minus]/2,
! and d rho/dr and d(rho omega^2)/dr likewise.
! Zone 1: the inner neighbour is the centre, r = 0, where Phi = 0, so it
! drops out (sub_diag(1) = 0). There d rho/dr = (rho(2) - rho(1))/dr and
! d(rho omega^2)/dr is taken as omega^2 d rho/dr.
      dr = min(radius(1),radius(2)-radius(1))
      weight_plus = dr/(radius(2)-radius(1))
      weight_minus = dr/radius(1)
      inv_dr2 = 1.0D0/dr**2
      dr_inv_r = 1.0D0/dr/radius(1)
      sub_diag(1) = 0.0D0
      drho_dr = (exp(ln10*log_density(2))-exp(ln10*log_density(1)))/dr
      diag(1) = -inv_dr2*(weight_plus+weight_minus)+dr_inv_r* &
           (weight_minus-weight_plus)-6.0D0/radius(1)**2 &
           -four_pi_g*drho_dr/gravity(1)
      super_diag(1) = weight_plus*(inv_dr2+dr_inv_r)
      rhs(1) = -four_pi_g*drho_dr*(radius(1)*angular_velocity(1))**2/3.0D0/ &
           gravity(1)
! The interior zones. The Phi'' difference divides by min(dr_minus,
! dr_plus) rather than their mean, so it is exact only on a uniform grid.
      do zone_index = 2,num_zones-1
         dr_plus = radius(zone_index+1)-radius(zone_index)
         dr_minus = radius(zone_index)-radius(zone_index-1)
         dr = min(dr_minus,dr_plus)
         weight_plus = dr/dr_plus
         weight_minus = dr/dr_minus
         inv_dr2 = 1.0D0/dr**2
         dr_inv_r = 1.0D0/dr/radius(zone_index)
         drho_dr = 0.5D0*(weight_plus*exp(ln10*log_density(zone_index+1))+ &
              (weight_minus-weight_plus)* &
              exp(ln10*log_density(zone_index))-weight_minus* &
              exp(ln10*log_density(zone_index-1)))/dr
         sub_diag(zone_index) = weight_minus*(inv_dr2-dr_inv_r)
         diag(zone_index) = -inv_dr2*(weight_plus+weight_minus)+dr_inv_r* &
              (weight_minus-weight_plus)-6.0D0/radius(zone_index)**2 &
              -four_pi_g*drho_dr/gravity(zone_index)
         super_diag(zone_index) = weight_plus*(inv_dr2+dr_inv_r)
         d_density_omega2 = 0.5D0*(weight_plus*density_omega2(zone_index+1) + &
              (weight_minus-weight_plus)*density_omega2(zone_index) - &
              weight_minus*density_omega2(zone_index-1))/dr
         rhs(zone_index) = -d_density_omega2*four_pi_g* &
              radius(zone_index)**2/3.0D0/gravity(zone_index)
      end do
! The surface: Phi' + 3 Phi/r = 0, and outside the star Phi(r) =
! Phi(R) (R/r)^3. Phi'' uses a ghost point at R + dr with that value;
! (2/r) Phi' - 6 Phi/r^2 becomes -12 Phi/r^2; the density and forcing
! terms are dropped (rhs = 0).
      dr = radius(num_zones) - radius(num_zones-1)
      radius_plus = radius(num_zones)+dr
      outer_boundary_ratio = (radius(num_zones)/radius_plus)**3
      inv_dr2 = 1.0D0/dr**2
      sub_diag(num_zones) = inv_dr2
      diag(num_zones) = (outer_boundary_ratio-2.0D0)*inv_dr2-1.2D1/ &
           radius(num_zones)**2
      super_diag(num_zones) = 0.0D0
      rhs(num_zones) = 0.0D0
! Solve by forward elimination and back substitution (the Thomas
! algorithm, tridag of Numerical Recipes).
      if (diag(1).eq.0.0D0) then
         write(*,*) 'compute_quadrupole: singular first pivot'
         ierr = 1
         return
      end if
      pivot = diag(1)
      solution(1) = rhs(1)/pivot
      do matrix_row = 2,num_zones
         gamma_elim(matrix_row) = super_diag(matrix_row-1)/pivot
         pivot = diag(matrix_row) - sub_diag(matrix_row)*gamma_elim(matrix_row)
         if (pivot.eq.0.0D0) then
            write(*,*) 'compute_quadrupole: singular pivot in elimination'
            ierr = 1
            return
         end if
         solution(matrix_row) = (rhs(matrix_row) - &
              sub_diag(matrix_row)*solution(matrix_row-1))/pivot
      end do
! Back substitution.
      do matrix_row = num_zones-1,1,-1
         solution(matrix_row) = solution(matrix_row) - &
              gamma_elim(matrix_row+1)*solution(matrix_row+1)
      end do
      do zone_index = 1,num_zones
         rot_scr%quadrupole_moment(zone_index) = solution(zone_index)
         rot_scr%local_gravity(zone_index) = gravity(zone_index)
      end do

      return
end subroutine compute_quadrupole

end module circulation

!----------------------------------------------------------------------
! circulation
!----------------------------------------------------------------------
! The velocities that drive rotational mixing: the circulation and
! instability velocities (circulation_velocities), their scales
! (diffusion_velocity_scales), the Zahn coupling factor and the
! quadrupole moment of the distorted star.
module circulation
      implicit none
      private
      public :: circulation_velocities, diffusion_velocity_scales, zahn_coupling_factor, compute_quadrupole

contains

! ---------------------------------------------------------------
!   SKIPPING CONVECTIVE REGIONS, CHECK FOR INSTABILITY AGAINST EACH
!   SECULAR MECHANISM.IF UNSTABLE,COMPUTE A CIRCULATION VELOCITY.
!   THE PROCEDURE FOR CALCULATING STABILITY IS AS FOLLOWS:
!      STABILITY IS CALCULATED FOR THE I/I-1 INTERFACE.
!      EACH STABILITY CRITERION IS EXPRESSED AS FOLLOWS:
!       A COMBINATION OF VARIABLES(RHO,P,ETC.)STORED AT EACH SHELL
!       MULTIPLIED BY A GRADIENT BETWEEN THE SHELLS MUST BE LESS THAN
!       SOME CRITICAL NUMBER.
! Computes, at each unstable interface, the Eddington-Sweet meridional
! circulation velocity (star%es_circulation_velocity, Zahn 1991/1992, with
! a quadrupole correction), the GSF-instability circulation velocity
! (star%gsf_circulation_velocity, Kippenhahn 1980 estimate), and the
! diffusive/dynamical shear velocity (star%secular_shear_velocity), then
! combines them into total_circulation_velocity (HV).
subroutine circulation_velocities(log_radius, radius, zone_min, zone_max, iteration, &
     am_transport_convective_flag, num_zones, omega, any_transport_active, &
     dlnomega_dlnr, dlnomega_dlnr_max, total_circulation_velocity, &
     total_luminosity, timestep, log_pressure)
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
      double precision, intent(in) :: total_luminosity, timestep, &
           log_pressure(json)
! --- locals ---
      integer :: i, i0, i1
      double precision :: dr, gmid, phistd, phis2, phi2, qw, &
           q1, q2, q0, qmu, qp, ddel, ddtest, ddtest2, qqq, &
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
      if (star%ctrl%no_am_transport_in_core) call zero_transport_in_core()

contains

! ---------------------------------------------------------------
! The previous iteration's velocities, for the averaging.
subroutine store_previous_velocities()
!  STORE THE OLD VELOCITY ESTIMATES FOR LATER USE.
      do i = zone_min,zone_max
         circ_scr%gsf_circulation_velocity_prev(i) = star%gsf_circulation_velocity(i)
         circ_scr%es_circulation_velocity_prev(i) = star%es_circulation_velocity(i)
         circ_scr%secular_shear_velocity_prev(i) = star%secular_shear_velocity(i)
      end do
      if (star%ctrl%use_diffusion_advection_transport) then
         do i = zone_min,zone_max
            rot_scr%es_advective_velocity_prev(i) = rot_scr%es_advective_velocity(i)
            rot_scr%es_diffusive_velocity_prev(i) = rot_scr%es_diffusive_velocity(i)
         end do
      end if
end subroutine store_previous_velocities

! ---------------------------------------------------------------
subroutine interpolate_to_interfaces()
!  MEAN MOLECULAR WEIGHT (AS WELL AS QUANTITIES WHICH DEPEND ON IT) AND
!  THE ANGULAR VELOCITY DISTRIBUTION CHANGE DURING A DIFFUSION TIMESTEP.
!  FIND THE NEW VALUES AT THE MIDPOINT IN RADIUS BETWEEN MASS SHELLS
!  WHERE STABILITY IS EVALUATED.  4 POINT LAGRANGIAN INTERPOLATION IS USED
!  FOR THE KINEMATIC VISCOSITY AND THERMAL DIFFUSIVITY.
!  TREATMENT OF FIRST INTERFACE.
      if (zone_min.lt.3) then
         rot_scr%kinematic_viscosity_interface(2)=log_lagrange4(mid_scr%visc, 1, 2)
         rot_scr%thermal_diffusivity_interface(2)=log_lagrange4(mid_scr%thdif, 1, 2)
         i0=3
      else
         i0=zone_min
      end if
!  TREATMENT OF LAST INTERFACE.
      if (zone_max.eq.num_zones) then
         rot_scr%kinematic_viscosity_interface(num_zones)= &
              log_lagrange4(mid_scr%visc, num_zones-3, num_zones)
         rot_scr%thermal_diffusivity_interface(num_zones)= &
              log_lagrange4(mid_scr%thdif, num_zones-3, num_zones)
         i1=num_zones-1
      else
         i1=zone_max
      end if
!  GENERAL CASE.
      do i = i0,i1
         rot_scr%kinematic_viscosity_interface(i)=log_lagrange4(mid_scr%visc, i-2, i)
         rot_scr%thermal_diffusivity_interface(i)=log_lagrange4(mid_scr%thdif, i-2, i)
      end do
! USE LINEAR INTERPOLATION FOR OMEGA AND MU.
      do i = 2,num_zones
         rot_scr%mean_molecular_weight_interface(i) = 0.5d0*(mid_scr%mu(i)+mid_scr%mu(i-1))
         rot_scr%omega_interface(i) = 0.5d0*(omega(i)+omega(i-1))
      end do
end subroutine interpolate_to_interfaces

! ---------------------------------------------------------------
subroutine eddington_sweet_velocity()
! MHP 8/03 OMITTED OLD KM1974 MERIDIONAL CIRCULATION VELOCITY ESTIMATE.
! EDDINGTON CIRCULATION, ZAHN (1991) VELOCITY ESTIMATE (COEFFICIENTS
! DEFINED IN rotation_stability_setup).
      do i = zone_min,zone_max
         if (am_transport_convective_flag(i).and. &
              am_transport_convective_flag(i-1)) cycle
! ORIGINAL ESTIMATE,USING DG/G = W**2 R**3 / GM.
         star%es_circulation_velocity(i) = rot_scr%omega_interface(i)**2* &
              (rot_scr%es_velocity_coeff1(i)+rot_scr%omega_interface(i)**2*rot_scr%es_velocity_coeff2(i))
! QUADRUPOLE TERM ADDED, DG/G COMPUTED AS PER ZAHN 1992.
! THE VARIATION IN GRAVITY ON A LEVEL SURFACE CAN BE EXPRESSED IN
! GENERAL AS 1/3 OMEGA**2 *D/DR((R**2 - PHI)/G), WHERE G IS THE
! LOCAL AVERAGE GRAVITY (~GM/R**2) AND PHI IS THE QUADRUPOLE MOMENT.
! THE CLASSICAL EXPRESSION NEGLECTS THE QUADRUPOLE AND THE SPATIAL
! DERIVATIVE OF THE MASS (PHISTD).  THE VELOCITY ESTIMATE FOR
! MERIDIONAL CIRCULATION IS THEREFORE MULTIPLED BY THE RATIO OF
! THE CORRECTED EXPRESSION TO THE CLASSICAL ONE (RAT).
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
         if (.not.star%ctrl%use_diffusion_advection_transport) then
            star%es_circulation_velocity(i) = &
                 abs(rot_scr%circulation_correction_ratio(i)*star%es_circulation_velocity(i)- &
                 rot_scr%es_shear_coeff(i)*qw)
         else
! MHP 05/02 ADD FACTOR OF 1/5 HERE
! MHP 06/02 ADDED D THETA/DT TERM
            rot_scr%theta_new(i) = rot_scr%es_relaxation_factor(i)*(rot_scr%theta_mean(i)*qw- &
                 rot_scr%theta_prev(i))/timestep
            rot_scr%es_advective_velocity(i) = 0.2d0* &
                 (rot_scr%circulation_correction_ratio(i)*star%es_circulation_velocity(i)+ &
                 rot_scr%theta_new(i))
            q1 = rot_scr%difad_shear_coeff1(i)*rot_scr%omega_interface(i)**2
            q2 = rot_scr%difad_shear_coeff2(i)*qw
            q0 = rot_scr%es_shear_coeff(i)*rot_scr%omega_interface(i)**2
            rot_scr%es_diffusive_velocity(i) = 0.2d0*(q0+q1+q2)
            star%es_circulation_velocity(i) = &
                 rot_scr%circulation_correction_ratio(i)*star%es_circulation_velocity(i)+ &
                 rot_scr%theta_new(i)+ &
                 (rot_scr%es_shear_coeff(i)+rot_scr%difad_shear_coeff1(i))*qw+ &
                 rot_scr%difad_shear_coeff2(i)*((omega(i)-omega(i-1))/dr)**2
            rot_scr%second_deriv_geom_factor_eqgrid(i) = &
                 rot_scr%second_deriv_geom_factor(i)*rot_scr%omega_interface(i)
         end if
      end do
end subroutine eddington_sweet_velocity

! ---------------------------------------------------------------
subroutine inhibit_by_mu_gradient()
! INHIBITION CAUSED BY GRADIENTS IN MEAN MOLECULAR WEIGHT.
! INCLUDE ONLY THE ZAHN & MAEDER 1998 TREATMENT, WHERE (DEL - DEL AD) IS
! REPLACED BY (DEL + DEL MU - DEL AD).
       do i = zone_min,zone_max
          qmu = log(mid_scr%mu(i))-log(mid_scr%mu(i-1))
          qp = ln10*(log_pressure(i)-log_pressure(i-1))
          ddel=rot_scr%del_grad_diff_interface(i)+ qmu/qp
          ddtest = max(rot_scr%del_grad_diff_interface(i),1.0d-3)
          ddtest2 = max(ddel,1.0d-3)
          qqq = ddtest/ddtest2
          star%es_circulation_velocity(i) = star%es_circulation_velocity(i)*qqq
          if (star%ctrl%use_diffusion_advection_transport) then
             rot_scr%es_advective_velocity(i) = rot_scr%es_advective_velocity(i)*qqq
             rot_scr%es_diffusive_velocity(i) = rot_scr%es_diffusive_velocity(i)*qqq
          end if
! ADD MU GRADIENTS TO VELOCITY ESTIMATES
          rot_scr%es_relaxation_factor(i) = rot_scr%es_relaxation_factor(i)*qqq
          rot_scr%velocity_coeff0(i) = rot_scr%velocity_coeff0(i)*qqq
       end do
end subroutine inhibit_by_mu_gradient

! ---------------------------------------------------------------
subroutine gsf_velocity()
!  GSF INSTABILITY.  VELOCITY DEFINED IN rotation_stability_setup.
! OMIT ALL BUT THE KIPPENHAHN 1980 ESTIMATE.
! IGSF IS RETAINED, AND THE VALUE DETERMINES WHAT IS PERMITTED TO
! INHIBIT IT.
! IGSF = 0: INHIBITED BY VISCOSITY AND MU GRADIENTS
! IGSF = 1: NEITHER VISCOSITY NOR MU GRADIENTS INHIBIT
! IGSF = 2: INHIBITED BY VISCOSITY ONLY
! NOTE THAT THE ABCD INSTABILITY - RELATED TO THE GSF INSTABILITY -
! WOULD NOT BE INHIBITED BY VISCOSITY/MU GRADIENTS, SO INCLUDING IT
! IS EQUIVALENT TO SETTING IGSF=1 (THE ORIGINAL K1980 ESTIMATE).
! THE INHIBITION FACTORS COME FROM SUBSEQUENT WORK BY SPRUIT.
! KIPPENHAHN (1980) ESTIMATE
      do i = zone_min,zone_max
         if (am_transport_convective_flag(i).and. &
              am_transport_convective_flag(i-1)) cycle
! MHP 8/93 STABILITY CONDITION ADDED, NEGLECTING THE EFFECTS OF
! MU GRADIENTS.
         if (star%ctrl%gsf_inhibition_mode.eq.2 .or. star%ctrl%gsf_inhibition_mode.eq.0) then
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
         if (star%ctrl%gsf_inhibition_mode.eq.0) then
            qwrmx=2.0d0*sqrt(rot_scr%gravity_factor_interface(i)* &
                 abs(mid_scr%mu(i)-mid_scr%mu(i-1)) &
                 /dr/rot_scr%mean_molecular_weight_interface(i))
            if (abs(dlnomega_dlnr(i)).lt.qwrmx) then
               star%gsf_circulation_velocity(i) = 0.0d0
               cycle
            else
! MHP 05/02 ADDED TESTS TO AVOID DIVIDE BY ZERO
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
! GSF IS TRIGGERED BY D OMEGA/DZ NONZERO (I.E. ROTATION NOT
! ON CYLINDERS,WHICH IS TRUE IN GENERAL IN OUR MODELS), OR
! BY D LN(J/M)/DR < 0, WHICH SOMETIMES OCCURS.
! IF D LN(J/M)/DR < 0, CHECK SECOND CRITERIA AND USE THE LARGEST
! VELOCITY; OTHERWISE,USE ONLY D OMEGA/DZ NON-ZERO CRITERION.
         if (omega(i)*radius(i)**2 .lt. omega(i-1)*radius(i-1)**2) then
            dlnjmdr = abs(2.0d0/rmid+(log(omega(i))- &
                      log(omega(i-1)))/dr)
            fx = max(2.0d0*dlnjmdr,0.25d0*dlnwdr)
! FACTOR OF R IN THE DENOMINATOR OCCURS BECAUSE ALL THE CIRCULATION
! VELOCITIES ARE LATER MULTIPLIED BY R (THE "LENGTH SCALE").
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
subroutine shear_velocities()
!  DIFFUSIVE AND DYNAMICAL SHEAR INSTABILITIES - REF. ENDAL&SOFIA PAPER II.
      do i = zone_min,zone_max
!  CHECK FOR OPERATION OF DYNAMICAL SHEAR.
!  IF DYNAMICAL SHEAR IS OPERATING,SET SECULAR SHEAR VELOCITY TO MAXIMUM
!  VALUE AND COMPUTE (LARGE) DYNAMICAL SHEAR VELOCITY.
         if (abs(dlnomega_dlnr(i)).gt.dlnomega_dlnr_max(i)) then
            qwr = abs(dlnomega_dlnr(i))
            star%secular_shear_velocity(i)=8.0d0/4.5d1* &
                 rot_scr%thermal_diffusivity_interface(i)* &
                 (qwr/dlnomega_dlnr_max(i))**2/rot_scr%radius_interface(i)
!  CORRECT GSF VELOCITY AS WELL.
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
! KIPPENHAHN (1980) ESTIMATE
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
!   FIND MAXIMUM GRADIENT IN OMEGA ALLOWED BY SECULAR SHEAR.
!   THE RUN OF QWRMAX INPUT IS THAT ALLOWED BY THE DYNAMICAL SHEAR;
!   THE SECULAR SHEAR RICHARDSON # IS RELATED BY
!   RICHNO(SECULAR) = PRANDTL# * CRITICAL REYNOLDS#/8 * RICHNO(DYNAMICAL)
!   PRANDTL # = KINEMATIC VISCOSITY/THERMOMETRIC DIFFUSIVITY.
! MHP 3/92 SQUARE ROOT OF PR# NEEDED, NOT PR # - ERROR CORRECTED!
! THE VELOCITY ESTIMATE HERE IS FROM ZAHN 1991.
         qwrmx = sqrt(rot_scr%kinematic_viscosity_interface(i)/ &
              rot_scr%thermal_diffusivity_interface(i)*1.25d-1*star%ctrl%critical_reynolds)* &
              dlnomega_dlnr_max(i)
         if (abs(dlnomega_dlnr(i)).gt.qwrmx) then
!  UNSTABLE; CHECK FOR MU GRADIENTS.
            if (abs((mid_scr%mu(i)-mid_scr%mu(i-1))/rot_scr%mean_molecular_weight_interface(i)) &
                 .lt.1.0d-10) then
               qwrmx2 = 0.0d0
               qwr = abs(dlnomega_dlnr(i)) - qwrmx
               star%secular_shear_velocity(i)=8.0d0/4.5d1* &
                    rot_scr%thermal_diffusivity_interface(i)* &
                    (qwr/dlnomega_dlnr_max(i))**2/rot_scr%radius_interface(i)
            else
!  CHECK FOR EFFECTS OF MU GRADIENT.
!  RICHNO = RHO/P*(-d lnT/d lnMU)*(del MU)/(del P)*(G/QWLNR)**2
!  WHERE -d lnT/d lnMU = (CON-1)/(1+3CON),CON=(a/3)T**4/P
!  GIVEN (1 - CON)P = CGAS*RHO*T/MU
               qwrmx2 = 2.0d0*sqrt(max(1.0d-20,rot_scr%mu_gradient_richardson_coeff(i)* &
                        abs((mid_scr%mu(i)-mid_scr%mu(i-1))/ &
                        rot_scr%mean_molecular_weight_interface(i))))
               if (abs(dlnomega_dlnr(i)).gt.qwrmx2) then
!  INTERFACE UNSTABLE WITH RESPECT TO BOTH CONDITIONS; CHOOSE THE
!  MAXIMUM GRADIENT IMPLIED BY THE SECOND CONDITION IF IT'S LARGER
!  THAN THE FIRST(I.E. IF A MU GRADIENT IS SLOWING J TRANSPORT).
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
subroutine apply_velocity_scales()
!  MULTIPLY THE VELOCITY ESTIMATES BY THE USER-DEFINED SCALE FACTORS
!  (MU-GRADIENT INHIBITION IS ALREADY INCLUDED ABOVE).
      do i = zone_min,zone_max
         star%es_circulation_velocity(i)=abs(star%ctrl%es_velocity_scale* &
              star%es_circulation_velocity(i))
         star%gsf_circulation_velocity(i) = star%ctrl%gsf_velocity_scale* &
              star%gsf_circulation_velocity(i)
         star%secular_shear_velocity(i)= star%ctrl%secular_shear_velocity_scale* &
              star%secular_shear_velocity(i)
      end do
! MHP 11/94
! REPEAT FOR DIF+AD
      if (star%ctrl%use_diffusion_advection_transport) then
         do i = zone_min,zone_max
               rot_scr%es_diffusive_velocity(i)=abs(star%ctrl%es_velocity_scale* &
                    rot_scr%es_diffusive_velocity(i))
               rot_scr%es_advective_velocity(i)=star%ctrl%es_velocity_scale* &
                    rot_scr%es_advective_velocity(i)
         end do
      end if
      if (star%ctrl%use_diffusion_advection_transport) then
! MHP 05/02
! CHANGED TO REFLECT THE DIFFERENT TREATMENT OF THE GSF INSTABILITY.
! THE ENDAL AND SOFIA VARIANT HAD V=THDIF*DLNWDR*W**2/R WHICH WAS
! THEN ADDED TO A DIFFUSION EQUATION D/DR(CON*V*R*DW/DR).
! IN THE ZAHN FORMULISM, WE ARE SOLVING AN EQUATION OF THE FORM
! D/DR(CON*V*W).  BECAUSE THE VELOCITY ESTIMATE ITSELF INCLUDES
! DW/DR THIS IS IN EFFECT A DIFFUSION TERM.  SO...WE NEED TO
! MULTIPLY THE ORIGINAL VELOCITY BY W AND DIVIDE BY DW/DR TO CAST
! IT AS AN ENTRY IN THE DIFFUSION EQUATION.
         do i = zone_min,zone_max
            rot_scr%es_diffusive_velocity(i) = rot_scr%es_diffusive_velocity(i)+ &
                 rot_scr%radius_interface(i)*(abs(star%gsf_circulation_velocity(i))+ &
                 abs(star%secular_shear_velocity(i)))
            rot_scr%shear_diffusion_coeff(i) = rot_scr%radius_interface(i)* &
                 abs(star%secular_shear_velocity(i))
            rot_scr%gsf_diffusion_coeff(i) = rot_scr%radius_interface(i)* &
                 abs(star%gsf_circulation_velocity(i))
         end do
      end if
end subroutine apply_velocity_scales

! ---------------------------------------------------------------
subroutine average_with_previous()
!  AVERAGE PREVIOUS AND NEW VELOCITY ESTIMATES AFTER THE FIRST ITERATION.
      do i = zone_min,zone_max
         star%gsf_circulation_velocity(i) = 0.5d0*(star%gsf_circulation_velocity(i) &
              + circ_scr%gsf_circulation_velocity_prev(i))
         star%es_circulation_velocity(i) = 0.5d0*(star%es_circulation_velocity(i) &
              + circ_scr%es_circulation_velocity_prev(i))
         star%secular_shear_velocity(i) = 0.5d0*(star%secular_shear_velocity(i) &
              + circ_scr%secular_shear_velocity_prev(i))
      end do
! MHP 11/94
      if (star%ctrl%use_diffusion_advection_transport) then
         do i = zone_min,zone_max
            rot_scr%es_advective_velocity(i) = 0.5d0*(rot_scr%es_advective_velocity(i) &
                 + rot_scr%es_advective_velocity_prev(i))
            rot_scr%es_diffusive_velocity(i) = 0.5d0*(rot_scr%es_diffusive_velocity(i) &
                 + rot_scr%es_diffusive_velocity_prev(i))
         end do
      end if
end subroutine average_with_previous

! ---------------------------------------------------------------
! The total velocity at each interface (below 1e-20, zero) and whether any
! transport acts.
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
subroutine zero_transport_in_core()
! 9/93 MIXING WITHOUT TRANSPORT ADDED.
! ZERO OUT COEFFICIENTS IN CORE TO AVOID NUMERICAL PROBLEMS IN
! THE H-BURNING SHELL.
      do i = zone_min,zone_max
         dell = rot_scr%luminosity_interface(i)/total_luminosity
         if (dell.lt.9.9d-1) then
            total_circulation_velocity(i) = 0.0d0
            star%gsf_circulation_velocity(i) = 0.0d0
            star%es_circulation_velocity(i) = 0.0d0
            star%secular_shear_velocity(i) = 0.0d0
! MHP 11/94
            rot_scr%es_advective_velocity(i) = 0.0d0
            rot_scr%es_diffusive_velocity(i) = 0.0d0
         else
            exit
         end if
      end do
end subroutine zero_transport_in_core


! 4-point Lagrangian interpolation in log(a) of a(k..k+3) to interface
! i, using the weights rot_scr%lagrange_interp_weights(1:4,i) set up
! by am_transport_grid. Token-identical to the former inline forms.
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
! 2/91 SUBROUTINE ALTERED TO ALLOW DIFFERENT FC'S FOR DIFFERENT
!      MECHANISMS (star%ctrl%es_mixing_scale etc.)
!
!  THIS SR DETERMINES THE RUN OF DIFFUSION COEFFICIENTS (VELOCITY
!  TIMES LENGTH SCALE) FOR THE DIFFUSION EQUATIONS.
!
!  INPUT VARIABLES
!
!       DECLARED:
!  radius_mid_prev (HRU) : RUN OF RADIUS (UNLOGGED) OF THE MIDPOINTS OF
!     THE SHELLS.
!  num_zones (M) : NUMBER OF SHELLS.
!  FC,FO,FW : USER PARAMETERS.  THE DIFFUSION COEFFICIENTS
!  CAN BE MODIFIED BY THESE THREE PARAMETERS.
!  VELM : RUN OF CONVECTIVE VELOCITIES.
!
!  OUTPUT VARIABLES
!
!       DECLARED
!  am_diffusion_coeff (COD) : RUN OF DIFFUSION COEFFICENTS FOR ANGULAR
!     MOMENTUM TRANSPORT.
!  mixing_diffusion_coeff (COD2) : SAME, WITH RESPECT TO COMPOSITION
!     TRANSPORT.
!  radius_mid (RMID) : radius_mid(I) IS DEFINED AS THE AVERAGE OF
!     radius_mid_prev(I) AND radius_mid_prev(I-1).
!
!  NOTE ON STORAGE: COD AND V ARE CALCULATED AT THE MIDPOINT
!  BETWEEN MASS POINTS. ELEMENT I CONTAINS THE INFORMATION FOR THE
!  (I,I-1) INTERFACE.
!
! radius_mid_prev is the run of (unlogged) shell-midpoint radii,
! matching circulation_velocities.f90's radius(json).
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
! MHP 9/14 ADDED LOOP TO ALLOW A CONSTANT BACKGROUND DIFFUSION COEFFICIENT
      con1 = c4pi*star%ctrl%difad_velocity_scale
      if (.not.star%ctrl%use_constant_background_diffusion) then
         if (.not.star%ctrl%lvfc) then
            con2 = con1*star%ctrl%mixing_velocity_scale
            do i = 2,num_zones
               am_diffusion_coeff(i)=con1*(star%es_circulation_velocity(i)+ &
                    star%gsf_circulation_velocity(i)+star%secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i)=con2*(star%ctrl%es_mixing_scale* &
                    star%es_circulation_velocity(i)+star%ctrl%gsf_mixing_scale* &
                    star%gsf_circulation_velocity(i)+star%ctrl%secular_shear_mixing_scale* &
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
         if (.not.star%ctrl%lvfc) then
            con2 = con1*star%ctrl%mixing_velocity_scale
            do i = 2,num_zones
               am_diffusion_coeff(i)=con1*(star%ctrl%constant_background_diffusion_coeff+ &
                    (star%es_circulation_velocity(i)+star%gsf_circulation_velocity(i)+ &
                    star%secular_shear_velocity(i))*radius_mid(i))
               mixing_diffusion_coeff(i)=con2*(star%ctrl%es_mixing_scale* &
                    star%es_circulation_velocity(i)+star%ctrl%gsf_mixing_scale* &
                    star%gsf_circulation_velocity(i)+star%ctrl%secular_shear_mixing_scale* &
                    star%secular_shear_velocity(i))*radius_mid(i)
            end do
         else
            do i = 2,num_zones
               am_diffusion_coeff(i) = con1*(star%es_circulation_velocity(i)+ &
                    star%gsf_circulation_velocity(i)+star%secular_shear_velocity(i))* &
                    radius_mid(i)
               mixing_diffusion_coeff(i) = am_diffusion_coeff(i)*star%vfc(i)
! MHP 8/13 ADD D.C. AFTER SCALE FACTOR FOR MIXING APPLIED
               am_diffusion_coeff(i) = am_diffusion_coeff(i)+ &
                    con1*star%ctrl%constant_background_diffusion_coeff
            end do
         end if
      end if
! OPTION TO SUPPRESS ANGULAR MOMENTUM TRANSPORT (BUT PERMIT MIXING).
      if (star%ctrl%no_am_transport_in_core) then
         do i = 2,num_zones
            am_diffusion_coeff(i) = 0.0d0
         end do
      end if
      return
end subroutine diffusion_velocity_scales

! ---------------------------------------------------------------
! This subroutine solves for the variable FC using the prescription
! of Zahn (1992). FC = MIN(1, C/30*U/(2V-ALPHA*U)), where
! ALPHA = 1/2 D(LN J/M)/D LN R and
! V = 1/(6*RHO*R) *D/DR(RHO* R**2 * U) (U = diffusion velocity).
! If FC = 0, skip calculation (no mixing).
subroutine zahn_coupling_factor(log_density, radius, diffusion_velocity, zone_min, &
     zone_max, angular_velocity)
      use rotation_scratch_lib

      use star_info_lib, only: star, json
      use physical_constants
      use math_lib
      implicit none

! INPUT VARIABLES :
! log_density - LOG10 OF DENSITY(GM/CM**3) AT ZONE CENTERS.
! radius - UNLOGGED RUN OF RADII AT ZONE CENTERS.
! diffusion_velocity - UNLOGGED RUN OF DIFFUSION VELOCITIES AT ZONE EDGES.
! zone_min,zone_max - STARTING AND ENDING SHELLS.
! angular_velocity - RUN OF ANGULAR VELOCITIES AT ZONE CENTERS.
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

! IF FC = 0, SKIP CALCULATION (NO MIXING).

      if (star%ctrl%mixing_velocity_scale.le.0.0D0) then
         do zone_index = zone_min, zone_max
            star%vfc(zone_index) = 0.0D0
         end do
         return
      endif
! DETERMINE ALPHA.
      do zone_index = zone_min, zone_max
         omega_mid = 0.5D0*(angular_velocity(zone_index)+ &
              angular_velocity(zone_index-1))
! (Fixed in the 2026 bug sweep:) the header defines alpha = 1/2 dln(r^2 Omega)/
! dln r (Zahn 1992) and the consumer forms |2V - alpha U| with
! circ_velocity already holding 2V; the inherited expression lacked
! the 1/2.
         half_dlnj_dlnr(zone_index) = 0.5D0*(angular_velocity(zone_index)* &
              radius(zone_index)**2-angular_velocity(zone_index-1)* &
              radius(zone_index-1)**2)/ &
              (radius(zone_index)-radius(zone_index-1))/ &
              (omega_mid*rot_scr%radius_interface(zone_index))
      end do
! DETERMINE 2V.
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
! FOR FIRST AND LAST SHELLS, CALCULATE U*D/DR(RHO*R**2) AND ASSUME
! DU/DR = VALUE DERIVED FROM THE LAST NEIGHBOR SHELL WITH NONZERO V.
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
! NOW COMPUTE RUN OF FC; THIS ASSUMES THAT FC = SAME FOR ALL
! MECHANISMS AND IS LIMITED TO A MAXIMUM OF 1.
      do zone_index = zone_min, zone_max
         if (diffusion_velocity(zone_index).le.0.0D0) then
            star%vfc(zone_index) = 0.0D0
         else
            denom_test = max(abs(circ_velocity(zone_index)- &
                 half_dlnj_dlnr(zone_index)*diffusion_velocity(zone_index))/ &
                 diffusion_velocity(zone_index),star%ctrl%mixing_velocity_scale)
            star%vfc(zone_index) = star%ctrl%mixing_velocity_scale/denom_test
         endif
      end do

      return
end subroutine zahn_coupling_factor

! ---------------------------------------------------------------
! Solves for the quadrupole moment of a rotating star, using the
! Sweet (1950) formulation as given in Zahn (1992). We are solving
! the equation
! 1/R D**2/DR2 (R*QUAD) -6*QUAD/R**2 - 4*PI*G*RHO*D(RHO)/DP * QUAD =
!  -4/3*PI*G*R**2/(GM/R**2) * D/DR (RHO*OMEGA**2)
! This can be expressed as D**2/DR2(QUAD) + 2/R DQUAD/DR -6*QUAD/R**2
! +4*PI*G*RHO*QUAD = 4*PI*G*RHO*R**2 *(OMEGA**2 - 2*OMEGA/G *DOMEGA/DR)
! using hydrostatic equilibrium and the ideal gas law
! D RHO/DR = DP/DR = -RHO*G.
! This defines a tridiagonal matrix system of the form
! A*QUAD(I-1)+B*QUAD(I)+C*QUAD(I+1) = D, which can be solved in
! the standard way for the run of QUAD(I) with appropriate boundary
! conditions. The central B.C. is that QUAD(0)=0; the surface B.C.
! is that there is no contribution to the quadrupole from outside the
! final shell, such that QUAD(I) varies as 1/R**3 (the exterior l=2
! solution).
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

! Tridiagonal-solve work arrays (Thomas algorithm) for the inline
! solve below.
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
! CENTER NUMERICAL DERIVATIVES, USING THE SMALLER OF THE LEFT AND
! RIGHT HANDED INTERVALS.
      dr = min(radius(1),radius(2)-radius(1))
      weight_plus = dr/(radius(2)-radius(1))
      weight_minus = dr/radius(1)
      inv_dr2 = 1.0D0/dr**2
      dr_inv_r = 1.0D0/dr/radius(1)
! CENTRAL BOUNDARY CONDITION : QUADRUPOLE GOES TO ZERO.
      sub_diag(1) = 0.0D0
      drho_dr = (exp(ln10*log_density(2))-exp(ln10*log_density(1)))/dr
      diag(1) = -inv_dr2*(weight_plus+weight_minus)+dr_inv_r* &
           (weight_minus-weight_plus)-6.0D0/radius(1)**2 &
           -four_pi_g*drho_dr/gravity(1)
      super_diag(1) = weight_plus*(inv_dr2+dr_inv_r)
      rhs(1) = -four_pi_g*drho_dr*(radius(1)*angular_velocity(1))**2/3.0D0/ &
           gravity(1)
! GENERAL CASE : SECOND DERIVATIVE BY A THREE-POINT DIFFERENCE ON THE
! (UNEQUAL) DR_MINUS/DR_PLUS SPACING, WEIGHTED BY THE SMALLER OF THE TWO.
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
! SURFACE B.C. D PHI/DR + 3 PHI/R = 0
! THE SECOND AND THIRD TERMS IN THE EXPRESSION CANCEL IN THIS CASE.
! OUTSIDE THE STAR PHI(R) = PHI(RSTAR)*(RSTAR/R)**3.
      dr = radius(num_zones) - radius(num_zones-1)
      radius_plus = radius(num_zones)+dr
      outer_boundary_ratio = (radius(num_zones)/radius_plus)**3
      inv_dr2 = 1.0D0/dr**2
      sub_diag(num_zones) = inv_dr2
      diag(num_zones) = (outer_boundary_ratio-2.0D0)*inv_dr2-1.2D1/ &
           radius(num_zones)**2
      super_diag(num_zones) = 0.0D0
      rhs(num_zones) = 0.0D0
! SOLUTION OF SYSTEM FROM NUMERICAL RECIPES.
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
! BACKSUBSTITUTION.
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

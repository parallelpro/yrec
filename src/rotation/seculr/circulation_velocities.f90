!----------------------------------------------------------------------
! vcirc
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original vcirc.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
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
!
!       SUBROUTINE VCIRC(HJM,HR,HRU,IMIN,IMAX,IT,LCZ,M,OMEGA,LDO,  ! KC 2025-05-31
!                        QWLNR,QWRMAX,HV,CLTOT,DT,HP)
subroutine circulation_velocities(log_radius, radius, zone_min, zone_max, iteration, &
     am_transport_convective_flag, num_zones, omega, any_transport_active, &
     dlnomega_dlnr, dlnomega_dlnr_max, total_circulation_velocity, &
     total_luminosity, timestep, log_pressure)
      use rotation_scratch_lib

      use star_info_lib, only: star, json
      use luout_lib
      use phys_const_lib
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
























! SPEED OF LIGHT.
!       DATA CC/2.99792458D10/
      integer :: i, i0, i1
      double precision :: dr, gmid, phistd, phis2, phi2, qw, &
           q1, q2, q0, qmu, qp, ddel, ddtest, ddtest2, qqq, v2, &
           qwrmx, fxx, qwrmx2, qwr, rmid, dlnwdr, dlnjmdr, fx, qwrmx_dyn, &
           wmin, wmid, dlnwdr0, dlnjmdr0, dell

! SKIPPING CONVECTIVE REGIONS, CHECK FOR INSTABILITY AGAINST EACH
! SECULAR MECHANISM.IF UNSTABLE,COMPUTE A CIRCULATION VELOCITY.
! THE PROCEDURE FOR CALCULATING STABILITY IS AS FOLLOWS:
!    STABILITY IS CALCULATED FOR THE I/I-1 INTERFACE.
!    EACH STABILITY CRITERION IS EXPRESSED AS FOLLOWS:
!     A COMBINATION OF VARIABLES(RHO,P,ETC.)STORED AT EACH SHELL
!     MULTIPLIED BY A GRADIENT BETWEEN THE SHELLS MUST BE LESS THAN
!     SOME CRITICAL NUMBER.

      any_transport_active = .false.
!  STORE THE OLD VELOCITY ESTIMATES FOR LATER USE.
      if (iteration.gt.1) then
         do i = zone_min,zone_max
            rot_scr%mu_gradient_velocity_prev(i) = circ_scr%mu_gradient_velocity(i)
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
      end if
      do i = 1,num_zones
         total_circulation_velocity(i) = 0.0d0
         star%es_circulation_velocity(i) = 0.0d0
         star%gsf_circulation_velocity(i) = 0.0d0
         star%secular_shear_velocity(i) = 0.0d0
         circ_scr%mu_gradient_velocity(i) = 0.0d0
      end do
!  MEAN MOLECULAR WEIGHT (AS WELL AS QUANTITIES WHICH DEPEND ON IT) AND
!  THE ANGULAR VELOCITY DISTRIBUTION CHANGE DURING A DIFFUSION TIMESTEP.
!  FIND THE NEW VALUES AT THE MIDPOINT IN RADIUS BETWEEN MASS SHELLS
!  WHERE STABILITY IS EVALUATED.  4 POINT LAGRANGIAN INTERPOLATION IS USED
!  FOR THE KINEMATIC VISCOSITY AND THERMAL DIFFUSIVITY.
!  TREATMENT OF FIRST INTERFACE.
      if (zone_min.lt.3) then
         rot_scr%kinematic_viscosity_interface(2)=exp(log(mix_scr%viscm(1))* &
              rot_scr%lagrange_interp_weights(1,2)+log(mix_scr%viscm(2))* &
              rot_scr%lagrange_interp_weights(2,2)+ &
              log(mix_scr%viscm(3))*rot_scr%lagrange_interp_weights(3,2)+log(mix_scr%viscm(4))* &
              rot_scr%lagrange_interp_weights(4,2))
         rot_scr%thermal_diffusivity_interface(2)=exp(log(mix_scr%thdifm(1))* &
              rot_scr%lagrange_interp_weights(1,2)+log(mix_scr%thdifm(2))* &
              rot_scr%lagrange_interp_weights(2,2) &
              +log(mix_scr%thdifm(3))*rot_scr%lagrange_interp_weights(3,2)+log(mix_scr%thdifm(4))* &
              rot_scr%lagrange_interp_weights(4,2))
         i0=3
      else
         i0=zone_min
      end if
!  TREATMENT OF LAST INTERFACE.
      if (zone_max.eq.num_zones) then
         rot_scr%kinematic_viscosity_interface(num_zones)= &
              exp(log(mix_scr%viscm(num_zones-3))*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              log(mix_scr%viscm(num_zones-2))*rot_scr%lagrange_interp_weights(2,num_zones) &
              +log(mix_scr%viscm(num_zones-1))*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              log(mix_scr%viscm(num_zones))*rot_scr%lagrange_interp_weights(4,num_zones))
         rot_scr%thermal_diffusivity_interface(num_zones)= &
              exp(log(mix_scr%thdifm(num_zones-3))* &
              rot_scr%lagrange_interp_weights(1,num_zones)+log(mix_scr%thdifm(num_zones-2))* &
              rot_scr%lagrange_interp_weights(2,num_zones)+log(mix_scr%thdifm(num_zones-1))* &
              rot_scr%lagrange_interp_weights(3,num_zones)+ &
              log(mix_scr%thdifm(num_zones))*rot_scr%lagrange_interp_weights(4,num_zones))
         i1=num_zones-1
      else
         i1=zone_max
      end if
!  GENERAL CASE.
      do i = i0,i1
         rot_scr%kinematic_viscosity_interface(i)=exp(log(mix_scr%viscm(i-2))* &
              rot_scr%lagrange_interp_weights(1,i)+log(mix_scr%viscm(i-1)) &
              *rot_scr%lagrange_interp_weights(2,i)+log(mix_scr%viscm(i))* &
              rot_scr%lagrange_interp_weights(3,i)+log(mix_scr%viscm(i+1)) &
              *rot_scr%lagrange_interp_weights(4,i))
         rot_scr%thermal_diffusivity_interface(i)=exp(log(mix_scr%thdifm(i-2))* &
              rot_scr%lagrange_interp_weights(1,i)+log(mix_scr%thdifm(i-1)) &
              *rot_scr%lagrange_interp_weights(2,i)+log(mix_scr%thdifm(i))* &
              rot_scr%lagrange_interp_weights(3,i)+log(mix_scr%thdifm(i+1)) &
              *rot_scr%lagrange_interp_weights(4,i))
      end do
! USE LINEAR INTERPOLATION FOR OMEGA AND MU.
      do i = 2,num_zones
         rot_scr%mean_molecular_weight_interface(i) = 0.5d0*(mix_scr%amum(i)+mix_scr%amum(i-1))
         rot_scr%omega_interface(i) = 0.5d0*(omega(i)+omega(i-1))
      end do
! MHP 8/03 OMITTED OLD KM1974 MERIDIONAL CIRCULATION VELOCITY ESTIMATE.
! THE IES FLAG IS THEREFORE NO LONGER IN USE.
!
! KIPPENHAHN AND MOLLENHOF(1974) MERIDIONAL CIRCULATION VELOCITY.
!      IF(IES.EQ.2)THEN
!         DO 30 I = IMIN,IMAX
!  SKIP CONVECTIVE INTERFACES.
!            IF(LCZ(I).AND.LCZ(I-1))GOTO 30
!  EDDINGTON CIRCULATION.  VELOCITY DEFINED IN SR SETUPV.
!            DR = HRU(I)-HRU(I-1)
!            DV = ((OMEGA(I)*HRU(I))**2-(OMEGA(I-1)*HRU(I-1))**2)/DR
! QUADRUPOLE TERM, AS PER ZAHN (1992), ADDED.
!            GMID = 0.5D0*(GG(I)+GG(I-1))
! THE VARIATION IN GRAVITY ON A LEVEL SURFACE CAN BE EXPRESSED IN
! GENERAL AS 1/3 OMEGA**2 *D/DR((R**2 - PHI)/G), WHERE G IS THE
! LOCAL AVERAGE GRAVITY (~GM/R**2) AND PHI IS THE QUADRUPOLE MOMENT.
! THE CLASSICAL EXPRESSION NEGLECTS THE QUADRUPOLE AND THE SPATIAL
! DERIVATIVE OF THE MASS (PHISTD).  THE VELOCITY ESTIMATE FOR
! MERIDIONAL CIRCULATION IS THEREFORE MULTIPLED BY THE RATIO OF
! THE CORRECTED EXPRESSION TO THE CLASSICAL ONE (RAT).
!            PHISTD = 2.0D0*CC23*WM(I)**2*RM(I)/GMID
!            PHIS2 = CC13*WM(I)**2*(HRU(I)**2/GG(I)-
!     *              HRU(I-1)**2/GG(I-1))/DR
!     *               +FACT6(I)*DV))
!   30    CONTINUE
!      ELSE
! ZAHN (1991) VELOCITY ESTIMATE.
      do i = zone_min,zone_max
         if (am_transport_convective_flag(i).and. &
              am_transport_convective_flag(i-1)) cycle
! ORIGINAL ESTIMATE,USING DG/G = W**2 R**3 / GM.
         star%es_circulation_velocity(i) = rot_scr%omega_interface(i)**2* &
              (rot_scr%es_velocity_coeff1(i)+rot_scr%omega_interface(i)**2*rot_scr%es_velocity_coeff2(i))
! QUADRUPOLE TERM ADDED, DG/G COMPUTED AS PER ZAHN 1992.
         dr = radius(i) - radius(i-1)
         gmid = 0.5d0*(rot_scr%local_gravity(i)+rot_scr%local_gravity(i-1))
         phistd = 2.0d0*cc23*rot_scr%omega_interface(i)**2*rot_scr%interface_radius(i)/gmid
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
!               VESA(I) = RAT(I)*VES(I)
!               VESD(I) = ABS(FES3(I)*WM(I)**2)
! MHP 06/02 ADDED D THETA/DT TERM
            rot_scr%theta_new(i) = rot_scr%es_relaxation_factor(i)*(rot_scr%theta_mean(i)*qw- &
                 rot_scr%theta_prev(i))/timestep
            rot_scr%es_advective_velocity(i) = 0.2d0* &
                 (rot_scr%circulation_correction_ratio(i)*star%es_circulation_velocity(i)+ &
                 rot_scr%theta_new(i))
            q1 = rot_scr%difad_shear_coeff1(i)*rot_scr%omega_interface(i)**2
            q2 = rot_scr%difad_shear_coeff2(i)*qw
            q0 = rot_scr%es_shear_coeff(i)*rot_scr%omega_interface(i)**2
!               VESD(I) = 0.2D0*ABS(FES3(I)*WM(I)**2)
!               VESD(I) = 0.2D0*ABS(Q0+Q1+Q2)
            rot_scr%es_diffusive_velocity(i) = 0.2d0*(q0+q1+q2)
! SECOND ORDER TERM
            rot_scr%vesd2(i) = 0.2d0*rot_scr%facd2(i)*rot_scr%omega_interface(i)
! THIRD ORDER TERM
            rot_scr%vesd3(i) = 0.2d0*rot_scr%facd3(i)*rot_scr%omega_interface(i)
!               VES(I) = ABS(RAT(I)*VES(I)+THN(I)-FES3(I)*QW)
            star%es_circulation_velocity(i) = &
                 rot_scr%circulation_correction_ratio(i)*star%es_circulation_velocity(i)+ &
                 rot_scr%theta_new(i)+ &
                 (rot_scr%es_shear_coeff(i)+rot_scr%difad_shear_coeff1(i))*qw+ &
                 rot_scr%difad_shear_coeff2(i)*((omega(i)-omega(i-1))/dr)**2
            rot_scr%second_deriv_geom_factor_eqgrid(i) = &
                 rot_scr%second_deriv_geom_factor(i)*rot_scr%omega_interface(i)
         end if
      end do
!      ENDIF
! INHIBITION CAUSED BY GRADIENTS IN MEAN MOLECULAR WEIGHT.
! INCLUDE ONLY THE ZAHN & MAEDER 1998 TREATMENT, WHERE (DEL - DEL AD) IS
! REPLACED BY (DEL + DEL MU - DEL AD).  RETAIN IMU FLAG FOR LATER USE.
!      IF(IMU.EQ.3)THEN
!  FICTITIOUS MU CURRENTS THAT OPPOSE ES CIRCULATION CALCULATED HERE.
!         DO 32 I = IMIN,IMAX
!            IF(LCZ(I).AND.LCZ(I-1))GOTO 32
! SQUARE ROOT OF TKH*,KIPPENHAHN,IAU#66,P.23,USING EQ.(12)FOR D ON P.25.
! NOTE FACTOR OF G IS SUBSUMED IN FACT4.
! LOCAL TIMESCALE ESTIMATE FOR MU INHIBITION.
! ASSUMES V/L = OMEGA.
!      ELSE IF(IMU.EQ.2)THEN
! MHP 06/02 REPLACE WITH THE ZAHN&MAEDER 1998 PRESCRIPTION
       do i = zone_min,zone_max
          qmu = log(mix_scr%amum(i))-log(mix_scr%amum(i-1))
          qp = ln10*(log_pressure(i)-log_pressure(i-1))
          ddel=rot_scr%del_grad_diff_interface(i)+ qmu/qp
          ddtest = max(rot_scr%del_grad_diff_interface(i),1.0d-3)
          ddtest2 = max(ddel,1.0d-3)
          qqq = ddtest/ddtest2
          star%es_circulation_velocity(i) = star%es_circulation_velocity(i)*qqq
          if (star%ctrl%use_diffusion_advection_transport) then
             rot_scr%es_advective_velocity(i) = rot_scr%es_advective_velocity(i)*qqq
             rot_scr%es_diffusive_velocity(i) = rot_scr%es_diffusive_velocity(i)*qqq
             rot_scr%vesd2(i) = rot_scr%vesd2(i)*qqq
             rot_scr%vesd3(i) = rot_scr%vesd3(i)*qqq
          end if
          q1 = rot_scr%difad_shear_coeff1(i)*rot_scr%omega_interface(i)**2*qqq
          q0 = rot_scr%es_shear_coeff(i)*rot_scr%omega_interface(i)**2*qqq
          dr = radius(i) - radius(i-1)
          v2 = 0.2d0*(q0+q1)*(omega(i)-omega(i-1))/dr/rot_scr%omega_interface(i)
! ADD MU GRADIENTS TO VELOCITY ESTIMATES
          rot_scr%es_relaxation_factor(i) = rot_scr%es_relaxation_factor(i)*qqq
          rot_scr%velocity_coeff0(i) = rot_scr%velocity_coeff0(i)*qqq
!         WRITE(*,911)I,WM(I),VESA(I),V2,Q0,Q1,THN(I),VES(I),QQQ
! 911           FORMAT(I5,1P8E12.3)
       end do
!         DO 33 I = IMIN,IMAX
!     *             /WM(I)**2
!   33    CONTINUE
! ALTERNATE EXPRESSION : D = D0/(1+R*DEL MU/MU)
! THIS IS ACTUALLY SIMILAR TO THE ZM98 FORMULA, THEREFORE OBSOLETE.
! SHOULD HAVE PUBLISHED IT - OH WELL.
!      ELSE
!         DO 34 I = IMIN,IMAX
!  GSF INSTABILITY.  VELOCITY DEFINED IN SR SETUPV.
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
!      IF(IGSF.EQ.0 .OR. IGSF.EQ.2 .OR. IGSF.EQ.1)THEN
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
            qwrmx=2.0d0*sqrt(rot_scr%interface_gravity_factor(i)* &
                 abs(mix_scr%amum(i)-mix_scr%amum(i-1)) &
                 /dr/rot_scr%mean_molecular_weight_interface(i))
            if (abs(dlnomega_dlnr(i)).lt.qwrmx) then
               star%gsf_circulation_velocity(i) = 0.0d0
               cycle
            else
! MHP 05/02 ADDED TESTS TO AVOID DIVIDE BY ZERO
!                  FXX2 = SQRT((ABS(QWLNR(I))-QWRMX)/QWRMX)
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
         if (omega(i)*radius(i)**2.lt.omega(i-1)*radius(i-1)**2) then
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
! OMIT JAMES AND KAHN ESTIMATE
!      ELSE
! JAMES AND KAHN (1971) ESTIMATE.
! NOTE FACTOR OF OMEGA**4, WHICH IS PART OF THE MERIDIONAL
! CIRCULATION VELOCITY ESTIMATE.
!         DO 41 I = IMIN,IMAX
!            IF(LCZ(I).AND.LCZ(I-1))GOTO 41
! MHP 8/93 STABILITY CONDITION ADDED, NEGLECTING THE EFFECTS OF
! MU GRADIENTS.
!     *                LOG(OMEGA(I-1)))/DR)
!            ELSE IF(IES.EQ.2)THEN
!  DIFFUSIVE AND DYNAMICAL SHEAR INSTABILITIES - REF. ENDAL&SOFIA PAPER II.
      do i = zone_min,zone_max
!  CHECK FOR OPERATION OF DYNAMICAL SHEAR.
!  IF DYNAMICAL SHEAR IS OPERATING,SET SECULAR SHEAR VELOCITY TO MAXIMUM
!  VALUE AND COMPUTE (LARGE) DYNAMICAL SHEAR VELOCITY.
         if (abs(dlnomega_dlnr(i)).gt.dlnomega_dlnr_max(i)) then
            qwr = abs(dlnomega_dlnr(i))
            star%secular_shear_velocity(i)=8.0d0/4.5d1* &
                 rot_scr%thermal_diffusivity_interface(i)* &
                 (qwr/dlnomega_dlnr_max(i))**2/rot_scr%interface_radius(i)
!  CORRECT GSF VELOCITY AS WELL.
            dr = radius(i)-radius(i-1)
            rmid = 0.5d0*(radius(i)+radius(i-1))
            if (omega(i).lt.omega(i-1)) then
               wmin = omega(i) + dlnomega_dlnr_max(i)*ln10* &
                    (log_radius(i)-log_radius(i-1))
               wmid = omega(i) + dlnomega_dlnr_max(i)*ln10* &
                    (log_radius(i)-log10(rmid))
            else
               wmin = omega(i) - dlnomega_dlnr_max(i)*ln10* &
                    (log_radius(i)-log_radius(i-1))
               wmid = omega(i) - dlnomega_dlnr_max(i)*ln10* &
                    (log_radius(i)-log10(rmid))
            end if
            dlnwdr0 = abs(log(omega(i))-log(omega(i-1)))/dr
            dlnjmdr0 = abs(2.0d0/rmid+(log(omega(i))- &
                         log(omega(i-1)))/dr)
            dlnjmdr = abs(2.0d0/rmid+(log(omega(i))-log(wmin))/dr)
            dlnwdr = abs(log(omega(i))-log(wmin))/dr
! MHP 05/02 ADDED FIRST BRANCH FOR IGSF=0 - SOMEHOW
! MISSED IN CODE CHANGES EARLIER!
! kippenhahn (1980) estimate
! OMIT BRANCH FOR JAMES AND KAHN.
!            IF(IGSF.EQ.2 .OR. IGSF.EQ.1 .OR. IGSF.EQ.0)THEN
            if (omega(i)*radius(i)**2.lt.wmin*radius(i-1)**2) then
               fx = max(2.0d0*dlnjmdr,0.25d0*dlnwdr)
               qwrmx_dyn = max(2.0d0*dlnjmdr0,0.25d0*dlnwdr0)
               star%gsf_circulation_velocity(i) = star%gsf_circulation_velocity(i)* &
                    fx/qwrmx_dyn
            else
               star%gsf_circulation_velocity(i)=star%gsf_circulation_velocity(i)* &
                    dlnwdr/dlnwdr0
            end if
!            ELSE IF(IGSF.EQ.3)THEN
!            ELSE IF(IES.EQ.2)THEN
            write(6,9911) i,omega(i),omega(i-1),wmin
 9911       format(1x,'DYNAMICAL SHEAR-SHELL',i5,1p,' WTOP',e11.3, &
                 ' WBOT',e11.3,' LIMIT',e11.3)
            cycle
! *** END OF CHANGED SECTION
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
            if (abs((mix_scr%amum(i)-mix_scr%amum(i-1))/rot_scr%mean_molecular_weight_interface(i)) &
                 .lt.1.0d-10) then
               qwrmx2 = 0.0d0
               qwr = abs(dlnomega_dlnr(i)) - qwrmx
               star%secular_shear_velocity(i)=8.0d0/4.5d1* &
                    rot_scr%thermal_diffusivity_interface(i)* &
                    (qwr/dlnomega_dlnr_max(i))**2/rot_scr%interface_radius(i)
            else
!  CHECK FOR EFFECTS OF MU GRADIENT.
!  RICHNO = RHO/P*(-d lnT/d lnMU)*(del MU)/(del P)*(G/QWLNR)**2
!  WHERE -d lnT/d lnMU = (CON-1)/(1+3CON),CON=(a/3)T**4/P
!  GIVEN (1 - CON)P = CGAS*RHO*T/MU
!      FACT = (RHOM/PM)*QTMU*DMU/AMUMI/DP*HGM**2
               qwrmx2 = 2.0d0*sqrt(max(1.0d-20,rot_scr%mu_gradient_richardson_coeff(i)* &
                        abs((mix_scr%amum(i)-mix_scr%amum(i-1))/ &
                        rot_scr%mean_molecular_weight_interface(i))))
               if (abs(dlnomega_dlnr(i)).gt.qwrmx2) then
!  INTERFACE UNSTABLE WITH RESPECT TO BOTH CONDITIONS; CHOOSE THE
!  MAXIMUM GRADIENT IMPLIED BY THE SECOND CONDITION IF IT'S LARGER
!  THAN THE FIRST(I.E. IF A MU GRADIENT IS SLOWING J TRANSPORT).
                  qwrmx = max(qwrmx2,qwrmx)
                  qwr = abs(dlnomega_dlnr(i)) - qwrmx
            star%secular_shear_velocity(i)=8.0d0/4.5d1* &
                 rot_scr%thermal_diffusivity_interface(i)* &
                 (qwr/dlnomega_dlnr_max(i))**2/rot_scr%interface_radius(i)
               end if
            end if
         end if
      end do
!  NOW DETERMINE WHETHER OR NOT MU GRADIENTS ARE STEEP ENOUGH TO
!  INHIBIT TRANSPORT.  MULTIPLY THE RESULTING VELOCITY ESTIMATES
!  BY THE USER DEFINED PARAMETERS FES AND FGSF.
!  IMU=3 KIPPENHAHN AND MOLLENHOF(1974)METHOD;IMU=2 LOCAL DAMPING
!  FACTOR METHOD.
! AGAIN, OMIT OBSOLETE MU GRADIENT TREATMENTS.
!     *               *FCC*VMU(I))
!     *               *FCC*VMU(I))
! MHP 05/02 ONLY DO THIS IF MU GRADIENTS NOT
! ALREADY ACCOUNTED FOR
! ALREADY INCLUDED - ONLY USE SCALE FACTOR
! MHP 8/03 MULTIPLY VELOCITY ESTIMATES BY USER PARAMETER
! SCALE FACTORS
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
                 rot_scr%interface_radius(i)*(abs(star%gsf_circulation_velocity(i))+ &
                 abs(star%secular_shear_velocity(i)))
            rot_scr%shear_diffusion_coeff(i) = rot_scr%interface_radius(i)* &
                 abs(star%secular_shear_velocity(i))
            rot_scr%gsf_diffusion_coeff(i) = rot_scr%interface_radius(i)* &
                 abs(star%gsf_circulation_velocity(i))
!            IF(VGSF(I).GT.0.0D0)THEN
! D LN W/DR
!               DR = HRU(I) - HRU(I-1)
!               QLNWR = ABS(LOG(OMEGA(I))-LOG(OMEGA(I-1))/DR)
! CEILING SET BY DYNAMICAL SHEAR
!               QLNWRMAX = ABS(QWRMAX(I)/(WM(I)*RM(I)))
! TAKE THE SMALLER OF THE TWO
         end do
      end if
!  AVERAGE PREVIOUS AND NEW VELOCITY ESTIMATES AFTER THE FIRST ITERATION.
      if (iteration.gt.1) then
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
      end if
      do i =zone_min,zone_max
         total_circulation_velocity(i) = star%gsf_circulation_velocity(i) + &
              star%es_circulation_velocity(i) + star%secular_shear_velocity(i)
         if (total_circulation_velocity(i).lt.1.0d-20) &
              total_circulation_velocity(i)=0.0d0
         if (total_circulation_velocity(i).gt.0.0d0) any_transport_active=.true.
      end do
! 9/93 MIXING WITHOUT TRANSPORT ADDED.
! ZERO OUT COEFFICIENTS IN CORE TO AVOID NUMERICAL PROBLEMS IN
! THE H-BURNING SHELL.
      if (star%ctrl%no_am_transport_in_core) then
         do i = zone_min,zone_max
            dell = rot_scr%interface_luminosity(i)/total_luminosity
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
      end if

      return
end subroutine circulation_velocities

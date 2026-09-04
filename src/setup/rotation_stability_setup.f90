!----------------------------------------------------------------------
! setupv
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setupv.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
!   THE PROCEDURE FOR CALCULATING STABILITY AGAINST ROTATIONALLY INDUCED
!   MIXING IS AS FOLLOWS:
!      STABILITY IS CALCULATED FOR THE I/I-1 INTERFACE; QUANTITIES RELATED
!      TO STABILITY AT THIS INTERFACE ARE STORED IN ARRAY ELEMENT(S) I.
!      EACH STABILITY CRITERION IS EXPRESSED AS FOLLOWS:
!         A COMBINATION OF VARIABLES(RHO,P,ETC.)STORED AT EACH SHELL
!         MULTIPLIED BY A GRADIENT BETWEEN THE SHELLS MUST BE LESS THAN
!         SOME CRITICAL NUMBER.
!  SETUPV CALCULATES THE FACTORS WHICH ARE BASED ON
!  QUANTITIES WHICH DO NOT CHANGE DURING A DIFFUSION TIMESTEP.
subroutine rotation_stability_setup(log_density, local_gravity, luminosity, log_pressure, &
     log_radius, mass_unlogged, log_temperature, transport_zone_begin, &
     transport_zone_end, num_zones, radius_unlogged, &
     dynamical_shear_omega_limit)
      use rotation_scratch_lib

      use star_info_lib, only: star, json
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: log_density(json), local_gravity(json), &
           luminosity(json), log_pressure(json), log_radius(json), &
           mass_unlogged(json), log_temperature(json)
      integer, intent(in) :: transport_zone_begin, transport_zone_end, &
           num_zones
      double precision, intent(out) :: radius_unlogged(json)
      double precision, intent(out) :: dynamical_shear_omega_limit(json)
      double precision :: specific_heat_interface(json)
! --- locals ---
      integer :: zone_idx, interior_begin, interior_end
      double precision :: grav_const
      double precision :: fac_local
      double precision :: dlnmu_dlnp, ddel_floor, ddel
      double precision :: pressure_scale_factor, temp_scale_factor, &
           eta_factor, ff_factor, specific_luminosity, cpigi_const
      double precision :: local_flux_factor
      double precision :: ht_temp_scale_prev, ht_temp_scale
      double precision :: mean_dlneps_dlnt, mean_neutrino_fraction
      double precision :: c1_factor, qchit, f1_local

!  FIND UNLOGGED RADII OF THE MODEL POINTS.
      do zone_idx = 1,num_zones
         radius_unlogged(zone_idx) = exp(ln10*log_radius(zone_idx))
      end do
!  FIND LOCATION (IN RADIUS) OF THE MIDPOINTS OF THE INTERFACES.
      do zone_idx = 2,num_zones
         rot_scr%interface_radius(zone_idx) = 0.5d0*(radius_unlogged(zone_idx) + &
              radius_unlogged(zone_idx-1))
      end do
!  FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS FOR ALL INTERFACES.
!  INTERFACE BETWEEN FIRST AND SECOND POINT USES THE FIRST 4 POINTS FOR
!  INTERPOLATION; BETWEEN LAST AND SECOND TO LAST USES THE LAST FOUR POINTS
!  AND BETWEEN I AND I-1 USES I-2,I-1,I,I+1.
!  TREATMENT OF FIRST INTERFACE.
      if(transport_zone_begin.lt.3)then
         call lagrange_weights(1, 2)
         interior_begin = 3
         call interpolate_to_interface(1, 2)
      else
         interior_begin = transport_zone_begin
      endif
!  TREATMENT OF LAST INTERFACE.
      if(transport_zone_end.eq.num_zones)then
         call lagrange_weights(num_zones-3, num_zones)
         interior_end = num_zones-1
         call interpolate_to_interface(num_zones-3, num_zones)
      else
         interior_end = transport_zone_end
      endif
!  COMPUTE INTERPOLATION FACTORS FOR ALL OTHER POINTS.
      do zone_idx = interior_begin,interior_end
         call lagrange_weights(zone_idx-2, zone_idx)
      end do
      grav_const = exp(ln10*cgl)
      do zone_idx = interior_begin,interior_end
!  USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND PHYSICAL VARIABLES
!  AT THE INTERFACES.
         call interpolate_to_interface(zone_idx-2, zone_idx)
      end do
      do zone_idx = transport_zone_begin,transport_zone_end
         rot_scr%delmi(zone_idx) = min(rot_scr%delmi(zone_idx),rot_scr%delami(zone_idx))
      end do
      do zone_idx = 1,num_zones
         rot_scr%es_velocity_coeff1(zone_idx) = 0.0d0
         rot_scr%es_velocity_coeff2(zone_idx) = 0.0d0
         rot_scr%es_shear_coeff(zone_idx) = 0.0d0
         rot_scr%fgsfj(zone_idx) = 0.0d0
         rot_scr%gsf_kippenhahn_coeff(zone_idx) = 0.0d0
         rot_scr%fact1(zone_idx) = 0.0d0
         rot_scr%fact2(zone_idx) = 0.0d0
         rot_scr%fact3(zone_idx) = 0.0d0
         rot_scr%fact4(zone_idx) = 0.0d0
         rot_scr%mu_gradient_richardson_coeff(zone_idx) = 0.0d0
! MHP 06/02
         rot_scr%es_relaxation_factor(zone_idx) = 0.0d0
         rot_scr%theta_mean(zone_idx) = 0.0d0
         rot_scr%difad_shear_coeff1(zone_idx) = 0.0d0
         rot_scr%difad_shear_coeff2(zone_idx) = 0.0d0
      end do
      do zone_idx = 2,num_zones
! ADDED FOR D THETA/DT TERM FROM ZAHN&MAEDER 1998
         fac_local = 2.0d0*(log_radius(zone_idx)+log_radius(zone_idx-1))-0.5d0* &
              (log10(mass_unlogged(zone_idx))+log10(mass_unlogged(zone_idx-1)))-cgl
         rot_scr%theta_mean(zone_idx) = exp(ln10*fac_local)
      end do
!  NOW COMPUTE STRUCTURAL QUANTITIES NEEDED TO EVALUATE VELOCITIES AT
!  ALL INTERFACES.
      cpigi_const = 4.0d0/c4pi/grav_const
      do zone_idx = transport_zone_begin,transport_zone_end
         if(.not.star%ctrl%use_diffusion_advection_transport)then
            dlnmu_dlnp = (log10(mix_scr%amum(zone_idx))-log10(mix_scr%amum(zone_idx-1)))/ &
                 (log_pressure(zone_idx)-log_pressure(zone_idx-1))
         else
            dlnmu_dlnp = 0.0d0
         endif
         if(zone_idx.eq.transport_zone_begin)then
            ddel_floor=max(1.0d-3,0.5d0*(mix_scr%del_adiabatic_mix(transport_zone_begin+1)- &
                 mix_scr%delm(transport_zone_begin+1)) &
                 +dlnmu_dlnp)
         else if(zone_idx.eq.transport_zone_end)then
            ddel_floor=max(1.0d-3,0.5d0*(mix_scr%del_adiabatic_mix(transport_zone_end-1)- &
                 mix_scr%delm(transport_zone_end-1)) &
                 +dlnmu_dlnp)
         else
            ddel_floor = 1.0d-3
         endif
         ddel = max(rot_scr%delami(zone_idx)-rot_scr%delmi(zone_idx),ddel_floor)
! MHP 06/02
! ADDED FOR ALTERNATE TREATMENT OF MU GRADIENTS
!         DDELM(I) = 0.5D0*(DELAMI(I)+DELAMI(I-1)-
!     *                     DELMI(I)-DELMI(I-1))
!         DDELM(I) = DDEL
         rot_scr%del_grad_diff_interface(zone_idx) = max(rot_scr%delami(zone_idx)-rot_scr%delmi(zone_idx)+dlnmu_dlnp,ddel_floor)
!         FESTIME(I) = PM(I)/(HGM(I)*DDEL*DM(I)*TM(I))
         rot_scr%es_relaxation_factor(zone_idx) = rot_scr%pm(zone_idx)/(rot_scr%interface_gravity_factor(zone_idx)*ddel*rot_scr%dm(zone_idx))
         pressure_scale_factor = rot_scr%pm(zone_idx)*rot_scr%interface_radius(zone_idx)**2/rot_scr%dm(zone_idx)/rot_scr%hs3(zone_idx)/grav_const
         temp_scale_factor = pressure_scale_factor/rot_scr%delmi(zone_idx)
         eta_factor = 2.0d0*cc23*rot_scr%interface_radius(zone_idx)**3/grav_const/rot_scr%hs3(zone_idx)
!        TTHERM = 8.0d0*CC23*CSIG*TM(I)**3/OPM(I)/DM(I)**2/CPM(I)
         rot_scr%gsf_kippenhahn_coeff(zone_idx)= 8.0d0*pressure_scale_factor*eta_factor/ddel
         ff_factor = rot_scr%pm(zone_idx)/(rot_scr%interface_gravity_factor(zone_idx)*ddel*specific_heat_interface(zone_idx)* &
              rot_scr%dm(zone_idx)*rot_scr%tm(zone_idx))
         fac_local = 2.0d0*eta_factor*ff_factor
         specific_luminosity = rot_scr%interface_luminosity(zone_idx)/rot_scr%hs3(zone_idx)
         rot_scr%es_velocity_coeff1(zone_idx) = fac_local*(specific_luminosity-rot_scr%epsilm(zone_idx))
         rot_scr%es_velocity_coeff2(zone_idx)= -0.5d0*fac_local*specific_luminosity*cpigi_const/rot_scr%dm(zone_idx)
         rot_scr%es_shear_coeff(zone_idx) = specific_luminosity*ff_factor*cc13*cpigi_const*rot_scr%interface_radius(zone_idx)/rot_scr%dm(zone_idx)
!         FES3(I) = EM*FF*CPIGI*(3.0D0*HTSC-RM(I))/DM(I)
         rot_scr%velocity_coeff0(zone_idx) = ff_factor
         rot_scr%velocity_coeff1a(zone_idx) = rot_scr%es_velocity_coeff1(zone_idx)/ff_factor
         rot_scr%velocity_coeff1b(zone_idx) = rot_scr%es_velocity_coeff2(zone_idx)/ff_factor
         rot_scr%velocity_coeff2a(zone_idx) = specific_luminosity*cpigi_const*3.0d0*temp_scale_factor/rot_scr%dm(zone_idx)
         rot_scr%velocity_coeff2b(zone_idx) = -specific_luminosity*cpigi_const*rot_scr%interface_radius(zone_idx)/rot_scr%dm(zone_idx)
         rot_scr%fgsfj(zone_idx) = abs(fac_local*specific_luminosity*temp_scale_factor)*rot_scr%interface_radius(zone_idx)
!  EDDINGTON CIRCULATION VELOCITY IS DEFINED AS
!  VES = FACT1*FACT2*OMEGA**2 (ENDAL AND SOFIA PAPER II).
!  NOTE THAT TO AVOID OVERFLOW 1/(DEL(AD)-DEL)IS SET TO A MAXIMUM OF 10^6.
!  THIS SHOULD ONLY BE AN ISSUE IF YOU HAVE ZERO OVERSHOOT AT THE BOUNDARY
!  OF A CONVECTION ZONE.  CAVEAT EMPTOR.
         rot_scr%fact1(zone_idx) = 1.0d0/rot_scr%qdtmi(zone_idx)/max(rot_scr%delami(zone_idx)-rot_scr%delmi(zone_idx),1.0d-3)
!         FACT2(I) = DELAMI(I)*ALM(I)*(RM(I)**3/CG2/HS3(I)**2)*
!     *              (2.0D0*RM(I)**2*(EPSILM(I)/ALM(I) - 1.0D0/HS3(I))
!    *              - 3.0D0/(C4PI*DM(I)*RM(I)))
         fac_local = rot_scr%delami(zone_idx)/rot_scr%interface_gravity_factor(zone_idx)
         rot_scr%fact2(zone_idx) = fac_local*eta_factor*(rot_scr%epsilm(zone_idx)-specific_luminosity)
         rot_scr%fact6(zone_idx) = -fac_local*specific_luminosity/c4pi/grav_const/rot_scr%dm(zone_idx)/rot_scr%interface_radius(zone_idx)
! GSF VELOCITY IS DEFINED AS
! VGSF = VES*FACT3/Hj**2/(2OMEGA/(D OMEGA/D LNR) + 1)
! WHERE HJ IS THE SPECIFIC ANGULAR MOMENTUM SCALE HEIGHT.
! AND FACT3 IS THE TEMPERATURE SCALE HEIGHT TIMES THE RADIUS.
!  HT = (DEL*D(LN P)/DR)**-1 = P*R**2/(DEL*RHO*GM)
!        FACT3(I) = PM(I)*RM(I)**3/(DELMI(I)*DM(I)*CG*HS3(I))
         rot_scr%fact3(zone_idx) = abs(temp_scale_factor*rot_scr%interface_radius(zone_idx))
!  FICTITIOUS MU CURRENTS THAT OPPOSE ES AND GSF CIRCULATION CALCULATED HERE.
!  VMU = FMU*FACT1*FACT4*mu*del(mu), FACT4=(D RHO/D MU)*PSCALH/TKH/MU
!  (KIPPENHAHN, IAU#66,P.22) USING TKH* = THE LOCAL KELVIN-HELMHOLTZ TIMESCALE.
!  FMU IS A USER PARAMETER, DMU IS DELTA MU BETWEEN THE SHELLS, AND QDMU IS
!  THE DERIVATIVE OF THE DENSITY WITH RESPECT TO MU.
!  IF VMU EXCEEDS VES AND VGSF THE INTERFACE IS STABLE.
!  D RHO/DMU = (1-(a/3)*T**4/P)*P*MU/(RHO*T*CGAS)
!  THE PRESSURE SCALE HEIGHT IS CGAS*T*R**2/GM (EXCLUDING A FACTOR MU
!  WHICH CANCELS WITH ANOTHER FACTOR OF MU WHICH ARISES ELSEWHERE) AND
!  THE KELVIN-HELMHOLTZ TIMESCALE IS G*M**2/(R*L).
!         FACT4(I) = ALM(I)*(RM(I)/HS3(I))**3/CG2/DM(I)*
!     *              (PM(I)-CA3*TM(I)**4)
         rot_scr%fact4(zone_idx) = rot_scr%interface_radius(zone_idx)**2/rot_scr%hs3(zone_idx)/grav_const/rot_scr%dm(zone_idx)* &
              (rot_scr%pm(zone_idx)-radiation_constant_over_3*rot_scr%tm(zone_idx)**4)/sqrt(rot_scr%interface_gravity_factor(zone_idx))
!  DYNAMICAL SHEAR INSTABILITY.
!  AN INTERFACE IS STABLE AGAINST THIS SHEAR WHENEVER THE RICHARDSON NUMBER
!  RICHNO = (RHO/P)*(DEL(AD)-DEL)*g**2/(D OMEGA/D LN R)**2>.25.
         dynamical_shear_omega_limit(zone_idx)=sqrt(rot_scr%dm(zone_idx)/rot_scr%pm(zone_idx)* &
              max(rot_scr%delami(zone_idx) - rot_scr%delmi(zone_idx),1.0d-3)) &
              *2.0d0*rot_scr%interface_gravity_factor(zone_idx)
!  THE DIFFUSIVE SHEAR IS INHIBITED BY MU GRADIENTS.
!  AN INTERFACE IS STABLE AGAINST THE DIFFUSIVE SHEAR IF
!  RICHARDSON NUMBER = FACT5*DMU/MU/(D OMEGA/D LNR)**2 > .25, WHERE
!  FACT5 = RHO/P*(-d lnT/d lnMU)/(del LN P)*G**2 AND
!   -d lnT/d lnMU = (CON-1)/(1+3CON),CON=(a/3)T**4/P
!  GIVEN (1 - CON)P = CGAS*RHO*T/MU
!  1 - CON IS BETA(CORRECTION FOR RADIATION PRESSURE).
         rot_scr%mu_gradient_richardson_coeff(zone_idx) = rot_scr%dm(zone_idx)*rot_scr%interface_gravity_factor(zone_idx)**2/ &
              rot_scr%pm(zone_idx)/ln10/(log_pressure(zone_idx)-log_pressure(zone_idx-1))* &
              (radiation_constant_over_3*rot_scr%tm(zone_idx)**4/rot_scr%pm(zone_idx)-1.0d0)/ &
              (1.0d0 + 3.0d0*radiation_constant_over_3*rot_scr%tm(zone_idx)**4/rot_scr%pm(zone_idx))
! MHP 3/92 ADD VECTOR FOR LOCAL KELVIN-HELMHOLTZ TIME SCALE.
!         TKH(I) = CG*HS3(I)**2/ABS(ALM(I))/RM(I)
      end do
! MHP 06/02 ADDED TERMS OF ORDER DW/DR FROM ZAH&MAEDER 1998
      if(star%ctrl%use_diffusion_advection_transport)then
         ht_temp_scale_prev = exp(ln10*(log_pressure(transport_zone_begin-1)+ &
              2.0d0*log_radius(transport_zone_begin-1)-log_density(transport_zone_begin-1))) &
              /mass_unlogged(transport_zone_begin-1)/grav_const/mix_scr%del_radiative_mix(transport_zone_begin-1)
         rot_scr%third_deriv_geom_factor(transport_zone_begin-1) = ht_temp_scale_prev
         do zone_idx = transport_zone_begin,transport_zone_end
!         DDEL = MAX(DELAMI(I)-DELMI(I),1.0D-3)
         pressure_scale_factor = rot_scr%pm(zone_idx)*rot_scr%interface_radius(zone_idx)**2/rot_scr%dm(zone_idx)/rot_scr%hs3(zone_idx)/grav_const
         temp_scale_factor = pressure_scale_factor/rot_scr%delmi(zone_idx)
         f1_local = rot_scr%pm(zone_idx)/(rot_scr%interface_gravity_factor(zone_idx)*rot_scr%del_grad_diff_interface(zone_idx)* &
              specific_heat_interface(zone_idx)*rot_scr%dm(zone_idx)*rot_scr%tm(zone_idx))
         c1_factor = cc23*rot_scr%interface_radius(zone_idx)**4/grav_const/rot_scr%hs3(zone_idx)
         rot_scr%second_deriv_geom_factor(zone_idx) = c1_factor
! D LN CHI/D LN T = 3 - D LN CAPPA/D LN T
         qchit = 3.0d0 - 0.5d0*(rot_scr%dlnkappa_dlnt(zone_idx)+rot_scr%dlnkappa_dlnt(zone_idx-1))
         ht_temp_scale = exp(ln10*(log_pressure(zone_idx)+2.0d0*log_radius(zone_idx)-log_density(zone_idx)))/ &
              mass_unlogged(zone_idx)/grav_const/mix_scr%delm(zone_idx)
         rot_scr%third_deriv_geom_factor(zone_idx) = ht_temp_scale
         mean_dlneps_dlnt = 0.5d0*(rot_scr%dlnepsilon_dlnt(zone_idx)+rot_scr%dlnepsilon_dlnt(zone_idx-1))
         mean_neutrino_fraction = 0.5d0*(rot_scr%neutrino_loss_fraction(zone_idx)+rot_scr%neutrino_loss_fraction(zone_idx-1))
         rot_scr%difad_shear_coeff1(zone_idx)= &
         - f1_local*c1_factor*(mean_neutrino_fraction*mean_dlneps_dlnt + &
              rot_scr%epsilm(zone_idx)*(1.0d0 - mean_neutrino_fraction - qchit))
         rot_scr%velocity_coeff2b(zone_idx) = rot_scr%velocity_coeff2b(zone_idx)+rot_scr%difad_shear_coeff1(zone_idx)/f1_local
         rot_scr%difad_shear_coeff2(zone_idx) = 0.0d0
         local_flux_factor = rot_scr%interface_luminosity(zone_idx)/c4pi/rot_scr%dm(zone_idx)/rot_scr%interface_radius(zone_idx)**2
         rot_scr%facd2(zone_idx) = f1_local*(local_flux_factor*qchit - rot_scr%epsilm(zone_idx)*temp_scale_factor)
         rot_scr%facd3(zone_idx) = -f1_local*local_flux_factor
         end do
      endif
      return

contains

!----------------------------------------------------------------------
! lagrange_weights
!----------------------------------------------------------------------
!  4-point Lagrangian interpolation weights at the interface
!  interface_idx (radius rot_scr%interface_radius(interface_idx)) from
!  the model points first_point..first_point+3, stored in
!  rot_scr%lagrange_interp_weights(1:4,interface_idx).
subroutine lagrange_weights(first_point, interface_idx)
      integer, intent(in) :: first_point, interface_idx
      double precision :: dr43, dr42, dr41, dr32, dr31, dr21
      double precision :: lag_denom1, lag_denom2, lag_denom3, lag_denom4
      double precision :: lag_x1, lag_x2, lag_x3, lag_x4

      dr43 = radius_unlogged(first_point+3) - radius_unlogged(first_point+2)
      dr42 = radius_unlogged(first_point+3) - radius_unlogged(first_point+1)
      dr41 = radius_unlogged(first_point+3) - radius_unlogged(first_point)
      dr32 = radius_unlogged(first_point+2) - radius_unlogged(first_point+1)
      dr31 = radius_unlogged(first_point+2) - radius_unlogged(first_point)
      dr21 = radius_unlogged(first_point+1) - radius_unlogged(first_point)
      lag_denom1 = -dr21*dr31*dr41
      lag_denom2 = dr21*dr32*dr42
      lag_denom3 = -dr31*dr32*dr43
      lag_denom4 = dr41*dr42*dr43
      lag_x1 = rot_scr%interface_radius(interface_idx) - radius_unlogged(first_point)
      lag_x2 = rot_scr%interface_radius(interface_idx) - radius_unlogged(first_point+1)
      lag_x3 = rot_scr%interface_radius(interface_idx) - radius_unlogged(first_point+2)
      lag_x4 = rot_scr%interface_radius(interface_idx) - radius_unlogged(first_point+3)
      rot_scr%lagrange_interp_weights(1,interface_idx) = (lag_x2*lag_x3*lag_x4)/lag_denom1
      rot_scr%lagrange_interp_weights(2,interface_idx) = (lag_x1*lag_x3*lag_x4)/lag_denom2
      rot_scr%lagrange_interp_weights(3,interface_idx) = (lag_x1*lag_x2*lag_x4)/lag_denom3
      rot_scr%lagrange_interp_weights(4,interface_idx) = (lag_x1*lag_x2*lag_x3)/lag_denom4
end subroutine lagrange_weights

!----------------------------------------------------------------------
! interp4
!----------------------------------------------------------------------
!  field interpolated to interface interface_idx with the weights
!  computed by lagrange_weights from points first_point..first_point+3.
double precision function interp4(field, first_point, interface_idx)
      double precision, intent(in) :: field(json)
      integer, intent(in) :: first_point, interface_idx
      interp4 = field(first_point)*rot_scr%lagrange_interp_weights(1,interface_idx)+ &
           field(first_point+1)*rot_scr%lagrange_interp_weights(2,interface_idx)+ &
           field(first_point+2)*rot_scr%lagrange_interp_weights(3,interface_idx)+ &
           field(first_point+3)*rot_scr%lagrange_interp_weights(4,interface_idx)
end function interp4

!----------------------------------------------------------------------
! interpolate_to_interface
!----------------------------------------------------------------------
!  Interpolate the structure variables to interface interface_idx from
!  model points first_point..first_point+3 (weights from
!  lagrange_weights). P, T and rho are interpolated in log10 and
!  unlogged; the luminosity is converted to cgs.
subroutine interpolate_to_interface(first_point, interface_idx)
      integer, intent(in) :: first_point, interface_idx

!  PRESSURE.
      rot_scr%pm(interface_idx)=exp(ln10*interp4(log_pressure, first_point, interface_idx))
!  TEMPERATURE.
      rot_scr%tm(interface_idx)=exp(ln10*interp4(log_temperature, first_point, interface_idx))
!  DENSITY.
      rot_scr%dm(interface_idx)=exp(ln10*interp4(log_density, first_point, interface_idx))
!  DEL (ACTUAL).
!  DEL (RADIATIVE) IS INTERPOLATED, AND DEL IS THE MIN OF DELA,DELR
!  (taken by the caller).
      rot_scr%delmi(interface_idx)=interp4(mix_scr%del_radiative_mix, first_point, interface_idx)
!  DEL(ADIABATIC).
      rot_scr%delami(interface_idx)=interp4(mix_scr%del_adiabatic_mix, first_point, interface_idx)
!  D LN RHO/D LN T.
      rot_scr%qdtmi(interface_idx)=interp4(mix_scr%qdtm, first_point, interface_idx)
!  UNLOGGED MASS INTERIOR TO THE INTERFACE.
      rot_scr%hs3(interface_idx)=interp4(mass_unlogged, first_point, interface_idx)
!  SPECIFIC ENERGY GENERATION RATE.
      rot_scr%epsilm(interface_idx)=interp4(mix_scr%esumm, first_point, interface_idx)
!  LUMINOSITY.
      rot_scr%interface_luminosity(interface_idx)=star%solar_luminosity_cgs*(interp4(luminosity, first_point, interface_idx))
!  LOCAL AVERAGE FORCE OF GRAVITY.
      rot_scr%interface_gravity_factor(interface_idx)=interp4(local_gravity, first_point, interface_idx)
!  specific heat
      specific_heat_interface(interface_idx)=interp4(mix_scr%cpm, first_point, interface_idx)
end subroutine interpolate_to_interface

end subroutine rotation_stability_setup

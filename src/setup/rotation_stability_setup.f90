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














! JvS 09/25 CHANGED CPM --> CPMI TO AVOID CONFLICT IN MDPHY
      double precision :: specific_heat_interface(json), opacity_interface(json)
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
!
! --- locals ---
      integer :: zone_idx, interior_begin, interior_end
      double precision :: dr43, dr42, dr41, dr32, dr31, dr21
      double precision :: lag_denom1, lag_denom2, lag_denom3, lag_denom4
      double precision :: lag_x1, lag_x2, lag_x3, lag_x4
      double precision :: grav_const, grav_const_sq
      double precision :: fac_local
      double precision :: dlnmu_dlnp, ddel_floor, ddel
      double precision :: pressure_scale_factor, temp_scale_factor, &
           eta_factor, ff_factor, specific_luminosity, cpigi_const
      double precision :: local_flux_factor
      double precision :: ht_temp_scale_prev, ht_temp_scale, ht_temp_scale2
      double precision :: dhtscale_dr, mean_dlneps_dlnt, mean_neutrino_fraction
      double precision :: c1_factor, qc1r, qqc1rr, qchit, qqchitr, dr_local, &
           qdr_local, f1_local, f2_local, f3_local, v0_local

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
         dr43 = radius_unlogged(4) - radius_unlogged(3)
         dr42 = radius_unlogged(4) - radius_unlogged(2)
         dr41 = radius_unlogged(4) - radius_unlogged(1)
         dr32 = radius_unlogged(3) - radius_unlogged(2)
         dr31 = radius_unlogged(3) - radius_unlogged(1)
         dr21 = radius_unlogged(2) - radius_unlogged(1)
         lag_denom1 = -dr21*dr31*dr41
         lag_denom2 = dr21*dr32*dr42
         lag_denom3 = -dr31*dr32*dr43
         lag_denom4 = dr41*dr42*dr43
         lag_x1 = rot_scr%interface_radius(2) - radius_unlogged(1)
         lag_x2 = rot_scr%interface_radius(2) - radius_unlogged(2)
         lag_x3 = rot_scr%interface_radius(2) - radius_unlogged(3)
         lag_x4 = rot_scr%interface_radius(2) - radius_unlogged(4)
         rot_scr%lagrange_interp_weights(1,2) = (lag_x2*lag_x3*lag_x4)/lag_denom1
         rot_scr%lagrange_interp_weights(2,2) = (lag_x1*lag_x3*lag_x4)/lag_denom2
         rot_scr%lagrange_interp_weights(3,2) = (lag_x1*lag_x2*lag_x4)/lag_denom3
         rot_scr%lagrange_interp_weights(4,2) = (lag_x1*lag_x2*lag_x3)/lag_denom4
         interior_begin = 3
         rot_scr%pm(2)=exp(ln10*(log_pressure(1)*rot_scr%lagrange_interp_weights(1,2)+ &
              log_pressure(2)*rot_scr%lagrange_interp_weights(2,2) &
              +log_pressure(3)*rot_scr%lagrange_interp_weights(3,2)+ &
              log_pressure(4)*rot_scr%lagrange_interp_weights(4,2)))
         rot_scr%tm(2)=exp(ln10*(log_temperature(1)*rot_scr%lagrange_interp_weights(1,2)+ &
              log_temperature(2)*rot_scr%lagrange_interp_weights(2,2) &
              +log_temperature(3)*rot_scr%lagrange_interp_weights(3,2)+ &
              log_temperature(4)*rot_scr%lagrange_interp_weights(4,2)))
         rot_scr%dm(2)=exp(ln10*(log_density(1)*rot_scr%lagrange_interp_weights(1,2)+ &
              log_density(2)*rot_scr%lagrange_interp_weights(2,2) &
              +log_density(3)*rot_scr%lagrange_interp_weights(3,2)+ &
              log_density(4)*rot_scr%lagrange_interp_weights(4,2)))
         rot_scr%delmi(2)=mix_scr%del_radiative_mix(1)*rot_scr%lagrange_interp_weights(1,2)+ &
              mix_scr%del_radiative_mix(2)*rot_scr%lagrange_interp_weights(2,2)+ &
              mix_scr%del_radiative_mix(3)*rot_scr%lagrange_interp_weights(3,2)+ &
              mix_scr%del_radiative_mix(4)*rot_scr%lagrange_interp_weights(4,2)
         rot_scr%delami(2)=mix_scr%del_adiabatic_mix(1)*rot_scr%lagrange_interp_weights(1,2)+ &
              mix_scr%del_adiabatic_mix(2)*rot_scr%lagrange_interp_weights(2,2)+ &
              mix_scr%del_adiabatic_mix(3)*rot_scr%lagrange_interp_weights(3,2)+ &
              mix_scr%del_adiabatic_mix(4)*rot_scr%lagrange_interp_weights(4,2)
         rot_scr%qdtmi(2)=mix_scr%qdtm(1)*rot_scr%lagrange_interp_weights(1,2)+mix_scr%qdtm(2)*rot_scr%lagrange_interp_weights(2,2)+ &
              mix_scr%qdtm(3)*rot_scr%lagrange_interp_weights(3,2)+mix_scr%qdtm(4)*rot_scr%lagrange_interp_weights(4,2)
         rot_scr%hs3(2)=mass_unlogged(1)*rot_scr%lagrange_interp_weights(1,2)+mass_unlogged(2)*rot_scr%lagrange_interp_weights(2,2)+ &
              mass_unlogged(3)*rot_scr%lagrange_interp_weights(3,2)+mass_unlogged(4)*rot_scr%lagrange_interp_weights(4,2)
         rot_scr%epsilm(2)=mix_scr%esumm(1)*rot_scr%lagrange_interp_weights(1,2)+mix_scr%esumm(2)*rot_scr%lagrange_interp_weights(2,2)+ &
              mix_scr%esumm(3)*rot_scr%lagrange_interp_weights(3,2)+mix_scr%esumm(4)*rot_scr%lagrange_interp_weights(4,2)
         rot_scr%interface_luminosity(2)=star%solar_luminosity_cgs*(luminosity(1)*rot_scr%lagrange_interp_weights(1,2)+ &
              luminosity(2)*rot_scr%lagrange_interp_weights(2,2)+ &
              luminosity(3)*rot_scr%lagrange_interp_weights(3,2)+luminosity(4)*rot_scr%lagrange_interp_weights(4,2))
         rot_scr%interface_gravity_factor(2)=local_gravity(1)*rot_scr%lagrange_interp_weights(1,2)+ &
              local_gravity(2)*rot_scr%lagrange_interp_weights(2,2)+ &
              local_gravity(3)*rot_scr%lagrange_interp_weights(3,2)+local_gravity(4)*rot_scr%lagrange_interp_weights(4,2)
!  opacity.
         opacity_interface(2)=mix_scr%om(1)*rot_scr%lagrange_interp_weights(1,2)+mix_scr%om(2)*rot_scr%lagrange_interp_weights(2,2)+ &
              mix_scr%om(3)*rot_scr%lagrange_interp_weights(3,2)+mix_scr%om(4)*rot_scr%lagrange_interp_weights(4,2)
!  specific heat
         specific_heat_interface(2)=mix_scr%cpm(1)*rot_scr%lagrange_interp_weights(1,2)+mix_scr%cpm(2)*rot_scr%lagrange_interp_weights(2,2)+ &
              mix_scr%cpm(3)*rot_scr%lagrange_interp_weights(3,2)+mix_scr%cpm(4)*rot_scr%lagrange_interp_weights(4,2)
      else
         interior_begin = transport_zone_begin
      endif
!  TREATMENT OF LAST INTERFACE.
      if(transport_zone_end.eq.num_zones)then
         dr43 = radius_unlogged(num_zones) - radius_unlogged(num_zones-1)
         dr42 = radius_unlogged(num_zones) - radius_unlogged(num_zones-2)
         dr41 = radius_unlogged(num_zones) - radius_unlogged(num_zones-3)
         dr32 = radius_unlogged(num_zones-1) - radius_unlogged(num_zones-2)
         dr31 = radius_unlogged(num_zones-1) - radius_unlogged(num_zones-3)
         dr21 = radius_unlogged(num_zones-2) - radius_unlogged(num_zones-3)
         lag_denom1 = -dr21*dr31*dr41
         lag_denom2 = dr21*dr32*dr42
         lag_denom3 = -dr31*dr32*dr43
         lag_denom4 = dr41*dr42*dr43
         lag_x1 = rot_scr%interface_radius(num_zones) - radius_unlogged(num_zones-3)
         lag_x2 = rot_scr%interface_radius(num_zones) - radius_unlogged(num_zones-2)
         lag_x3 = rot_scr%interface_radius(num_zones) - radius_unlogged(num_zones-1)
         lag_x4 = rot_scr%interface_radius(num_zones) - radius_unlogged(num_zones)
         rot_scr%lagrange_interp_weights(1,num_zones) = (lag_x2*lag_x3*lag_x4)/lag_denom1
         rot_scr%lagrange_interp_weights(2,num_zones) = (lag_x1*lag_x3*lag_x4)/lag_denom2
         rot_scr%lagrange_interp_weights(3,num_zones) = (lag_x1*lag_x2*lag_x4)/lag_denom3
         rot_scr%lagrange_interp_weights(4,num_zones) = (lag_x1*lag_x2*lag_x3)/lag_denom4
         interior_end = num_zones-1
         rot_scr%pm(num_zones)=exp(ln10*(log_pressure(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              log_pressure(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones) &
              +log_pressure(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              log_pressure(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones)))
         rot_scr%tm(num_zones)=exp(ln10*(log_temperature(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              log_temperature(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones) &
              +log_temperature(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              log_temperature(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones)))
         rot_scr%dm(num_zones)=exp(ln10*(log_density(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              log_density(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones) &
              +log_density(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              log_density(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones)))
         rot_scr%delmi(num_zones)=mix_scr%del_radiative_mix(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              mix_scr%del_radiative_mix(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones)+ &
              mix_scr%del_radiative_mix(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              mix_scr%del_radiative_mix(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones)
         rot_scr%delami(num_zones)=mix_scr%del_adiabatic_mix(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              mix_scr%del_adiabatic_mix(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones)+ &
              mix_scr%del_adiabatic_mix(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              mix_scr%del_adiabatic_mix(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones)
         rot_scr%qdtmi(num_zones)=mix_scr%qdtm(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              mix_scr%qdtm(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones)+ &
              mix_scr%qdtm(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              mix_scr%qdtm(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones)
         rot_scr%hs3(num_zones)=mass_unlogged(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              mass_unlogged(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones)+ &
              mass_unlogged(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              mass_unlogged(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones)
         rot_scr%epsilm(num_zones)=mix_scr%esumm(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              mix_scr%esumm(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones)+ &
              mix_scr%esumm(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              mix_scr%esumm(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones)
         rot_scr%interface_luminosity(num_zones)=star%solar_luminosity_cgs*(luminosity(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              luminosity(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones)+ &
              luminosity(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              luminosity(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones))
         rot_scr%interface_gravity_factor(num_zones)=local_gravity(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              local_gravity(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones)+ &
              local_gravity(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              local_gravity(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones)
!  opacity.
         opacity_interface(num_zones)=mix_scr%om(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              mix_scr%om(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones)+ &
              mix_scr%om(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              mix_scr%om(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones)
!  specific heat
         specific_heat_interface(num_zones)=mix_scr%cpm(num_zones-3)*rot_scr%lagrange_interp_weights(1,num_zones)+ &
              mix_scr%cpm(num_zones-2)*rot_scr%lagrange_interp_weights(2,num_zones)+ &
              mix_scr%cpm(num_zones-1)*rot_scr%lagrange_interp_weights(3,num_zones)+ &
              mix_scr%cpm(num_zones)*rot_scr%lagrange_interp_weights(4,num_zones)
      else
         interior_end = transport_zone_end
      endif
!  COMPUTE INTERPOLATION FACTORS FOR ALL OTHER POINTS.
      do zone_idx = interior_begin,interior_end
         dr43 = radius_unlogged(zone_idx+1) - radius_unlogged(zone_idx)
         dr42 = radius_unlogged(zone_idx+1) - radius_unlogged(zone_idx-1)
         dr41 = radius_unlogged(zone_idx+1) - radius_unlogged(zone_idx-2)
         dr32 = radius_unlogged(zone_idx) - radius_unlogged(zone_idx-1)
         dr31 = radius_unlogged(zone_idx) - radius_unlogged(zone_idx-2)
         dr21 = radius_unlogged(zone_idx-1) - radius_unlogged(zone_idx-2)
         lag_denom1 = -dr21*dr31*dr41
         lag_denom2 = dr21*dr32*dr42
         lag_denom3 = -dr31*dr32*dr43
         lag_denom4 = dr41*dr42*dr43
         lag_x1 = rot_scr%interface_radius(zone_idx) - radius_unlogged(zone_idx-2)
         lag_x2 = rot_scr%interface_radius(zone_idx) - radius_unlogged(zone_idx-1)
         lag_x3 = rot_scr%interface_radius(zone_idx) - radius_unlogged(zone_idx)
         lag_x4 = rot_scr%interface_radius(zone_idx) - radius_unlogged(zone_idx+1)
         rot_scr%lagrange_interp_weights(1,zone_idx) = (lag_x2*lag_x3*lag_x4)/lag_denom1
         rot_scr%lagrange_interp_weights(2,zone_idx) = (lag_x1*lag_x3*lag_x4)/lag_denom2
         rot_scr%lagrange_interp_weights(3,zone_idx) = (lag_x1*lag_x2*lag_x4)/lag_denom3
         rot_scr%lagrange_interp_weights(4,zone_idx) = (lag_x1*lag_x2*lag_x3)/lag_denom4
      end do
      grav_const = exp(ln10*cgl)
      grav_const_sq = grav_const**2
      do zone_idx = interior_begin,interior_end
!  USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND PHYSICAL VARIABLES
!  AT THE INTERFACES.
!  PRESSURE.
         rot_scr%pm(zone_idx)=exp(ln10*(log_pressure(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              log_pressure(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx) &
              +log_pressure(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              log_pressure(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx)))
!  TEMPERATURE.
         rot_scr%tm(zone_idx)=exp(ln10*(log_temperature(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              log_temperature(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx) &
              +log_temperature(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              log_temperature(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx)))
!  DENSITY.
         rot_scr%dm(zone_idx)=exp(ln10*(log_density(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              log_density(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx) &
              +log_density(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              log_density(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx)))
!  DEL (ACTUAL).
!  DEL (RADIATIVE) IS INTERPOLATED, AND DEL IS THE MIN OF DELA,DELR.
         rot_scr%delmi(zone_idx)=mix_scr%del_radiative_mix(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              mix_scr%del_radiative_mix(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx)+ &
              mix_scr%del_radiative_mix(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              mix_scr%del_radiative_mix(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx)
!  DEL(ADIABATIC).
         rot_scr%delami(zone_idx)=mix_scr%del_adiabatic_mix(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              mix_scr%del_adiabatic_mix(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx)+ &
              mix_scr%del_adiabatic_mix(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              mix_scr%del_adiabatic_mix(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx)
!  D LN RHO/D LN T.
         rot_scr%qdtmi(zone_idx)=mix_scr%qdtm(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              mix_scr%qdtm(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx)+ &
              mix_scr%qdtm(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              mix_scr%qdtm(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx)
!  UNLOGGED MASS INTERIOR TO THE INTERFACE.
         rot_scr%hs3(zone_idx)=mass_unlogged(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              mass_unlogged(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx)+ &
              mass_unlogged(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              mass_unlogged(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx)
!  SPECIFIC ENERGY GENERATION RATE.
         rot_scr%epsilm(zone_idx)=mix_scr%esumm(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              mix_scr%esumm(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx)+ &
              mix_scr%esumm(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              mix_scr%esumm(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx)
!  LUMINOSITY.
         rot_scr%interface_luminosity(zone_idx)=star%solar_luminosity_cgs*(luminosity(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              luminosity(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx)+ &
              luminosity(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              luminosity(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx))
!  LOCAL AVERAGE FORCE OF GRAVITY.
         rot_scr%interface_gravity_factor(zone_idx)=local_gravity(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              local_gravity(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx)+ &
              local_gravity(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              local_gravity(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx)
!  opacity.
         opacity_interface(zone_idx)=mix_scr%om(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              mix_scr%om(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx)+ &
              mix_scr%om(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              mix_scr%om(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx)
!  specific heat
         specific_heat_interface(zone_idx)=mix_scr%cpm(zone_idx-2)*rot_scr%lagrange_interp_weights(1,zone_idx)+ &
              mix_scr%cpm(zone_idx-1)*rot_scr%lagrange_interp_weights(2,zone_idx)+ &
              mix_scr%cpm(zone_idx)*rot_scr%lagrange_interp_weights(3,zone_idx)+ &
              mix_scr%cpm(zone_idx+1)*rot_scr%lagrange_interp_weights(4,zone_idx)
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
         f2_local = rot_scr%interface_luminosity(zone_idx)*rot_scr%interface_radius(zone_idx)/rot_scr%hs3(zone_idx)/3.0d0
         f3_local = 0.75d0*rot_scr%hs3(zone_idx)/(cpi*rot_scr%dm(zone_idx)*rot_scr%interface_radius(zone_idx)**3)
         v0_local = -f1_local*f2_local*f3_local
!         V0 = -PM(I)*ALM(I)/(C4PI*CG**2*HS3(I)*DDEL*CPM(I)*
!     *         DM(I)**2*TM(I))
         ff_factor = rot_scr%interface_radius(zone_idx)**3/rot_scr%hs3(zone_idx)
         c1_factor = cc23*rot_scr%interface_radius(zone_idx)**4/grav_const/rot_scr%hs3(zone_idx)
         rot_scr%second_deriv_geom_factor(zone_idx) = c1_factor
         qc1r = 4.0d0*cc23*ff_factor/grav_const* &
              (1.0d0-cpi*rot_scr%dm(zone_idx)*ff_factor)
         dr_local = exp10(log_radius(zone_idx))-exp10(log_radius(zone_idx-1))
         qdr_local = (exp10(log_density(zone_idx))-exp10(log_density(zone_idx-1)))/dr_local
         qqc1rr = 8.0d0*ff_factor/grav_const/rot_scr%interface_radius(zone_idx)* &
         (1.0d0-cc13*cpi*ff_factor*(1.0d1*rot_scr%dm(zone_idx)-rot_scr%interface_radius(zone_idx)*qdr_local) &
          + cc13*8.0d0*(cpi*rot_scr%dm(zone_idx)*ff_factor)**2)
! D LN CHI/D LN T = 3 - D LN CAPPA/D LN T
         qchit = 3.0d0 - 0.5d0*(rot_scr%dlnkappa_dlnt(zone_idx)+rot_scr%dlnkappa_dlnt(zone_idx-1))
         qqchitr = (rot_scr%dlnkappa_dlnt(zone_idx-1)-rot_scr%dlnkappa_dlnt(zone_idx))/dr_local
         ht_temp_scale = exp(ln10*(log_pressure(zone_idx)+2.0d0*log_radius(zone_idx)-log_density(zone_idx)))/ &
              mass_unlogged(zone_idx)/grav_const/mix_scr%delm(zone_idx)
         rot_scr%third_deriv_geom_factor(zone_idx) = ht_temp_scale
         ht_temp_scale2 = dr_local/ln10/(log_temperature(zone_idx-1)-log_temperature(zone_idx))
         dhtscale_dr = (abs(ht_temp_scale)-abs(ht_temp_scale_prev))/dr_local
         ht_temp_scale_prev = ht_temp_scale
         mean_dlneps_dlnt = 0.5d0*(rot_scr%dlnepsilon_dlnt(zone_idx)+rot_scr%dlnepsilon_dlnt(zone_idx-1))
         mean_neutrino_fraction = 0.5d0*(rot_scr%neutrino_loss_fraction(zone_idx)+rot_scr%neutrino_loss_fraction(zone_idx-1))
!         FACT7(I)= -V0*(QQCHITR*C1+QCHIT*QC1R)
!         FACT7(I)= -V0*(QQCHITR*C1)
         rot_scr%difad_shear_coeff1(zone_idx)= &
         - f1_local*c1_factor*(mean_neutrino_fraction*mean_dlneps_dlnt + &
              rot_scr%epsilm(zone_idx)*(1.0d0 - mean_neutrino_fraction - qchit))
         rot_scr%velocity_coeff2b(zone_idx) = rot_scr%velocity_coeff2b(zone_idx)+rot_scr%difad_shear_coeff1(zone_idx)/f1_local
!         FACT7(I)=V0*((QHTR-QCHIT)*QC1R-QQCHITR*C1+HTSC*QQC1RR)
!     *            - F1*EPSILM(I)*(HTSC*QC1R + C1*(QETM - QCHIT))
         rot_scr%difad_shear_coeff2(zone_idx) = 0.0d0
!         FACT8(I)=V0*((QHTR-QCHIT)*C1+HTSC*2.0D0*QC1R)
!     *            - F1*EPSILM(I)*HTSC*C1
! Q variables not used
         local_flux_factor = rot_scr%interface_luminosity(zone_idx)/c4pi/rot_scr%dm(zone_idx)/rot_scr%interface_radius(zone_idx)**2
         rot_scr%facd2(zone_idx) = f1_local*(local_flux_factor*qchit - rot_scr%epsilm(zone_idx)*temp_scale_factor)
         rot_scr%facd3(zone_idx) = -f1_local*local_flux_factor
         end do
      endif
      return
end subroutine rotation_stability_setup

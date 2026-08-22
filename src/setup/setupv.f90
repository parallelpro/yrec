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
subroutine setupv(log_density, local_gravity, luminosity, log_pressure, &
     log_radius, mass_unlogged, log_temperature, transport_zone_begin, &
     transport_zone_end, num_zones, radius_unlogged, &
     dynamical_shear_omega_limit)

      use star_info_lib, only: star
      use star_info_lib, only: star
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: log_density(json), local_gravity(json), &
           luminosity(json), log_pressure(json), log_radius(json), &
           mass_unlogged(json), log_temperature(json)
      integer, intent(in) :: transport_zone_begin, transport_zone_end, &
           num_zones
      double precision, intent(out) :: radius_unlogged(json)
      double precision, intent(out) :: dynamical_shear_omega_limit(json)














! JvS 09/25 CHANGED CPM --> CPMI TO AVOID CONFLICT IN MDPHY
      double precision :: specific_heat_interface(json), opacity_interface(json)
      save
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
    5 continue
      end do
!  FIND LOCATION (IN RADIUS) OF THE MIDPOINTS OF THE INTERFACES.
      do zone_idx = 2,num_zones
         star%rot%interface_radius(zone_idx) = 0.5d0*(radius_unlogged(zone_idx) + &
              radius_unlogged(zone_idx-1))
   10 continue
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
         lag_x1 = star%rot%interface_radius(2) - radius_unlogged(1)
         lag_x2 = star%rot%interface_radius(2) - radius_unlogged(2)
         lag_x3 = star%rot%interface_radius(2) - radius_unlogged(3)
         lag_x4 = star%rot%interface_radius(2) - radius_unlogged(4)
         star%rot%lagrange_interp_weights(1,2) = (lag_x2*lag_x3*lag_x4)/lag_denom1
         star%rot%lagrange_interp_weights(2,2) = (lag_x1*lag_x3*lag_x4)/lag_denom2
         star%rot%lagrange_interp_weights(3,2) = (lag_x1*lag_x2*lag_x4)/lag_denom3
         star%rot%lagrange_interp_weights(4,2) = (lag_x1*lag_x2*lag_x3)/lag_denom4
         interior_begin = 3
         star%rot%pm(2)=exp(ln10*(log_pressure(1)*star%rot%lagrange_interp_weights(1,2)+ &
              log_pressure(2)*star%rot%lagrange_interp_weights(2,2) &
              +log_pressure(3)*star%rot%lagrange_interp_weights(3,2)+ &
              log_pressure(4)*star%rot%lagrange_interp_weights(4,2)))
         star%rot%tm(2)=exp(ln10*(log_temperature(1)*star%rot%lagrange_interp_weights(1,2)+ &
              log_temperature(2)*star%rot%lagrange_interp_weights(2,2) &
              +log_temperature(3)*star%rot%lagrange_interp_weights(3,2)+ &
              log_temperature(4)*star%rot%lagrange_interp_weights(4,2)))
         star%rot%dm(2)=exp(ln10*(log_density(1)*star%rot%lagrange_interp_weights(1,2)+ &
              log_density(2)*star%rot%lagrange_interp_weights(2,2) &
              +log_density(3)*star%rot%lagrange_interp_weights(3,2)+ &
              log_density(4)*star%rot%lagrange_interp_weights(4,2)))
         star%rot%delmi(2)=star%mix_phys%del_radiative_mix(1)*star%rot%lagrange_interp_weights(1,2)+ &
              star%mix_phys%del_radiative_mix(2)*star%rot%lagrange_interp_weights(2,2)+ &
              star%mix_phys%del_radiative_mix(3)*star%rot%lagrange_interp_weights(3,2)+ &
              star%mix_phys%del_radiative_mix(4)*star%rot%lagrange_interp_weights(4,2)
         star%rot%delami(2)=star%mix_phys%del_adiabatic_mix(1)*star%rot%lagrange_interp_weights(1,2)+ &
              star%mix_phys%del_adiabatic_mix(2)*star%rot%lagrange_interp_weights(2,2)+ &
              star%mix_phys%del_adiabatic_mix(3)*star%rot%lagrange_interp_weights(3,2)+ &
              star%mix_phys%del_adiabatic_mix(4)*star%rot%lagrange_interp_weights(4,2)
         star%rot%qdtmi(2)=star%mix_phys%qdtm(1)*star%rot%lagrange_interp_weights(1,2)+star%mix_phys%qdtm(2)*star%rot%lagrange_interp_weights(2,2)+ &
              star%mix_phys%qdtm(3)*star%rot%lagrange_interp_weights(3,2)+star%mix_phys%qdtm(4)*star%rot%lagrange_interp_weights(4,2)
         star%rot%hs3(2)=mass_unlogged(1)*star%rot%lagrange_interp_weights(1,2)+mass_unlogged(2)*star%rot%lagrange_interp_weights(2,2)+ &
              mass_unlogged(3)*star%rot%lagrange_interp_weights(3,2)+mass_unlogged(4)*star%rot%lagrange_interp_weights(4,2)
         star%rot%epsilm(2)=star%mix_phys%esumm(1)*star%rot%lagrange_interp_weights(1,2)+star%mix_phys%esumm(2)*star%rot%lagrange_interp_weights(2,2)+ &
              star%mix_phys%esumm(3)*star%rot%lagrange_interp_weights(3,2)+star%mix_phys%esumm(4)*star%rot%lagrange_interp_weights(4,2)
         star%rot%interface_luminosity(2)=solar_luminosity_cgs*(luminosity(1)*star%rot%lagrange_interp_weights(1,2)+ &
              luminosity(2)*star%rot%lagrange_interp_weights(2,2)+ &
              luminosity(3)*star%rot%lagrange_interp_weights(3,2)+luminosity(4)*star%rot%lagrange_interp_weights(4,2))
         star%rot%interface_gravity_factor(2)=local_gravity(1)*star%rot%lagrange_interp_weights(1,2)+ &
              local_gravity(2)*star%rot%lagrange_interp_weights(2,2)+ &
              local_gravity(3)*star%rot%lagrange_interp_weights(3,2)+local_gravity(4)*star%rot%lagrange_interp_weights(4,2)
!  opacity.
         opacity_interface(2)=star%mix_phys%om(1)*star%rot%lagrange_interp_weights(1,2)+star%mix_phys%om(2)*star%rot%lagrange_interp_weights(2,2)+ &
              star%mix_phys%om(3)*star%rot%lagrange_interp_weights(3,2)+star%mix_phys%om(4)*star%rot%lagrange_interp_weights(4,2)
!  specific heat
         specific_heat_interface(2)=star%mix_phys%cpm(1)*star%rot%lagrange_interp_weights(1,2)+star%mix_phys%cpm(2)*star%rot%lagrange_interp_weights(2,2)+ &
              star%mix_phys%cpm(3)*star%rot%lagrange_interp_weights(3,2)+star%mix_phys%cpm(4)*star%rot%lagrange_interp_weights(4,2)
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
         lag_x1 = star%rot%interface_radius(num_zones) - radius_unlogged(num_zones-3)
         lag_x2 = star%rot%interface_radius(num_zones) - radius_unlogged(num_zones-2)
         lag_x3 = star%rot%interface_radius(num_zones) - radius_unlogged(num_zones-1)
         lag_x4 = star%rot%interface_radius(num_zones) - radius_unlogged(num_zones)
         star%rot%lagrange_interp_weights(1,num_zones) = (lag_x2*lag_x3*lag_x4)/lag_denom1
         star%rot%lagrange_interp_weights(2,num_zones) = (lag_x1*lag_x3*lag_x4)/lag_denom2
         star%rot%lagrange_interp_weights(3,num_zones) = (lag_x1*lag_x2*lag_x4)/lag_denom3
         star%rot%lagrange_interp_weights(4,num_zones) = (lag_x1*lag_x2*lag_x3)/lag_denom4
         interior_end = num_zones-1
         star%rot%pm(num_zones)=exp(ln10*(log_pressure(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              log_pressure(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones) &
              +log_pressure(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              log_pressure(num_zones)*star%rot%lagrange_interp_weights(4,num_zones)))
         star%rot%tm(num_zones)=exp(ln10*(log_temperature(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              log_temperature(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones) &
              +log_temperature(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              log_temperature(num_zones)*star%rot%lagrange_interp_weights(4,num_zones)))
         star%rot%dm(num_zones)=exp(ln10*(log_density(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              log_density(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones) &
              +log_density(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              log_density(num_zones)*star%rot%lagrange_interp_weights(4,num_zones)))
         star%rot%delmi(num_zones)=star%mix_phys%del_radiative_mix(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              star%mix_phys%del_radiative_mix(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones)+ &
              star%mix_phys%del_radiative_mix(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              star%mix_phys%del_radiative_mix(num_zones)*star%rot%lagrange_interp_weights(4,num_zones)
         star%rot%delami(num_zones)=star%mix_phys%del_adiabatic_mix(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              star%mix_phys%del_adiabatic_mix(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones)+ &
              star%mix_phys%del_adiabatic_mix(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              star%mix_phys%del_adiabatic_mix(num_zones)*star%rot%lagrange_interp_weights(4,num_zones)
         star%rot%qdtmi(num_zones)=star%mix_phys%qdtm(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              star%mix_phys%qdtm(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones)+ &
              star%mix_phys%qdtm(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              star%mix_phys%qdtm(num_zones)*star%rot%lagrange_interp_weights(4,num_zones)
         star%rot%hs3(num_zones)=mass_unlogged(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              mass_unlogged(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones)+ &
              mass_unlogged(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              mass_unlogged(num_zones)*star%rot%lagrange_interp_weights(4,num_zones)
         star%rot%epsilm(num_zones)=star%mix_phys%esumm(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              star%mix_phys%esumm(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones)+ &
              star%mix_phys%esumm(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              star%mix_phys%esumm(num_zones)*star%rot%lagrange_interp_weights(4,num_zones)
         star%rot%interface_luminosity(num_zones)=solar_luminosity_cgs*(luminosity(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              luminosity(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones)+ &
              luminosity(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              luminosity(num_zones)*star%rot%lagrange_interp_weights(4,num_zones))
         star%rot%interface_gravity_factor(num_zones)=local_gravity(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              local_gravity(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones)+ &
              local_gravity(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              local_gravity(num_zones)*star%rot%lagrange_interp_weights(4,num_zones)
!  opacity.
         opacity_interface(num_zones)=star%mix_phys%om(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              star%mix_phys%om(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones)+ &
              star%mix_phys%om(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              star%mix_phys%om(num_zones)*star%rot%lagrange_interp_weights(4,num_zones)
!  specific heat
         specific_heat_interface(num_zones)=star%mix_phys%cpm(num_zones-3)*star%rot%lagrange_interp_weights(1,num_zones)+ &
              star%mix_phys%cpm(num_zones-2)*star%rot%lagrange_interp_weights(2,num_zones)+ &
              star%mix_phys%cpm(num_zones-1)*star%rot%lagrange_interp_weights(3,num_zones)+ &
              star%mix_phys%cpm(num_zones)*star%rot%lagrange_interp_weights(4,num_zones)
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
         lag_x1 = star%rot%interface_radius(zone_idx) - radius_unlogged(zone_idx-2)
         lag_x2 = star%rot%interface_radius(zone_idx) - radius_unlogged(zone_idx-1)
         lag_x3 = star%rot%interface_radius(zone_idx) - radius_unlogged(zone_idx)
         lag_x4 = star%rot%interface_radius(zone_idx) - radius_unlogged(zone_idx+1)
         star%rot%lagrange_interp_weights(1,zone_idx) = (lag_x2*lag_x3*lag_x4)/lag_denom1
         star%rot%lagrange_interp_weights(2,zone_idx) = (lag_x1*lag_x3*lag_x4)/lag_denom2
         star%rot%lagrange_interp_weights(3,zone_idx) = (lag_x1*lag_x2*lag_x4)/lag_denom3
         star%rot%lagrange_interp_weights(4,zone_idx) = (lag_x1*lag_x2*lag_x3)/lag_denom4
   20 continue
      end do
      grav_const = exp(ln10*cgl)
      grav_const_sq = grav_const**2
      do zone_idx = interior_begin,interior_end
!  USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND PHYSICAL VARIABLES
!  AT THE INTERFACES.
!  PRESSURE.
         star%rot%pm(zone_idx)=exp(ln10*(log_pressure(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              log_pressure(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx) &
              +log_pressure(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              log_pressure(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx)))
!  TEMPERATURE.
         star%rot%tm(zone_idx)=exp(ln10*(log_temperature(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              log_temperature(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx) &
              +log_temperature(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              log_temperature(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx)))
!  DENSITY.
         star%rot%dm(zone_idx)=exp(ln10*(log_density(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              log_density(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx) &
              +log_density(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              log_density(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx)))
!  DEL (ACTUAL).
!  DEL (RADIATIVE) IS INTERPOLATED, AND DEL IS THE MIN OF DELA,DELR.
         star%rot%delmi(zone_idx)=star%mix_phys%del_radiative_mix(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              star%mix_phys%del_radiative_mix(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx)+ &
              star%mix_phys%del_radiative_mix(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              star%mix_phys%del_radiative_mix(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx)
!  DEL(ADIABATIC).
         star%rot%delami(zone_idx)=star%mix_phys%del_adiabatic_mix(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              star%mix_phys%del_adiabatic_mix(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx)+ &
              star%mix_phys%del_adiabatic_mix(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              star%mix_phys%del_adiabatic_mix(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx)
!  D LN RHO/D LN T.
         star%rot%qdtmi(zone_idx)=star%mix_phys%qdtm(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              star%mix_phys%qdtm(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx)+ &
              star%mix_phys%qdtm(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              star%mix_phys%qdtm(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx)
!  UNLOGGED MASS INTERIOR TO THE INTERFACE.
         star%rot%hs3(zone_idx)=mass_unlogged(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              mass_unlogged(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx)+ &
              mass_unlogged(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              mass_unlogged(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx)
!  SPECIFIC ENERGY GENERATION RATE.
         star%rot%epsilm(zone_idx)=star%mix_phys%esumm(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              star%mix_phys%esumm(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx)+ &
              star%mix_phys%esumm(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              star%mix_phys%esumm(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx)
!  LUMINOSITY.
         star%rot%interface_luminosity(zone_idx)=solar_luminosity_cgs*(luminosity(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              luminosity(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx)+ &
              luminosity(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              luminosity(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx))
!  LOCAL AVERAGE FORCE OF GRAVITY.
         star%rot%interface_gravity_factor(zone_idx)=local_gravity(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              local_gravity(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx)+ &
              local_gravity(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              local_gravity(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx)
!  opacity.
         opacity_interface(zone_idx)=star%mix_phys%om(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              star%mix_phys%om(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx)+ &
              star%mix_phys%om(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              star%mix_phys%om(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx)
!  specific heat
         specific_heat_interface(zone_idx)=star%mix_phys%cpm(zone_idx-2)*star%rot%lagrange_interp_weights(1,zone_idx)+ &
              star%mix_phys%cpm(zone_idx-1)*star%rot%lagrange_interp_weights(2,zone_idx)+ &
              star%mix_phys%cpm(zone_idx)*star%rot%lagrange_interp_weights(3,zone_idx)+ &
              star%mix_phys%cpm(zone_idx+1)*star%rot%lagrange_interp_weights(4,zone_idx)
   30 continue
      end do
      do zone_idx = transport_zone_begin,transport_zone_end
         star%rot%delmi(zone_idx) = min(star%rot%delmi(zone_idx),star%rot%delami(zone_idx))
   35 continue
      end do
      do zone_idx = 1,num_zones
         star%rot%es_velocity_coeff1(zone_idx) = 0.0d0
         star%rot%es_velocity_coeff2(zone_idx) = 0.0d0
         star%rot%es_shear_coeff(zone_idx) = 0.0d0
         star%rot%fgsfj(zone_idx) = 0.0d0
         star%rot%gsf_kippenhahn_coeff(zone_idx) = 0.0d0
         star%rot%fact1(zone_idx) = 0.0d0
         star%rot%fact2(zone_idx) = 0.0d0
         star%rot%fact3(zone_idx) = 0.0d0
         star%rot%fact4(zone_idx) = 0.0d0
         star%rot%mu_gradient_richardson_coeff(zone_idx) = 0.0d0
! MHP 06/02
         star%rot%es_relaxation_factor(zone_idx) = 0.0d0
         star%rot%theta_mean(zone_idx) = 0.0d0
         star%rot%difad_shear_coeff1(zone_idx) = 0.0d0
         star%rot%difad_shear_coeff2(zone_idx) = 0.0d0
      end do
      do zone_idx = 2,num_zones
! ADDED FOR D THETA/DT TERM FROM ZAHN&MAEDER 1998
         fac_local = 2.0d0*(log_radius(zone_idx)+log_radius(zone_idx-1))-0.5d0* &
              (log10(mass_unlogged(zone_idx))+log10(mass_unlogged(zone_idx-1)))-cgl
         star%rot%theta_mean(zone_idx) = exp(ln10*fac_local)
      end do
!  NOW COMPUTE STRUCTURAL QUANTITIES NEEDED TO EVALUATE VELOCITIES AT
!  ALL INTERFACES.
      cpigi_const = 4.0d0/c4pi/grav_const
      do zone_idx = transport_zone_begin,transport_zone_end
         if(.not.use_diffusion_advection_transport)then
            dlnmu_dlnp = (log10(star%mix_phys%amum(zone_idx))-log10(star%mix_phys%amum(zone_idx-1)))/ &
                 (log_pressure(zone_idx)-log_pressure(zone_idx-1))
         else
            dlnmu_dlnp = 0.0d0
         endif
         if(zone_idx.eq.transport_zone_begin)then
            ddel_floor=max(1.0d-3,0.5d0*(star%mix_phys%del_adiabatic_mix(transport_zone_begin+1)- &
                 star%mix_phys%delm(transport_zone_begin+1)) &
                 +dlnmu_dlnp)
         else if(zone_idx.eq.transport_zone_end)then
            ddel_floor=max(1.0d-3,0.5d0*(star%mix_phys%del_adiabatic_mix(transport_zone_end-1)- &
                 star%mix_phys%delm(transport_zone_end-1)) &
                 +dlnmu_dlnp)
         else
            ddel_floor = 1.0d-3
         endif
         ddel = max(star%rot%delami(zone_idx)-star%rot%delmi(zone_idx),ddel_floor)
! MHP 06/02
! ADDED FOR ALTERNATE TREATMENT OF MU GRADIENTS
!         DDELM(I) = 0.5D0*(DELAMI(I)+DELAMI(I-1)-
!     *                     DELMI(I)-DELMI(I-1))
!         DDELM(I) = DDEL
         star%rot%del_grad_diff_interface(zone_idx) = max(star%rot%delami(zone_idx)-star%rot%delmi(zone_idx)+dlnmu_dlnp,ddel_floor)
!         FESTIME(I) = PM(I)/(HGM(I)*DDEL*DM(I)*TM(I))
         star%rot%es_relaxation_factor(zone_idx) = star%rot%pm(zone_idx)/(star%rot%interface_gravity_factor(zone_idx)*ddel*star%rot%dm(zone_idx))
         pressure_scale_factor = star%rot%pm(zone_idx)*star%rot%interface_radius(zone_idx)**2/star%rot%dm(zone_idx)/star%rot%hs3(zone_idx)/grav_const
         temp_scale_factor = pressure_scale_factor/star%rot%delmi(zone_idx)
         eta_factor = 2.0d0*cc23*star%rot%interface_radius(zone_idx)**3/grav_const/star%rot%hs3(zone_idx)
!        TTHERM = 8.0d0*CC23*CSIG*TM(I)**3/OPM(I)/DM(I)**2/CPM(I)
         star%rot%gsf_kippenhahn_coeff(zone_idx)= 8.0d0*pressure_scale_factor*eta_factor/ddel
         ff_factor = star%rot%pm(zone_idx)/(star%rot%interface_gravity_factor(zone_idx)*ddel*specific_heat_interface(zone_idx)* &
              star%rot%dm(zone_idx)*star%rot%tm(zone_idx))
         fac_local = 2.0d0*eta_factor*ff_factor
         specific_luminosity = star%rot%interface_luminosity(zone_idx)/star%rot%hs3(zone_idx)
         star%rot%es_velocity_coeff1(zone_idx) = fac_local*(specific_luminosity-star%rot%epsilm(zone_idx))
         star%rot%es_velocity_coeff2(zone_idx)= -0.5d0*fac_local*specific_luminosity*cpigi_const/star%rot%dm(zone_idx)
         star%rot%es_shear_coeff(zone_idx) = specific_luminosity*ff_factor*cc13*cpigi_const*star%rot%interface_radius(zone_idx)/star%rot%dm(zone_idx)
!         FES3(I) = EM*FF*CPIGI*(3.0D0*HTSC-RM(I))/DM(I)
         star%rot%velocity_coeff0(zone_idx) = ff_factor
         star%rot%velocity_coeff1a(zone_idx) = star%rot%es_velocity_coeff1(zone_idx)/ff_factor
         star%rot%velocity_coeff1b(zone_idx) = star%rot%es_velocity_coeff2(zone_idx)/ff_factor
         star%rot%velocity_coeff2a(zone_idx) = specific_luminosity*cpigi_const*3.0d0*temp_scale_factor/star%rot%dm(zone_idx)
         star%rot%velocity_coeff2b(zone_idx) = -specific_luminosity*cpigi_const*star%rot%interface_radius(zone_idx)/star%rot%dm(zone_idx)
         star%rot%fgsfj(zone_idx) = abs(fac_local*specific_luminosity*temp_scale_factor)*star%rot%interface_radius(zone_idx)
!  EDDINGTON CIRCULATION VELOCITY IS DEFINED AS
!  VES = FACT1*FACT2*OMEGA**2 (ENDAL AND SOFIA PAPER II).
!  NOTE THAT TO AVOID OVERFLOW 1/(DEL(AD)-DEL)IS SET TO A MAXIMUM OF 10^6.
!  THIS SHOULD ONLY BE AN ISSUE IF YOU HAVE ZERO OVERSHOOT AT THE BOUNDARY
!  OF A CONVECTION ZONE.  CAVEAT EMPTOR.
         star%rot%fact1(zone_idx) = 1.0d0/star%rot%qdtmi(zone_idx)/max(star%rot%delami(zone_idx)-star%rot%delmi(zone_idx),1.0d-3)
!         FACT2(I) = DELAMI(I)*ALM(I)*(RM(I)**3/CG2/HS3(I)**2)*
!     *              (2.0D0*RM(I)**2*(EPSILM(I)/ALM(I) - 1.0D0/HS3(I))
!    *              - 3.0D0/(C4PI*DM(I)*RM(I)))
         fac_local = star%rot%delami(zone_idx)/star%rot%interface_gravity_factor(zone_idx)
         star%rot%fact2(zone_idx) = fac_local*eta_factor*(star%rot%epsilm(zone_idx)-specific_luminosity)
         star%rot%fact6(zone_idx) = -fac_local*specific_luminosity/c4pi/grav_const/star%rot%dm(zone_idx)/star%rot%interface_radius(zone_idx)
! GSF VELOCITY IS DEFINED AS
! VGSF = VES*FACT3/Hj**2/(2OMEGA/(D OMEGA/D LNR) + 1)
! WHERE HJ IS THE SPECIFIC ANGULAR MOMENTUM SCALE HEIGHT.
! AND FACT3 IS THE TEMPERATURE SCALE HEIGHT TIMES THE RADIUS.
!  HT = (DEL*D(LN P)/DR)**-1 = P*R**2/(DEL*RHO*GM)
!        FACT3(I) = PM(I)*RM(I)**3/(DELMI(I)*DM(I)*CG*HS3(I))
         star%rot%fact3(zone_idx) = abs(temp_scale_factor*star%rot%interface_radius(zone_idx))
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
         star%rot%fact4(zone_idx) = star%rot%interface_radius(zone_idx)**2/star%rot%hs3(zone_idx)/grav_const/star%rot%dm(zone_idx)* &
              (star%rot%pm(zone_idx)-radiation_constant_over_3*star%rot%tm(zone_idx)**4)/sqrt(star%rot%interface_gravity_factor(zone_idx))
!  DYNAMICAL SHEAR INSTABILITY.
!  AN INTERFACE IS STABLE AGAINST THIS SHEAR WHENEVER THE RICHARDSON NUMBER
!  RICHNO = (RHO/P)*(DEL(AD)-DEL)*g**2/(D OMEGA/D LN R)**2>.25.
         dynamical_shear_omega_limit(zone_idx)=sqrt(star%rot%dm(zone_idx)/star%rot%pm(zone_idx)* &
              max(star%rot%delami(zone_idx) - star%rot%delmi(zone_idx),1.0d-3)) &
              *2.0d0*star%rot%interface_gravity_factor(zone_idx)
!  THE DIFFUSIVE SHEAR IS INHIBITED BY MU GRADIENTS.
!  AN INTERFACE IS STABLE AGAINST THE DIFFUSIVE SHEAR IF
!  RICHARDSON NUMBER = FACT5*DMU/MU/(D OMEGA/D LNR)**2 > .25, WHERE
!  FACT5 = RHO/P*(-d lnT/d lnMU)/(del LN P)*G**2 AND
!   -d lnT/d lnMU = (CON-1)/(1+3CON),CON=(a/3)T**4/P
!  GIVEN (1 - CON)P = CGAS*RHO*T/MU
!  1 - CON IS BETA(CORRECTION FOR RADIATION PRESSURE).
         star%rot%mu_gradient_richardson_coeff(zone_idx) = star%rot%dm(zone_idx)*star%rot%interface_gravity_factor(zone_idx)**2/ &
              star%rot%pm(zone_idx)/ln10/(log_pressure(zone_idx)-log_pressure(zone_idx-1))* &
              (radiation_constant_over_3*star%rot%tm(zone_idx)**4/star%rot%pm(zone_idx)-1.0d0)/ &
              (1.0d0 + 3.0d0*radiation_constant_over_3*star%rot%tm(zone_idx)**4/star%rot%pm(zone_idx))
! MHP 3/92 ADD VECTOR FOR LOCAL KELVIN-HELMHOLTZ TIME SCALE.
!         TKH(I) = CG*HS3(I)**2/ABS(ALM(I))/RM(I)
   40 continue
      end do
! MHP 06/02 ADDED TERMS OF ORDER DW/DR FROM ZAH&MAEDER 1998
      if(use_diffusion_advection_transport)then
         ht_temp_scale_prev = exp(ln10*(log_pressure(transport_zone_begin-1)+ &
              2.0d0*log_radius(transport_zone_begin-1)-log_density(transport_zone_begin-1))) &
              /mass_unlogged(transport_zone_begin-1)/grav_const/star%mix_phys%del_radiative_mix(transport_zone_begin-1)
         star%rot%third_deriv_geom_factor(transport_zone_begin-1) = ht_temp_scale_prev
         do zone_idx = transport_zone_begin,transport_zone_end
!         DDEL = MAX(DELAMI(I)-DELMI(I),1.0D-3)
         pressure_scale_factor = star%rot%pm(zone_idx)*star%rot%interface_radius(zone_idx)**2/star%rot%dm(zone_idx)/star%rot%hs3(zone_idx)/grav_const
         temp_scale_factor = pressure_scale_factor/star%rot%delmi(zone_idx)
         f1_local = star%rot%pm(zone_idx)/(star%rot%interface_gravity_factor(zone_idx)*star%rot%del_grad_diff_interface(zone_idx)* &
              specific_heat_interface(zone_idx)*star%rot%dm(zone_idx)*star%rot%tm(zone_idx))
         f2_local = star%rot%interface_luminosity(zone_idx)*star%rot%interface_radius(zone_idx)/star%rot%hs3(zone_idx)/3.0d0
         f3_local = 0.75d0*star%rot%hs3(zone_idx)/(cpi*star%rot%dm(zone_idx)*star%rot%interface_radius(zone_idx)**3)
         v0_local = -f1_local*f2_local*f3_local
!         V0 = -PM(I)*ALM(I)/(C4PI*CG**2*HS3(I)*DDEL*CPM(I)*
!     *         DM(I)**2*TM(I))
         ff_factor = star%rot%interface_radius(zone_idx)**3/star%rot%hs3(zone_idx)
         c1_factor = cc23*star%rot%interface_radius(zone_idx)**4/grav_const/star%rot%hs3(zone_idx)
         star%rot%second_deriv_geom_factor(zone_idx) = c1_factor
         qc1r = 4.0d0*cc23*ff_factor/grav_const* &
              (1.0d0-cpi*star%rot%dm(zone_idx)*ff_factor)
         dr_local = 10.0d0**log_radius(zone_idx)-10.0d0**log_radius(zone_idx-1)
         qdr_local = (10.0d0**log_density(zone_idx)-10.0d0**log_density(zone_idx-1))/dr_local
         qqc1rr = 8.0d0*ff_factor/grav_const/star%rot%interface_radius(zone_idx)* &
         (1.0d0-cc13*cpi*ff_factor*(1.0d1*star%rot%dm(zone_idx)-star%rot%interface_radius(zone_idx)*qdr_local) &
          + cc13*8.0d0*(cpi*star%rot%dm(zone_idx)*ff_factor)**2)
! D LN CHI/D LN T = 3 - D LN CAPPA/D LN T
         qchit = 3.0d0 - 0.5d0*(star%rot%dlnkappa_dlnt(zone_idx)+star%rot%dlnkappa_dlnt(zone_idx-1))
         qqchitr = (star%rot%dlnkappa_dlnt(zone_idx-1)-star%rot%dlnkappa_dlnt(zone_idx))/dr_local
         ht_temp_scale = exp(ln10*(log_pressure(zone_idx)+2.0d0*log_radius(zone_idx)-log_density(zone_idx)))/ &
              mass_unlogged(zone_idx)/grav_const/star%mix_phys%delm(zone_idx)
         star%rot%third_deriv_geom_factor(zone_idx) = ht_temp_scale
         ht_temp_scale2 = dr_local/ln10/(log_temperature(zone_idx-1)-log_temperature(zone_idx))
         dhtscale_dr = (abs(ht_temp_scale)-abs(ht_temp_scale_prev))/dr_local
         ht_temp_scale_prev = ht_temp_scale
         mean_dlneps_dlnt = 0.5d0*(star%rot%dlnepsilon_dlnt(zone_idx)+star%rot%dlnepsilon_dlnt(zone_idx-1))
         mean_neutrino_fraction = 0.5d0*(star%rot%neutrino_loss_fraction(zone_idx)+star%rot%neutrino_loss_fraction(zone_idx-1))
!         FACT7(I)= -V0*(QQCHITR*C1+QCHIT*QC1R)
!         FACT7(I)= -V0*(QQCHITR*C1)
         star%rot%difad_shear_coeff1(zone_idx)= &
         - f1_local*c1_factor*(mean_neutrino_fraction*mean_dlneps_dlnt + &
              star%rot%epsilm(zone_idx)*(1.0d0 - mean_neutrino_fraction - qchit))
         star%rot%velocity_coeff2b(zone_idx) = star%rot%velocity_coeff2b(zone_idx)+star%rot%difad_shear_coeff1(zone_idx)/f1_local
!         FACT7(I)=V0*((QHTR-QCHIT)*QC1R-QQCHITR*C1+HTSC*QQC1RR)
!     *            - F1*EPSILM(I)*(HTSC*QC1R + C1*(QETM - QCHIT))
         star%rot%difad_shear_coeff2(zone_idx) = 0.0d0
!         FACT8(I)=V0*((QHTR-QCHIT)*C1+HTSC*2.0D0*QC1R)
!     *            - F1*EPSILM(I)*HTSC*C1
! Q variables not used
!         Q1 = -V0*(QQCHITR*C1)
!         Q2 = -F1*C1*FNU*QETM
!         Q3 = -F1*C1*EPSILM(I)*(1.0D0 - FNU - QCHIT)
!         Q4 = F1*FL*QCHIT
!         Q5 = -F1*EPSILM(I)*HTSC
         local_flux_factor = star%rot%interface_luminosity(zone_idx)/c4pi/star%rot%dm(zone_idx)/star%rot%interface_radius(zone_idx)**2
         star%rot%facd2(zone_idx) = f1_local*(local_flux_factor*qchit - star%rot%epsilm(zone_idx)*temp_scale_factor)
         star%rot%facd3(zone_idx) = -f1_local*local_flux_factor
!         WRITE(*,911)I,V0,FES3(I),FACT7(I),FACD2(I),FACD3(I),DDELM(I)
!         WRITE(*,912)Q1,Q2,Q3,Q4,Q5
! 911     FORMAT(I5,1P6E12.3)
! 912     FORMAT(5X,1P5E12.3)
         end do
      endif
      return
end subroutine setupv

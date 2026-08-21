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

      use mdphy_lib
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





! MHP 06/02 ADDED FACT7 AND FACT8 FOR DIF/AD TREATMENT
! common/vfact/: fact1-fact5/fact7/fact8, all used here (originally
! FACT1-FACT5,FACT7,FACT8). Naming matches vcirc.f90.
      double precision :: fact1(json), fact2(json), fact3(json), fact4(json), &
           mu_gradient_richardson_coeff(json), difad_shear_coeff1(json), &
           difad_shear_coeff2(json)
      common/vfact/ fact1, fact2, fact3, fact4, mu_gradient_richardson_coeff, &
           difad_shear_coeff1, difad_shear_coeff2

! common/intfac/: lagrange_interp_weights (originally FACI), the
! 4-point Lagrangian interpolation weights used to interpolate model
! quantities onto zone interfaces. Naming matches vcirc.f90.
      double precision :: lagrange_interp_weights(4,json)
      common/intfac/ lagrange_interp_weights

! MHP 3/92 COMMON BLOCK ADDED FOR SHEAR VELOCITY
! common/taukh/: fact6/es_velocity_coeff1/es_velocity_coeff2/fgsfj/
! gsf_kippenhahn_coeff/es_shear_coeff, all used here. Naming matches
! vcirc.f90.
      double precision :: fact6(json), es_velocity_coeff1(json), &
           es_velocity_coeff2(json), fgsfj(json), gsf_kippenhahn_coeff(json), &
           es_shear_coeff(json)
      common/taukh/ fact6, es_velocity_coeff1, es_velocity_coeff2, fgsfj, &
           gsf_kippenhahn_coeff, es_shear_coeff

! common/intvar/: all members used here (interface-interpolated model
! quantities). Naming matches vcirc.f90.
      double precision :: interface_luminosity(json), delami(json), &
           delmi(json), dm(json), epsilm(json), interface_gravity_factor(json), &
           hs3(json), pm(json), qdtmi(json), interface_radius(json), tm(json)
      common/intvar/ interface_luminosity, delami, delmi, dm, epsilm, &
           interface_gravity_factor, hs3, pm, qdtmi, interface_radius, tm

! MHP 06/02
! Time change of theta
! common/oldrot2/: del_grad_diff_interface/es_relaxation_factor/
! theta_mean are used here; tho/theta_new/theta_prev/qwrst/wmst/qwrmst
! are unused placeholders. Naming matches vcirc.f90.
      double precision :: tho(json), theta_new(json), theta_mean(json), &
           del_grad_diff_interface(json), es_relaxation_factor(json), &
           theta_prev(json), qwrst(json), wmst(json), qwrmst(json)
      common/oldrot2/ tho, theta_new, theta_mean, del_grad_diff_interface, &
           es_relaxation_factor, theta_prev, qwrst, wmst, qwrmst

! CAPPA AND EPSILON DERIVATIVES
! common/rotder/: this batch's own block (no reuse precedent found
! elsewhere in the already-converted sources). dlnkappa_dlnrho/
! dlnkappa_dlnt are the opacity's logarithmic derivatives w.r.t.
! density/temperature; dlnepsilon_dlnrho/dlnepsilon_dlnt are the
! energy-generation rate's logarithmic derivatives; neutrino_loss_
! fraction is the fraction of energy generation lost to neutrinos.
! All used here.
      double precision :: dlnkappa_dlnrho(json), dlnkappa_dlnt(json), &
           dlnepsilon_dlnrho(json), dlnepsilon_dlnt(json), &
           neutrino_loss_fraction(json)
      common/rotder/ dlnkappa_dlnrho, dlnkappa_dlnt, dlnepsilon_dlnrho, &
           dlnepsilon_dlnt, neutrino_loss_fraction

! DEFINITION TERMS FOR THE SECOND AND THIRD DERIVATIVE
! TERMS IN THE ZAHN FORMULATION
! common/egridder/: all members used here. Naming matches vcirc.f90/
! dadcoeft.f90.
      double precision :: second_deriv_geom_factor_eqgrid(json), &
           third_deriv_geom_factor_eqgrid(json), second_deriv_geom_factor(json), &
           third_deriv_geom_factor(json)
      common/egridder/ second_deriv_geom_factor_eqgrid, &
           third_deriv_geom_factor_eqgrid, second_deriv_geom_factor, &
           third_deriv_geom_factor

! MHP 06/02
! HIGHER ORDER TERMS IN M.C. VELOCITY
! common/difad3/: facd2/facd3/fv2b(velocity_coeff2b)/eq_velocity_
! coeff2b are used here; the remaining members are unused
! placeholders. Naming matches vcirc.f90/dadcoeft.f90.
      double precision :: facd2(json), facd3(json), vesd2(json), &
           vesd3(json), am_2nd_deriv_coeff(json), am_3rd_deriv_coeff(json), &
           geometric_factor(json), velocity_coeff0(json), &
           velocity_coeff1a(json), velocity_coeff1b(json), &
           velocity_coeff2a(json), velocity_coeff2b(json), &
           eq_velocity_coeff0(json), eq_velocity_coeff1a(json), &
           eq_velocity_coeff1b(json), eq_velocity_coeff2a(json), &
           eq_velocity_coeff2b(json), shear_diffusion_coeff(json), &
           gsf_diffusion_coeff(json), shear_diffusion_coeff_eqgrid(json), &
           gsf_diffusion_coeff_eqgrid(json)
      common/difad3/ facd2, facd3, vesd2, vesd3, am_2nd_deriv_coeff, &
           am_3rd_deriv_coeff, geometric_factor, velocity_coeff0, &
           velocity_coeff1a, velocity_coeff1b, velocity_coeff2a, &
           velocity_coeff2b, eq_velocity_coeff0, eq_velocity_coeff1a, &
           eq_velocity_coeff1b, eq_velocity_coeff2a, eq_velocity_coeff2b, &
           shear_diffusion_coeff, gsf_diffusion_coeff, &
           shear_diffusion_coeff_eqgrid, gsf_diffusion_coeff_eqgrid


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
      do 5 zone_idx = 1,num_zones
         radius_unlogged(zone_idx) = exp(ln10*log_radius(zone_idx))
    5 continue
!  FIND LOCATION (IN RADIUS) OF THE MIDPOINTS OF THE INTERFACES.
      do 10 zone_idx = 2,num_zones
         interface_radius(zone_idx) = 0.5d0*(radius_unlogged(zone_idx) + &
              radius_unlogged(zone_idx-1))
   10 continue
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
         lag_x1 = interface_radius(2) - radius_unlogged(1)
         lag_x2 = interface_radius(2) - radius_unlogged(2)
         lag_x3 = interface_radius(2) - radius_unlogged(3)
         lag_x4 = interface_radius(2) - radius_unlogged(4)
         lagrange_interp_weights(1,2) = (lag_x2*lag_x3*lag_x4)/lag_denom1
         lagrange_interp_weights(2,2) = (lag_x1*lag_x3*lag_x4)/lag_denom2
         lagrange_interp_weights(3,2) = (lag_x1*lag_x2*lag_x4)/lag_denom3
         lagrange_interp_weights(4,2) = (lag_x1*lag_x2*lag_x3)/lag_denom4
         interior_begin = 3
         pm(2)=exp(ln10*(log_pressure(1)*lagrange_interp_weights(1,2)+ &
              log_pressure(2)*lagrange_interp_weights(2,2) &
              +log_pressure(3)*lagrange_interp_weights(3,2)+ &
              log_pressure(4)*lagrange_interp_weights(4,2)))
         tm(2)=exp(ln10*(log_temperature(1)*lagrange_interp_weights(1,2)+ &
              log_temperature(2)*lagrange_interp_weights(2,2) &
              +log_temperature(3)*lagrange_interp_weights(3,2)+ &
              log_temperature(4)*lagrange_interp_weights(4,2)))
         dm(2)=exp(ln10*(log_density(1)*lagrange_interp_weights(1,2)+ &
              log_density(2)*lagrange_interp_weights(2,2) &
              +log_density(3)*lagrange_interp_weights(3,2)+ &
              log_density(4)*lagrange_interp_weights(4,2)))
         delmi(2)=mix_phys%del_radiative_mix(1)*lagrange_interp_weights(1,2)+ &
              mix_phys%del_radiative_mix(2)*lagrange_interp_weights(2,2)+ &
              mix_phys%del_radiative_mix(3)*lagrange_interp_weights(3,2)+ &
              mix_phys%del_radiative_mix(4)*lagrange_interp_weights(4,2)
         delami(2)=mix_phys%del_adiabatic_mix(1)*lagrange_interp_weights(1,2)+ &
              mix_phys%del_adiabatic_mix(2)*lagrange_interp_weights(2,2)+ &
              mix_phys%del_adiabatic_mix(3)*lagrange_interp_weights(3,2)+ &
              mix_phys%del_adiabatic_mix(4)*lagrange_interp_weights(4,2)
         qdtmi(2)=mix_phys%qdtm(1)*lagrange_interp_weights(1,2)+mix_phys%qdtm(2)*lagrange_interp_weights(2,2)+ &
              mix_phys%qdtm(3)*lagrange_interp_weights(3,2)+mix_phys%qdtm(4)*lagrange_interp_weights(4,2)
         hs3(2)=mass_unlogged(1)*lagrange_interp_weights(1,2)+mass_unlogged(2)*lagrange_interp_weights(2,2)+ &
              mass_unlogged(3)*lagrange_interp_weights(3,2)+mass_unlogged(4)*lagrange_interp_weights(4,2)
         epsilm(2)=mix_phys%esumm(1)*lagrange_interp_weights(1,2)+mix_phys%esumm(2)*lagrange_interp_weights(2,2)+ &
              mix_phys%esumm(3)*lagrange_interp_weights(3,2)+mix_phys%esumm(4)*lagrange_interp_weights(4,2)
         interface_luminosity(2)=solar_luminosity_cgs*(luminosity(1)*lagrange_interp_weights(1,2)+ &
              luminosity(2)*lagrange_interp_weights(2,2)+ &
              luminosity(3)*lagrange_interp_weights(3,2)+luminosity(4)*lagrange_interp_weights(4,2))
         interface_gravity_factor(2)=local_gravity(1)*lagrange_interp_weights(1,2)+ &
              local_gravity(2)*lagrange_interp_weights(2,2)+ &
              local_gravity(3)*lagrange_interp_weights(3,2)+local_gravity(4)*lagrange_interp_weights(4,2)
!  opacity.
         opacity_interface(2)=mix_phys%om(1)*lagrange_interp_weights(1,2)+mix_phys%om(2)*lagrange_interp_weights(2,2)+ &
              mix_phys%om(3)*lagrange_interp_weights(3,2)+mix_phys%om(4)*lagrange_interp_weights(4,2)
!  specific heat
         specific_heat_interface(2)=mix_phys%cpm(1)*lagrange_interp_weights(1,2)+mix_phys%cpm(2)*lagrange_interp_weights(2,2)+ &
              mix_phys%cpm(3)*lagrange_interp_weights(3,2)+mix_phys%cpm(4)*lagrange_interp_weights(4,2)
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
         lag_x1 = interface_radius(num_zones) - radius_unlogged(num_zones-3)
         lag_x2 = interface_radius(num_zones) - radius_unlogged(num_zones-2)
         lag_x3 = interface_radius(num_zones) - radius_unlogged(num_zones-1)
         lag_x4 = interface_radius(num_zones) - radius_unlogged(num_zones)
         lagrange_interp_weights(1,num_zones) = (lag_x2*lag_x3*lag_x4)/lag_denom1
         lagrange_interp_weights(2,num_zones) = (lag_x1*lag_x3*lag_x4)/lag_denom2
         lagrange_interp_weights(3,num_zones) = (lag_x1*lag_x2*lag_x4)/lag_denom3
         lagrange_interp_weights(4,num_zones) = (lag_x1*lag_x2*lag_x3)/lag_denom4
         interior_end = num_zones-1
         pm(num_zones)=exp(ln10*(log_pressure(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              log_pressure(num_zones-2)*lagrange_interp_weights(2,num_zones) &
              +log_pressure(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              log_pressure(num_zones)*lagrange_interp_weights(4,num_zones)))
         tm(num_zones)=exp(ln10*(log_temperature(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              log_temperature(num_zones-2)*lagrange_interp_weights(2,num_zones) &
              +log_temperature(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              log_temperature(num_zones)*lagrange_interp_weights(4,num_zones)))
         dm(num_zones)=exp(ln10*(log_density(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              log_density(num_zones-2)*lagrange_interp_weights(2,num_zones) &
              +log_density(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              log_density(num_zones)*lagrange_interp_weights(4,num_zones)))
         delmi(num_zones)=mix_phys%del_radiative_mix(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              mix_phys%del_radiative_mix(num_zones-2)*lagrange_interp_weights(2,num_zones)+ &
              mix_phys%del_radiative_mix(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              mix_phys%del_radiative_mix(num_zones)*lagrange_interp_weights(4,num_zones)
         delami(num_zones)=mix_phys%del_adiabatic_mix(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              mix_phys%del_adiabatic_mix(num_zones-2)*lagrange_interp_weights(2,num_zones)+ &
              mix_phys%del_adiabatic_mix(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              mix_phys%del_adiabatic_mix(num_zones)*lagrange_interp_weights(4,num_zones)
         qdtmi(num_zones)=mix_phys%qdtm(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              mix_phys%qdtm(num_zones-2)*lagrange_interp_weights(2,num_zones)+ &
              mix_phys%qdtm(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              mix_phys%qdtm(num_zones)*lagrange_interp_weights(4,num_zones)
         hs3(num_zones)=mass_unlogged(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              mass_unlogged(num_zones-2)*lagrange_interp_weights(2,num_zones)+ &
              mass_unlogged(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              mass_unlogged(num_zones)*lagrange_interp_weights(4,num_zones)
         epsilm(num_zones)=mix_phys%esumm(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              mix_phys%esumm(num_zones-2)*lagrange_interp_weights(2,num_zones)+ &
              mix_phys%esumm(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              mix_phys%esumm(num_zones)*lagrange_interp_weights(4,num_zones)
         interface_luminosity(num_zones)=solar_luminosity_cgs*(luminosity(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              luminosity(num_zones-2)*lagrange_interp_weights(2,num_zones)+ &
              luminosity(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              luminosity(num_zones)*lagrange_interp_weights(4,num_zones))
         interface_gravity_factor(num_zones)=local_gravity(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              local_gravity(num_zones-2)*lagrange_interp_weights(2,num_zones)+ &
              local_gravity(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              local_gravity(num_zones)*lagrange_interp_weights(4,num_zones)
!  opacity.
         opacity_interface(num_zones)=mix_phys%om(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              mix_phys%om(num_zones-2)*lagrange_interp_weights(2,num_zones)+ &
              mix_phys%om(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              mix_phys%om(num_zones)*lagrange_interp_weights(4,num_zones)
!  specific heat
         specific_heat_interface(num_zones)=mix_phys%cpm(num_zones-3)*lagrange_interp_weights(1,num_zones)+ &
              mix_phys%cpm(num_zones-2)*lagrange_interp_weights(2,num_zones)+ &
              mix_phys%cpm(num_zones-1)*lagrange_interp_weights(3,num_zones)+ &
              mix_phys%cpm(num_zones)*lagrange_interp_weights(4,num_zones)
      else
         interior_end = transport_zone_end
      endif
!  COMPUTE INTERPOLATION FACTORS FOR ALL OTHER POINTS.
      do 20 zone_idx = interior_begin,interior_end
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
         lag_x1 = interface_radius(zone_idx) - radius_unlogged(zone_idx-2)
         lag_x2 = interface_radius(zone_idx) - radius_unlogged(zone_idx-1)
         lag_x3 = interface_radius(zone_idx) - radius_unlogged(zone_idx)
         lag_x4 = interface_radius(zone_idx) - radius_unlogged(zone_idx+1)
         lagrange_interp_weights(1,zone_idx) = (lag_x2*lag_x3*lag_x4)/lag_denom1
         lagrange_interp_weights(2,zone_idx) = (lag_x1*lag_x3*lag_x4)/lag_denom2
         lagrange_interp_weights(3,zone_idx) = (lag_x1*lag_x2*lag_x4)/lag_denom3
         lagrange_interp_weights(4,zone_idx) = (lag_x1*lag_x2*lag_x3)/lag_denom4
   20 continue
      grav_const = exp(ln10*cgl)
      grav_const_sq = grav_const**2
      do 30 zone_idx = interior_begin,interior_end
!  USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND PHYSICAL VARIABLES
!  AT THE INTERFACES.
!  PRESSURE.
         pm(zone_idx)=exp(ln10*(log_pressure(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              log_pressure(zone_idx-1)*lagrange_interp_weights(2,zone_idx) &
              +log_pressure(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              log_pressure(zone_idx+1)*lagrange_interp_weights(4,zone_idx)))
!  TEMPERATURE.
         tm(zone_idx)=exp(ln10*(log_temperature(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              log_temperature(zone_idx-1)*lagrange_interp_weights(2,zone_idx) &
              +log_temperature(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              log_temperature(zone_idx+1)*lagrange_interp_weights(4,zone_idx)))
!  DENSITY.
         dm(zone_idx)=exp(ln10*(log_density(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              log_density(zone_idx-1)*lagrange_interp_weights(2,zone_idx) &
              +log_density(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              log_density(zone_idx+1)*lagrange_interp_weights(4,zone_idx)))
!  DEL (ACTUAL).
!  DEL (RADIATIVE) IS INTERPOLATED, AND DEL IS THE MIN OF DELA,DELR.
         delmi(zone_idx)=mix_phys%del_radiative_mix(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              mix_phys%del_radiative_mix(zone_idx-1)*lagrange_interp_weights(2,zone_idx)+ &
              mix_phys%del_radiative_mix(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              mix_phys%del_radiative_mix(zone_idx+1)*lagrange_interp_weights(4,zone_idx)
!  DEL(ADIABATIC).
         delami(zone_idx)=mix_phys%del_adiabatic_mix(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              mix_phys%del_adiabatic_mix(zone_idx-1)*lagrange_interp_weights(2,zone_idx)+ &
              mix_phys%del_adiabatic_mix(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              mix_phys%del_adiabatic_mix(zone_idx+1)*lagrange_interp_weights(4,zone_idx)
!  D LN RHO/D LN T.
         qdtmi(zone_idx)=mix_phys%qdtm(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              mix_phys%qdtm(zone_idx-1)*lagrange_interp_weights(2,zone_idx)+ &
              mix_phys%qdtm(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              mix_phys%qdtm(zone_idx+1)*lagrange_interp_weights(4,zone_idx)
!  UNLOGGED MASS INTERIOR TO THE INTERFACE.
         hs3(zone_idx)=mass_unlogged(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              mass_unlogged(zone_idx-1)*lagrange_interp_weights(2,zone_idx)+ &
              mass_unlogged(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              mass_unlogged(zone_idx+1)*lagrange_interp_weights(4,zone_idx)
!  SPECIFIC ENERGY GENERATION RATE.
         epsilm(zone_idx)=mix_phys%esumm(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              mix_phys%esumm(zone_idx-1)*lagrange_interp_weights(2,zone_idx)+ &
              mix_phys%esumm(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              mix_phys%esumm(zone_idx+1)*lagrange_interp_weights(4,zone_idx)
!  LUMINOSITY.
         interface_luminosity(zone_idx)=solar_luminosity_cgs*(luminosity(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              luminosity(zone_idx-1)*lagrange_interp_weights(2,zone_idx)+ &
              luminosity(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              luminosity(zone_idx+1)*lagrange_interp_weights(4,zone_idx))
!  LOCAL AVERAGE FORCE OF GRAVITY.
         interface_gravity_factor(zone_idx)=local_gravity(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              local_gravity(zone_idx-1)*lagrange_interp_weights(2,zone_idx)+ &
              local_gravity(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              local_gravity(zone_idx+1)*lagrange_interp_weights(4,zone_idx)
!  opacity.
         opacity_interface(zone_idx)=mix_phys%om(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              mix_phys%om(zone_idx-1)*lagrange_interp_weights(2,zone_idx)+ &
              mix_phys%om(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              mix_phys%om(zone_idx+1)*lagrange_interp_weights(4,zone_idx)
!  specific heat
         specific_heat_interface(zone_idx)=mix_phys%cpm(zone_idx-2)*lagrange_interp_weights(1,zone_idx)+ &
              mix_phys%cpm(zone_idx-1)*lagrange_interp_weights(2,zone_idx)+ &
              mix_phys%cpm(zone_idx)*lagrange_interp_weights(3,zone_idx)+ &
              mix_phys%cpm(zone_idx+1)*lagrange_interp_weights(4,zone_idx)
   30 continue
      do 35 zone_idx = transport_zone_begin,transport_zone_end
         delmi(zone_idx) = min(delmi(zone_idx),delami(zone_idx))
   35 continue
      do zone_idx = 1,num_zones
         es_velocity_coeff1(zone_idx) = 0.0d0
         es_velocity_coeff2(zone_idx) = 0.0d0
         es_shear_coeff(zone_idx) = 0.0d0
         fgsfj(zone_idx) = 0.0d0
         gsf_kippenhahn_coeff(zone_idx) = 0.0d0
         fact1(zone_idx) = 0.0d0
         fact2(zone_idx) = 0.0d0
         fact3(zone_idx) = 0.0d0
         fact4(zone_idx) = 0.0d0
         mu_gradient_richardson_coeff(zone_idx) = 0.0d0
! MHP 06/02
         es_relaxation_factor(zone_idx) = 0.0d0
         theta_mean(zone_idx) = 0.0d0
         difad_shear_coeff1(zone_idx) = 0.0d0
         difad_shear_coeff2(zone_idx) = 0.0d0
      end do
      do zone_idx = 2,num_zones
! ADDED FOR D THETA/DT TERM FROM ZAHN&MAEDER 1998
         fac_local = 2.0d0*(log_radius(zone_idx)+log_radius(zone_idx-1))-0.5d0* &
              (log10(mass_unlogged(zone_idx))+log10(mass_unlogged(zone_idx-1)))-cgl
         theta_mean(zone_idx) = exp(ln10*fac_local)
      end do
!  NOW COMPUTE STRUCTURAL QUANTITIES NEEDED TO EVALUATE VELOCITIES AT
!  ALL INTERFACES.
      cpigi_const = 4.0d0/c4pi/grav_const
      do 40 zone_idx = transport_zone_begin,transport_zone_end
         if(.not.use_diffusion_advection_transport)then
            dlnmu_dlnp = (log10(mix_phys%amum(zone_idx))-log10(mix_phys%amum(zone_idx-1)))/ &
                 (log_pressure(zone_idx)-log_pressure(zone_idx-1))
         else
            dlnmu_dlnp = 0.0d0
         endif
         if(zone_idx.eq.transport_zone_begin)then
            ddel_floor=max(1.0d-3,0.5d0*(mix_phys%del_adiabatic_mix(transport_zone_begin+1)- &
                 mix_phys%delm(transport_zone_begin+1)) &
                 +dlnmu_dlnp)
         else if(zone_idx.eq.transport_zone_end)then
            ddel_floor=max(1.0d-3,0.5d0*(mix_phys%del_adiabatic_mix(transport_zone_end-1)- &
                 mix_phys%delm(transport_zone_end-1)) &
                 +dlnmu_dlnp)
         else
            ddel_floor = 1.0d-3
         endif
         ddel = max(delami(zone_idx)-delmi(zone_idx),ddel_floor)
! MHP 06/02
! ADDED FOR ALTERNATE TREATMENT OF MU GRADIENTS
!         DDELM(I) = 0.5D0*(DELAMI(I)+DELAMI(I-1)-
!     *                     DELMI(I)-DELMI(I-1))
!         DDELM(I) = DDEL
         del_grad_diff_interface(zone_idx) = max(delami(zone_idx)-delmi(zone_idx)+dlnmu_dlnp,ddel_floor)
!         FESTIME(I) = PM(I)/(HGM(I)*DDEL*DM(I)*TM(I))
         es_relaxation_factor(zone_idx) = pm(zone_idx)/(interface_gravity_factor(zone_idx)*ddel*dm(zone_idx))
         pressure_scale_factor = pm(zone_idx)*interface_radius(zone_idx)**2/dm(zone_idx)/hs3(zone_idx)/grav_const
         temp_scale_factor = pressure_scale_factor/delmi(zone_idx)
         eta_factor = 2.0d0*cc23*interface_radius(zone_idx)**3/grav_const/hs3(zone_idx)
!        TTHERM = 8.0d0*CC23*CSIG*TM(I)**3/OPM(I)/DM(I)**2/CPM(I)
         gsf_kippenhahn_coeff(zone_idx)= 8.0d0*pressure_scale_factor*eta_factor/ddel
         ff_factor = pm(zone_idx)/(interface_gravity_factor(zone_idx)*ddel*specific_heat_interface(zone_idx)* &
              dm(zone_idx)*tm(zone_idx))
         fac_local = 2.0d0*eta_factor*ff_factor
         specific_luminosity = interface_luminosity(zone_idx)/hs3(zone_idx)
         es_velocity_coeff1(zone_idx) = fac_local*(specific_luminosity-epsilm(zone_idx))
         es_velocity_coeff2(zone_idx)= -0.5d0*fac_local*specific_luminosity*cpigi_const/dm(zone_idx)
         es_shear_coeff(zone_idx) = specific_luminosity*ff_factor*cc13*cpigi_const*interface_radius(zone_idx)/dm(zone_idx)
!         FES3(I) = EM*FF*CPIGI*(3.0D0*HTSC-RM(I))/DM(I)
         velocity_coeff0(zone_idx) = ff_factor
         velocity_coeff1a(zone_idx) = es_velocity_coeff1(zone_idx)/ff_factor
         velocity_coeff1b(zone_idx) = es_velocity_coeff2(zone_idx)/ff_factor
         velocity_coeff2a(zone_idx) = specific_luminosity*cpigi_const*3.0d0*temp_scale_factor/dm(zone_idx)
         velocity_coeff2b(zone_idx) = -specific_luminosity*cpigi_const*interface_radius(zone_idx)/dm(zone_idx)
         fgsfj(zone_idx) = abs(fac_local*specific_luminosity*temp_scale_factor)*interface_radius(zone_idx)
!  EDDINGTON CIRCULATION VELOCITY IS DEFINED AS
!  VES = FACT1*FACT2*OMEGA**2 (ENDAL AND SOFIA PAPER II).
!  NOTE THAT TO AVOID OVERFLOW 1/(DEL(AD)-DEL)IS SET TO A MAXIMUM OF 10^6.
!  THIS SHOULD ONLY BE AN ISSUE IF YOU HAVE ZERO OVERSHOOT AT THE BOUNDARY
!  OF A CONVECTION ZONE.  CAVEAT EMPTOR.
         fact1(zone_idx) = 1.0d0/qdtmi(zone_idx)/max(delami(zone_idx)-delmi(zone_idx),1.0d-3)
!         FACT2(I) = DELAMI(I)*ALM(I)*(RM(I)**3/CG2/HS3(I)**2)*
!     *              (2.0D0*RM(I)**2*(EPSILM(I)/ALM(I) - 1.0D0/HS3(I))
!    *              - 3.0D0/(C4PI*DM(I)*RM(I)))
         fac_local = delami(zone_idx)/interface_gravity_factor(zone_idx)
         fact2(zone_idx) = fac_local*eta_factor*(epsilm(zone_idx)-specific_luminosity)
         fact6(zone_idx) = -fac_local*specific_luminosity/c4pi/grav_const/dm(zone_idx)/interface_radius(zone_idx)
! GSF VELOCITY IS DEFINED AS
! VGSF = VES*FACT3/Hj**2/(2OMEGA/(D OMEGA/D LNR) + 1)
! WHERE HJ IS THE SPECIFIC ANGULAR MOMENTUM SCALE HEIGHT.
! AND FACT3 IS THE TEMPERATURE SCALE HEIGHT TIMES THE RADIUS.
!  HT = (DEL*D(LN P)/DR)**-1 = P*R**2/(DEL*RHO*GM)
!        FACT3(I) = PM(I)*RM(I)**3/(DELMI(I)*DM(I)*CG*HS3(I))
         fact3(zone_idx) = abs(temp_scale_factor*interface_radius(zone_idx))
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
         fact4(zone_idx) = interface_radius(zone_idx)**2/hs3(zone_idx)/grav_const/dm(zone_idx)* &
              (pm(zone_idx)-radiation_constant_over_3*tm(zone_idx)**4)/sqrt(interface_gravity_factor(zone_idx))
!  DYNAMICAL SHEAR INSTABILITY.
!  AN INTERFACE IS STABLE AGAINST THIS SHEAR WHENEVER THE RICHARDSON NUMBER
!  RICHNO = (RHO/P)*(DEL(AD)-DEL)*g**2/(D OMEGA/D LN R)**2>.25.
         dynamical_shear_omega_limit(zone_idx)=sqrt(dm(zone_idx)/pm(zone_idx)* &
              max(delami(zone_idx) - delmi(zone_idx),1.0d-3)) &
              *2.0d0*interface_gravity_factor(zone_idx)
!  THE DIFFUSIVE SHEAR IS INHIBITED BY MU GRADIENTS.
!  AN INTERFACE IS STABLE AGAINST THE DIFFUSIVE SHEAR IF
!  RICHARDSON NUMBER = FACT5*DMU/MU/(D OMEGA/D LNR)**2 > .25, WHERE
!  FACT5 = RHO/P*(-d lnT/d lnMU)/(del LN P)*G**2 AND
!   -d lnT/d lnMU = (CON-1)/(1+3CON),CON=(a/3)T**4/P
!  GIVEN (1 - CON)P = CGAS*RHO*T/MU
!  1 - CON IS BETA(CORRECTION FOR RADIATION PRESSURE).
         mu_gradient_richardson_coeff(zone_idx) = dm(zone_idx)*interface_gravity_factor(zone_idx)**2/ &
              pm(zone_idx)/ln10/(log_pressure(zone_idx)-log_pressure(zone_idx-1))* &
              (radiation_constant_over_3*tm(zone_idx)**4/pm(zone_idx)-1.0d0)/ &
              (1.0d0 + 3.0d0*radiation_constant_over_3*tm(zone_idx)**4/pm(zone_idx))
! MHP 3/92 ADD VECTOR FOR LOCAL KELVIN-HELMHOLTZ TIME SCALE.
!         TKH(I) = CG*HS3(I)**2/ABS(ALM(I))/RM(I)
   40 continue
! MHP 06/02 ADDED TERMS OF ORDER DW/DR FROM ZAH&MAEDER 1998
      if(use_diffusion_advection_transport)then
         ht_temp_scale_prev = exp(ln10*(log_pressure(transport_zone_begin-1)+ &
              2.0d0*log_radius(transport_zone_begin-1)-log_density(transport_zone_begin-1))) &
              /mass_unlogged(transport_zone_begin-1)/grav_const/mix_phys%del_radiative_mix(transport_zone_begin-1)
         third_deriv_geom_factor(transport_zone_begin-1) = ht_temp_scale_prev
         do zone_idx = transport_zone_begin,transport_zone_end
!         DDEL = MAX(DELAMI(I)-DELMI(I),1.0D-3)
         pressure_scale_factor = pm(zone_idx)*interface_radius(zone_idx)**2/dm(zone_idx)/hs3(zone_idx)/grav_const
         temp_scale_factor = pressure_scale_factor/delmi(zone_idx)
         f1_local = pm(zone_idx)/(interface_gravity_factor(zone_idx)*del_grad_diff_interface(zone_idx)* &
              specific_heat_interface(zone_idx)*dm(zone_idx)*tm(zone_idx))
         f2_local = interface_luminosity(zone_idx)*interface_radius(zone_idx)/hs3(zone_idx)/3.0d0
         f3_local = 0.75d0*hs3(zone_idx)/(cpi*dm(zone_idx)*interface_radius(zone_idx)**3)
         v0_local = -f1_local*f2_local*f3_local
!         V0 = -PM(I)*ALM(I)/(C4PI*CG**2*HS3(I)*DDEL*CPM(I)*
!     *         DM(I)**2*TM(I))
         ff_factor = interface_radius(zone_idx)**3/hs3(zone_idx)
         c1_factor = cc23*interface_radius(zone_idx)**4/grav_const/hs3(zone_idx)
         second_deriv_geom_factor(zone_idx) = c1_factor
         qc1r = 4.0d0*cc23*ff_factor/grav_const* &
              (1.0d0-cpi*dm(zone_idx)*ff_factor)
         dr_local = 10.0d0**log_radius(zone_idx)-10.0d0**log_radius(zone_idx-1)
         qdr_local = (10.0d0**log_density(zone_idx)-10.0d0**log_density(zone_idx-1))/dr_local
         qqc1rr = 8.0d0*ff_factor/grav_const/interface_radius(zone_idx)* &
         (1.0d0-cc13*cpi*ff_factor*(1.0d1*dm(zone_idx)-interface_radius(zone_idx)*qdr_local) &
          + cc13*8.0d0*(cpi*dm(zone_idx)*ff_factor)**2)
! D LN CHI/D LN T = 3 - D LN CAPPA/D LN T
         qchit = 3.0d0 - 0.5d0*(dlnkappa_dlnt(zone_idx)+dlnkappa_dlnt(zone_idx-1))
         qqchitr = (dlnkappa_dlnt(zone_idx-1)-dlnkappa_dlnt(zone_idx))/dr_local
         ht_temp_scale = exp(ln10*(log_pressure(zone_idx)+2.0d0*log_radius(zone_idx)-log_density(zone_idx)))/ &
              mass_unlogged(zone_idx)/grav_const/mix_phys%delm(zone_idx)
         third_deriv_geom_factor(zone_idx) = ht_temp_scale
         ht_temp_scale2 = dr_local/ln10/(log_temperature(zone_idx-1)-log_temperature(zone_idx))
         dhtscale_dr = (abs(ht_temp_scale)-abs(ht_temp_scale_prev))/dr_local
         ht_temp_scale_prev = ht_temp_scale
         mean_dlneps_dlnt = 0.5d0*(dlnepsilon_dlnt(zone_idx)+dlnepsilon_dlnt(zone_idx-1))
         mean_neutrino_fraction = 0.5d0*(neutrino_loss_fraction(zone_idx)+neutrino_loss_fraction(zone_idx-1))
!         FACT7(I)= -V0*(QQCHITR*C1+QCHIT*QC1R)
!         FACT7(I)= -V0*(QQCHITR*C1)
         difad_shear_coeff1(zone_idx)= &
         - f1_local*c1_factor*(mean_neutrino_fraction*mean_dlneps_dlnt + &
              epsilm(zone_idx)*(1.0d0 - mean_neutrino_fraction - qchit))
         velocity_coeff2b(zone_idx) = velocity_coeff2b(zone_idx)+difad_shear_coeff1(zone_idx)/f1_local
!         FACT7(I)=V0*((QHTR-QCHIT)*QC1R-QQCHITR*C1+HTSC*QQC1RR)
!     *            - F1*EPSILM(I)*(HTSC*QC1R + C1*(QETM - QCHIT))
         difad_shear_coeff2(zone_idx) = 0.0d0
!         FACT8(I)=V0*((QHTR-QCHIT)*C1+HTSC*2.0D0*QC1R)
!     *            - F1*EPSILM(I)*HTSC*C1
! Q variables not used
!         Q1 = -V0*(QQCHITR*C1)
!         Q2 = -F1*C1*FNU*QETM
!         Q3 = -F1*C1*EPSILM(I)*(1.0D0 - FNU - QCHIT)
!         Q4 = F1*FL*QCHIT
!         Q5 = -F1*EPSILM(I)*HTSC
         local_flux_factor = interface_luminosity(zone_idx)/c4pi/dm(zone_idx)/interface_radius(zone_idx)**2
         facd2(zone_idx) = f1_local*(local_flux_factor*qchit - epsilm(zone_idx)*temp_scale_factor)
         facd3(zone_idx) = -f1_local*local_flux_factor
!         WRITE(*,911)I,V0,FES3(I),FACT7(I),FACD2(I),FACD3(I),DDELM(I)
!         WRITE(*,912)Q1,Q2,Q3,Q4,Q5
! 911     FORMAT(I5,1P6E12.3)
! 912     FORMAT(5X,1P5E12.3)
         end do
      endif
      return
end subroutine setupv

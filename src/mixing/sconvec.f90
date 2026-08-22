!----------------------------------------------------------------------
! sconvec
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original sconvec.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! SEMI-CONVECTION IS TESTED FOR ALL CONVECTIVE REGIONS IN THIS SUBROUTINE.
! THE REFERENCE FOR THIS SECTION IS :
! CASTELLANI ET AL.(1971)ASTROPHYSICS AND SPACE SCIENCE, VOL. 10, PP.340-349.
!
! AT THE EDGE OF ALL CONVECTIVE REGIONS, THE RADIATIVE TEMPERATURE GRADIENT
! IS COMPUTED WITH BOTH THE LOCAL ABUNDANCES AND THOSE IN THE NEIGHBORING
! CONVECTION ZONE.  IF DEL(RAD)>DEL(AD) WHEN THE ABUNDANCES ARE PERTURBED,
! THEN THE REGION IS CONSIDERED UNSTABLE TO THE GROWTH OF A SEMI-CONVECTIVE
! REGION AND THE CODE COMPUTES ITS EXTENT.
!
! INPUT VARIABLES :
!
! timestep : model timestep (seconds)
! log_density,log_luminosity,log_pressure,log_radius,log_mass,
!   log_temperature : model run of density, luminosity, pressure, radius,
!   mass, and temperature.
! composition : model run of mass fractions of different chemical species;
!   composition(1,...) is hydrogen and composition(3,...) is the metals.
! mixed_zone_bounds, num_zones_mixed : the locations of the edges of
!   convective regions, stored in pairs, and the number of distinct
!   convection zones. E.g. mixed_zone_bounds(3,1) is the id of the
!   bottom shell in the 3rd cz from the center and
!   mixed_zone_bounds(3,2) is the id of the top shell in the same cz.
!
! OUTPUT VARIABLES :
!
! mixed_zone_bounds, num_zones_mixed can be altered by the subroutine.
subroutine sconvec(timestep, composition, log_density, log_luminosity, &
     log_pressure, log_radius, log_mass, log_temperature, num_zones, &
     mixed_zone_bounds, num_zones_mixed, log_teff, ierr)

      use luout_lib
      use const_lib
      use eos_lib
      use kap_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: timestep
      double precision, intent(in) :: composition(15,json)
      double precision, intent(inout) :: log_density(json)
      double precision, intent(in) :: log_luminosity(json), &
           log_pressure(json), log_radius(json), log_mass(json), &
           log_temperature(json)
      integer, intent(in) :: num_zones
      integer, intent(inout) :: mixed_zone_bounds(12,2)
      integer, intent(inout) :: num_zones_mixed
      double precision, intent(in) :: log_teff

! MHP 8/25 Removed unused variables
!      CHARACTER*256 FLAOL, FPUREZ
! MHP 8/25 Removed character file names from common block
! former common/nwlaol/: not used in this file.


      double precision :: ion_fraction(3)
      logical :: only_check_core
      integer :: loop_upper_bound, zone_idx, edge_side
      logical :: up_semiconv_flag, down_semiconv_flag
      integer :: cz_edge_idx, adjacent_radiative_idx
      logical :: want_derivatives, use_conductive_opacity_flag, &
           in_atmosphere, is_convective
      integer :: saha_state
      double precision :: log_luminosity_zone, log_mass_zone, &
           log_pressure_zone, log_temperature_zone, log_density_zone, &
           log_radius_zone
      double precision :: hydrogen_fraction, metal_fraction
      double precision :: pressure_rotation_factor, &
           temperature_rotation_factor
      integer :: current_zone_idx
      double precision :: temperature_k, pressure_cgs, density_cgs, beta, &
           beta_inverse, beta14, specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      double precision :: opacity, log_opacity, dlnkap_dlnrho, dlnkap_dlnt
      double precision :: actual_gradient, radiative_gradient, &
           dgrad_dt_component, dgrad_dp_component, dgrad_dr_component, &
           convective_velocity
      double precision :: boundary_mean_molecular_weight, &
           perturbed_radiative_gradient, max_overshoot_radius
      double precision :: gradient_ratio, local_radiative_gradient
      double precision :: radius_sum
      logical :: reached_max_extent
      double precision :: radius_prev, radius_curr, radius_diff
      integer :: search_begin, search_end, search_step, search_zone_idx, &
           new_edge_idx
      integer :: k_idx, pair_idx

! RUN THROUGH ALL THE CONVECTION ZONES.
!$$$      DO 210 I = 1,NZONE

      integer, intent(out) :: ierr

      ierr = 0

      only_check_core = .true.
      if (only_check_core) then
         loop_upper_bound = 1
      else
         loop_upper_bound = num_zones_mixed
      end if
      do zone_idx = 1, loop_upper_bound

! DETERMINE IF THIS REGION IS A CORE CONVECTION ZONE, SURFACE CZ,
! OR INTERMEDIATE CZ.
! LFLAGU IS T IF SEMICONVECTION ABOVE THE CZ NEEDS TO BE COMPUTED;
! LFLAGD IS T IF SEMICONVECTION BELOW THE CZ NEEDS TO BE COMPUTED.
         if (mixed_zone_bounds(zone_idx,1).eq.1) then
! CONVECTIVE CORE
! CHECK FOR A FULLY CONVECTIVE STAR; SKIP THIS SUBROUTINE IF THERE IS ONE.
            if (mixed_zone_bounds(zone_idx,2).eq.num_zones) return
            up_semiconv_flag = .true.
            down_semiconv_flag = .false.
         else if (mixed_zone_bounds(zone_idx,2).eq.num_zones) then
! CONVECTIVE ENVELOPE
            up_semiconv_flag = .false.
            down_semiconv_flag = .true.
         else
! INTERMEDIATE CONVECTION ZONE (NOT INCLUDING CENTRAL OR SURFACE POINT).
            up_semiconv_flag = .true.
            down_semiconv_flag = .true.
         end if
! CHECK SEMICONVECTION BELOW (K=1) AND ABOVE (K=2) THE CZ.
         do edge_side = 1, 2
! SKIP SEMI-CONVECTION BELOW A CENTRAL CZ AND ABOVE A SURFACE ONE.
            if (edge_side.eq.1.and..not.down_semiconv_flag) goto 200
            if (edge_side.eq.2.and..not.up_semiconv_flag) goto 200
! JMC AND JMR ARE THE LOCATIONS OF THE EDGE OF THE CONVECTIVE ZONE
! AND THE FIRST RADIATIVE POINT OUTSIDE IT RESPECTIVELY.
            cz_edge_idx = mixed_zone_bounds(zone_idx,edge_side)
            if (edge_side.eq.1) adjacent_radiative_idx = cz_edge_idx - 1
            if (edge_side.eq.2) adjacent_radiative_idx = cz_edge_idx + 1
! SET UP FLAGS FOR CALLS TO THE BASIC PHYSICS ROUTINES.
            want_derivatives = .false.
            use_conductive_opacity_flag = .true.
            in_atmosphere = .true.
            saha_state = 0
! USE THE STRUCTURE VARIABLES FOR THE FIRST POINT OUTSIDE THE CZ
! AND THE ABUNDANCES OF THE CZ TO DETERMINE THE MODIFIED MU AND
! RADIATIVE AND ADIABATIC TEMPERATURE GRADIENTS.
            log_luminosity_zone = log_luminosity(adjacent_radiative_idx)
            log_mass_zone = log_mass(adjacent_radiative_idx)
            log_pressure_zone = log_pressure(adjacent_radiative_idx)
            log_temperature_zone = log_temperature(adjacent_radiative_idx)
            log_density_zone = log_density(adjacent_radiative_idx)
! MHP 10/02 added definition of RL
            log_radius_zone = log_radius(adjacent_radiative_idx)
            hydrogen_fraction = composition(1,cz_edge_idx)
            metal_fraction = composition(3,cz_edge_idx)
!*** ADD ROTATION VECTORS ***
            pressure_rotation_factor = 1.0d0
            temperature_rotation_factor = 1.0d0
            current_zone_idx = cz_edge_idx
            call eos_get(log_temperature_zone, temperature_k, &
                 log_pressure_zone, pressure_cgs, log_density_zone, &
                 density_cgs, hydrogen_fraction, metal_fraction, beta, &
                 beta_inverse, beta14, ion_fraction, specific_gas_constant, &
                 ion_mean_weight_inverse, electron_mean_weight_inverse, &
                 electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
                 specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
                 dlnrho_dlnp_dt, adiabatic_gradient_dt, &
                 adiabatic_gradient_dp, specific_heat_cp_dt, &
                 specific_heat_cp_dp, want_derivatives, in_atmosphere, &
                 saha_state, composition_at_zone=composition(:,cz_edge_idx))
! DBG 12/95 GET OPACITY
            call kap_get(log_density_zone, log_temperature_zone, &
                 hydrogen_fraction, metal_fraction, opacity, log_opacity, &
                 dlnkap_dlnrho, dlnkap_dlnt, ion_fraction)
            call tpgrad(log_temperature_zone, temperature_k, &
                 log_pressure_zone, pressure_cgs, density_cgs, &
                 log_radius_zone, log_mass_zone, log_luminosity_zone, &
                 opacity, dlnrho_dlnt, dlnrho_dlnp, dlnkap_dlnt, &
                 dlnkap_dlnrho, specific_heat_cp, actual_gradient, &
                 radiative_gradient, adiabatic_gradient, dlnrho_dlnt_dt, &
                 dlnrho_dlnp_dt, adiabatic_gradient_dt, &
                 adiabatic_gradient_dp, dgrad_dt_component, &
                 dgrad_dp_component, dgrad_dr_component, specific_heat_cp_dt, &
                 specific_heat_cp_dp, convective_velocity, want_derivatives, &
                 is_convective, pressure_rotation_factor, &
                 temperature_rotation_factor, log_teff, ierr)
            if (ierr /= 0) return
! SKIP IF ZONE IS STABLE WITH THE CORE COMPOSITION.
            if (radiative_gradient.lt.adiabatic_gradient) goto 200
! STORE MEAN MOLECULAR WEIGHT, ADJUSTED RADIATIVE TEMPERATURE GRADIENT,
! AND THE QUANTITY (DELR - DELA)/DELR.
            boundary_mean_molecular_weight = ion_mean_weight_inverse + &
                 electron_mean_weight_inverse
            perturbed_radiative_gradient = radiative_gradient
            max_overshoot_radius = 1.0d0 - adiabatic_gradient/ &
                 radiative_gradient
! REPEAT CALL WITH THE LOCAL COMPOSITION.
            hydrogen_fraction = composition(1,adjacent_radiative_idx)
            metal_fraction = composition(3,adjacent_radiative_idx)
            current_zone_idx = adjacent_radiative_idx
            call eos_get(log_temperature_zone, temperature_k, &
                 log_pressure_zone, pressure_cgs, log_density_zone, &
                 density_cgs, hydrogen_fraction, metal_fraction, beta, &
                 beta_inverse, beta14, ion_fraction, specific_gas_constant, &
                 ion_mean_weight_inverse, electron_mean_weight_inverse, &
                 electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
                 specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
                 dlnrho_dlnp_dt, adiabatic_gradient_dt, &
                 adiabatic_gradient_dp, specific_heat_cp_dt, &
                 specific_heat_cp_dp, want_derivatives, in_atmosphere, &
                 saha_state, &
                 composition_at_zone=composition(:,adjacent_radiative_idx))
! DBG 12/95 GET OPACITY
            call kap_get(log_density_zone, log_temperature_zone, &
                 hydrogen_fraction, metal_fraction, opacity, log_opacity, &
                 dlnkap_dlnrho, dlnkap_dlnt, ion_fraction)
            call tpgrad(log_temperature_zone, temperature_k, &
                 log_pressure_zone, pressure_cgs, density_cgs, &
                 log_radius_zone, log_mass_zone, log_luminosity_zone, &
                 opacity, dlnrho_dlnt, dlnrho_dlnp, dlnkap_dlnt, &
                 dlnkap_dlnrho, specific_heat_cp, actual_gradient, &
                 radiative_gradient, adiabatic_gradient, dlnrho_dlnt_dt, &
                 dlnrho_dlnp_dt, adiabatic_gradient_dt, &
                 adiabatic_gradient_dp, dgrad_dt_component, &
                 dgrad_dp_component, dgrad_dr_component, specific_heat_cp_dt, &
                 specific_heat_cp_dp, convective_velocity, want_derivatives, &
                 is_convective, pressure_rotation_factor, &
                 temperature_rotation_factor, log_teff, ierr)
            if (ierr /= 0) return
            log_density(adjacent_radiative_idx) = log_density_zone
! FDEL IS THE RATIO OF THE GRADIENTS WITH THE OLD COMP AND NEW ONE.
            gradient_ratio = perturbed_radiative_gradient/radiative_gradient
            local_radiative_gradient = radiative_gradient
! RMAX IS THE MAXIMUM RADIAL DISTANCE THAT OVERSHOOT MAY PENETRATE, AND
! IS THE TIMESTEP(DELTS) MULTIPLIED BY THE VELOCITY OF PROPAGATION OF
! THE INSTABILITY (VP IN EQUATION 5PRIME, P. 347).
            max_overshoot_radius = max_overshoot_radius*(log_luminosity( &
                 cz_edge_idx)*solar_luminosity_cgs/(10.0d0*c4pi* &
                 dexp(ln10*log_pressure(cz_edge_idx))))* &
                 (timestep/dexp(ln10*(log_radius(cz_edge_idx)+ &
                 log_radius(cz_edge_idx))))
! INITIALIZE SUM.
            radius_sum = 0.0d0
! LMAX IS SET TO TRUE IF THE OVERSHOOT REACHES ITS MAXIMUM EXTENT.
            reached_max_extent = .false.
            if (edge_side.eq.1) then
               radius_prev = exp(ln10*log_radius(adjacent_radiative_idx+1))
               search_begin = adjacent_radiative_idx - 1
               search_end = 1
               search_step = -1
            else if (edge_side.eq.2) then
               radius_prev = exp(ln10*log_radius(adjacent_radiative_idx-1))
               search_begin = adjacent_radiative_idx + 1
               search_end = num_zones
               search_step = 1
            end if
            do search_zone_idx = search_begin, search_end, search_step
! TEST ON MAXIMUM OVERSHOOTING LIMIT IN RADIUS
               if (edge_side.eq.1) radius_curr = &
                    exp(ln10*log_radius(search_zone_idx+1))
               if (edge_side.eq.2) radius_curr = &
                    exp(ln10*log_radius(search_zone_idx-1))
               radius_diff = radius_curr - radius_prev
               radius_sum = radius_sum + (1.0d0 - &
                    boundary_mean_molecular_weight/( &
                    ion_mean_weight_inverse+electron_mean_weight_inverse))* &
                    radius_diff
               radius_prev = radius_curr
               if (radius_sum.gt.max_overshoot_radius) goto 33
! TEST ON RESCALED RAD. GRADIENT FOR CONVECTION
               log_luminosity_zone = log_luminosity(search_zone_idx)
               log_mass_zone = log_mass(search_zone_idx)
               log_pressure_zone = log_pressure(search_zone_idx)
               log_temperature_zone = log_temperature(search_zone_idx)
               log_density_zone = log_density(search_zone_idx)
               hydrogen_fraction = composition(1,search_zone_idx)
               metal_fraction = composition(3,search_zone_idx)
               current_zone_idx = search_zone_idx
               call eos_get(log_temperature_zone, temperature_k, &
                    log_pressure_zone, pressure_cgs, log_density_zone, &
                    density_cgs, hydrogen_fraction, metal_fraction, beta, &
                    beta_inverse, beta14, ion_fraction, &
                    specific_gas_constant, ion_mean_weight_inverse, &
                    electron_mean_weight_inverse, &
                    electron_degeneracy_parameter, dlnrho_dlnt, &
                    dlnrho_dlnp, specific_heat_cp, adiabatic_gradient, &
                    dlnrho_dlnt_dt, dlnrho_dlnp_dt, adiabatic_gradient_dt, &
                    adiabatic_gradient_dp, specific_heat_cp_dt, &
                    specific_heat_cp_dp, want_derivatives, in_atmosphere, &
                    saha_state, &
                    composition_at_zone=composition(:,search_zone_idx))
! DBG 12/95 GET OPACITY
               call kap_get(log_density_zone, log_temperature_zone, &
                    hydrogen_fraction, metal_fraction, opacity, &
                    log_opacity, dlnkap_dlnrho, dlnkap_dlnt, ion_fraction)
               call tpgrad(log_temperature_zone, temperature_k, &
                    log_pressure_zone, pressure_cgs, density_cgs, &
                    log_radius_zone, log_mass_zone, log_luminosity_zone, &
                    opacity, dlnrho_dlnt, dlnrho_dlnp, dlnkap_dlnt, &
                    dlnkap_dlnrho, specific_heat_cp, actual_gradient, &
                    radiative_gradient, adiabatic_gradient, dlnrho_dlnt_dt, &
                    dlnrho_dlnp_dt, adiabatic_gradient_dt, &
                    adiabatic_gradient_dp, dgrad_dt_component, &
                    dgrad_dp_component, dgrad_dr_component, &
                    specific_heat_cp_dt, specific_heat_cp_dp, &
                    convective_velocity, want_derivatives, is_convective, &
                    pressure_rotation_factor, temperature_rotation_factor, &
                    log_teff, ierr)
               if (ierr /= 0) return
               log_density(search_zone_idx) = log_density_zone
! EXIT IF ZONE IS RADIATIVELY STABLE.
               if (gradient_ratio*radiative_gradient.lt.adiabatic_gradient) &
                    goto 34
   32       continue
            end do
            if (edge_side.eq.1) search_zone_idx = 0
            if (edge_side.eq.2) search_zone_idx = num_zones + 1
   33       reached_max_extent = .true.
   34       continue
! ASSIGN THE NEW EDGE LOCATION TO MXZONE.
            if (edge_side.eq.1) then
               mixed_zone_bounds(zone_idx,edge_side) = search_zone_idx + 1
               new_edge_idx = search_zone_idx + 1
            else if (edge_side.eq.2) then
               mixed_zone_bounds(zone_idx,edge_side) = search_zone_idx - 1
               new_edge_idx = search_zone_idx - 1
            end if
            write(short_file_unit,601) cz_edge_idx, new_edge_idx, &
                 reached_max_extent, perturbed_radiative_gradient, &
                 local_radiative_gradient
  601       format(1x,'CZ OLD EDGE ',i3,' EXTENDED TO-', &
                 i3,' LIMIT=',l1,' RAD.GRADS-IN/OUT',2f8.4)
  200    continue
         end do
  210 continue
      end do
!  CHECK FOR MERGERS OF NEARBY CONVECTION ZONES CAUSED BY SEMI-CONVECTION.
      if (num_zones_mixed.le.1) return
      k_idx = 1
   85 continue
!  CHECK IF 'TOP' OF ONE REGION IS ABOVE 'BOTTOM' OF THE NEXT ONE.
      if (mixed_zone_bounds(k_idx,2).gt.mixed_zone_bounds(k_idx+1,1)) then
!  IF THIS OCCURS, TWO CONVECTION ZONES HAVE MERGED.
         write(short_file_unit,93) ((mixed_zone_bounds(zone_idx,pair_idx), &
              pair_idx=1,2),zone_idx=k_idx,k_idx+1), &
              mixed_zone_bounds(k_idx,1), mixed_zone_bounds(k_idx+1,2)
   93    format(2x,'MIXED ZONES MERGED DUE TO SEMICONVECTION' &
              /2x,'OLD',2('[',i3,'-',i3,']'), &
              ' NEW','[',i3,'-',i3,']')
         mixed_zone_bounds(k_idx+1,1) = mixed_zone_bounds(k_idx,1)
         do zone_idx = k_idx, num_zones_mixed-1
            do pair_idx = 1, 2
               mixed_zone_bounds(zone_idx,pair_idx) = &
                    mixed_zone_bounds(zone_idx+1,pair_idx)
   95       continue
            end do
   90    continue
         end do
         num_zones_mixed = num_zones_mixed - 1
         if (k_idx.le.num_zones_mixed-1) then
            goto 85
         else
            return
         end if
      end if
      k_idx = k_idx + 1
      if (k_idx.le.num_zones_mixed-1) goto 85

      return
end subroutine sconvec

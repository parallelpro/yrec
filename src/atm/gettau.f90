!----------------------------------------------------------------------
! gettau
!----------------------------------------------------------------------
! G Somers, 3/17
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original gettau.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This routine devises a new, systematic method for determining the
! convective overturn timescale. It will calculate the overturn
! timescale one pressure scale height above the BCZ for whichever of
! the following three conditions holds.
!
! 1) 1 PSCA above the BCZ is in the interior model.
! 2) 1 PSCA above the BCZ is in the envelope model.
! 3) The BCZ is in the interior, but 1 PSCA above is in the envelope.
!
! If none of these are satisfied, TAUCZ is set to 0.0.
!
! If 1): Call ENVINT to determine PPHOT, TAUINT(NEW) to determine TAUCZ
! If 2): Call ENVINT to determine PPHOT, which calls TAUCAL for TAUCZ.
! If 3): Stitch together interior and envelope.
subroutine gettau(composition, log_radius, log_pressure, log_density, &
     enclosed_mass, log_temperature, fp, ft, log_teff, log_total_mass, &
     log_luminosity_lsun, num_zones, convective_flag, radius_at_bcz)
      use turnover_lib
      use scrtch_lib
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: log_radius(json), log_pressure(json), &
           log_density(json)
      double precision, intent(in) :: enclosed_mass(json)
      double precision, intent(inout) :: log_temperature(json)
      double precision, intent(in) :: fp(json), ft(json)
      double precision, intent(in) :: log_teff, log_total_mass, &
           log_luminosity_lsun
      integer, intent(in) :: num_zones
      logical, intent(in) :: convective_flag(json)
      double precision, intent(out) :: radius_at_bcz


! common/envgen/: not used in this file. Naming matches wrtmod.f90.
      double precision :: atm_step_size, envelope_step_size
      logical :: envelope_generation_flag
      common/envgen/ atm_step_size, envelope_step_size, envelope_generation_flag

! common/intatm/: not used in this file. Naming matches envint.f90.
      double precision :: atm_error_tol, atm_step_initial, atm_step_begin, &
           atm_step_min, atm_step_max
      common/intatm/ atm_error_tol, atm_step_initial, atm_step_begin, &
           atm_step_min, atm_step_max

! common/intenv/: not used in this file. Naming matches envint.f90.
      double precision :: env_error_tol, env_step_begin, env_step_min, &
           env_step_max
      common/intenv/ env_error_tol, env_step_begin, env_step_min, env_step_max




! common/envstruct/: env_log10_radius/env_log10_pressure/
! env_log10_density/env_log10_mass/env_log10_temperature/
! env_convective_velocity/env_gradients/num_env_points are used here.
! Naming matches envint.f90.
      double precision :: env_log10_pressure(json), env_log10_temperature(json), &
           env_log10_mass(json), env_log10_density(json), env_log10_radius(json), &
           env_hydrogen_fraction(json), env_metal_fraction(json)
      logical :: env_convective_flag(json)
      double precision :: env_gradients(3,json), env_convective_velocity(json), &
           env_beta(json)
      double precision :: env_gamma1(json), env_specific_heat_cp(json), &
           env_ion_fraction(3,json)
      double precision :: env_opacity(json), env_luminosity(json), &
           env_dlnrho_dlnt(json)
      integer :: num_env_points
      common/envstruct/ env_log10_pressure, env_log10_temperature, &
           env_log10_mass, env_log10_density, env_log10_radius, &
           env_hydrogen_fraction, env_metal_fraction, env_convective_flag, &
           env_gradients, env_convective_velocity, env_beta, &
           env_gamma1, env_specific_heat_cp, env_ion_fraction, &
           env_opacity, env_luminosity, env_dlnrho_dlnt, num_env_points


      double precision :: combined_composition(15,json)
      double precision :: combined_radius(json), combined_pressure(json), &
           combined_density(json), combined_mass(json), combined_temperature(json), &
           combined_gravity(json), combined_velocity(json)
      double precision :: combined_grad1(json), combined_grad2(json)
      logical :: combined_convective_flag(json)
      double precision :: envint_dummy1(4), envint_dummy2(3), &
           envint_dummy3(3), envint_dummy4(3)
      integer :: katm, kenv, ksaha

! common/spots/: spot_filling_factor/spot_temp_contrast are used here.
! Naming matches wrtmod.f90.
      double precision :: spot_filling_factor, spot_temp_contrast
      logical :: spot_depth_varies
      common/spots/ spot_filling_factor, spot_temp_contrast, spot_depth_varies

      save

      integer :: zone_index, species_index, bcz_top_zone
      logical :: fully_convective_flag
      double precision :: pressure_diff_check
      double precision :: fp_surface, ft_surface
      integer :: envint_unused_flag
      double precision :: luminosity_linear
      integer :: ixx_flag
      logical :: print_flag, surface_bc_flag, pulsation_output_flag
      double precision :: hydrogen_fraction, metal_fraction
      double precision :: log_radius_surface, log_gravity_surface
      double precision :: pressure_limit
      double precision :: spot_adjusted_log_teff
      integer :: combined_num_points

! TAUCZ = 0.0
      turnover%convective_turnover_timescale = 0.0
      turnover%pphot = 0.0

! Check if 1 PSCA above BCZ is within envelope. If shell_diag%so, only the interior
! model should be considered for TAUCZ. Set LCALCENV = .FALSE.. If not,
! have ENVINT calculate the full structure, and stitch the envelope
! and interior together.
!
! Determine where the BCZ is.
      do 71 zone_index = num_zones-1,1,-1
         if (.not.convective_flag(zone_index)) goto 81
   71 continue
      fully_convective_flag = .true.
      zone_index = 0
   81 bcz_top_zone = zone_index + 1
! Compare the pressure there to the surface pressure.
      pressure_diff_check = log_pressure(bcz_top_zone)-log_pressure(num_zones)
! IF CHKPRS > 1.0, AT LEAST 3 PSCAS, SO LOOK INTERIOR.
! IF CHKPRS < 1.0, 1 PSCA UP MIGHT BE IN ENV. STITCH TOGETHER.
! FINALLY, IF LNEWTCZ = .FALSE., MAKE SURE TO CALC AN ENV.
      if (pressure_diff_check.lt.1.0.or..not.use_new_turnover_timescale) then
         calc_envelope_flag = .true.
      else
         calc_envelope_flag = .false.
      endif
!
! CALL ENVINT
!
      fp_surface = fp(num_zones)
      ft_surface = ft(num_zones)
      envint_unused_flag = 0
      luminosity_linear = dexp(ln10*log_luminosity_lsun)
      katm = 0
      kenv = 0
      ksaha = 0
      ixx_flag=0
      print_flag=.true.
      surface_bc_flag = .false.
      pulsation_output_flag = .false.
      hydrogen_fraction = composition(1,num_zones)
      metal_fraction = composition(3,num_zones)
      log_radius_surface = 0.5D0*(log_luminosity_lsun + log10_solar_luminosity &
           - 4.0D0*log_teff - c4pil - csigl)
      log_gravity_surface = cgl + log_total_mass - log_radius_surface - &
           log_radius_surface
      pressure_limit = log_pressure(num_zones)
! G Somers 10/14, FOR SPOTTED RUNS, FIND THE
! PRESSURE AT THE AMBIENT TEMPERATURE ATEFFL
      if (convective_flag(num_zones).and.spot_filling_factor.ne.0.0.and. &
           spot_temp_contrast.ne.1.0) then
         spot_adjusted_log_teff = log_teff - 0.25*log10(spot_filling_factor* &
              spot_temp_contrast**4.0+1.0-spot_filling_factor)
      else
         spot_adjusted_log_teff = log_teff
      endif
      call envint(luminosity_linear,fp_surface,ft_surface,log_gravity_surface, &
           log_total_mass,ixx_flag,print_flag,surface_bc_flag, &
           pressure_limit,log_radius_surface,spot_adjusted_log_teff, &
           hydrogen_fraction,metal_fraction,envint_dummy1,envint_unused_flag, &
           katm,kenv,ksaha, &
           envint_dummy2,envint_dummy3,envint_dummy4,pulsation_output_flag)
! PPHOT is now set, and structure variables are caluclated if
! LCALCENV = .TRUE..
!
! IF LNEWTCZ IS FALSE, THEN TAUCZ WILL HAVE BEEN SET IN ENVINT, IF THE
! TOP OF THE INTERIOR MODEL IS RADIATIVE. CHECK IF TAUCZ = 0.0. IF NOT,
! THEN GO TO THE END.
!
      if (turnover%convective_turnover_timescale.ne.0.0) goto 100
!
! COLLECT THE NECESSARY STRUCTURE VARIABLES INTO DUMMY VECTORS.
      combined_num_points = num_zones
      do zone_index=1,combined_num_points
         do species_index=1,15
            combined_composition(species_index,zone_index) = &
                 composition(species_index,zone_index)
         enddo
         combined_radius(zone_index) = log_radius(zone_index)
         combined_pressure(zone_index) = log_pressure(zone_index)
         combined_density(zone_index) = log_density(zone_index)
         combined_mass(zone_index) = enclosed_mass(zone_index)
         combined_gravity(zone_index) = enclosed_mass(zone_index)* &
              exp(ln10*(cgl-2.0D0*log_radius(zone_index)))
         combined_temperature(zone_index) = log_temperature(zone_index)

         combined_velocity(zone_index) = shell_diag%svel(zone_index)
         combined_grad1(zone_index) = shell_diag%del_grad(1,zone_index)
         combined_grad2(zone_index) = shell_diag%del_grad(3,zone_index)
         combined_convective_flag(zone_index) = convective_flag(zone_index)
      enddo
!
      if (calc_envelope_flag.and.use_new_turnover_timescale) then
! IF CHKPRS < 1, THEN STITCH THE ENVELOPE ONTO THE INTERIOR.
! ENVELOPE WAS JUST INTEGRATED IN ENVINT ABOVE, SO USE THAT RUN.
! THIS CODE BORROWED FROM STITCH.F.
         combined_num_points = num_zones+num_env_points-1
         do zone_index=num_zones+1,combined_num_points
            do species_index=1,15
               combined_composition(species_index,zone_index) = &
                    composition(species_index,num_zones+1)
            enddo
            combined_radius(zone_index) = env_log10_radius(zone_index-num_zones+1)
            combined_pressure(zone_index) = env_log10_pressure(zone_index-num_zones+1)
            combined_density(zone_index) = env_log10_density(zone_index-num_zones+1)
            combined_mass(zone_index) = exp(ln10*(log_total_mass+ &
                 env_log10_mass(zone_index-num_zones+1)))
            combined_gravity(zone_index) = combined_mass(zone_index)* &
                 exp(ln10*(cgl-2.0D0*combined_radius(zone_index)))
            combined_temperature(zone_index) = env_log10_temperature(zone_index-num_zones+1)

            combined_velocity(zone_index) = env_convective_velocity(zone_index-num_zones+1)
            combined_grad1(zone_index) = env_gradients(1,zone_index-num_zones+1)
            combined_grad2(zone_index) = env_gradients(2,zone_index-num_zones+1)
            combined_convective_flag(zone_index) = env_convective_flag(zone_index-num_zones+1)
         enddo
      endif
! CALL TAUINT
      if (use_new_turnover_timescale) then
!          CALL TAUINTNEW(HCOMPF,HS2,HSF,LCF,HRF,HPF,HDF,HGF,MM,M,HVF,
!      *                  DELF1,DELF2,HSTOT,RBCZ)  ! KC 2025-05-31
         call tauintnew(combined_mass,combined_convective_flag,combined_radius, &
              combined_pressure,combined_density,combined_gravity, &
              combined_num_points,num_zones,combined_velocity, &
              combined_grad1,combined_grad2,radius_at_bcz)
      else
!          CALL TAUINT(HCOMPF,HS2,HSF,LCF,HRF,HPF,HDF,HGF,MM,HVF,
!      *               DELF1,DELF2,HSTOT)  ! KC 2025-05-31
         call tauint(combined_mass,combined_convective_flag,combined_radius, &
              combined_pressure,combined_density,combined_gravity, &
              combined_num_points,combined_velocity,combined_grad1,combined_grad2)
      endif
! RETURN FULL FUNCTIONALITY TO ENVINT
 100  continue
      calc_envelope_flag = .true.
      return
end subroutine gettau

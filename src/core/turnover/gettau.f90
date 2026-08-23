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
      use atm_lib
      use envint_lib, only: atm_get
      use envstruct_lib
      use star_info_lib, only: star, i_grad_ad, i_grad_rad
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










      double precision :: combined_composition(15,json)
      double precision :: combined_radius(json), combined_pressure(json), &
           combined_density(json), combined_mass(json), combined_temperature(json), &
           combined_gravity(json), combined_velocity(json)
      double precision :: combined_grad1(json), combined_grad2(json)
      logical :: combined_convective_flag(json)
      double precision :: atm_get_dummy1(4), atm_get_dummy2(3), &
           atm_get_dummy3(3), atm_get_dummy4(3)
      integer :: katm, kenv, ksaha
      integer :: zone_index, species_index, bcz_top_zone
      logical :: fully_convective_flag
      double precision :: pressure_diff_check
      double precision :: fp_surface, ft_surface
      integer :: atm_get_unused_flag
      double precision :: luminosity_linear
      integer :: ixx_flag
      logical :: print_flag, surface_bc_flag, pulsation_output_flag
      double precision :: hydrogen_fraction, metal_fraction
      double precision :: log_radius_surface, log_gravity_surface
      double precision :: pressure_limit
      double precision :: spot_adjusted_log_teff
      integer :: combined_num_points

! TAUCZ = 0.0
      star%turnover%convective_turnover_timescale = 0.0
      star%turnover%pphot = 0.0

! Check if 1 PSCA above BCZ is within envelope. If so, only the interior
! model should be considered for TAUCZ. Set LCALCENV = .FALSE.. If not,
! have ENVINT calculate the full structure, and stitch the envelope
! and interior together.
!
! Determine where the BCZ is.
      do zone_index = num_zones-1,1,-1
         if (.not.convective_flag(zone_index)) exit
      end do
      if (zone_index < (1)) then
      fully_convective_flag = .true.
      zone_index = 0
      end if
      bcz_top_zone = zone_index + 1
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
      atm_get_unused_flag = 0
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
      call atm_get(luminosity_linear,fp_surface,ft_surface,log_gravity_surface, &
           log_total_mass,ixx_flag,print_flag,surface_bc_flag, &
           pressure_limit,log_radius_surface,spot_adjusted_log_teff, &
           hydrogen_fraction,metal_fraction,atm_get_dummy1,atm_get_unused_flag, &
           katm,kenv,ksaha, &
           atm_get_dummy2,atm_get_dummy3,atm_get_dummy4,pulsation_output_flag)
! PPHOT is now set, and structure variables are caluclated if
! LCALCENV = .TRUE..
!
! IF LNEWTCZ IS FALSE, THEN TAUCZ WILL HAVE BEEN SET IN ENVINT, IF THE
! TOP OF THE INTERIOR MODEL IS RADIATIVE. CHECK IF TAUCZ = 0.0. IF NOT,
! THEN GO TO THE END.
!
      if (star%turnover%convective_turnover_timescale.eq.0.0) then
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

         combined_velocity(zone_index) = star%diag%svel(zone_index)
         combined_grad1(zone_index) = star%diag%del_grad(i_grad_rad,zone_index)
         combined_grad2(zone_index) = star%diag%del_grad(i_grad_ad,zone_index)
         combined_convective_flag(zone_index) = convective_flag(zone_index)
      enddo
!
      if (calc_envelope_flag.and.use_new_turnover_timescale) then
! IF CHKPRS < 1, THEN STITCH THE ENVELOPE ONTO THE INTERIOR.
! ENVELOPE WAS JUST INTEGRATED IN ENVINT ABOVE, SO USE THAT RUN.
! THIS CODE BORROWED FROM STITCH.F.
         combined_num_points = num_zones+env_struct%num_env_points-1
         do zone_index=num_zones+1,combined_num_points
            do species_index=1,15
               combined_composition(species_index,zone_index) = &
                    composition(species_index,num_zones+1)
            enddo
            combined_radius(zone_index) = env_struct%env_log10_radius(zone_index-num_zones+1)
            combined_pressure(zone_index) = env_struct%env_log10_pressure(zone_index-num_zones+1)
            combined_density(zone_index) = env_struct%env_log10_density(zone_index-num_zones+1)
            combined_mass(zone_index) = exp(ln10*(log_total_mass+ &
                 env_struct%env_log10_mass(zone_index-num_zones+1)))
            combined_gravity(zone_index) = combined_mass(zone_index)* &
                 exp(ln10*(cgl-2.0D0*combined_radius(zone_index)))
            combined_temperature(zone_index) = env_struct%env_log10_temperature(zone_index-num_zones+1)

            combined_velocity(zone_index) = env_struct%env_convective_velocity(zone_index-num_zones+1)
            combined_grad1(zone_index) = env_struct%env_gradients(1,zone_index-num_zones+1)
            combined_grad2(zone_index) = env_struct%env_gradients(2,zone_index-num_zones+1)
            combined_convective_flag(zone_index) = env_struct%env_convective_flag(zone_index-num_zones+1)
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
      end if
      calc_envelope_flag = .true.
      return
end subroutine gettau

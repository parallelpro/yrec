!----------------------------------------------------------------------
! calcad
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original calcad.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! CALCAD: links together and integrates sound speed profiles, for
! acoustic depth to convection zone, acoustic depth of atmosphere.
! INPUTS:
!      log_radius: array of radii for the interiors solution
!      envelope_cz_log_radius: physical radius of the convection zone
!           (logged), calculated in wrtout
!      num_shells: last interior point
!      log_density: interior's density
!      log_pressure: interior pressure
!      log_temperature: interior temperature
!      log_luminosity_lsun: logged luminosity
!      shape_factor_fp: rotational distortion, needed for envelope
!           integration
!      shape_factor_ft: rotational distortion, needed for envelope
!           integration
!      log_total_mass: needed for ENVINT
!      log_teff: logged effective temperature
!      composition: vector of compositions for interior
!      age_gyr: age
!      envelope_cz_bottom_index: index of the bottom of the surface
!           convection zone
!
! OUTPUTS:
!       Adjusts the values in the common block for acoustic depth:
!            acoustic_depth_cz_fraction = acoustic depth to CZ / sound
!                 travel time from surface-center
!            acoustic_depth_cz_seconds = acoustic depth to CZ (seconds)
!            acoustic_crossing_time_seconds = sound travel time from surface
!                 to center (seconds)
!            eos_adiabatic_gradient = adiabatic gradient from OPAL 2006
!                 EOS or SCV EOS
!
! COMMON BLOCKED THINGS:
!
!        calcad_ageout_output_active = Logical flag for output of CALCAD
!             files at AGEOUT ages
!        AGEOUT(5) = Ages for which you want output. Hard coded in
!             parmin
!        IACAT =
!        IJLAST = Again, for output
!        LJLAST = Again, for output
!        LJWRT = ouput on/off toggle
!        LADON = toggle on calcad imposed 4d interpolation in opacity
subroutine calcad(log_radius, envelope_cz_log_radius, num_shells, &
     log_density, log_pressure, log_temperature, log_luminosity_lsun, &
     shape_factor_fp, shape_factor_ft, log_total_mass, &
!      *                  LPRT, TEFFL, HCOMP, NKK, DAGE, DDAGE, JENV)  ! KC 2025-05-31
     log_teff, composition, age_gyr, envelope_cz_bottom_index)
      use star_info_lib, only: star
      use eos_lib
      use atm_lib
      use envint_lib, only: atm_get
      use atm_table_lib
      use star_info_lib
      use envstruct_lib
      use luout_lib
      use const_lib
      use numerics_lib
      implicit none
      integer, parameter :: nts = 63, nps = 76

      double precision, intent(in) :: log_radius(json)
      double precision, intent(in) :: envelope_cz_log_radius
      integer, intent(in) :: num_shells
      double precision, intent(in) :: log_density(json), log_pressure(json)
      double precision, intent(in) :: log_temperature(json)
      double precision, intent(in) :: log_luminosity_lsun
      double precision, intent(in) :: shape_factor_fp(json), &
           shape_factor_ft(json)
      double precision, intent(in) :: log_total_mass
      double precision, intent(in) :: log_teff
      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: age_gyr
      integer, intent(in) :: envelope_cz_bottom_index












! former common/llot95/: this file declared only a single unused
! scalar for a block whose canonical layout (established in
! getopal95.f90) is a much larger set of OPAL95 pure-Z-table arrays
! (opacity_table%opal95_fixed_z_opacity/opal95_fixed_z, now
! use-associated) -- the single-scalar declaration was preserved
! exactly from the original (COMMON/LLOT95/ ZTAB) while this block was
! still raw COMMON, but now that it's a module, the mismatched local
! layout can't be carried forward: it maps onto (an unused corner of)
! the same canonical type instead. Never read or set in this file.








! G Somers END
! --- locals ---
      double precision :: spline_radius_neighborhood(7), &
           spline_speed_neighborhood(7), spline_second_deriv(7), &
           spline_interp_value(1)
      double precision :: atm_get_dummy1(4), atm_get_dummy2(3), &
           atm_get_dummy3(3), atm_get_dummy4(3)
      double precision :: hydrogen_fraction, metal_fraction, log10_gravity, &
           log10_radius_local, log10_pressure_limit, luminosity_linear
      double precision :: atm_step_begin_saved, atm_step_min_saved, &
           atm_step_max_saved, env_step_begin_saved, env_step_min_saved, &
           env_step_max_saved
      double precision :: star_radius_cm(json), star_temperature_1e6k(json), &
           star_density_cgs(json), star_hydrogen_fraction(json), &
           star_pressure_cgs(json), star_metal_fraction(json), &
           star_inverse_sound_speed(json)
      double precision :: star_inverse_sound_speed_to_cz(json), &
           star_radius_to_cz(json)
!      REAL*8 ATMOSX, ATMOST(JSON),ATMOSD(JSON),ATMOSP(JSON),ADELAD(JSON),ATMOSC(JSON),
!     *       ATGAM1(JSON), ATMOSR(JSON),GM1(JSON)
      double precision :: atmosphere_hydrogen_fraction, local_gamma1(json)
      integer :: zone_idx, integration_count, neighborhood_idx, remainder, &
           cz_zone_index, grid_count, iendj_unused, first_match_flag, &
           cz_segment_count
      integer :: atm_get_unused_flag, katm, kenv, ksaha, ixx_flag
      integer :: klo, khi
      double precision :: cz_radius_cm(1)
! boole's output argument is a length-1 array (see numerics/boole.f90);
! these hold its result before copying into the plain scalars
! acoustic_crossing_time_seconds/acoustic_depth_cz_seconds (common/acdpth/), matching
! this file's existing cz_radius_cm/spline_interp_value convention.
      double precision :: acoustic_crossing_time_arr(1), acoustic_depth_to_cz_arr(1)
      double precision :: pressure_rotation_factor, temperature_rotation_factor
      logical :: print_flag, surface_bc_flag, pulsation_output_flag
      double precision :: spot_adjusted_teff

! (The per-zone eos scratch locals that used to live here -- the
! star_*_local copies, the eqstat2 output block, and the Cox & Giuli
! chi/cv intermediates -- moved into eos_lib.f90's eos_get_gamma1
! along with the dispatch itself; 2026, ROADMAP.md stage 1.)

! Initialize values:
      star%acoustic_depth_cz_fraction = 0.0d0
      star%acoustic_depth_cz_seconds = 0.0d0
      star%acoustic_crossing_time_seconds = 0.0d0
      star%heii_zone_acoustic_width = 0.0d0
      star%acoustic_depth_heii = 0.0d0
      atmosphere_hydrogen_fraction = 0.0d0
      star%atmosphere_sound_travel_time = 0.0d0
!      KTSAV = 3

! Need to calculate the sound speed in the envelope: first, stitch
! envelope and interior together (so that T, P, R all agree at fitting
! point. Then calculate sound speed in the envelope and interior
! and integrate.

! First, do the envelope integration to make sure the values in ENVSTRUCT
! are up-to-date (this section of code grabbed from wrtmod). First, integrate the atm and store values with grey atm
      atm_step_begin_saved = star%job%atm_step_begin
      atm_step_min_saved = star%job%atm_step_min
      atm_step_max_saved = star%job%atm_step_max
      env_step_begin_saved = star%job%env_step_begin
      env_step_min_saved = star%job%env_step_min
      env_step_max_saved = star%job%env_step_max
      star%job%atm_step_begin = star%ctrl%atm_step_size
      star%job%atm_step_min = star%ctrl%atm_step_size
      star%job%atm_step_max = star%ctrl%atm_step_size
      star%job%env_step_begin = star%ctrl%envelope_step_size
      star%job%env_step_min = star%ctrl%envelope_step_size
      star%job%env_step_max = star%ctrl%envelope_step_size
      atm_get_unused_flag = 0
      luminosity_linear = dexp(ln10*log_luminosity_lsun)
      pressure_rotation_factor = shape_factor_fp(num_shells)
      temperature_rotation_factor = shape_factor_ft(num_shells)
      katm = 0
      kenv = 0
      ksaha = 0
      pulsation_output_flag = .true.
      ixx_flag = 0
      print_flag = .true.
      surface_bc_flag = .false.
      hydrogen_fraction = composition(1,num_shells)
      metal_fraction = composition(3,num_shells)
      log10_radius_local = 0.5d0*(log_luminosity_lsun+star%log10_solar_luminosity- &
           4.0d0*log_teff-c4pil-csigl)
      log10_gravity = cgl+star%stotal-log10_radius_local-log10_radius_local
      log10_pressure_limit = log_pressure(num_shells)
!      IF (KTTAU .EQ. 0) LAOLY = .TRUE.
!      IF (KTTAU .EQ. 3) LAOLY = .FALSE.      ! for grey atm intergration: stores values in common block
!      G Somers 10/14, FOR SPOTTED RUNS, FIND THE PRESSURE AT
!      THE AMBIENT TEMPERATURE ATEFFL
      if(envelope_cz_bottom_index.eq.num_shells.and.star%ctrl%spot_filling_factor.ne.0.0.and. &
           star%ctrl%spot_temp_contrast.ne.1.0)then
            spot_adjusted_teff = log_teff-0.25*log10(star%ctrl%spot_filling_factor* &
                 star%ctrl%spot_temp_contrast**4.0+1.0-star%ctrl%spot_filling_factor)
      else
            spot_adjusted_teff = log_teff
      endif
      call atm_get(luminosity_linear,pressure_rotation_factor, &
           temperature_rotation_factor,log10_gravity,log_total_mass,ixx_flag, &
           print_flag,surface_bc_flag, &
           log10_pressure_limit,log10_radius_local,spot_adjusted_teff, &
           hydrogen_fraction,metal_fraction,atm_get_dummy1,atm_get_unused_flag, &
           katm,kenv,ksaha, &
           atm_get_dummy2,atm_get_dummy3,atm_get_dummy4,pulsation_output_flag)
!      G Somers END



! Stitch interior and envelope together, convert into ESAC06 units
      do zone_idx=1,num_shells
            star_radius_cm(zone_idx)=dexp(ln10*log_radius(zone_idx))            ! get unlogged units
            star_temperature_1e6k(zone_idx)=dexp(ln10*log_temperature(zone_idx))/1.0d6      ! ESAC06 takes T in units of 10^6K
            star_density_cgs(zone_idx)=dexp(ln10*log_density(zone_idx))
            star_pressure_cgs(zone_idx)=dexp(ln10*log_pressure(zone_idx))
            star_hydrogen_fraction(zone_idx)=composition(1,zone_idx)
            star_metal_fraction(zone_idx)=composition(3,num_shells)
      end do
      do zone_idx=1,env_struct%num_env_points-1
            star_radius_cm(num_shells+zone_idx)=dexp(ln10*env_struct%env_log10_radius(zone_idx+1))
            star_temperature_1e6k(num_shells+zone_idx)=dexp(ln10*env_struct%env_log10_temperature(zone_idx+1))/1.0d6
            star_density_cgs(num_shells+zone_idx)=dexp(ln10*env_struct%env_log10_density(zone_idx+1))
            star_pressure_cgs(num_shells+zone_idx)=dexp(ln10*env_struct%env_log10_pressure(zone_idx+1))
            star_hydrogen_fraction(num_shells+zone_idx)=composition(1,num_shells)
            star_metal_fraction(num_shells+zone_idx)=composition(3,num_shells)
      end do


! Call EOS interpolator via eos_lib's eos_get_gamma1, which performs
! the OPAL2006-vs-Yale/SCV dispatch (and the Cox & Giuli gamma1
! construction on the Yale/SCV path) that historically lived here as
! direct esac06/eqstat2 calls -- see eos_lib.f90 (2026, phase three,
! ROADMAP.md stage 1).
      do zone_idx=1,num_shells+env_struct%num_env_points-1
            call eos_get_gamma1(star_hydrogen_fraction(zone_idx), &
                 star_metal_fraction(zone_idx), &
                 star_temperature_1e6k(zone_idx), &
                 star_density_cgs(zone_idx),star_pressure_cgs(zone_idx), &
                 local_gamma1(zone_idx),star%eos_adiabatic_gradient(zone_idx), &
                 ksaha)
            star_inverse_sound_speed(zone_idx)=1.0d0/ &
                 sqrt(local_gamma1(zone_idx)*star_pressure_cgs(zone_idx)/star_density_cgs(zone_idx))
      end do




! Integrate for the full acoustic depth for normalization
      integration_count=env_struct%num_env_points-1+num_shells-1
      remainder=mod((env_struct%num_env_points-1+num_shells-1),4)
      if (remainder.ne.0) then
            grid_count=integration_count+4-remainder+1
      else
            grid_count=integration_count+1
      endif
      if (integration_count.le.3) then
            integration_count=4
            grid_count=1+1
      endif

      call boole(star_radius_cm,star_inverse_sound_speed,integration_count, &
           grid_count,acoustic_crossing_time_arr)
      star%acoustic_crossing_time_seconds = acoustic_crossing_time_arr(1)




! Find the location of the CZ, but only if one exists
      if (envelope_cz_bottom_index.lt.num_shells) then
            cz_radius_cm=dexp(ln10*(envelope_cz_log_radius+star%log10_solar_radius))


!      get the location of the convection zone, and count nonzero entries:
                  iendj_unused=0
                  first_match_flag=0
                  do zone_idx=1,integration_count
                        if (star_radius_cm(zone_idx).le.cz_radius_cm(1) .and. &
                             star_radius_cm(zone_idx+1).ge.cz_radius_cm(1)) then
                              cz_zone_index = zone_idx
                        endif

!                  Interpolate to find Cs at the exact CZ radius
                        if (zone_idx.eq.cz_zone_index .and. first_match_flag.eq.0 .and. &
                             cz_zone_index.gt.3) then
                              first_match_flag=1
                              do neighborhood_idx=1,7
                                    spline_radius_neighborhood(neighborhood_idx)= &
                                         star_radius_cm(cz_zone_index-4+neighborhood_idx)
                                    spline_speed_neighborhood(neighborhood_idx)= &
                                         star_inverse_sound_speed(cz_zone_index-4+neighborhood_idx)
                              end do
                              call splinj(spline_radius_neighborhood, &
                                   spline_speed_neighborhood,spline_second_deriv,7)
                              call splint(spline_radius_neighborhood, &
                                   spline_speed_neighborhood,7,spline_second_deriv, &
                                   cz_radius_cm(1),spline_interp_value(1),klo,khi)
                              star_radius_cm(cz_zone_index)=cz_radius_cm(1)
                              star_inverse_sound_speed(cz_zone_index)=spline_interp_value(1)
                        endif
                        if (zone_idx.eq.cz_zone_index .and. first_match_flag.eq.0 .and. &
                             cz_zone_index.le.3) then
                              first_match_flag=1
                              if (cz_zone_index .eq. 1) then
                                    cz_radius_cm(1) = star_radius_cm(cz_zone_index)
                                    spline_interp_value(1) = star_inverse_sound_speed(cz_zone_index)
                              else
                                    do neighborhood_idx=1,7
                                          spline_radius_neighborhood(neighborhood_idx)= &
                                               star_radius_cm(cz_zone_index-2+neighborhood_idx)
                                          spline_speed_neighborhood(neighborhood_idx)= &
                                               star_inverse_sound_speed(cz_zone_index-2+neighborhood_idx)
                                    end do
                                    call splinj(spline_radius_neighborhood, &
                                         spline_speed_neighborhood,spline_second_deriv,7)
                                    call splint(spline_radius_neighborhood, &
                                         spline_speed_neighborhood,7,spline_second_deriv, &
                                         cz_radius_cm(1),spline_interp_value(1),klo,khi)
                                    star_radius_cm(cz_zone_index)=cz_radius_cm(1)
                                    star_inverse_sound_speed(cz_zone_index)=spline_interp_value(1)
                              endif
                        endif

                  end do
                  do zone_idx=1,integration_count-cz_zone_index+1
                        star_radius_to_cz(zone_idx)=star_radius_cm(cz_zone_index-1+zone_idx)
                        star_inverse_sound_speed_to_cz(zone_idx)= &
                             star_inverse_sound_speed(cz_zone_index-1+zone_idx)
                  end do
! Then call Boole from surface to cz:
            remainder=mod((integration_count-cz_zone_index+1),4)
            cz_segment_count=integration_count-cz_zone_index+1
            if (remainder.ne.0) then
                  grid_count=integration_count-cz_zone_index+1+4-remainder+1
            else
                  grid_count=integration_count-cz_zone_index+1+1
            endif
            call boole(star_radius_to_cz,star_inverse_sound_speed_to_cz, &
                 cz_segment_count,grid_count,acoustic_depth_to_cz_arr)
            star%acoustic_depth_cz_seconds = acoustic_depth_to_cz_arr(1)
      else
            star%acoustic_depth_cz_seconds = 0.0d0
      endif

! Deal with the atmosphere seperately
!      Invert vectors and shift units
!      Adopt the atmosphere metallicity to be the same as the last envelope pt.
!      ATMOSX=STARX(M+NUMENV-1)
! Determine length of Atmosphere vector:
!      IATCNT=0
!      DO 56, I=1,JSON-1


! Put the atmosphere in actual units, not differential
!      DO 39, I=1,IATCNT-1
!            ACATMR(IATCNT-I)=ACATMR(IATCNT-I+1)+ACATMR(IATCNT-I)
!39      CONTINUE
!      DO 40, I=1,IATCNT
!            ACATMR(I)=ACATMR(I)+STARR(M+NUMENV-1)
!40      CONTINUE

!      DO 37,I=1,IATCNT

!      Pass atm vectors to ESAC06
!      IORDER=9
!      IRAD=1
!      DO 38, I=1,IATCNT


!       Then call Boole to do the integration:
!
!      CALL BOOLE(ATMOSR,ATMOSC,V,NN,TATMOS)

! Output normalized acoustic depth
      star%acoustic_depth_cz_fraction=0.0d0
      star%acoustic_depth_cz_fraction=star%acoustic_depth_cz_seconds/star%acoustic_crossing_time_seconds

! Output acoustic depth info to ISHORT
      write(short_file_unit,67)star%acoustic_depth_cz_seconds,star%acoustic_crossing_time_seconds, &
           star%acoustic_depth_cz_fraction,star%atmosphere_sound_travel_time
67       format(1X,'Acoustic depth to CZ:',F14.8,2X,'Acoustic depth to center', &
     F13.7,2X,'Normalized taucz:',F11.9, 'Acoustic depth of atmopshere:',F16.8,2X)

! 555            CONTINUE
!--------------------------------------------------------------
! Save all vectors of interest when the end of a kind card is reached.
      if(star%job%calcad_ageout_output_active)then

            write(unit=star%ctrl%calcad_file_unit,fmt=1506) env_struct%num_env_points,num_shells,cz_zone_index
1506            format(1X, 'Number of points in envelope:',I5,2X, &
     'Number of points in interior:',I5,2X,'Index near Rcz:' &
     ,I5,2X)

                  do zone_idx=1,num_shells+env_struct%num_env_points-1
                        if (zone_idx .le. integration_count-cz_zone_index+1) then
!                         WRITE(UNIT=ICLCD,FMT=1504),DAGE, STARR(I), STARC(I),
                        write(unit=star%ctrl%calcad_file_unit,fmt=1504) age_gyr, &
                             star_radius_cm(zone_idx),star_inverse_sound_speed(zone_idx), &
     star_radius_to_cz(zone_idx), star_inverse_sound_speed_to_cz(zone_idx), &
     star%eos_adiabatic_gradient(zone_idx), local_gamma1(zone_idx), star_pressure_cgs(zone_idx), &
     star_temperature_1e6k(zone_idx),star_density_cgs(zone_idx),star_hydrogen_fraction(zone_idx)

                        else
!                         WRITE(UNIT=ICLCD,FMT=1504) DAGE, STARR(I), STARC(I),
                        write(unit=star%ctrl%calcad_file_unit,fmt=1504) age_gyr, &
                             star_radius_cm(zone_idx),star_inverse_sound_speed(zone_idx), &
     0.0d0, 0.0d0, star%eos_adiabatic_gradient(zone_idx), local_gamma1(zone_idx), star_pressure_cgs(zone_idx), &
     star_temperature_1e6k(zone_idx),star_density_cgs(zone_idx),star_hydrogen_fraction(zone_idx)

                        endif


                  end do
1504                  format(1X,11E16.8)
!                  DO 1520 I=1,IATCNT
!                        WRITE(UNIT=IACAT,FMT=1521),DAGE, ATMOSR(I), ATMOSC(I),
!    *                     DELADJ(I), ATGAM1(I), ATMOSP(I), ATMOST(I), ATMOSX
!
!1520                  CONTINUE
!1521                  FORMAT(1X,8E16.8)


      endif



!      DO 444, I=1,JSON






! 50      CONTINUE
      return

end subroutine calcad

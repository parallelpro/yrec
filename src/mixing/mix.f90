!----------------------------------------------------------------------
! mix
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mix.f; only variable names, source form, and comment style
! were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! The energy generation routine (engeb) is used for calculations of
! epsilon and its derivatives in the Henyey portion of the code
! (coefft) prior to calling mix the first time, store original
! variables in vector hcompp. rezone hcompp along with hcomp in
! hpoint. mix is called prior to rezoning and on the third level of
! iteration in the Henyey portion.
!
! OUTLINE OF MIX :
!
! 1) call convec to determine location of cz edges and edges of mixed
!    regions with overshoot.
! 2) determine run of rates for given run of comp, rho, t.
!    if desired, the sr spits out solar neutrino fluxes to i/o unit isnu.
! 3) use kemcom to implicitly solve for abundances in radiative regions.
! 4) solve for burning in each convection zone.
! 5) compute semi-convection if desired.
! 6) compute microscopic diffusion if desired.
!
! The rate subroutine called by mix is a stripped-down version of engeb
! (called nrate) which has no derivatives or energy yields. it does
! contain neutrino fluxes for solar neutrino calculations.
subroutine mix(timestep, iteration_level, timestep_years, core_cz_edge, &
     envelope_cz_edge, mixed_zone_bounds_no_overshoot, ierr)
      use star_info_lib, only: star, i_be9, i_c12, i_c13, i_h1, i_h2, i_he3, i_he4, i_li6, i_li7, i_metals, i_n14, i_n15, i_o16, i_o17, i_o18, json
      use luout_lib
      use const_lib
      use net_lib
      use burn_lib

      implicit none

      double precision, intent(in) :: timestep
      integer, intent(in) :: iteration_level
      double precision, intent(out) :: timestep_years
      integer, intent(out) :: core_cz_edge, envelope_cz_edge
      integer, intent(inout) :: mixed_zone_bounds_no_overshoot(12,2)


















! rate_pp..frac_be7_electron: per-zone reaction rates/branching
! fractions (originally HR1-HR13,HF1,HF2). Naming and ordering match
! the rates.f90 dummy-argument list (rate_pp,...,rate_zero13,
! frac_c12_alpha,frac_be7_electron) that these arrays are filled from
! below.
      double precision :: rate_pp(json), rate_he3_he3(json), &
           rate_he3_he4(json), rate_c12_p(json), rate_c13_p(json), &
           rate_n14_p(json), rate_o16_p(json), rate_c13_alpha(json), &
           rate_zero9(json), rate_c12_alpha(json), rate_n14_alpha(json), &
           rate_triple_alpha(json), rate_zero13(json), &
           frac_c12_alpha(json), frac_be7_electron(json)
      double precision :: species_sum(15)
      integer :: radiative_zone_bounds(13,2)
      logical :: deep_mix_flag(json)
      integer :: num_species
      double precision :: max_mixed_mass
      integer :: zone_idx, copy_idx, clear_idx
      integer :: num_radiative_zones, num_mixed_zones, &
           num_mixed_zones_no_overshoot
      double precision :: deuterium_test
      double precision :: log_density_zone, log_temperature_zone, &
           hydrogen_fraction, helium_fraction, metal_fraction, &
           he3_fraction, c12_fraction, c13_fraction, n14_fraction, &
           n15_fraction, o16_fraction, o17_fraction, o18_fraction, &
           deuterium_fraction, li6_fraction, li7_fraction, be9_fraction
      integer :: radiative_region_idx, inner_zone_idx, mixed_zone_idx
      integer :: zone_begin, zone_end
      double precision :: dt_gyr, dc_dt, do_dt, dx_dt, dy_dt
      integer :: mix_start, mix_end, shell_idx, species_idx
      double precision :: weight_sum
      double precision :: total_mass_unlogged, dlnp_dr(json)
      integer :: search_idx
      double precision :: mass_fraction_above, settling_timescale, &
           max_settling_dt, settling_dt
      integer :: num_settling_steps, substep_idx

! NSPEC IS THE NUMBER OF SPECIES BEING TRACKED.
      integer, intent(out) :: ierr

      ierr = 0

      if (use_extended_composition) then
         num_species = 15
      else
         num_species = 11
      end if
!  DDAGE IS THE TIMESTEP IN YEARS.
      timestep_years = timestep/seconds_per_year
!
! DETERMINE THE LOCATION OF CONVECTION ZONES WITH AND WITHOUT OVERSHOOT.
!
! ******* INTRODUCE DPENV PARAMETER FOR DEEP MIXING *******************
!  DPENV MIXES THE STAR FROM THE SURFACE TO MASS FRACTION DPENV
!  E.G. DPENV = 0.7 MEANS THE OUTER .3 OF THE STAR IS MIXED
      if (star%ctrl%dpenv.lt.1.0d0 .and. iteration_level.gt.1) then
! MIX FROM CENTER TO A FIXED MASS FRACTION
         max_mixed_mass = star%ctrl%dpenv*exp(ln10*star%log_total_mass)
         do zone_idx = 1, star%nz
            if (star%m(zone_idx).gt.max_mixed_mass) then
               do copy_idx = zone_idx, star%nz
                  deep_mix_flag(copy_idx) = star%convective_flag(copy_idx)
               end do
               exit
            else
               deep_mix_flag(zone_idx) = .true.
            end if
         end do
! MHP 1/95 CHANGE - DPENV MIXES TO A FIXED FRACTION OF THE MAXIMUM
! LUMINOSITY.
      else
         do copy_idx = 1, star%nz
            deep_mix_flag(copy_idx) = star%convective_flag(copy_idx)
         end do
      end if

      call convec(star%xa, star%logRho, star%logP, star%logR, &
           star%log_mass, star%logT, deep_mix_flag, star%nz, &
           radiative_zone_bounds, star%mixed_zone_bounds, &
           mixed_zone_bounds_no_overshoot, core_cz_edge, envelope_cz_edge, &
           num_radiative_zones, num_mixed_zones, &
           num_mixed_zones_no_overshoot)

! FIND BURNING RATES (HR1- HR13,HF1,HF2).
      if (use_mass_accretion .and. star%ctrl%mass_accretion_rate.gt.0.0d0) then
         deuterium_test = max(star%xa(i_h2,star%nz), &
              star%ctrl%accreted_composition(12))
      else
         deuterium_test = star%xa(i_h2,star%nz)
      end if
      do zone_idx = 1, star%nz
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
         if (star%logT(zone_idx).le.star%ctrl%tcut(1)) exit
! SCALAR VARIABLES ARE USED IN THE CALLS TO THE ENERGY GENERATION ROUTINES.
! SET SCALARS EQUAL TO THE GLOBAL ARRAYS FOR THE VARIABLES OF INTEREST.
! DL-LOG(DENSITY),TL-LOG TEMPERATURE,X***-MASS FRACTION OF SPECIES ***,
! WITH HYDROGEN,HELIUM,AND METALS DENOTED AS USUAL BY X,Y,Z.
         log_density_zone = star%logRho(zone_idx)
         log_temperature_zone = star%logT(zone_idx)
         hydrogen_fraction = star%xa(i_h1,zone_idx)
         helium_fraction = star%xa(i_he4,zone_idx)
         metal_fraction = star%xa(i_metals,zone_idx)
         he3_fraction = star%xa(i_he3,zone_idx)
         c12_fraction = star%xa(i_c12,zone_idx)
         c13_fraction = star%xa(i_c13,zone_idx)
         n14_fraction = star%xa(i_n14,zone_idx)
         n15_fraction = star%xa(i_n15,zone_idx)
         o16_fraction = star%xa(i_o16,zone_idx)
         o17_fraction = star%xa(i_o17,zone_idx)
         o18_fraction = star%xa(i_o18,zone_idx)
         deuterium_fraction = star%xa(i_h2,zone_idx)
         li6_fraction = star%xa(i_li6,zone_idx)
         li7_fraction = star%xa(i_li7,zone_idx)
         be9_fraction = star%xa(i_be9,zone_idx)
! SETUP NUCLEAR ENERGY TERMS
         call rates(log_density_zone, log_temperature_zone, &
              hydrogen_fraction, helium_fraction, he3_fraction, &
              c12_fraction, c13_fraction, n14_fraction, o16_fraction, &
              o18_fraction, zone_idx, rate_pp, rate_he3_he3, rate_he3_he4, &
              rate_c12_p, rate_c13_p, rate_n14_p, rate_o16_p, &
              rate_c13_alpha, rate_zero9, rate_c12_alpha, rate_n14_alpha, &
              rate_triple_alpha, rate_zero13, frac_c12_alpha, &
              frac_be7_electron)
! MHP 5/02 COMPUTE RATE FOR DEUTERIUM BURNING
         if (deuterium_test.gt.1.0d-11) then
            call deutrate(log_density_zone, log_temperature_zone, &
                 hydrogen_fraction, zone_idx, iteration_level)
         else
            star%light_burn%deuterium_burning_rate(zone_idx) = 0.0d0
         end if
      end do
      if (zone_idx > (star%nz)) then
      zone_idx = star%nz + 1
      end if
      do clear_idx = zone_idx, star%nz
         rate_pp(clear_idx) = 0.0d0
         rate_he3_he3(clear_idx) = 0.0d0
         rate_he3_he4(clear_idx) = 0.0d0
         rate_c12_p(clear_idx) = 0.0d0
         rate_c13_p(clear_idx) = 0.0d0
         rate_n14_p(clear_idx) = 0.0d0
         rate_o16_p(clear_idx) = 0.0d0
         rate_c13_alpha(clear_idx) = 0.0d0
         rate_zero9(clear_idx) = 0.0d0
         rate_c12_alpha(clear_idx) = 0.0d0
         rate_n14_alpha(clear_idx) = 0.0d0
         rate_triple_alpha(clear_idx) = 0.0d0
         rate_zero13(clear_idx) = 0.0d0
         frac_c12_alpha(clear_idx) = 0.0d0
         frac_be7_electron(clear_idx) = 0.0d0
! MHP 5/02 ZERO OUT DEUTERIUM BURNING RATE
         star%light_burn%deuterium_burning_rate(clear_idx) = 0.0d0
      end do
!
!  NOW IMPLICITLY SOLVE FOR THE NEW ABUNDANCES AT THE END OF THE
!  TIMESTEP.  THIS IS DONE SHELL BY SHELL FOR RADIATIVE REGIONS,
!  AND FOR EACH CONVECTION ZONE AS A UNIT.
!
! RADIATIVE ZONES.
!
      do radiative_region_idx = 1, num_radiative_zones
         do inner_zone_idx = radiative_zone_bounds(radiative_region_idx,1), &
              radiative_zone_bounds(radiative_region_idx,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
            if (star%logT(inner_zone_idx).le.star%ctrl%tcut(1)) exit
            zone_begin = inner_zone_idx
            zone_end = inner_zone_idx
            call kemcom(star%logT, zone_begin, zone_end, rate_pp, &
                 rate_he3_he3, rate_he3_he4, rate_c12_p, rate_c13_p, &
                 rate_n14_p, rate_o16_p, rate_c13_alpha, rate_c12_alpha, &
                 rate_n14_alpha, rate_triple_alpha, frac_c12_alpha, &
                 star%dm, star%xa, timestep_years, ierr)
            if (ierr /= 0) return
         end do
      end do
      if (radiative_region_idx > num_radiative_zones) then
      end if
!
! CONVECTION ZONES.
! NOTE KEMCOM ALSO AUTOMATICALLY HOMOGENIZE CONVECTION ZONES.
!
      do mixed_zone_idx = 1, num_mixed_zones
         zone_begin = star%mixed_zone_bounds(mixed_zone_idx,1)
         zone_end = star%mixed_zone_bounds(mixed_zone_idx,2)
         call kemcom(star%logT, zone_begin, zone_end, rate_pp, &
              rate_he3_he3, rate_he3_he4, rate_c12_p, rate_c13_p, &
              rate_n14_p, rate_o16_p, rate_c13_alpha, rate_c12_alpha, &
              rate_n14_alpha, rate_triple_alpha, frac_c12_alpha, &
              star%dm, star%xa, timestep_years, ierr)
         if (ierr /= 0) return
      end do
      do zone_idx = 1, star%nz
         star%rot%reaction_rate_by_zone(1,zone_idx) = rate_pp(zone_idx)
         star%rot%reaction_rate_by_zone(2,zone_idx) = rate_he3_he3(zone_idx)
         star%rot%reaction_rate_by_zone(3,zone_idx) = rate_he3_he4(zone_idx)
         star%rot%reaction_rate_by_zone(4,zone_idx) = rate_c12_p(zone_idx)
         star%rot%reaction_rate_by_zone(5,zone_idx) = rate_c13_p(zone_idx)
         star%rot%reaction_rate_by_zone(6,zone_idx) = rate_n14_p(zone_idx)
         star%rot%reaction_rate_by_zone(7,zone_idx) = rate_o16_p(zone_idx)
         star%rot%reaction_rate_by_zone(8,zone_idx) = rate_c13_alpha(zone_idx)
         star%rot%reaction_rate_by_zone(9,zone_idx) = rate_zero9(zone_idx)
         star%rot%reaction_rate_by_zone(10,zone_idx) = rate_c12_alpha(zone_idx)
         star%rot%reaction_rate_by_zone(11,zone_idx) = rate_n14_alpha(zone_idx)
         star%rot%reaction_rate_by_zone(12,zone_idx) = rate_triple_alpha(zone_idx)
         star%rot%reaction_rate_by_zone(13,zone_idx) = rate_zero13(zone_idx)
         star%rot%reaction_rate_by_zone(14,zone_idx) = frac_c12_alpha(zone_idx)
         star%rot%reaction_rate_by_zone(15,zone_idx) = frac_be7_electron(zone_idx)
      end do
!
! IF ITLVL=1 THEN THE RATES OF HYDROGEN AND HELIUM BURNING ARE
!   COMPUTED EXPLICITLY ASSUMING EQUILIBRIUM HE3 AND CN.
!   IF THIS IS THE CASE, CALL OLD ENERGY GENERATION ROUTINE AND
!   OVERWRITE THE HYDROGEN ABUNDANCES (H BURNING) AND Z,C12,O16
!   ABUNDANCES (HE BURNING).
      if (iteration_level.eq.1) then
!
!  EXPLICITLY SOLVE FOR THE NEW HYDROGEN AND HELIUM BURNING RATES
!  AT THE BEGINNING OF THE TIMESTEP, ASSUMING EQUILIBRIUM HE3 AND CN
!  CYCLE ABUNDANCES.  THIS IS DONE SHELL BY SHELL FOR RADIATIVE REGIONS,
!  AND FOR EACH CONVECTION ZONE AS A UNIT.
!
! RADIATIVE ZONES.
!
         dt_gyr = timestep_years*1.0d-9
         do radiative_region_idx = 1, num_radiative_zones
            do inner_zone_idx = &
                 radiative_zone_bounds(radiative_region_idx,1), &
                 radiative_zone_bounds(radiative_region_idx,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
               if (star%logT(inner_zone_idx).le.star%ctrl%tcut(1)) exit
               zone_begin = inner_zone_idx
               zone_end = inner_zone_idx
               call eqburn(rate_pp, rate_he3_he3, rate_he3_he4, &
                    rate_c12_p, rate_c13_p, rate_n14_p, rate_o16_p, &
                    rate_c12_alpha, rate_triple_alpha, star%dm, &
                    star%logT, zone_begin, zone_end, dc_dt, do_dt, &
                    dx_dt, dy_dt, c12_fraction, o16_fraction, &
                    hydrogen_fraction, metal_fraction)
!  USE THE EXPLICIT HYDROGEN BURNING RATE.
               if (dx_dt.ne.0.0d0) then
                  star%xa(i_h1,inner_zone_idx) = hydrogen_fraction + &
                       dx_dt*dt_gyr
               end if
!  USE THE HELIUM BURNING RATE FROM EQBURN AND THE CARBON,ALPHA
!  BURNING RATE
               if (dy_dt.ne.0.0d0) then
                  star%xa(i_metals,inner_zone_idx) = metal_fraction - &
                       dy_dt*dt_gyr
                  star%xa(i_c12,inner_zone_idx) = c12_fraction + &
                       dc_dt*dt_gyr
                  star%xa(i_o16,inner_zone_idx) = o16_fraction + &
                       do_dt*dt_gyr
               end if
            end do
         end do
         if (radiative_region_idx > num_radiative_zones) then
         end if
!
! CONVECTION ZONES.
! NOTE KEMCOM ALSO AUTOMATICALLY HOMOGENIZE CONVECTION ZONES.
!
         do mixed_zone_idx = 1, num_mixed_zones
            zone_begin = star%mixed_zone_bounds(mixed_zone_idx,1)
            zone_end = star%mixed_zone_bounds(mixed_zone_idx,2)
            call eqburn(rate_pp, rate_he3_he3, rate_he3_he4, rate_c12_p, &
                 rate_c13_p, rate_n14_p, rate_o16_p, rate_c12_alpha, &
                 rate_triple_alpha, star%dm, star%logT, &
                 zone_begin, zone_end, dc_dt, do_dt, dx_dt, dy_dt, &
                 c12_fraction, o16_fraction, hydrogen_fraction, &
                 metal_fraction)
!  USE THE EXPLICIT HYDROGEN BURNING RATE.
            if (dx_dt.ne.0.0d0) then
               hydrogen_fraction = hydrogen_fraction + dx_dt*dt_gyr
               do zone_idx = zone_begin, zone_end
                  star%xa(i_h1,zone_idx) = hydrogen_fraction
               end do
            end if
!  USE THE HELIUM BURNING RATE FROM EQBURN AND THE CARBON,ALPHA
!  BURNING RATE
            if (dy_dt.ne.0.0d0) then
               metal_fraction = metal_fraction - dy_dt*dt_gyr
               c12_fraction = c12_fraction + dc_dt*dt_gyr
               o16_fraction = o16_fraction + do_dt*dt_gyr
               do zone_idx = zone_begin, zone_end
                  star%xa(i_metals,zone_idx) = metal_fraction
                  star%xa(i_c12,zone_idx) = c12_fraction
                  star%xa(i_o16,zone_idx) = o16_fraction
               end do
            end if
         end do
      end if
!
! DETERMINE EXTENT OF SEMI-CONVECTION IF APPLICABLE.
!
      if (lsemic) then
         if (iteration_level.gt.1) &
              call sconvec(timestep, star%xa, star%logRho, &
              star%luminosity_lsun, star%logP, star%logR, star%log_mass, &
              star%logT, star%nz, star%mixed_zone_bounds, &
              num_mixed_zones, star%log_Teff, ierr)
              if (ierr /= 0) return
      end if
      if (lsemic .or. (iteration_level .eq. 1)) then
!
!    MIX CONVECTIVE REGIONS IN ORDER.
!    THIS NEEDS TO BE DONE IF SEMI-CONVECTION IS BEING CHECKED, OR
!    IF EXPLICIT H AND HE BURNING IS BEING DONE.
!
! NZONE IS THE NUMBER OF DISTINCT CONVECTION ZONES.
! MXZONE(J,1) AND MXZONE(J,2) ARE THE FIRST AND LAST SHELLS CONTAINED
! IN THE JTH CONVECTION ZONE.
         do mixed_zone_idx = 1, num_mixed_zones
! I1 AND I2 ARE THE FIRST AND LAST CONVECTIVE SHELLS IN THE GIVEN REGION.
            mix_start = star%mixed_zone_bounds(mixed_zone_idx,1)
            mix_end = star%mixed_zone_bounds(mixed_zone_idx,2)
            if (mix_start.ne.1 .and. mix_start.ge.mix_end) cycle
! INITIALIZE SUMS.
            weight_sum = 0.0d0
            do species_idx = 1, num_species
               species_sum(species_idx) = 0.0d0
            end do
!  ADD UP THE TOTAL MASS OF EACH SPECIES IN THE CONVECTIVE REGION.
!  (HS2 IS THE MASS CONTAINED WITHIN A SHELL IN GRAMS).
            do shell_idx = mix_start, mix_end
               weight_sum = weight_sum + star%dm(shell_idx)
               do species_idx = 1, num_species
                  species_sum(species_idx) = species_sum(species_idx) + &
                       star%xa(species_idx,shell_idx)*star%dm(shell_idx)
               end do
            end do
!  DIVIDE BY THE TOTAL MASS TO FIND THE MEAN MASS FRACTION IN THE REGION.
            do species_idx = 1, num_species
               species_sum(species_idx) = species_sum(species_idx)/weight_sum
            end do
!  APPLY THE MEAN MASS FRACTION OF ALL SPECIES THROUGHOUT THE CZ.
            do shell_idx = mix_start, mix_end
               do species_idx = 1, num_species
                  star%xa(species_idx,shell_idx) = species_sum(species_idx)
               end do
            end do
         end do
      end if
!  WRITE OUT THE LOCATIONS OF MIXED REGIONS.
      if (num_mixed_zones.ge.1) then
         write(short_file_unit,110) ((star%mixed_zone_bounds(zone_idx,inner_zone_idx), &
              inner_zone_idx=1,2),zone_idx=1,num_mixed_zones)
  110    format(' ZONES MIXED IN ORDER--',12('(',i5,',',i5,') ') )
      end if
!
! MICROSCOPIC DIFFUSION OF HELIUM.
!
!***MHP 12/89 MIX MODIFIED TO INCLUDE CALL TO GRAVITATIONAL SETTLING
!      ROUTINE USING THE BAHCALL AND LOEB METHOD.
! FIRST DEFINE VARIABLES NEEDED FOR SETTLING -
! HQPR=VECTOR OF D LN P/DR.
! STOT=TOTAL STELLAR MASS(UNLOGGED).
      if (diffuse_helium_active) then
      settling: do
         if (star%xa(i_h1,1).lt.star%ctrl%hydrogen_diffusion_floor) then
            diffuse_helium_active = .false.
            exit settling
         end if
         total_mass_unlogged = exp(ln10*star%log_total_mass)
         do zone_idx = 1, star%nz
            dlnp_dr(zone_idx) = -exp(ln10*(star%logRho(zone_idx)+cgl+ &
                 star%log_mass(zone_idx)-2.0d0*star%logR(zone_idx)- &
                 star%logP(zone_idx)))
         end do
! MHP 6/90 CHANGE ADDED : THE TIMESTEP FOR SETTLING IS RESTRICTED TO
!   A FRACTION OF THE TIMESCALE FOR SETTLING AT THE OUTER BOUNDARY.
!   THE OUTER BOUNDARY IS EITHER THE SURFACE CONVECTION ZONE OR THE
!   FIRST POINT WHERE THE HELIUM ABUNDANCE RISES ABOVE A USER-SPECIFIED
!   MINIMUM VALUE (YMIN).
!
!   LOCATE OUTER BOUNDARY.
         if (.not.star%convective_flag(star%nz)) then
            write(short_file_unit,911)
! DBG 2/92 CHANGED STOP TO JUST A WARNING MESSAGE, EXECUTION CONTINUES
  911       format(1x,'NO SURFACE CZ - DIFFUSION NOT MEANINGFUL')
            exit settling
         end if
         do search_idx = envelope_cz_edge, 1, -1
            if (star%xa(i_he4,search_idx).gt.star%ctrl%helium_diffusion_min) exit
         end do
         if (search_idx < (1)) then
!   Y<YMIN FOR THE WHOLE STAR IF THE CODE GETS HERE.
         exit settling
         end if
!  FM IS THE MASS FRACTION ABOVE THE OUTER POINT.
         mass_fraction_above = (total_mass_unlogged- &
              star%m(search_idx))/total_mass_unlogged
!  TSCALE IS THE TIMESCALE FOR SETTLING OF HELIUM AT THE OUTER
!  BOUNDARY (MICHAUD ET AL 1984, APJ V.282,P.206)
         settling_timescale = 4.348d21*seconds_per_year*mass_fraction_above/ &
              exp(ln10*1.5d0*star%logT(search_idx))
!  RESTRICT TIMESTEP TO THE MINIMUM OF THE MODEL TIMESTEP AND
!  A USER SPECIFIED FRACTION (DT_GS) OF THE SETTLING TIMESCALE.
         max_settling_dt = star%ctrl%settling_timestep_fraction*settling_timescale
         num_settling_steps = int(timestep/max_settling_dt)
         if (mod(max_settling_dt,timestep).ne.0.0d0) num_settling_steps = &
              num_settling_steps + 1
         settling_dt = timestep/dfloat(num_settling_steps)
         do substep_idx = 1, num_settling_steps
! PERFORM GRAVITATIONAL SETTLING. IF LNEWDIF = TRUE, USE THE NEW ROUTINES
! IN MICRODIFF. ELSE, USE THE OLD ROUTINES IN GRSETT.
            if (star%ctrl%use_new_diffusion_routines) then
               call microdiff(settling_dt, star%xa, dlnp_dr, &
                    star%logR, star%logRho, star%m, &
                    star%logT, deep_mix_flag, star%nz, &
                    total_mass_unlogged)
            else
               call grsett(settling_dt, star%xa, dlnp_dr, star%logR, &
                    star%logRho, star%m, star%logT, &
                    deep_mix_flag, star%nz, total_mass_unlogged)
            end if
         end do
      exit settling
      end do settling
      end if
! RENORMALIZE COMPOSITION TO GUARD AGAINST ANOMALIES (I.E. SMALL NEGATIVE
! ABUNDANCES...).
      do zone_idx = 1, star%nz
         do species_idx = 1, num_species
            star%xa(species_idx,zone_idx) = &
                 max(star%xa(species_idx,zone_idx),0.0d0)
            star%xa(species_idx,zone_idx) = &
                 min(star%xa(species_idx,zone_idx),1.0d0)
         end do
         star%xa(i_metals,zone_idx) = min(star%xa(i_metals,zone_idx), &
              1.0d0-star%xa(i_h1,zone_idx)-star%xa(i_he3,zone_idx))
         star%xa(i_he4,zone_idx) = 1.0d0 - star%xa(i_h1,zone_idx) - &
              star%xa(i_metals,zone_idx) - star%xa(i_he3,zone_idx)
      end do
! MHP 1/95 ADDED CALL TO RESET JENV,JCORE FOR DEEP MIXING.
      if (star%ctrl%dpenv.lt.1.0d0 .and. iteration_level.gt.1) then
         call convec(star%xa, star%logRho, star%logP, star%logR, &
              star%log_mass, star%logT, star%convective_flag, star%nz, &
              radiative_zone_bounds, star%mixed_zone_bounds, &
              mixed_zone_bounds_no_overshoot, core_cz_edge, &
              envelope_cz_edge, num_radiative_zones, num_mixed_zones, &
              num_mixed_zones_no_overshoot)
      end if
! MHP 5/02 ADDED DEUTERIUM BURNING
      if (use_extended_composition) then
! RADIATIVE ZONES.
!
         dt_gyr = timestep_years*1.0d-9
         do radiative_region_idx = 1, num_radiative_zones
            do inner_zone_idx = &
                 radiative_zone_bounds(radiative_region_idx,1), &
                 radiative_zone_bounds(radiative_region_idx,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
               if (star%logT(inner_zone_idx).le.star%ctrl%tcut(1)) exit
               zone_begin = inner_zone_idx
               zone_end = inner_zone_idx
               call dburn(zone_begin, zone_end, star%nz, star%dm, &
                    star%xa, dt_gyr)
            end do
         end do
         if (radiative_region_idx > num_radiative_zones) then
         end if
!
! CONVECTION ZONES.
! NOTE KEMCOM ALSO AUTOMATICALLY HOMOGENIZE CONVECTION ZONES.
!
         do mixed_zone_idx = 1, num_mixed_zones
            zone_begin = star%mixed_zone_bounds(mixed_zone_idx,1)
            zone_end = star%mixed_zone_bounds(mixed_zone_idx,2)
            call dburn(zone_begin, zone_end, star%nz, star%dm, &
                 star%xa, dt_gyr)
         end do
      end if
      return
end subroutine mix

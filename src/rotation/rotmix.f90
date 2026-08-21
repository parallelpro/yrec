!----------------------------------------------------------------------
! rotmix
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rotmix.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! 1) USE KEMCOM TO IMPLICITLY SOLVE FOR ABUNDANCES IN RADIATIVE REGIONS.
! 2) SOLVE FOR BURNING IN EACH CONVECTION ZONE.
! Then (if diffuse_helium_active) performs gravitational settling of
! helium (Bahcall & Loeb method), and finally renormalizes the
! composition to guard against small negative/overflowing abundances.
!
!       SUBROUTINE ROTMIX(DELTS,HCOMP,HS2,HT,ITLVL,M,MRZONE,MXZONE,  ! KC 2025-05-31
!                  NRZONE,NZONE
!                  ,HSTOT,HD,HS,HR,HP,LC,HS1)
subroutine rotmix(timestep, composition, shell_mass, log_temperature, &
     num_zones, radiative_zone_bounds, convective_zone_bounds, &
     num_radiative_zones, num_convective_zones, log_total_mass, &
     log_density, log_mass, log_radius, log_pressure, convective_flag, &
     enclosed_mass)

      use mdphy_lib
      use scrtch_lib
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: shell_mass(json), log_temperature(json)
      integer, intent(in) :: num_zones
      integer, intent(in) :: radiative_zone_bounds(13,2), &
           convective_zone_bounds(12,2)
      integer, intent(in) :: num_radiative_zones, num_convective_zones
      double precision, intent(in) :: log_total_mass
      double precision, intent(in) :: log_density(json), log_mass(json), &
           log_radius(json), log_pressure(json)
      logical, intent(in) :: convective_flag(json)
      double precision, intent(in) :: enclosed_mass(json)

! common/burn/: reaction_rate_by_zone (originally HCOMPM) stores the
! per-zone reaction rates/branching fractions computed in mix.f90; read
! here and copied into the rate_*/frac_* locals below before calling
! KEMCOM. Naming matches mix.f90.
      double precision :: reaction_rate_by_zone(15,json)
      common/burn/ reaction_rate_by_zone












      integer :: num_species
      double precision :: timestep_years
      double precision :: rate_pp(json), rate_he3_he3(json), &
           rate_he3_he4(json), rate_c12_p(json), rate_c13_p(json), &
           rate_n14_p(json), rate_o16_p(json), rate_c13_alpha(json), &
           rate_zero9(json), rate_c12_alpha(json), rate_n14_alpha(json), &
           rate_triple_alpha(json), rate_zero13(json), &
           frac_c12_alpha(json), frac_be7_electron(json)
      double precision :: dlnp_dr_settling(json), del_grad2_save(json)
      logical :: am_transport_convective_flag(json)
      double precision :: total_mass
      save

      integer :: zone_idx, region_idx, species_idx, burn_zone_start, &
           burn_zone_end, outer_boundary_zone, num_settling_substeps, &
           substep_idx
      double precision :: mass_fraction_above, settling_timescale, &
           max_settling_dt, settling_dt

! NSPEC IS THE NUMBER OF SPECIES BEING TRACKED.
      if (use_extended_composition) then
         num_species = 15
      else
         num_species = 11
      end if
!  DDAGE IS THE TIMESTEP IN YEARS.
      timestep_years = timestep/seconds_per_year
      do zone_idx = 1,num_zones
         rate_pp(zone_idx) = reaction_rate_by_zone(1,zone_idx)
         rate_he3_he3(zone_idx) = reaction_rate_by_zone(2,zone_idx)
         rate_he3_he4(zone_idx) = reaction_rate_by_zone(3,zone_idx)
         rate_c12_p(zone_idx) = reaction_rate_by_zone(4,zone_idx)
         rate_c13_p(zone_idx) = reaction_rate_by_zone(5,zone_idx)
         rate_n14_p(zone_idx) = reaction_rate_by_zone(6,zone_idx)
         rate_o16_p(zone_idx) = reaction_rate_by_zone(7,zone_idx)
         rate_c13_alpha(zone_idx) = reaction_rate_by_zone(8,zone_idx)
         rate_zero9(zone_idx) = reaction_rate_by_zone(9,zone_idx)
         rate_c12_alpha(zone_idx) = reaction_rate_by_zone(10,zone_idx)
         rate_n14_alpha(zone_idx) = reaction_rate_by_zone(11,zone_idx)
         rate_triple_alpha(zone_idx) = reaction_rate_by_zone(12,zone_idx)
         rate_zero13(zone_idx) = reaction_rate_by_zone(13,zone_idx)
         frac_c12_alpha(zone_idx) = reaction_rate_by_zone(14,zone_idx)
         frac_be7_electron(zone_idx) = reaction_rate_by_zone(15,zone_idx)
      end do
!
!  NOW IMPLICITLY SOLVE FOR THE NEW ABUNDANCES AT THE END OF THE
!  TIMESTEP.  THIS IS DONE SHELL BY SHELL FOR RADIATIVE REGIONS,
!  AND FOR EACH CONVECTION ZONE AS A UNIT.
!
! RADIATIVE ZONES.
!
      do 40 region_idx = 1,num_radiative_zones
         do 30 zone_idx = radiative_zone_bounds(region_idx,1), &
              radiative_zone_bounds(region_idx,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
            if (log_temperature(zone_idx).le.tcut(1)) goto 45
            burn_zone_start = zone_idx
            burn_zone_end = zone_idx
            call kemcom(log_temperature,burn_zone_start,burn_zone_end, &
                 rate_pp,rate_he3_he3,rate_he3_he4,rate_c12_p,rate_c13_p, &
                 rate_n14_p,rate_o16_p, &
!      *                   HR8,HR9,HR10,HR11,HR12,HR13,HF1,HS2,HCOMP,
!      *                   DDAGE,ITLVL)  ! KC 2025-05-31
                 rate_c13_alpha,rate_c12_alpha,rate_n14_alpha, &
                 rate_triple_alpha,frac_c12_alpha,shell_mass,composition, &
                 timestep_years)
   30    continue
   40 continue
   45 continue
!
! CONVECTION ZONES.
! NOTE KEMCOM ALSO AUTOMATICALLY HOMOGENIZE CONVECTION ZONES.
!
      do 50 region_idx = 1,num_convective_zones
         burn_zone_start = convective_zone_bounds(region_idx,1)
         burn_zone_end = convective_zone_bounds(region_idx,2)
         call kemcom(log_temperature,burn_zone_start,burn_zone_end, &
              rate_pp,rate_he3_he3,rate_he3_he4,rate_c12_p,rate_c13_p, &
              rate_n14_p,rate_o16_p, &
!      *                HR8,HR9,HR10,HR11,HR12,HR13,HF1,HS2,HCOMP,
!      *                DDAGE,ITLVL)  ! KC 2025-05-31
              rate_c13_alpha,rate_c12_alpha,rate_n14_alpha, &
              rate_triple_alpha,frac_c12_alpha,shell_mass,composition, &
              timestep_years)
   50 continue
!
! MICROSCOPIC DIFFUSION OF HELIUM.
!
!***BC 5/92 ROTMIX MODIFIED TO INCLUDE CALL TO GRAVITATIONAL SETTLING
!      ROUTINE USING THE BAHCALL AND LOEB METHOD.
! FIRST DEFINE VARIABLES NEEDED FOR SETTLING -
! HQPR=VECTOR OF D LN P/DR.
! STOT=TOTAL STELLAR MASS(UNLOGGED).
      if (diffuse_helium_active) then
         if (composition(1,1).lt.hydrogen_diffusion_floor) then
            diffuse_helium_active=.false.
            goto 170
         end if
! MHP 6/90 CHANGE ADDED : THE TIMESTEP FOR SETTLING IS RESTRICTED TO
!   A FRACTION OF THE TIMESCALE FOR SETTLING AT THE OUTER BOUNDARY.
!   THE OUTER BOUNDARY IS EITHER THE SURFACE CONVECTION ZONE OR THE
!   FIRST POINT WHERE THE HELIUM ABUNDANCE RISES ABOVE A USER-SPECIFIED
!   MINIMUM VALUE (YMIN).
!
!   LOCATE OUTER BOUNDARY.
         if (.not.convective_flag(num_zones)) then
            write(*,911)
            write(short_file_unit,911)
  911       format(1x,'NO SURFACE CZ - DIFFUSION NOT MEANINGFUL'/ &
                 'STOPPED IN SUBROUTINE MIX')
            stop
         end if
! 7/92 MHP STATEMENT ADDED FOR EXIT IF OVERSHOOT CAUSES A FULLY CONVECTIVE
! CASE.  AVOIDS OCCASIONAL PRE-MS CRASH IN FIRST RADIATIVE MODEL.
         if (convective_zone_bounds(num_convective_zones,1).le.2 .and. &
              convective_zone_bounds(num_convective_zones,2).eq.num_zones) &
              goto 170
         outer_boundary_zone = radiative_zone_bounds(num_radiative_zones,2)
         do 140 zone_idx = outer_boundary_zone,1,-1
            if (composition(2,zone_idx).gt.helium_diffusion_min) goto 150
  140    continue
!   Y<YMIN FOR THE WHOLE STAR IF THE CODE GETS HERE.
         goto 170
  150    continue
         total_mass=exp(ln10*log_total_mass)
         do 130 zone_idx = 1,num_zones
            del_grad2_save(zone_idx) = shell_diag%del_grad(2,zone_idx)
            shell_diag%del_grad(2,zone_idx) = mix_phys%delm(zone_idx)
            dlnp_dr_settling(zone_idx)=-exp(ln10*(log_density(zone_idx)+ &
                 cgl+log_mass(zone_idx)-2.0d0*log_radius(zone_idx)- &
                 log_pressure(zone_idx)))
  130    continue
! ***BC 6/92 only check for timestep cutting if JMAX is large
!
         if (outer_boundary_zone.ge.2) then
!  FM IS THE MASS FRACTION ABOVE THE OUTER POINT.
            mass_fraction_above = (total_mass- &
                 enclosed_mass(outer_boundary_zone))/total_mass
!  TSCALE IS THE TIMESCALE FOR SETTLING OF HELIUM AT THE OUTER
!  BOUNDARY (MICHAUD ET AL 1984, APJ V.282,P.206)
            settling_timescale = 4.348d21*seconds_per_year* &
                 mass_fraction_above/ &
                 exp(ln10*1.5d0*log_temperature(outer_boundary_zone))
!  RESTRICT TIMESTEP TO THE MINIMUM OF THE MODEL TIMESTEP AND
!  A USER SPECIFIED FRACTION (DT_GS) OF THE SETTLING TIMESCALE.
            write(69,*) 'JMAX=',outer_boundary_zone,' FM= ', &
                 mass_fraction_above,' TSCALE=',settling_timescale
            max_settling_dt = settling_timestep_fraction*settling_timescale
            num_settling_substeps = int(timestep/max_settling_dt)
            if (mod(max_settling_dt,timestep).ne.0.0d0.or. &
                 num_settling_substeps.eq.0) &
                 num_settling_substeps=num_settling_substeps+1
            write(69,*)'DTMAX=',max_settling_dt,' DELTS=',timestep, &
                 ' NSTEP=',num_settling_substeps
            settling_dt = timestep/dfloat(num_settling_substeps)
         else
            num_settling_substeps=1
            settling_dt = timestep
         end if
         do region_idx = 1,num_radiative_zones
            do zone_idx = radiative_zone_bounds(region_idx,1), &
                 radiative_zone_bounds(region_idx,2)
               am_transport_convective_flag(zone_idx) = .false.
            end do
         end do
         do region_idx = 1,num_convective_zones
            do zone_idx = convective_zone_bounds(region_idx,1), &
                 convective_zone_bounds(region_idx,2)
               am_transport_convective_flag(zone_idx) = .true.
            end do
         end do
         do 160 substep_idx = 1,num_settling_substeps
! PERFORM GRAVITATIONAL SETTLING. IF LNEWDIF = TRUE, USE THE NEW ROUTINES
! IN MICRODIFF. ELSE, USE THE OLD ROUTINES IN GRSETT.
            if (use_new_diffusion_routines) then
               call microdiff(settling_dt,composition,dlnp_dr_settling, &
                    log_radius,log_density,enclosed_mass,log_temperature, &
                    am_transport_convective_flag,num_zones,total_mass)
            else
               call grsett(settling_dt,composition,dlnp_dr_settling, &
                    log_radius,log_density,enclosed_mass,log_temperature, &
                    am_transport_convective_flag,num_zones,total_mass)
            end if
  160    continue
         do zone_idx = 1,num_zones
            shell_diag%del_grad(2,zone_idx) = del_grad2_save(zone_idx)
         end do
  170    continue
      end if
! END OF GRAVITATIONAL SETTLING
!
! RENORMALIZE COMPOSITION TO GUARD AGAINST ANOMALIES (I.E. SMALL NEGATIVE
! ABUNDANCES...).
      do 180 zone_idx = 1,num_zones
         do 175 species_idx = 1,num_species
            composition(species_idx,zone_idx) = &
                 max(composition(species_idx,zone_idx),0.0d0)
            composition(species_idx,zone_idx) = &
                 min(composition(species_idx,zone_idx),1.0d0)
  175    continue
         composition(3,zone_idx) = min(composition(3,zone_idx), &
              1.0d0-composition(1,zone_idx)-composition(4,zone_idx))
         composition(2,zone_idx) = 1.0d0 - composition(1,zone_idx) - &
              composition(3,zone_idx) - composition(4,zone_idx)
  180 continue
      return
end subroutine rotmix

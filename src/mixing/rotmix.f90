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
! The reaction rates are taken from rot_scr%reaction_rate_by_zone, filled
! by mix.f90 earlier in the timestep.
subroutine rotmix(timestep, composition, shell_mass, log_temperature, &
     num_zones, radiative_zone_bounds, convective_zone_bounds, &
     num_radiative_zones, num_convective_zones, log_total_mass, &
     log_density, log_mass, log_radius, log_pressure, convective_flag, &
     enclosed_mass, ierr)
      use rotation_scratch_lib

      use star_info_lib, only: star, json
      use luout_lib
      use run_log_lib, only: solver_diagnostics
      use phys_const_lib
      use math_lib
      implicit none

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

      integer :: num_species
      double precision :: timestep_years
      double precision :: rate_pp(json), rate_he3_he3(json), &
           rate_he3_he4(json), rate_c12_p(json), rate_c13_p(json), &
           rate_n14_p(json), rate_o16_p(json), rate_c13_alpha(json), &
           rate_c12_alpha(json), rate_n14_alpha(json), &
           rate_triple_alpha(json), frac_c12_alpha(json)
      double precision :: dlnp_dr_settling(json), del_grad2_save(json)
      logical :: am_transport_convective_flag(json)
      double precision :: total_mass
      integer :: zone_idx, region_idx, species_idx, burn_zone_start, &
           burn_zone_end, outer_boundary_zone, num_settling_substeps, &
           substep_idx
      double precision :: mass_fraction_above, settling_timescale, &
           max_settling_dt, settling_dt

! NSPEC IS THE NUMBER OF SPECIES BEING TRACKED.
      integer, intent(out) :: ierr

      ierr = 0

      if (star%job%use_extended_composition) then
         num_species = 15
      else
         num_species = 11
      end if
!  DDAGE IS THE TIMESTEP IN YEARS.
      timestep_years = timestep/seconds_per_year
      do zone_idx = 1,num_zones
         rate_pp(zone_idx) = rot_scr%reaction_rate_by_zone(1,zone_idx)
         rate_he3_he3(zone_idx) = rot_scr%reaction_rate_by_zone(2,zone_idx)
         rate_he3_he4(zone_idx) = rot_scr%reaction_rate_by_zone(3,zone_idx)
         rate_c12_p(zone_idx) = rot_scr%reaction_rate_by_zone(4,zone_idx)
         rate_c13_p(zone_idx) = rot_scr%reaction_rate_by_zone(5,zone_idx)
         rate_n14_p(zone_idx) = rot_scr%reaction_rate_by_zone(6,zone_idx)
         rate_o16_p(zone_idx) = rot_scr%reaction_rate_by_zone(7,zone_idx)
         rate_c13_alpha(zone_idx) = rot_scr%reaction_rate_by_zone(8,zone_idx)
         rate_c12_alpha(zone_idx) = rot_scr%reaction_rate_by_zone(10,zone_idx)
         rate_n14_alpha(zone_idx) = rot_scr%reaction_rate_by_zone(11,zone_idx)
         rate_triple_alpha(zone_idx) = rot_scr%reaction_rate_by_zone(12,zone_idx)
         frac_c12_alpha(zone_idx) = rot_scr%reaction_rate_by_zone(14,zone_idx)
      end do
!
!  NOW IMPLICITLY SOLVE FOR THE NEW ABUNDANCES AT THE END OF THE
!  TIMESTEP.  THIS IS DONE SHELL BY SHELL FOR RADIATIVE REGIONS,
!  AND FOR EACH CONVECTION ZONE AS A UNIT.
!
! RADIATIVE ZONES.
!
      do region_idx = 1,num_radiative_zones
         do zone_idx = radiative_zone_bounds(region_idx,1), &
              radiative_zone_bounds(region_idx,2)
! EXIT LOOP ONCE T DROPS BELOW NUCLEAR REACTION T CUTOFF
            if (log_temperature(zone_idx).le.star%ctrl%nuclear_logT_cutoffs(1)) exit
            burn_zone_start = zone_idx
            burn_zone_end = zone_idx
            call solve_composition(log_temperature,burn_zone_start,burn_zone_end, &
                 rate_pp,rate_he3_he3,rate_he3_he4,rate_c12_p,rate_c13_p, &
                 rate_n14_p,rate_o16_p, &
                 rate_c13_alpha,rate_c12_alpha,rate_n14_alpha, &
                 rate_triple_alpha,frac_c12_alpha,shell_mass,composition, &
                 timestep_years, ierr)
            if (ierr /= 0) return
         end do
      end do
!
! CONVECTION ZONES.
! NOTE KEMCOM ALSO AUTOMATICALLY HOMOGENIZE CONVECTION ZONES.
!
      do region_idx = 1,num_convective_zones
         burn_zone_start = convective_zone_bounds(region_idx,1)
         burn_zone_end = convective_zone_bounds(region_idx,2)
         call solve_composition(log_temperature,burn_zone_start,burn_zone_end, &
              rate_pp,rate_he3_he3,rate_he3_he4,rate_c12_p,rate_c13_p, &
              rate_n14_p,rate_o16_p, &
              rate_c13_alpha,rate_c12_alpha,rate_n14_alpha, &
              rate_triple_alpha,frac_c12_alpha,shell_mass,composition, &
              timestep_years, ierr)
         if (ierr /= 0) return
      end do
!
! MICROSCOPIC DIFFUSION OF HELIUM.
!
!***BC 5/92 ROTMIX MODIFIED TO INCLUDE CALL TO GRAVITATIONAL SETTLING
!      ROUTINE USING THE BAHCALL AND LOEB METHOD.
! FIRST DEFINE VARIABLES NEEDED FOR SETTLING -
! HQPR=VECTOR OF D LN P/DR.
! STOT=TOTAL STELLAR MASS(UNLOGGED).
      if (star%job%diffuse_helium_active) then
      settling: do
         if (composition(1,1).lt.star%ctrl%hydrogen_diffusion_floor) then
            star%job%diffuse_helium_active=.false.
            exit settling
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
            write(run_log_unit,911)
  911       format(1x,'NO SURFACE CZ - DIFFUSION NOT MEANINGFUL'/ &
                 'STOPPED IN SUBROUTINE MIX')
            ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the driver-side
            ! call sites (core/main, core/crrect, core/starin, setup/hpoint)
            ! preserve the historical stop on a nonzero return.
            ierr = 1
            return
         end if
! 7/92 MHP STATEMENT ADDED FOR EXIT IF OVERSHOOT CAUSES A FULLY CONVECTIVE
! CASE.  AVOIDS OCCASIONAL PRE-MS CRASH IN FIRST RADIATIVE MODEL.
         if (convective_zone_bounds(num_convective_zones,1).gt.2 .or. convective_zone_bounds(num_convective_zones,2).ne.num_zones) then
         outer_boundary_zone = radiative_zone_bounds(num_radiative_zones,2)
         do zone_idx = outer_boundary_zone,1,-1
            if (composition(2,zone_idx).gt.star%ctrl%helium_diffusion_min) exit
         end do
         if (zone_idx < (1)) then
!   Y<YMIN FOR THE WHOLE STAR IF THE CODE GETS HERE.
         exit settling
         end if
         total_mass=exp(ln10*log_total_mass)
         do zone_idx = 1,num_zones
            del_grad2_save(zone_idx) = star%gradT(zone_idx)
            star%gradT(zone_idx) = mix_scr%delm(zone_idx)
            dlnp_dr_settling(zone_idx)=-exp(ln10*(log_density(zone_idx)+ &
                 cgl+log_mass(zone_idx)-2.0d0*log_radius(zone_idx)- &
                 log_pressure(zone_idx)))
         end do
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
! (2026: this settling diagnostic wrote to never-opened unit 69,
! leaving stray fort.69 files -- now solver forensics in the run
! log, like its peers)
            if (solver_diagnostics()) write(run_log_unit,*) 'JMAX=', &
                 outer_boundary_zone,' FM= ', &
                 mass_fraction_above,' TSCALE=',settling_timescale
            max_settling_dt = star%ctrl%settling_timestep_fraction*settling_timescale
            num_settling_substeps = int(timestep/max_settling_dt)
            if (mod(max_settling_dt,timestep).ne.0.0d0.or. &
                 num_settling_substeps.eq.0) &
                 num_settling_substeps=num_settling_substeps+1
            if (solver_diagnostics()) write(run_log_unit,*) &
                 'DTMAX=',max_settling_dt,' DELTS=',timestep, &
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
         do substep_idx = 1,num_settling_substeps
! PERFORM GRAVITATIONAL SETTLING. IF LNEWDIF = TRUE, USE THE NEW ROUTINES
! IN MICRODIFF. ELSE, USE THE OLD ROUTINES IN GRSETT.
            if (star%ctrl%use_new_diffusion_routines) then
               call microdiff(settling_dt,composition,dlnp_dr_settling, &
                    log_radius,log_density,enclosed_mass,log_temperature, &
                    am_transport_convective_flag,num_zones,total_mass, ierr)
               if (ierr /= 0) return
            else
               call gravitational_settling(settling_dt,composition,dlnp_dr_settling, &
                    log_radius,log_density,enclosed_mass,log_temperature, &
                    am_transport_convective_flag,num_zones,total_mass, ierr)
               if (ierr /= 0) return
            end if
         end do
         do zone_idx = 1,num_zones
            star%gradT(zone_idx) = del_grad2_save(zone_idx)
         end do
         end if
      exit settling
      end do settling
      end if
! END OF GRAVITATIONAL SETTLING
!
! RENORMALIZE COMPOSITION TO GUARD AGAINST ANOMALIES (I.E. SMALL NEGATIVE
! ABUNDANCES...).
      do zone_idx = 1,num_zones
         do species_idx = 1,num_species
            composition(species_idx,zone_idx) = &
                 max(composition(species_idx,zone_idx),0.0d0)
            composition(species_idx,zone_idx) = &
                 min(composition(species_idx,zone_idx),1.0d0)
         end do
         composition(3,zone_idx) = min(composition(3,zone_idx), &
              1.0d0-composition(1,zone_idx)-composition(4,zone_idx))
         composition(2,zone_idx) = 1.0d0 - composition(1,zone_idx) - &
              composition(3,zone_idx) - composition(4,zone_idx)
      end do
      return
end subroutine rotmix

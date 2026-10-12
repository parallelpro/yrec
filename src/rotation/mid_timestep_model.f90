!----------------------------------------------------------------------
! mid_timestep_model
!----------------------------------------------------------------------
! The rotation sub-step. A diffusion sub-step is shorter than the model
! step: advance_rotation_substep interpolates the structure linearly in
! time between the step's start and end, burns (burn_and_settle, which
! also settles helium), finds the omega the changed radii imply under
! conservation of angular momentum, and locates the edge of a central
! convection zone and the base of the surface one with that zone's mass
! and moment of inertia. The midpoint structure (eta_squared_mid ...
! am_transport_convective_flag_mid) goes into rot_scr, where the rest of
! the sub-step reads it.
module mid_timestep_model
      use am_regions, only: am_convective_regions
      use omega_profile, only: solid_body_omega, omega_from_j
      implicit none
      private :: burn_and_settle
contains

subroutine advance_rotation_substep(full_timestep, sub_timestep, time_fraction, first_call, &
     moment_of_inertia_cz, cz_mass_bottom, cz_mass_top, &
     core_boundary_zone, envelope_boundary_zone, fully_convective_flag, &
     surface_cz_active, &
     radiative_zone_bounds, convective_zone_bounds, num_radiative_zones, &
     num_convective_zones, ierr)
      use rotation_scratch_lib
      use star_info_lib, only: star, i_eps_grav, i_eps_neu, i_h2, json
      use net_lib, only: burn_deuterium_region, min_deuterium_substep, &
           num_rate_rows
      use physical_constants
      use burn_step, only: surface_accretion_now, burn_regions, &
           network_rates_of_model
      use math_lib
      implicit none

      double precision, intent(in) :: full_timestep, sub_timestep
      double precision, intent(in) :: time_fraction
      logical, intent(in) :: first_call
      double precision, intent(out) :: moment_of_inertia_cz
      double precision, intent(out) :: cz_mass_bottom, cz_mass_top
      integer, intent(out) :: core_boundary_zone, envelope_boundary_zone
      logical, intent(out) :: fully_convective_flag
      logical, intent(out) :: surface_cz_active
      integer, intent(out) :: radiative_zone_bounds(13,2), &
           convective_zone_bounds(12,2)
      integer, intent(out) :: num_radiative_zones, num_convective_zones
      integer, intent(out) :: ierr
! --- locals ---
! (deuterium_rate_mid/_start live in rot_scr: carried across sub-steps)
      integer :: num_species_tracked
      integer :: i, j, k
      double precision :: step_fraction_ratio
! the composition rate vectors of the midpoint model, for burn_and_settle
      double precision :: rate_by_zone(num_rate_rows, json)
      double precision :: total_epsilon
      integer, allocatable :: regions(:,:)
      integer :: solid_zone_start, solid_zone_end

      ierr = 0

      num_species_tracked = 11
      if (star%job%use_extended_composition) num_species_tracked = 15
!  INITIALIZE COMPOSITION ARRAY THE FIRST TIME THROUGH (RESET star%xa
!  TO THE START-OF-TIMESTEP COMPOSITION star%xa_start).
! MHP 06/02 SAVE THE CONVECTIVE STATE FROM THE PRIOR STEP.
      if (first_call) then
         do j = 1,star%nz
            rot_scr%convective_flag_prev(j) = star%convective_flag_start(j)
            mix_scr%amum(j) = rot_scr%old_amu(j)
            do i = 1,num_species_tracked
               star%xa(i,j) = star%xa_start(i,j)
            end do
         end do
      else
         do j = 1,star%nz
            rot_scr%convective_flag_prev(j) = rot_scr%convective_flag_mid(j)
         end do
      endif
!  INTERPOLATE IN THE MODEL VARIABLES AND AUXILLARY PHYSICS.
      step_fraction_ratio = sub_timestep/full_timestep
      do j = 1,star%nz
         rot_scr%log_density_mid(j) = star%logRho_start(j) + &
              time_fraction*(star%logRho(j)-star%logRho_start(j))
         rot_scr%luminosity_lsun_mid(j) = star%luminosity_lsun_start(j) + &
              time_fraction*(star%luminosity_lsun(j)-star%luminosity_lsun_start(j))
         rot_scr%log_pressure_mid(j) = star%logP_start(j) + &
              time_fraction*(star%logP(j)-star%logP_start(j))
         rot_scr%log_radius_mid(j) = star%logR_start(j) + &
              time_fraction*(star%logR(j)-star%logR_start(j))
         rot_scr%log_temperature_mid(j) = star%logT_start(j) + &
              time_fraction*(star%logT(j)-star%logT_start(j))
         rot_scr%hg_mid(j) = star%old_hg(j) + time_fraction*(star%mean_gravity(j) - star%old_hg(j))
         mix_scr%del_adiabatic_mix(j) = rot_scr%old_del_adiabatic_mix(j) + &
              time_fraction*(star%grada(j)-rot_scr%old_del_adiabatic_mix(j))
         mix_scr%delm(j) = rot_scr%old_delm(j) + time_fraction*(star%gradT(j) - rot_scr%old_delm(j))
         mix_scr%del_radiative_mix(j) = rot_scr%old_del_radiative_mix(j) + &
              time_fraction*(star%gradr(j) - rot_scr%old_del_radiative_mix(j))
         mix_scr%esumm(j) = rot_scr%old_esum(j) + time_fraction*(star%eps_total(j) - rot_scr%old_esum(j))
         mix_scr%viscm(j) = rot_scr%old_visc(j) + time_fraction*(star%visc(j) - rot_scr%old_visc(j))
         mix_scr%thdifm(j) = rot_scr%old_thdif(j) + time_fraction*(star%thdif(j) - rot_scr%old_thdif(j))
         mix_scr%cpm(j) = rot_scr%old_cp(j) + time_fraction*(star%cp(j) - rot_scr%old_cp(j))
         mix_scr%qdtm(j) = rot_scr%old_qdt(j) + time_fraction*(star%qdt(j) - rot_scr%old_qdt(j))
         mix_scr%om(j) = rot_scr%old_om(j) + time_fraction*(star%opacity_zone(j) - rot_scr%old_om(j))
         mix_scr%amum(j) = mix_scr%amum(j) + step_fraction_ratio*(star%mu(j) - rot_scr%old_amu(j))
! MHP 6/00 ADDED TOTAL ENERGY GENERATION
         total_epsilon = star%eps_total(j)+star%eps_channels(i_eps_neu,j)+star%eps_channels(i_eps_grav,j)
         mix_scr%epsm(j) = rot_scr%old_eps(j)+time_fraction*(total_epsilon-rot_scr%old_eps(j))
      end do
!  CHECK FOR ADVANCING OR RECEDING CONVECTIVE REGIONS.USE INTERPOLATED
!  RADIATIVE AND ADIABATIC TEMPERATURE GRADIENTS TO DETERMINE WHETHER
!  OR NOT A ZONE IS CONVECTIVE IF IT CHANGES STATE OVER THE COURSE OF A
!  TIMESTEP.
      do i = 1,star%nz
         if (star%convective_flag(i).eqv.rot_scr%convective_flag_prev(i)) then
            rot_scr%convective_flag_mid(i) = star%convective_flag(i)
            mix_scr%velm(i) = rot_scr%old_vel(i) + time_fraction*(star%conv_vel(i)-rot_scr%old_vel(i))
         else
            if (mix_scr%del_adiabatic_mix(i).lt.mix_scr%del_radiative_mix(i)) then
               rot_scr%convective_flag_mid(i) = .true.
               mix_scr%velm(i) = max(rot_scr%old_vel(i),star%conv_vel(i))
            else
               rot_scr%convective_flag_mid(i) = .false.
               mix_scr%velm(i) = 0.0D0
            endif
         endif
      end do
! (THE ORIGINAL'S MHP 06/02 RE-DISTRIBUTION OF J FOR SHELLS THAT
!  DROPPED OUT OF A RECEDING CONVECTION ZONE WAS HARDWIRED OFF AND IS
!  NOT COMPUTED.)
!  CONVECTIVE OVERSHOOT APPLIED TO NORMAL CONVECTION ZONES.
      call am_convective_regions(star%xa,rot_scr%log_density_mid,rot_scr%log_pressure_mid,rot_scr%log_radius_mid, &
           star%log_mass,rot_scr%log_temperature_mid,rot_scr%convective_flag_mid,star%nz, &
           rot_scr%am_transport_convective_flag_mid,radiative_zone_bounds, &
           convective_zone_bounds,num_radiative_zones,num_convective_zones, star%iov1, star%iov2, ierr)
      if (ierr /= 0) return
!  REBURN THE ORIGINAL MIXTURE FOR THE SMALL DIFFUSION TIME STEP.
! MHP 6/00 - FOR THERMODYNAMIC CONSISTENCY USE CURRENT T, NOT END OF
! TIMESTEP T (REPLACED HT IN CALL WITH HTM)
! the network's rates on the midpoint model
      call network_rates_of_model(rot_scr%log_density_mid, &
           rot_scr%log_temperature_mid, star%xa, star%nz, rate_by_zone)
      call burn_and_settle(sub_timestep,star%xa,star%dm,rot_scr%log_temperature_mid, &
           rate_by_zone, star%nz,radiative_zone_bounds,convective_zone_bounds, &
           num_radiative_zones,num_convective_zones &
           ,star%log_total_mass,rot_scr%log_density_mid,star%log_mass,rot_scr%log_radius_mid, &
           rot_scr%log_pressure_mid,rot_scr%am_transport_convective_flag_mid,star%m, ierr)
      if (ierr /= 0) return
      do i = 1,star%nz
         do j = 1,num_species_tracked
            rot_scr%composition_snapshot(j,i) = star%xa(j,i)
         end do
      end do
! MHP 05/02 ADDED DEUTERIUM BURNING
      if (star%xa_start(i_h2,star%nz).gt.1.0D-14) then
! INCREMENT THE TIMESTEP
         if (first_call) then
            do i = 1,star%nz
               rot_scr%deuterium_rate_mid(i) = star%deuterium_burning_rate_start(i)+ &
                    step_fraction_ratio*(star%deuterium_burning_rate(i)- &
                    star%deuterium_burning_rate_start(i))
               rot_scr%deuterium_rate_mid_start(i) = star%deuterium_burning_rate_start(i)
            end do
         else
            do i = 1,star%nz
               rot_scr%deuterium_rate_mid_start(i) = rot_scr%deuterium_rate_mid(i)
               rot_scr%deuterium_rate_mid(i) = rot_scr%deuterium_rate_mid(i)+ &
                    step_fraction_ratio*(star%deuterium_burning_rate(i)- &
                    star%deuterium_burning_rate_start(i))
            end do
         endif
! EACH BURN REGION OF THE MIDPOINT MODEL (burn_regions).
         regions = burn_regions(rot_scr%log_temperature_mid, &
              radiative_zone_bounds, num_radiative_zones, &
              convective_zone_bounds, num_convective_zones)
         do k = 1, size(regions,2)
            call burn_deuterium_region(regions(1,k), regions(2,k), star%nz, &
                 star%dm, star%xa, sub_timestep, star%xa_start, &
                 rot_scr%deuterium_rate_mid_start, rot_scr%deuterium_rate_mid, &
                 surface_accretion_now(step_fraction_ratio*star%accreted_mass), &
                 min_deuterium_substep)
         end do
      endif
!  DETERMINE EXTENT OF CENTRAL CONVECTION ZONE.
!  IMIN IS THE FIRST ZONE ABOVE A CENTRAL CONVECTION ZONE, AND THUS THE
!  FIRST ZONE CONSIDERED FOR STABILITY AGAINST ROTATIONAL INSTABILITIES.
      if (rot_scr%am_transport_convective_flag_mid(1)) then
         do i = 2,star%nz
            if (.not.rot_scr%am_transport_convective_flag_mid(i)) exit
         end do
         core_boundary_zone = max(2,i-1)
      else
         core_boundary_zone = 2
      endif
!  DETERMINE EXTENT OF SURFACE CONVECTION ZONE.
      fully_convective_flag = .false.
      if (rot_scr%am_transport_convective_flag_mid(star%nz)) then
!  SURFACE C.Z. EXISTS.  FIND LOWEST SHELL (IMAX), WHICH IS ALSO THE
!  UPPERMOST ZONE CONSIDERED FOR STABILITY AGAINST ROTATIONALLY INDUCED MIXING.
         do i = star%nz-1,1,-1
            if (.not.rot_scr%am_transport_convective_flag_mid(i)) exit
         end do
!  (LOOP COMPLETED WITHOUT EXIT, I.E. I = 0: FULLY CONVECTIVE.)
         if (i < 1) fully_convective_flag = .true.
         envelope_boundary_zone = i + 1
!  HSTOP IS THE MASS AT THE TOP OF THE C.Z.
!  HSBOT IS THE MASS AT THE BOTTOM OF THE C.Z.
         cz_mass_top = exp(ln10*star%log_total_mass)
         if (envelope_boundary_zone.gt.1) then
            cz_mass_bottom = 0.5D0*(star%m(envelope_boundary_zone)+ &
                 star%m(envelope_boundary_zone-1))
         else
            cz_mass_bottom = 0.0D0
         endif
!  LCZSUR=T IF A SURFACE C.Z.DEEP ENOUGH FOR ANGULAR MOMENTUM LOSS EXISTS
         if ((cz_mass_top-cz_mass_bottom)/star%solar_mass_cgs.gt.0.0D0) then
            surface_cz_active = .true.
         else
            surface_cz_active= .false.
         endif
      else
!  NO SURFACE C.Z.
         envelope_boundary_zone = star%nz
         surface_cz_active = .false.
      endif
! JNT 09/25 FOR 05/15 IMPJMOD=1 IS THE SAME AS LSOLID
      if (.not.star%ctrl%force_solid_body_rotation .and. (star%ctrl%solid_body_mode_flag.ne.1)) then
!  NOW FIND THE RUN OF ROTATION VARIABLES THAT ARE CONSISTENT WITH THE
!  INTERMEDIATE STRUCTURE AND THE RUN OF SPECIFIC ANGULAR MOMENTUM J/M.
!  J/M = I/M * OMEGA.
!  FIRST GUESS AT MOMENT OF INERTIA OF DISTORTED SPERICAL SHELLS:
!  I/M = 2/3 R**2.
         do i = 1,star%nz
            rot_scr%moment_of_inertia_mid(i) = cc23*exp(ln10*2.0D0*rot_scr%log_radius_mid(i))
            rot_scr%omega_mid(i) = star%j_rot(i)/rot_scr%moment_of_inertia_mid(i)
            rot_scr%moment_of_inertia_mid(i) = star%dm(i)*rot_scr%moment_of_inertia_mid(i)
         end do
!  SOLVE FOR THE ANGULAR VELOCITIES OF THE SHELLS GIVEN THE SPECIFIC
!  ANGULAR MOMENTUM (HJM) AND A FIRST GUESS AT THE ANGULAR VELOCITY(OMEGAM)
!  AND MOMENT OF INERTIA (HIM).
         call omega_from_j(rot_scr%log_density_mid,star%j_rot,rot_scr%log_radius_mid, &
              star%log_mass,star%dm,rot_scr%am_transport_convective_flag_mid,star%nz, &
              rot_scr%eta_squared_mid,rot_scr%moment_of_inertia_mid,rot_scr%omega_mid,rot_scr%qiw_mid, &
              rot_scr%mean_radius_mid)
!  FIND TOTAL MOMENT OF INERTIA OF THE SURFACE C.Z. IF APPLICABLE.
         if (surface_cz_active) then
            moment_of_inertia_cz = rot_scr%moment_of_inertia_mid(star%nz)
            do i = star%nz-1,1,-1
               if (.not.rot_scr%am_transport_convective_flag_mid(i)) exit
               moment_of_inertia_cz = moment_of_inertia_cz + rot_scr%moment_of_inertia_mid(i)
            end do
         endif
      else
         solid_zone_start = 1
         solid_zone_end = star%nz
         call solid_body_omega(rot_scr%log_density_mid,star%j_rot,rot_scr%log_radius_mid, &
              star%log_mass,star%dm,solid_zone_start,solid_zone_end, &
              rot_scr%eta_squared_mid,rot_scr%moment_of_inertia_mid,rot_scr%omega_mid,rot_scr%qiw_mid, &
              rot_scr%mean_radius_mid,star%nz)
         if (surface_cz_active) then
            moment_of_inertia_cz = rot_scr%moment_of_inertia_mid(1)
            do i = 2,star%nz
               moment_of_inertia_cz = moment_of_inertia_cz + rot_scr%moment_of_inertia_mid(i)
            end do
         endif
      endif
      return
end subroutine advance_rotation_substep

! ---------------------------------------------------------------
! Burns each burn region (burn_step's burn_regions: shell by shell in
! the radiative regions, each convection zone as a unit, homogenized) by
! burn_network_region. Then (if diffuse_helium_active) performs gravitational settling of
! helium (Bahcall & Loeb method), and finally renormalizes the
! composition to guard against small negative/overflowing abundances.
! rate_by_zone holds the composition rate vectors of the model being
! burned (network_rates_of_model: mid_timestep_model passes the rotation
! sub-step's midpoint).
subroutine burn_and_settle(timestep, composition, shell_mass, log_temperature, &
     rate_by_zone, num_zones, radiative_zone_bounds, convective_zone_bounds, &
     num_radiative_zones, num_convective_zones, log_total_mass, &
     log_density, log_mass, log_radius, log_pressure, convective_flag, &
     enclosed_mass, ierr)
      use rotation_scratch_lib

      use star_info_lib, only: star, json
      use io_units
      use run_log, only: solver_diagnostics
      use physical_constants
      use math_lib
      use burn_step, only: burn_network_region, burn_regions
      use element_diffusion_lib, only: diffuse_elements
      use net_lib, only: num_rate_rows
      implicit none

      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: shell_mass(json), log_temperature(json)
      double precision, intent(in) :: rate_by_zone(num_rate_rows, json)
      integer, intent(in) :: num_zones
      integer, intent(in) :: radiative_zone_bounds(13,2), &
           convective_zone_bounds(12,2)
      integer, intent(in) :: num_radiative_zones, num_convective_zones
      double precision, intent(in) :: log_total_mass
      double precision, intent(in) :: log_density(json), log_mass(json), &
           log_radius(json), log_pressure(json)
      logical, intent(in) :: convective_flag(json)
! inout: the settling routines convert it to Bahcall & Loeb units in
! place and back, which can change its last bit
      double precision, intent(inout) :: enclosed_mass(json)

      integer :: num_species
      double precision :: dlnp_dr_settling(json)
      logical :: am_transport_convective_flag(json)
      double precision :: total_mass
      integer, allocatable :: regions(:,:)
      integer :: zone_idx, region_idx, species_idx, outer_boundary_zone, num_settling_substeps, &
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
! THE NEW ABUNDANCES AT THE END OF THE (SUB-)STEP, timestep (s).
      regions = burn_regions(log_temperature, radiative_zone_bounds, &
           num_radiative_zones, convective_zone_bounds, num_convective_zones)
      do region_idx = 1, size(regions,2)
         call burn_network_region(log_temperature, regions(1,region_idx), &
              regions(2,region_idx), rate_by_zone, shell_mass, composition, &
              timestep, ierr)
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
            ! a nonzero ierr stops the run
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
! gravitational settling, with the mid-step temperature gradient
            call diffuse_elements(settling_dt,composition,dlnp_dr_settling, &
                 log_radius,log_density,enclosed_mass,log_temperature, &
                 mix_scr%delm,am_transport_convective_flag,num_zones, &
                 total_mass, ierr)
            if (ierr /= 0) return
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
end subroutine burn_and_settle

end module mid_timestep_model

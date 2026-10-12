!----------------------------------------------------------------------
! mid_timestep_model
!----------------------------------------------------------------------
! The model of one rotation sub-step. A sub-step of the angular momentum
! transport is shorter than the model step, so advance_rotation_substep
! builds the model at a point inside the step by linear interpolation in
! time between the step's start and end. It then burns and settles the
! composition over the sub-step (burn_and_settle), burns deuterium,
! locates the convection zones (the edge of a central one; the base,
! mass and moment of inertia of the surface one), and finds the omega
! that keeps each shell's specific angular momentum on the new radii.
! The result goes into rot_scr (log_density_mid ...
! am_transport_convective_flag_mid) and mid_scr, whose fields carry the
! names of the star% fields they interpolate; the rest of the sub-step
! reads it from there.
module mid_timestep_model
      use rotation_configuration, only: rotation_settings
      use am_regions, only: am_convective_regions
      use omega_profile, only: solid_body_omega, omega_from_j
      implicit none
      private :: burn_and_settle
contains

! ---------------------------------------------------------------
! The model of one rotation sub-step.
! Inputs: the model step full_timestep and the sub-step sub_timestep
! (s); time_fraction, the point of the step to interpolate to
! (evolve_angular_momentum passes the fraction of the step elapsed at
! the end of this sub-step); first_call on the step's first sub-step.
! Method: interpolate the structure and the auxiliary physics linearly
! in time (rot_scr%*_mid, mid_scr); decide the convective state of the
! zones that changed it during the step; find the convection zones with
! their overshoot (am_convective_regions); burn and settle over the
! sub-step (burn_and_settle) and keep the result in
! rot_scr%composition_snapshot; burn deuterium; locate the central and
! surface convection zones; find omega from j on the new structure.
! Outputs: core_boundary_zone and envelope_boundary_zone, the first and
! last interfaces the transport examines; fully_convective_flag; the
! surface convection zone (surface_cz_active, cz_mass_bottom and
! cz_mass_top in g, moment_of_inertia_cz in g cm^2); the radiative and
! convective zone tables. A nonzero ierr stops the run.
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
! (deuterium_rate_mid and deuterium_rate_mid_start live in rot_scr,
! carried from one sub-step to the next)
      integer :: num_species_tracked
      integer :: i, j, k
      double precision :: step_fraction_ratio
! the composition rate vectors of the interpolated model, for
! burn_and_settle
      double precision :: rate_by_zone(num_rate_rows, json)
      double precision :: total_epsilon
      integer, allocatable :: regions(:,:)
      integer :: solid_zone_start, solid_zone_end

      ierr = 0

      num_species_tracked = 11
      if (rotation_settings%use_extended_composition) num_species_tracked = 15
      call start_substep_composition()
      call interpolate_structure()
      call settle_convective_state()
! A shell left behind by a receding convection zone keeps its j/m. (The
! original's correction of j for such shells, MHP 6/02, was hardwired
! off and is not implemented.)
! The convection zones of the interpolated model with their overshoot,
! the zones convective for angular momentum transport
! (am_transport_convective_flag_mid), and the zone tables.
      call am_convective_regions(star%xa,rot_scr%log_density_mid,rot_scr%log_pressure_mid,rot_scr%log_radius_mid, &
           star%log_mass,rot_scr%log_temperature_mid,rot_scr%convective_flag_mid,star%nz, &
           rot_scr%am_transport_convective_flag_mid,radiative_zone_bounds, &
           convective_zone_bounds,num_radiative_zones,num_convective_zones, star%iov1, star%iov2, ierr)
      if (ierr /= 0) return
! Burn the start-of-sub-step composition (star%xa_start) over the
! sub-step. For thermodynamic consistency (MHP 6/00) the network rates
! and temperatures are those of the interpolated model, not those of
! the step's end.
      call network_rates_of_model(rot_scr%log_density_mid, &
           rot_scr%log_temperature_mid, star%xa, star%nz, rate_by_zone)
      call burn_and_settle(sub_timestep,star%xa,star%dm,rot_scr%log_temperature_mid, &
           rate_by_zone, star%nz,radiative_zone_bounds,convective_zone_bounds, &
           num_radiative_zones,num_convective_zones &
           ,star%log_total_mass,rot_scr%log_density_mid,star%log_mass,rot_scr%log_radius_mid, &
           rot_scr%log_pressure_mid,rot_scr%am_transport_convective_flag_mid,star%m, ierr)
      if (ierr /= 0) return
! The burned and settled composition, which secular_transport restores
! (slots 1-4) at the start of each of its later iterations.
      do i = 1,star%nz
         do j = 1,num_species_tracked
            rot_scr%composition_snapshot(j,i) = star%xa(j,i)
         end do
      end do
! Deuterium burns while the surface still has some at the start of the
! sub-step (X_D > 1e-14).
      if (star%xa_start(i_h2,star%nz).gt.1.0D-14) call burn_deuterium_over_substep()
      call locate_convection_zones()
      call rotate_midpoint_model()

contains

! ---------------------------------------------------------------
! The state the sub-step starts from. On the step's first sub-step: the
! start-of-step composition (star%xa = star%xa_start), the start-of-step
! convective flags and the previous model's mean molecular weight. On
! later sub-steps: the previous sub-step's convective flags (MHP 6/02).
subroutine start_substep_composition()
      if (first_call) then
         do j = 1,star%nz
            rot_scr%convective_flag_prev(j) = star%convective_flag_start(j)
            mid_scr%mu(j) = rot_scr%old_mu(j)
            do i = 1,num_species_tracked
               star%xa(i,j) = star%xa_start(i,j)
            end do
         end do
      else
         do j = 1,star%nz
            rot_scr%convective_flag_prev(j) = rot_scr%convective_flag_mid(j)
         end do
      endif
end subroutine start_substep_composition

! ---------------------------------------------------------------
! The structure and the auxiliary physics interpolated linearly in time
! between the start and the end of the step, at f = time_fraction:
!   q_mid = q_start + f (q_end - q_start)
! for log rho, L, log P, log r, log T, the mean gravity and mid_scr's
! gradients, energy generation, viscosity, thermal diffusivity, c_P,
! qdt and opacity (start values from star%*_start, star%old_hg and
! rot_scr%old_*). The mean molecular weight instead accumulates each
! sub-step's share of its change,
!   mu_mid = mu_mid + (sub_timestep/full_timestep) (mu_end - mu_old),
! which reaches the same point of the step (the caller takes the share
! back out on a redo).
subroutine interpolate_structure()
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
         rot_scr%mean_gravity_mid(j) = star%old_hg(j) + time_fraction*(star%mean_gravity(j) - star%old_hg(j))
         mid_scr%grada(j) = rot_scr%old_grada(j) + &
              time_fraction*(star%grada(j)-rot_scr%old_grada(j))
         mid_scr%gradT(j) = rot_scr%old_gradT(j) + time_fraction*(star%gradT(j) - rot_scr%old_gradT(j))
         mid_scr%gradr(j) = rot_scr%old_gradr(j) + &
              time_fraction*(star%gradr(j) - rot_scr%old_gradr(j))
         mid_scr%eps_total(j) = rot_scr%old_eps_total(j) + time_fraction*(star%eps_total(j) - rot_scr%old_eps_total(j))
         mid_scr%visc(j) = rot_scr%old_visc(j) + time_fraction*(star%visc(j) - rot_scr%old_visc(j))
         mid_scr%thdif(j) = rot_scr%old_thdif(j) + time_fraction*(star%thdif(j) - rot_scr%old_thdif(j))
         mid_scr%cp(j) = rot_scr%old_cp(j) + time_fraction*(star%cp(j) - rot_scr%old_cp(j))
         mid_scr%qdt(j) = rot_scr%old_qdt(j) + time_fraction*(star%qdt(j) - rot_scr%old_qdt(j))
         mid_scr%opacity(j) = rot_scr%old_opacity(j) + time_fraction*(star%opacity_zone(j) - rot_scr%old_opacity(j))
         mid_scr%mu(j) = mid_scr%mu(j) + step_fraction_ratio*(star%mu(j) - rot_scr%old_mu(j))
! The total energy generation, eps_total + eps_neu + eps_grav with the
! stored signs (MHP 6/00).
         total_epsilon = star%eps_total(j)+star%eps_channels(i_eps_neu,j)+star%eps_channels(i_eps_grav,j)
         mid_scr%eps_plus_neu_grav(j) = rot_scr%old_eps_plus_neu_grav(j)+time_fraction*(total_epsilon-rot_scr%old_eps_plus_neu_grav(j))
      end do
end subroutine interpolate_structure

! ---------------------------------------------------------------
! The convective state of each zone at the interpolated model, for
! advancing and receding convection zones. A zone convective (or not)
! both at the step's end and at the previous sub-step keeps that state,
! its convective velocity interpolated in time. A zone that changed
! state is convective if the interpolated gradients satisfy
! grad_ad < grad_rad (Schwarzschild), with the larger of the old and new
! convective velocities; otherwise it is radiative with zero velocity.
subroutine settle_convective_state()
      do i = 1,star%nz
         if (star%convective_flag(i).eqv.rot_scr%convective_flag_prev(i)) then
            rot_scr%convective_flag_mid(i) = star%convective_flag(i)
            mid_scr%conv_vel(i) = rot_scr%old_conv_vel(i) + time_fraction*(star%conv_vel(i)-rot_scr%old_conv_vel(i))
         else
            if (mid_scr%grada(i).lt.mid_scr%gradr(i)) then
               rot_scr%convective_flag_mid(i) = .true.
               mid_scr%conv_vel(i) = max(rot_scr%old_conv_vel(i),star%conv_vel(i))
            else
               rot_scr%convective_flag_mid(i) = .false.
               mid_scr%conv_vel(i) = 0.0D0
            endif
         endif
      end do
end subroutine settle_convective_state

! ---------------------------------------------------------------
! Deuterium burning over the sub-step (MHP 5/02). The burning rate is
! linear in time between its start- and end-of-step values r_0 and r_1:
! deuterium_rate_mid_start is its value at the start of the sub-step
! (the previous sub-step's end), and deuterium_rate_mid advances by the
! sub-step's share,
!   rate_mid = rate_mid + (sub_timestep/full_timestep) (r_1 - r_0),
! from r_0 on the first sub-step.
subroutine burn_deuterium_over_substep()
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
! Each burn region of the interpolated model (burn_regions: a radiative
! shell hotter than the nuclear cutoff, or a mixed region as a unit),
! with the sub-step's share of the accreted mass. burn_deuterium_region
! starts from star%xa_start, so in a region that still has deuterium
! its X and He3 replace those burn_and_settle left.
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
end subroutine burn_deuterium_over_substep

! ---------------------------------------------------------------
! The convection zones that bound the transport, from the flags of
! angular momentum transport (overshoot included).
subroutine locate_convection_zones()
! The central convection zone. core_boundary_zone (legacy IMIN) indexes
! the first interface the transport examines: max(2, i-1), with i the
! first zone that is not convective. That is the outermost zone of a
! central convection zone, so the first interface examined is the one
! below the zone's edge (interface i); 2 without a central convection
! zone.
      if (rot_scr%am_transport_convective_flag_mid(1)) then
         do i = 2,star%nz
            if (.not.rot_scr%am_transport_convective_flag_mid(i)) exit
         end do
         core_boundary_zone = max(2,i-1)
      else
         core_boundary_zone = 2
      endif
! The surface convection zone.
      fully_convective_flag = .false.
      if (rot_scr%am_transport_convective_flag_mid(star%nz)) then
! envelope_boundary_zone (legacy IMAX) is its lowest zone; the interface
! below that zone, the zone's base, is the last interface the transport
! examines.
         do i = star%nz-1,1,-1
            if (.not.rot_scr%am_transport_convective_flag_mid(i)) exit
         end do
! No radiative zone found (i = 0): the model is fully convective.
         if (i < 1) fully_convective_flag = .true.
         envelope_boundary_zone = i + 1
! The masses (g) at the top of the convection zone (the total mass) and
! at its base (between its lowest zone and the zone below; zero when it
! reaches the centre).
         cz_mass_top = exp(ln10*star%log_total_mass)
         if (envelope_boundary_zone.gt.1) then
            cz_mass_bottom = 0.5D0*(star%m(envelope_boundary_zone)+ &
                 star%m(envelope_boundary_zone-1))
         else
            cz_mass_bottom = 0.0D0
         endif
! surface_cz_active (legacy LCZSUR): the surface convection zone is deep
! enough for the wind's angular momentum loss; the threshold on its mass
! is zero, so any convection zone of positive mass qualifies.
         if ((cz_mass_top-cz_mass_bottom)/star%solar_mass_cgs.gt.0.0D0) then
            surface_cz_active = .true.
         else
            surface_cz_active= .false.
         endif
      else
! No surface convection zone: the transport runs to the last interface.
         envelope_boundary_zone = star%nz
         surface_cz_active = .false.
      endif
end subroutine locate_convection_zones

! ---------------------------------------------------------------
! The rotation of the interpolated model: omega, the moments of inertia
! and the shape (eta_squared_mid, mean_radius_mid, di_domega_mid) that
! carry each shell's specific angular momentum star%j_rot on the new
! radii, and the moment of inertia of the surface convection zone.
! solid_body_mode_flag = 1 acts like force_solid_body_rotation.
subroutine rotate_midpoint_model()
      if (.not.rotation_settings%force_solid_body_rotation .and. (rotation_settings%solid_body_mode_flag.ne.1)) then
! Differential rotation: j/m = (I/m) omega shell by shell. First guess,
! the thin spherical shell: I/m = (2/3) r^2, omega = j/(I/m), and
! I = dm (I/m).
         do i = 1,star%nz
            rot_scr%moment_of_inertia_mid(i) = cc23*exp(ln10*2.0D0*rot_scr%log_radius_mid(i))
            rot_scr%omega_mid(i) = star%j_rot(i)/rot_scr%moment_of_inertia_mid(i)
            rot_scr%moment_of_inertia_mid(i) = star%dm(i)*rot_scr%moment_of_inertia_mid(i)
         end do
! omega of each shell from star%j_rot and that first guess, iterating
! with the moments of inertia of the distorted shells; the convection
! zones rotate by their rotation law (omega_from_j).
         call omega_from_j(rot_scr%log_density_mid,star%j_rot,rot_scr%log_radius_mid, &
              star%log_mass,star%dm,rot_scr%am_transport_convective_flag_mid,star%nz, &
              rot_scr%eta_squared_mid,rot_scr%moment_of_inertia_mid,rot_scr%omega_mid,rot_scr%di_domega_mid, &
              rot_scr%mean_radius_mid)
! The moment of inertia of the surface convection zone, if there is one.
         if (surface_cz_active) then
            moment_of_inertia_cz = rot_scr%moment_of_inertia_mid(star%nz)
            do i = star%nz-1,1,-1
               if (.not.rot_scr%am_transport_convective_flag_mid(i)) exit
               moment_of_inertia_cz = moment_of_inertia_cz + rot_scr%moment_of_inertia_mid(i)
            end do
         endif
      else
! Solid-body rotation: one omega for the whole star, which the wind then
! spins down as a unit, so moment_of_inertia_cz is the star's total.
         solid_zone_start = 1
         solid_zone_end = star%nz
         call solid_body_omega(rot_scr%log_density_mid,star%j_rot,rot_scr%log_radius_mid, &
              star%log_mass,star%dm,solid_zone_start,solid_zone_end, &
              rot_scr%eta_squared_mid,rot_scr%moment_of_inertia_mid,rot_scr%omega_mid,rot_scr%di_domega_mid, &
              rot_scr%mean_radius_mid,star%nz)
         if (surface_cz_active) then
            moment_of_inertia_cz = rot_scr%moment_of_inertia_mid(1)
            do i = 2,star%nz
               moment_of_inertia_cz = moment_of_inertia_cz + rot_scr%moment_of_inertia_mid(i)
            end do
         endif
      endif
end subroutine rotate_midpoint_model

end subroutine advance_rotation_substep

! ---------------------------------------------------------------
! The composition over a rotation sub-step of timestep (s), in place in
! composition (star%xa layout).
! 1. Burn each burn region (burn_step's burn_regions: each radiative
!    shell hotter than the nuclear cutoff alone, each mixed region as a
!    unit, homogenized) from the start-of-sub-step composition
!    star%xa_start (burn_network_region). rate_by_zone holds the
!    composition rate vectors of the model being burned
!    (network_rates_of_model; advance_rotation_substep passes those of
!    its interpolated model) and log_temperature its temperatures.
! 2. While diffuse_helium_active, settle helium (settle_helium).
! 3. Clip and renormalize the abundances (renormalize_composition).
! A nonzero ierr stops the run.
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

! num_species (legacy NSPEC): the species tracked, 11, or 15 with the
! extended composition
      integer :: num_species
      double precision :: dlnp_dr_settling(json)
      logical :: am_transport_convective_flag(json)
      double precision :: total_mass
      integer, allocatable :: regions(:,:)
      integer :: zone_idx, region_idx, species_idx, outer_boundary_zone, num_settling_substeps, &
           substep_idx
      double precision :: mass_fraction_above, settling_timescale, &
           max_settling_dt, settling_dt

      integer, intent(out) :: ierr

      ierr = 0

      if (rotation_settings%use_extended_composition) then
         num_species = 15
      else
         num_species = 11
      end if
! The abundances at the end of the sub-step, region by region.
      regions = burn_regions(log_temperature, radiative_zone_bounds, &
           num_radiative_zones, convective_zone_bounds, num_convective_zones)
      do region_idx = 1, size(regions,2)
         call burn_network_region(log_temperature, regions(1,region_idx), &
              regions(2,region_idx), rate_by_zone, shell_mass, composition, &
              timestep, ierr)
         if (ierr /= 0) return
      end do
      if (star%job%diffuse_helium_active) then
         call settle_helium()
         if (ierr /= 0) return
      end if
      call renormalize_composition()

contains

! ---------------------------------------------------------------
! Gravitational settling of helium over the sub-step (BC 5/92), by
! element_diffusion_lib's diffuse_elements (the Bahcall & Loeb
! formalism), in sub-steps no longer than settling_timestep_fraction of
! the settling timescale at the outer boundary (Michaud et al. 1984;
! MHP 6/90). The outer boundary is the top of the outermost radiative
! region, just below the surface convection zone. Settling stops for
! good once the central X drops below hydrogen_diffusion_floor; with no
! surface convection zone the run stops (ierr = 1).
subroutine settle_helium()
      settling: do
         if (composition(1,1).lt.rotation_settings%hydrogen_diffusion_floor) then
            star%job%diffuse_helium_active=.false.
            exit settling
         end if
! Settling needs a surface convection zone.
         if (.not.convective_flag(num_zones)) then
            write(*,911)
            write(run_log_unit,911)
  911       format(1x,'NO SURFACE CZ - DIFFUSION NOT MEANINGFUL'/ &
                 'STOPPED IN SUBROUTINE MIX')
            ! a nonzero ierr stops the run
            ierr = 1
            return
         end if
! No settling when the outermost mixed region, with its overshoot, runs
! from zone 2 or below to the surface: the model is fully convective.
! This avoids an occasional pre-main-sequence crash in the first
! radiative model (MHP 7/92).
         if (convective_zone_bounds(num_convective_zones,1).gt.2 .or. convective_zone_bounds(num_convective_zones,2).ne.num_zones) then
! The outer boundary: the top zone of the outermost radiative region.
! Nothing settles if no zone at or below it has Y above
! helium_diffusion_min (legacy YMIN).
         outer_boundary_zone = radiative_zone_bounds(num_radiative_zones,2)
         do zone_idx = outer_boundary_zone,1,-1
            if (composition(2,zone_idx).gt.rotation_settings%helium_diffusion_min) exit
         end do
         if (zone_idx < (1)) then
         exit settling
         end if
! The total mass (g), and d ln P/dr = -rho G m/(r^2 P) from hydrostatic
! equilibrium.
         total_mass=exp(ln10*log_total_mass)
         do zone_idx = 1,num_zones
            dlnp_dr_settling(zone_idx)=-exp(ln10*(log_density(zone_idx)+ &
                 cgl+log_mass(zone_idx)-2.0d0*log_radius(zone_idx)- &
                 log_pressure(zone_idx)))
         end do
! The settling step is limited only when the outer boundary lies above
! zone 1 (BC 6/92).
         if (outer_boundary_zone.ge.2) then
! The mass fraction above the outer boundary, f_m = (M - m_b)/M.
            mass_fraction_above = (total_mass- &
                 enclosed_mass(outer_boundary_zone))/total_mass
! The settling timescale of helium at the outer boundary (Michaud et
! al. 1984, ApJ 282, 206), in s:
!   tau = 4.348e21 yr f_m/T_b^1.5.
            settling_timescale = 4.348d21*seconds_per_year* &
                 mass_fraction_above/ &
                 exp(ln10*1.5d0*log_temperature(outer_boundary_zone))
! Settling sub-steps no longer than settling_timestep_fraction (legacy
! DT_GS) times tau, and one if timestep is shorter. The count is
! int(timestep/max_settling_dt) + 1 unless max_settling_dt equals
! timestep: the mod test takes its arguments as (max_settling_dt,
! timestep), as in the original, so it adds one even when timestep is
! an exact multiple of max_settling_dt.
            if (solver_diagnostics()) write(run_log_unit,*) 'JMAX=', &
                 outer_boundary_zone,' FM= ', &
                 mass_fraction_above,' TSCALE=',settling_timescale
            max_settling_dt = rotation_settings%settling_timestep_fraction*settling_timescale
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
! The settling's convective flags: false over the radiative regions,
! true over the mixed ones (overshoot included).
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
! Settle, with the temperature gradient of the interpolated model
! (mid_scr%gradT).
            call diffuse_elements(settling_dt,composition,dlnp_dr_settling, &
                 log_radius,log_density,enclosed_mass,log_temperature, &
                 mid_scr%gradT,am_transport_convective_flag,num_zones, &
                 total_mass, ierr)
            if (ierr /= 0) return
         end do
         end if
      exit settling
      end do settling
end subroutine settle_helium

! ---------------------------------------------------------------
! Guards against anomalies such as small negative abundances: each
! abundance clipped to [0, 1], then Z = min(Z, 1 - X - He3) and
! Y = 1 - X - Z - He3, so that slots 1-4 sum to one.
subroutine renormalize_composition()
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
end subroutine renormalize_composition

end subroutine burn_and_settle

end module mid_timestep_model

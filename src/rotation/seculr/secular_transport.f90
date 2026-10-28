!----------------------------------------------------------------------
! secular_transport
!----------------------------------------------------------------------
! The transport of angular momentum, and the mixing of the composition
! it drives, by the secular rotational instabilities over one rotation
! sub-step (apply_secular_transport, legacy SECULR).
!
! apply_secular_transport is the iteration; each of its steps is a
! procedure of this module with its inputs and outputs in its argument
! list:
!   set_up_interfaces            once per sub-step: the structure at the
!                                interfaces, the stability coefficients
!                                and the quadrupole moment (into stab)
!   compute_velocities           each iteration: the mechanisms'
!                                velocities and their sum U
!   restart_iteration            from the second iteration: omega
!                                averaged with the sub-step's start,
!                                new velocities, j and composition back
!                                to the start
!   spin_down_uncoupled_surface  the wind on a surface convection zone
!                                the transport does not reach
!   transport_unstable_regions   each run of interfaces with U > 0 on
!                                its own grid: angular momentum
!                                (region_transport, with the wind's
!                                terms from magnetic_braking's
!                                outer_point_wind_terms), then the
!                                composition
!   check_iterate                the checks on j and the composition,
!                                and convergence
!   mix_composition              the rotational mixing of a range of
!                                species
!
! Each quantity has one name in rotation's drivers and their helpers
! (total_velocity is U, the sum of the mechanisms' velocities, legacy
! HV; the eq_ prefix marks a value on the transport region's equally
! spaced grid). The star-blind kernels keep generic names for what they
! are handed: am_diffusion_coeffs' diffusion_coeff and
! magnetic_braking's domega_start/domega_end for wind_loss_explicit/
! wind_loss_implicit (WIND1/WIND2).
! The sub-step's interpolated structure is read from sub_model
! (rotation_scratch_lib, filled by mid_timestep_model each sub-step),
! the interface structure and coefficients from stab, the transport
! grid from chi_grid, and the model's arrays and scalars (log_mass, m,
! dm, j_rot, xa, nz, log_L, log_total_mass, star_mass, log_Teff) from
! star.
module secular_transport
      use rotation_configuration, only: rotation_settings, disk_locking_engaged, &
           solid_body_everywhere, num_tracked_species, last_mu_species
      use transport_checks, only: check_mixed_species
      use rotation_scratch_lib, only: substep_model_state, substep_physics_state, &
           stability_state, transport_grid_state, circulation_velocity_state
      use star_info_lib, only: star_info, json
      use magnetic_braking, only: wind_conditions
      implicit none
      private
      public :: apply_secular_transport

! The surface's part in a sub-step's transport, which the spin-down of
! an uncoupled surface convection zone and the region that reaches the
! surface share: whether the model has a surface convection zone
! (cz_active) and disk locking is engaged (disk_locked), the surface
! omega at the sub-step's start, and the wind's terms for a region's
! outermost grid point (outer_point_wind_terms). surface_implicit is the
! surface region's implicit term from an earlier iteration of the
! sub-step, if surface_implicit_set, with which the next one is
! averaged; the per-region terms of a region below the surface (zero)
! do not touch it (docs/bug_sweep.md O27).
      type :: surface_wind_state
         logical :: cz_active = .false., disk_locked = .false.
         double precision :: omega_start = 0.0d0
         double precision :: loss_explicit = 0.0d0, loss_implicit = 0.0d0
         double precision :: surface_implicit = 0.0d0
         logical :: surface_implicit_set = .false.
      end type surface_wind_state

contains

! ---------------------------------------------------------------
! One sub-step (sub_timestep, s) of angular momentum transport and
! rotational mixing over the interfaces zone_min..zone_max (interface i
! lies between zones i-1 and i; mid_timestep_model's core_boundary_zone
! and envelope_boundary_zone).
! Method: the circulation and instability velocities U at the
! interfaces (total_velocity: Eddington-Sweet circulation, secular
! and dynamical shear, the GSF instability) give the diffusion
! coefficients D_omega and D_mix. Each run of interfaces with U > 0 is
! an unstable region, solved on its own grid, equally spaced in chi:
! the implicit diffusion equation for omega (with the wind when the
! region reaches the surface convection zone), mapped back to j/m, then
! the mixing of X, Y, Z and He3. U depends on omega, so the sub-step
! iterates, at most max_diffusion_iters times: each later iteration
! restarts from the saved j and composition with U from omega averaged
! between the sub-step's start and the previous iterate (and averaged
! with the previous U), until check_angular_momentum finds j converged.
! The other species mix once, at the end.
! Outputs: star%j_rot and star%xa; star%omega_mid and the shape from
! the new j; redo_flag, with sub_timestep halved and cut_count raised,
! when a check asks for a shorter sub-step; diffusion_solve_ok;
! omega_surface and cz_moment_of_inertia of the surface convection zone;
! mixing_diffusion_coeff and total_velocity; star%am_diffusion_coeff
! and star%mixing_diffusion_coeff for the profile. A nonzero ierr stops
! the run.
subroutine apply_secular_transport(s, sub_model, sub_physics, stab, chi_grid, circ_scr, &
     sub_timestep, specific_angular_momentum_saved, &
     zone_min, zone_max, redo_flag, cut_count, &
     cz_moment_of_inertia, cz_mass_bottom, cz_mass_top, omega_surface, &
     surface_cz_active, mixing_diffusion_coeff, total_velocity, &
     diffusion_solve_ok, ierr)
      type(star_info), intent(inout) :: s
      type(substep_model_state), intent(inout) :: sub_model
      type(substep_physics_state), intent(inout) :: sub_physics
      type(stability_state), intent(inout) :: stab
      type(transport_grid_state), intent(inout) :: chi_grid
      type(circulation_velocity_state), intent(inout) :: circ_scr
      double precision, intent(inout) :: sub_timestep
      double precision, intent(in) :: specific_angular_momentum_saved(json)
      integer, intent(in) :: zone_min, zone_max
      logical, intent(out) :: redo_flag
      integer, intent(inout) :: cut_count
! cz_moment_of_inertia: on entry the caller's moment of inertia of the
! surface convection zone (the whole star's, for a solid body), which
! the solid-body spin-down uses as it is; set anew where the zone's is
! summed here
      double precision, intent(inout) :: cz_moment_of_inertia
      double precision, intent(in) :: cz_mass_bottom
      double precision, intent(in) :: cz_mass_top
      double precision, intent(inout) :: omega_surface
      logical, intent(in) :: surface_cz_active
      double precision, intent(inout) :: mixing_diffusion_coeff(json)
      double precision, intent(inout) :: total_velocity(json)
      logical, intent(inout) :: diffusion_solve_ok
      integer, intent(out) :: ierr
! D_omega at the interfaces (diffusion_velocity_scales); j/m after the
! previous iterate (check_angular_momentum's reference); omega at the
! sub-step's start
      double precision :: am_diffusion_coeff(json), specific_angular_momentum_prev(json), &
           omega_start(json)
! the first and last zones' own convective flags (see
! flag_end_zones_convective)
      logical :: first_zone_convective, last_zone_convective
! some interface has U > 0
      logical :: transport_active
      integer :: iteration
      type(surface_wind_state) :: wind
! the region bounds the final mixing rescans for itself
      integer :: final_zone_begin, final_zone_end

      ierr = 0
! A sub-step that is stable from the start returns before any check
! sets redo_flag, so it starts out false (otherwise a previous redo's
! .true. would survive and rewind the clock again).
      redo_flag = .false.
      call refuse_advection_diffusion(ierr)
      if (ierr /= 0) return

! Disk locking (MHP 9/94): while it is engaged, a transport region that
! reaches the surface holds its surface omega fixed, and an uncoupled
! surface convection zone feels no wind.
      wind%cz_active = surface_cz_active
      wind%disk_locked = disk_locking_engaged(s%disk_gate_age_gyr)
      first_zone_convective = sub_model%am_transport_convective_flag(1)
      last_zone_convective = sub_model%am_transport_convective_flag(s%nz)
      call set_up_interfaces(s, sub_model, sub_physics, stab, zone_min, zone_max, ierr)
      if (ierr /= 0) return
! The velocities of the unperturbed model, the first iteration's. The
! Zahn coupling factor star%vfc is zero unless compute_velocities sets
! it.
      s%vfc(1:s%nz) = 0.0d0
      iteration = 1
      call compute_velocities(s, sub_model, sub_physics, stab, circ_scr, &
           zone_min, zone_max, iteration, total_velocity, transport_active)
! The state at the start of the sub-step: j/m (check_angular_momentum's
! first previous iterate), omega, and the surface omega for the wind
! torque. No iterate has been accepted yet.
      specific_angular_momentum_prev(1:s%nz) = specific_angular_momentum_saved(1:s%nz)
      omega_start(1:s%nz) = s%omega_mid(1:s%nz)
      wind%omega_start = sub_model%omega_start(s%nz)
      diffusion_solve_ok = .false.

! The iterations. The first uses the velocities of the unperturbed
! model. Each later one recomputes them from omega averaged between the
! sub-step's start and the previous iteration's result
! (restart_iteration), averaged in turn with the previous velocities,
! for a more accurate run of diffusion coefficients.
      do iteration = 1,rotation_settings%max_diffusion_iters
         omega_surface = s%omega_mid(s%nz)
         if (iteration > 1) call restart_iteration(s, sub_model, sub_physics, stab, &
              circ_scr, zone_min, zone_max, iteration, omega_start, &
              specific_angular_momentum_saved, total_velocity, transport_active)
         call spin_down_uncoupled_surface(s, sub_model, wind, zone_max, &
              total_velocity(zone_max) == 0.0d0, sub_timestep, cz_mass_bottom, &
              cz_mass_top, omega_surface, cz_moment_of_inertia, ierr)
         if (ierr /= 0) return
! Stable against every mechanism (no U > 0, legacy LDO false): nothing
! to transport.
         if (.not. transport_active) return
         call flag_end_zones_convective(sub_model, s%nz)
! The diffusion coefficients at every interface from the velocities:
! D_omega (legacy COD) and D_mix; on the star, they are the profile
! columns D_omega and D_mix_rot (the step's last sub-step).
         call set_diffusion_coefficients(s, stab, am_diffusion_coeff, mixing_diffusion_coeff)
         call transport_unstable_regions(s, sub_model, sub_physics, stab, chi_grid, wind, &
              zone_min, zone_max, iteration, sub_timestep, total_velocity, &
              am_diffusion_coeff, mixing_diffusion_coeff, omega_start, &
              specific_angular_momentum_saved, cz_moment_of_inertia, omega_surface, &
              diffusion_solve_ok, ierr)
         if (ierr /= 0) return
         call restore_end_zone_flags(sub_model, s%nz, first_zone_convective, &
              last_zone_convective)
         call check_iterate(s, sub_model, sub_physics, iteration, &
              specific_angular_momentum_prev, specific_angular_momentum_saved, &
              sub_timestep, cut_count, diffusion_solve_ok, redo_flag, ierr)
         if (ierr /= 0) return
! Converged.
         if (diffusion_solve_ok) exit
! A check halved the sub-step: the caller redoes it.
         if (redo_flag) return
      end do

! Converged, or at the iteration limit, where the last iterate is
! accepted (no redo, no error; check_angular_momentum only prints).
! diffusion_solve_ok = true also selects diffuse_composition_driver's
! final mode, which rescans the interfaces and builds each region's
! grid itself. Then the species the iterations leave alone (slots 5 to
! 11, or to 15 with the extended composition) mix once across every
! unstable region.
      diffusion_solve_ok = .true.
      call flag_end_zones_convective(sub_model, s%nz)
      call mix_composition(s, sub_model, sub_physics, chi_grid, sub_timestep, &
           mixing_diffusion_coeff, total_velocity, zone_min, zone_max, &
           final_zone_begin, final_zone_end, diffusion_solve_ok, &
           last_mu_species + 1, num_tracked_species(), ierr)
      if (ierr /= 0) return
! Those species get the checks the iterations give slots 1-4: an
! anomalous abundance halves the sub-step for the caller to redo
! (docs/bug_sweep.md O14).
      call check_mixed_species(s%xa, last_mu_species + 1, num_tracked_species(), s%nz, &
           sub_timestep, cut_count, redo_flag, ierr)
      if (ierr /= 0) return
      call restore_end_zone_flags(sub_model, s%nz, first_zone_convective, &
           last_zone_convective)
end subroutine apply_secular_transport

! ---------------------------------------------------------------
! The advection-diffusion transport (use_diffusion_advection_transport,
! legacy LDIFAD) is commented out (2026-10-12, the author's decision;
! see am_transport_coefficients): it never worked, as the equal-grid
! coefficients it read (eq_velocity_coeff*, geometric_factor,
! am_2nd/3rd_deriv_coeff, dchi_dr_edge, mixing_geometric_factor,
! equatorial_radius) were set nowhere, nor were their F77 originals.
! The flag acts only with rotation and instability transport on, so the
! many non-rotating decks that set it run unaffected; here, where it
! would act, the run stops (ierr = 1).
subroutine refuse_advection_diffusion(ierr)
      integer, intent(inout) :: ierr
      if (rotation_settings%use_diffusion_advection_transport) then
         write(*,'(a)') 'apply_secular_transport: use_diffusion_advection_transport' &
              //' (LDIFAD) is not available -- the advection-diffusion' &
              //' transport is commented out. Set it .false.'
         ierr = 1
      end if
end subroutine refuse_advection_diffusion

! ---------------------------------------------------------------
! The diffusion's boundary conditions treat the first and last zones as
! convective, a crude central and surface boundary condition. The
! caller keeps their own flags and gives them back after each pass
! (restore_end_zone_flags).
subroutine flag_end_zones_convective(sub_model, num_zones)
      type(substep_model_state), intent(inout) :: sub_model
      integer, intent(in) :: num_zones
      sub_model%am_transport_convective_flag(1) = .true.
      sub_model%am_transport_convective_flag(num_zones) = .true.
end subroutine flag_end_zones_convective

subroutine restore_end_zone_flags(sub_model, num_zones, first_zone_convective, &
     last_zone_convective)
      type(substep_model_state), intent(inout) :: sub_model
      integer, intent(in) :: num_zones
      logical, intent(in) :: first_zone_convective, last_zone_convective
      sub_model%am_transport_convective_flag(1) = first_zone_convective
      sub_model%am_transport_convective_flag(num_zones) = last_zone_convective
end subroutine restore_end_zone_flags

! ---------------------------------------------------------------
! Once per sub-step, what every iteration shares (into stab): the
! structure at the interfaces, the parts of the velocities that depend
! only on it, the zones' unlogged radii and the dynamical-shear limit on
! |d omega/d ln r| (rotation_stability_setup); d omega/d ln r of the
! sub-step's omega; and the quadrupole moment of the distorted star, for
! the correction of the Eddington-Sweet circulation (compute_quadrupole,
! which sets ierr on a zero pivot).
subroutine set_up_interfaces(s, sub_model, sub_physics, stab, zone_min, zone_max, ierr)
      use stability_coefficients, only: rotation_stability_setup
      use circulation, only: compute_quadrupole
      type(star_info), intent(in) :: s
      type(substep_model_state), intent(in) :: sub_model
      type(substep_physics_state), intent(in) :: sub_physics
      type(stability_state), intent(inout) :: stab
      integer, intent(in) :: zone_min, zone_max
      integer, intent(out) :: ierr
      call rotation_stability_setup(stab, sub_physics, &
           sub_model%log_density, sub_model%mean_gravity, sub_model%luminosity_lsun, &
           sub_model%log_pressure, sub_model%log_radius, s%m, sub_model%log_temperature, &
           zone_min, zone_max, s%nz)
      call compute_omega_gradients(zone_min, zone_max, stab%radius_interface, &
           sub_model%log_radius, s%omega_mid, stab%domega_dlnr)
      call compute_quadrupole(stab, &
           sub_model%log_density, sub_model%mean_gravity, s%omega_mid, s%nz, ierr)
end subroutine set_up_interfaces

! ---------------------------------------------------------------
! The circulation and instability velocities of this iteration
! (total_velocity, U; each mechanism's on the star) and whether any
! transport acts (transport_active), from the surface luminosity (erg/s)
! and stab. Under use_zahn_coupling_factor (legacy LVFC), also the Zahn
! (1992) coupling factor star%vfc, by which diffusion_velocity_scales
! scales D_mix (gated on LVFC as in the F77 seculr.f, not on LDIFAD).
subroutine compute_velocities(s, sub_model, sub_physics, stab, circ_scr, &
     zone_min, zone_max, iteration, total_velocity, transport_active)
      use circulation, only: circulation_velocities, zahn_coupling_factor
      type(star_info), intent(inout) :: s
      type(substep_model_state), intent(in) :: sub_model
      type(substep_physics_state), intent(in) :: sub_physics
      type(stability_state), intent(inout) :: stab
      type(circulation_velocity_state), intent(inout) :: circ_scr
      integer, intent(in) :: zone_min, zone_max, iteration
      double precision, intent(out) :: total_velocity(json)
      logical, intent(out) :: transport_active
      call circulation_velocities(stab, sub_physics, circ_scr, &
           sub_model%log_radius, zone_min, zone_max, iteration, &
           sub_model%am_transport_convective_flag, s%nz, s%omega_mid, transport_active, &
           total_velocity, s%solar_luminosity_cgs*sub_model%luminosity_lsun(s%nz), &
           sub_model%log_pressure, s%mechanism_velocity)
      if (rotation_settings%use_zahn_coupling_factor) &
           call zahn_coupling_factor(stab, sub_model%log_density, stab%radius, &
                total_velocity, zone_min, zone_max, s%omega_mid, s%vfc)
end subroutine compute_velocities

! ---------------------------------------------------------------
! The start of the second and later iterations. omega is averaged with
! its value at the start of the sub-step,
!   omega_mid = (omega_mid + omega_start)/2,
! and its gradient and the velocities are recomputed from it. Then j/m
! goes back to its start, composition slots 1-4 to
! sub_model%composition_snapshot, and the convection zones are remixed
! to their proper depth (MHP 10/91): otherwise a deepening convection
! zone takes the local abundance at its old base.
subroutine restart_iteration(s, sub_model, sub_physics, stab, circ_scr, &
     zone_min, zone_max, iteration, omega_start, specific_angular_momentum_saved, &
     total_velocity, transport_active)
      use convection_lib, only: homogenize_convection_zones
      type(star_info), intent(inout) :: s
      type(substep_model_state), intent(inout) :: sub_model
      type(substep_physics_state), intent(in) :: sub_physics
      type(stability_state), intent(inout) :: stab
      type(circulation_velocity_state), intent(inout) :: circ_scr
      integer, intent(in) :: zone_min, zone_max, iteration
      double precision, intent(in) :: omega_start(json), &
           specific_angular_momentum_saved(json)
      double precision, intent(out) :: total_velocity(json)
      logical, intent(out) :: transport_active
      integer :: i, j
      do i = 1,s%nz
         s%omega_mid(i) = 0.5D0*(s%omega_mid(i)+omega_start(i))
      end do
      call compute_omega_gradients(zone_min, zone_max, stab%radius_interface, &
           sub_model%log_radius, s%omega_mid, stab%domega_dlnr)
      call compute_velocities(s, sub_model, sub_physics, stab, circ_scr, &
           zone_min, zone_max, iteration, total_velocity, transport_active)
      do i = 1,s%nz
         s%j_rot(i) = specific_angular_momentum_saved(i)
         do j = 1,last_mu_species
            s%xa(j,i) = sub_model%composition_snapshot(j,i)
         end do
      end do
      call homogenize_convection_zones(s%xa, s%dm, sub_model%am_transport_convective_flag, &
           s%nz)
end subroutine restart_iteration

! ---------------------------------------------------------------
! The wind's spin-down (brake_convection_zone, which removes the lost
! j/m from every zone it brakes):
! 1. of a surface convection zone the diffusion does not couple to the
!    interior (cz_uncoupled: U = 0 at its base, interface zone_max),
!    when disk locking is not engaged: zones zone_max..nz, whose
!    moment of inertia becomes cz_moment_of_inertia;
! 2. otherwise, of the whole star when it rotates as a solid body
!    (solid_body_mode_flag = 1 acts like force_solid_body_rotation),
!    with the caller's cz_moment_of_inertia; omega and the shape are
!    then recomputed. This branch does not test disk locking.
subroutine spin_down_uncoupled_surface(s, sub_model, wind, zone_max, cz_uncoupled, &
     sub_timestep, cz_mass_bottom, cz_mass_top, omega_surface, cz_moment_of_inertia, ierr)
      use magnetic_braking, only: brake_convection_zone
      use omega_profile, only: solid_body_omega
      use physical_constants, only: ln10
      use math_lib
      type(star_info), intent(inout) :: s
      type(substep_model_state), intent(inout) :: sub_model
      type(surface_wind_state), intent(in) :: wind
      integer, intent(in) :: zone_max
      logical, intent(in) :: cz_uncoupled
      double precision, intent(in) :: sub_timestep, cz_mass_bottom, cz_mass_top
      double precision, intent(inout) :: omega_surface, cz_moment_of_inertia
      integer, intent(out) :: ierr
      integer :: i
      ierr = 0
      if (wind%cz_active .and. cz_uncoupled .and. .not. wind%disk_locked) then
         cz_moment_of_inertia = 0.0D0
         do i = zone_max,s%nz
            cz_moment_of_inertia = cz_moment_of_inertia + sub_model%moment_of_inertia(i)
         end do
         call brake_convection_zone(s%log_L, sub_timestep, cz_mass_bottom, &
              cz_mass_top, zone_max, s%nz, rotation_settings%use_wind_torque, omega_surface, &
              s%star_mass, s%log_Teff, cz_moment_of_inertia, &
              s%j_rot, wind_now(s), ierr)
      else if (wind%cz_active .and. solid_body_everywhere()) then
! The torque on the whole star: from mass 0 to the total mass, zones
! 1..nz; then one omega for the whole star from its new j.
         call brake_convection_zone(s%log_L, sub_timestep, 0.0D0, &
              exp(ln10*s%log_total_mass), 1, s%nz, rotation_settings%use_wind_torque, &
              omega_surface, s%star_mass, s%log_Teff, cz_moment_of_inertia, &
              s%j_rot, wind_now(s), ierr)
         if (ierr /= 0) return
         call solid_body_omega(sub_model%log_density, s%j_rot, sub_model%log_radius, &
              s%log_mass, s%dm, 1, s%nz, &
              sub_model%eta_squared, sub_model%moment_of_inertia, s%omega_mid, &
              sub_model%di_domega, sub_model%mean_radius, s%nz, s%disk_gate_age_gyr)
      endif
end subroutine spin_down_uncoupled_surface

! ---------------------------------------------------------------
! The diffusion coefficients at every interface from the velocities
! (diffusion_velocity_scales): D_omega (am_diffusion_coeff, legacy COD)
! and D_mix (mixing_diffusion_coeff), copied to the star for the
! profile.
subroutine set_diffusion_coefficients(s, stab, am_diffusion_coeff, mixing_diffusion_coeff)
      use circulation, only: diffusion_velocity_scales
      type(star_info), intent(inout) :: s
      type(stability_state), intent(in) :: stab
      double precision, intent(out) :: am_diffusion_coeff(json)
      double precision, intent(inout) :: mixing_diffusion_coeff(json)
! the interface radii, which diffusion_velocity_scales also returns
      double precision :: radius_mid(json)
      integer :: i
      call diffusion_velocity_scales(stab%radius, s%nz, radius_mid, &
           am_diffusion_coeff, mixing_diffusion_coeff, s%mechanism_velocity, s%vfc)
      do i = 1, s%nz
         s%am_diffusion_coeff(i) = am_diffusion_coeff(i)
         s%mixing_diffusion_coeff(i) = mixing_diffusion_coeff(i)
      end do
end subroutine set_diffusion_coefficients

! ---------------------------------------------------------------
! Each unstable region (a run of interfaces with U > 0,
! next_transport_region) is solved separately, scanning outward from
! zone_min. On the region's grid, equally spaced in chi, a normalized
! combination of the variables the Henyey mesh is placed by,
!   chi = log M/dm + L/(L_tot dl) - log P/dp
! (am_transport_grid; chi_grid%eq_chi_spacing its spacing): the
! angular momentum transport (diffuse_region_angular_momentum, with the
! wind's terms for its outermost point from outer_point_wind_terms when
! the region reaches the surface convection zone; none under solid-body
! rotation, force_solid_body_rotation or solid_body_mode_flag = 1, MHP
! 3/09), then the mixing of composition slots 1-4 (X, Y, Z, He3), the
! species that set the gradient of the mean molecular weight. A region
! of a single unstable interface (chi_grid%single_interface, legacy
! LDUM2) is skipped.
subroutine transport_unstable_regions(s, sub_model, sub_physics, stab, chi_grid, wind, &
     zone_min, zone_max, iteration, sub_timestep, total_velocity, &
     am_diffusion_coeff, mixing_diffusion_coeff, omega_start, &
     specific_angular_momentum_saved, cz_moment_of_inertia, omega_surface, &
     diffusion_solve_ok, ierr)
      use transport_regions, only: next_transport_region
      use am_grid, only: am_transport_grid
      use magnetic_braking, only: outer_point_wind_terms
      use region_transport, only: diffuse_region_angular_momentum
      type(star_info), intent(inout) :: s
      type(substep_model_state), intent(in) :: sub_model
      type(substep_physics_state), intent(in) :: sub_physics
      type(stability_state), intent(in) :: stab
      type(transport_grid_state), intent(inout) :: chi_grid
      type(surface_wind_state), intent(inout) :: wind
      integer, intent(in) :: zone_min, zone_max, iteration
      double precision, intent(in) :: sub_timestep
      double precision, intent(in) :: total_velocity(json), am_diffusion_coeff(json), &
           mixing_diffusion_coeff(json), omega_start(json), &
           specific_angular_momentum_saved(json)
      double precision, intent(inout) :: cz_moment_of_inertia, omega_surface
      logical, intent(in) :: diffusion_solve_ok
      integer, intent(out) :: ierr
      logical :: in_unstable_region, region_found, fix_omega_at_surface
      integer :: scan_start_zone, zone_begin, zone_end

      ierr = 0
      in_unstable_region = .false.
      scan_start_zone = zone_min
      region_loop: do
         call next_transport_region(total_velocity, zone_max, scan_start_zone, &
              in_unstable_region, zone_begin, zone_end, region_found)
! No further unstable region.
         if (.not. region_found) exit region_loop
         call am_transport_grid(chi_grid, sub_physics, &
              am_diffusion_coeff, mixing_diffusion_coeff, sub_model%log_density, &
              sub_model%moment_of_inertia, specific_angular_momentum_saved, &
              sub_model%luminosity_lsun, sub_model%log_pressure, sub_model%log_radius, &
              stab%radius, s%log_mass, s%m, s%dm, s%log_total_mass, zone_begin, zone_end, &
              sub_model%am_transport_convective_flag, s%nz, omega_start)
         if (chi_grid%single_interface) then
            if (scan_start_zone <= zone_max) then
               cycle region_loop
            else
               exit region_loop
            endif
         endif
         if (.not. solid_body_everywhere()) then
            call outer_point_wind_terms(wind%cz_active .and. zone_end == zone_max, &
                 wind%disk_locked, iteration, s%log_L, sub_timestep, &
                 chi_grid%eq_moment_of_inertia(chi_grid%num_eq_points), s%omega_mid(s%nz), &
                 s%star_mass, s%log_Teff, wind%omega_start, wind_now(s), &
                 wind%surface_implicit, wind%surface_implicit_set, &
                 wind%loss_explicit, wind%loss_implicit, fix_omega_at_surface, &
                 cz_moment_of_inertia, omega_surface, ierr)
            if (ierr /= 0) return
            call diffuse_region_angular_momentum(chi_grid, sub_timestep, &
                 wind%loss_explicit, wind%loss_implicit, fix_omega_at_surface, &
                 s%dm, zone_begin, zone_end, sub_model%am_transport_convective_flag, &
                 s%nz, s%j_rot, ierr)
            if (ierr /= 0) return
         endif
         call mix_composition(s, sub_model, sub_physics, chi_grid, sub_timestep, &
              mixing_diffusion_coeff, total_velocity, zone_min, zone_max, &
              zone_begin, zone_end, diffusion_solve_ok, 1, last_mu_species, ierr)
         if (ierr /= 0) return
! On to the next region, unless the scan has passed zone_max.
         if (scan_start_zone > zone_max) exit region_loop
      end do region_loop
end subroutine transport_unstable_regions

! ---------------------------------------------------------------
! The end of an iteration. check_angular_momentum stops the run on a
! negative j/m, tests convergence against the previous iterate, finds
! omega and the shape from the new j, and smooths any new reversal of
! the omega gradient. Its diffusion_converged argument, the
! commented-out advection-diffusion solve's own convergence, is false
! (docs/bug_sweep.md O28). Unless it asked for a redo,
! check_composition checks the mixed composition (an anomalous
! abundance halves sub_timestep and asks for a redo) and recomputes the
! mean molecular weight. Mixing without angular momentum transport
! (no_am_transport_in_core; diffusion_velocity_scales zeros D_omega)
! accepts the first iteration (MHP 9/93).
subroutine check_iterate(s, sub_model, sub_physics, iteration, &
     specific_angular_momentum_prev, specific_angular_momentum_saved, &
     sub_timestep, cut_count, diffusion_solve_ok, redo_flag, ierr)
      use transport_checks, only: check_angular_momentum, check_composition
      type(star_info), intent(inout) :: s
      type(substep_model_state), intent(inout) :: sub_model
      type(substep_physics_state), intent(inout) :: sub_physics
      integer, intent(in) :: iteration
      double precision, intent(inout) :: specific_angular_momentum_prev(json)
      double precision, intent(in) :: specific_angular_momentum_saved(json)
      double precision, intent(inout) :: sub_timestep
      integer, intent(inout) :: cut_count
      logical, intent(inout) :: diffusion_solve_ok
      logical, intent(out) :: redo_flag
      integer, intent(out) :: ierr
      call check_angular_momentum(s, sub_model%log_density, specific_angular_momentum_prev, &
           specific_angular_momentum_saved, sub_model%log_radius, s%log_mass, s%dm, &
           iteration, sub_model%am_transport_convective_flag, s%nz, sub_timestep, &
           sub_model%eta_squared, sub_model%moment_of_inertia, s%j_rot, cut_count, &
           diffusion_solve_ok, redo_flag, s%omega_mid, sub_model%di_domega, &
           sub_model%mean_radius, .false., ierr)
      if (ierr /= 0) return
      if (.not. redo_flag) call check_composition(s, sub_physics, s%xa, iteration, &
           s%nz, sub_timestep, cut_count, diffusion_solve_ok, redo_flag, ierr)
      if (ierr /= 0) return
      if (rotation_settings%no_am_transport_in_core) diffusion_solve_ok = .true.
end subroutine check_iterate

! ---------------------------------------------------------------
! The rotational mixing of species species_begin..species_end over the
! sub-step (diffuse_composition_driver): across the region zone_begin..
! zone_end, on its grid, during the iterations; across every region,
! which the driver rescans (setting zone_begin and zone_end), once
! diffusion_solve_ok is set.
subroutine mix_composition(s, sub_model, sub_physics, chi_grid, sub_timestep, &
     mixing_diffusion_coeff, total_velocity, zone_min, zone_max, &
     zone_begin, zone_end, diffusion_solve_ok, species_begin, species_end, ierr)
      use rotational_mixing, only: diffuse_composition_driver
      type(star_info), intent(inout) :: s
      type(substep_model_state), intent(in) :: sub_model
      type(substep_physics_state), intent(in) :: sub_physics
      type(transport_grid_state), intent(inout) :: chi_grid
      double precision, intent(in) :: sub_timestep
      double precision, intent(in) :: mixing_diffusion_coeff(json), total_velocity(json)
      integer, intent(in) :: zone_min, zone_max
      integer, intent(inout) :: zone_begin, zone_end
      logical, intent(in) :: diffusion_solve_ok
      integer, intent(in) :: species_begin, species_end
      integer, intent(out) :: ierr
      call diffuse_composition_driver(chi_grid, sub_physics, sub_timestep, mixing_diffusion_coeff, &
           sub_model%log_density, sub_model%luminosity_lsun, &
           sub_model%log_pressure, sub_model%log_radius, s%log_mass, s%m, s%dm, &
           s%log_total_mass, total_velocity, zone_begin, zone_end, zone_max, &
           zone_min, sub_model%am_transport_convective_flag, diffusion_solve_ok, &
           s%nz, s%xa, species_begin, species_end, ierr)
end subroutine mix_composition

! ---------------------------------------------------------------
! The step's state the wind laws read: tau_cz and log P_phot at the
! start and the end of the step, and the fraction of it elapsed.
function wind_now(s) result(conditions)
      type(star_info), intent(in) :: s
      type(wind_conditions) :: conditions
      conditions = wind_conditions( &
           turnover_timescale_start=s%convective_turnover_timescale_old, &
           turnover_timescale_end=s%convective_turnover_timescale, &
           log_pphot_start=s%pphot0, log_pphot_end=s%pphot, &
           step_fraction=s%fracstep)
end function wind_now

! ---------------------------------------------------------------
! d omega/d ln r at each interface i =
! zone_min..zone_max (between zones i-1 and i), centred on the interface
! radius r_c = radius_interface(i) (cm):
!   domega_dlnr(i) = (omega_i - omega_{i-1})/4
!                      * (1/(ln r_i - ln r_c) + 1/(ln r_c - ln r_{i-1})),
! with log_radius = log10 r and omega (rad/s) at the zone centres. This
! is (omega_i - omega_{i-1})/(ln r_i - ln r_{i-1}) when r_c lies midway
! in ln r. Entries outside zone_min..zone_max are left as they are.
pure subroutine compute_omega_gradients(zone_min, zone_max, radius_interface, log_radius, &
     omega, domega_dlnr)
      use physical_constants, only: ln10
      use math_lib
      integer, intent(in) :: zone_min, zone_max
      double precision, intent(in) :: radius_interface(:), log_radius(:), omega(:)
      double precision, intent(inout) :: domega_dlnr(:)
      integer :: i
      double precision :: log_radius_center
! dlnr_weight (legacy DR1): the sum of the two one-sided inverse ln r
! spacings.
      double precision :: dlnr_weight
      do i = zone_min,zone_max
         log_radius_center = log10(radius_interface(i))
         dlnr_weight = 1.0D0/ln10/(log_radius(i)-log_radius_center)+ &
              1.0D0/ln10/(log_radius_center-log_radius(i-1))
         domega_dlnr(i) = 0.25D0*(omega(i)-omega(i-1))*dlnr_weight
      end do
end subroutine compute_omega_gradients

end module secular_transport

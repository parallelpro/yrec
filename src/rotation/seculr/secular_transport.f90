!----------------------------------------------------------------------
! secular_transport
!----------------------------------------------------------------------
! The transport of angular momentum, and the mixing of the composition
! it drives, by the secular rotational instabilities over one rotation
! sub-step (apply_secular_transport, legacy SECULR), and the tridiagonal
! solve of its diffusion equation (solve_angular_momentum_transport).
!
! A few local and dummy quantities have different names in the callees
! this file dispatches to. The name used here is the one most callees
! share (or the more descriptive one); the others are listed here,
! following henyey_iteration's convention:
!  - dynamical_shear_omega_limit (QWRMAX): circulation_velocities calls
!    the same quantity dlnomega_dlnr_max.
!  - radius_unlogged (HRU): circulation_velocities,
!    zahn_coupling_factor, compute_quadrupole and am_transport_grid call
!    it radius; diffusion_velocity_scales calls it radius_mid_prev.
!  - mixing_diffusion_coeff (COD2, dummy arg): diffuse_composition_driver
!    calls it diffusion_coeff.
!  - diffusion_velocity (HV, dummy arg): circulation_velocities calls it
!    total_circulation_velocity.
!  - eq_mixing_diffusion_coeff (ECOD2): diffuse_composition_driver calls it
!    equally_spaced_diffusion_coeff.
!  - eq_angular_momentum (EJ): equal_grid_to_model calls it angular_momentum.
!  - eq_am_diffusion_coeff (ECOD): am_diffusion_coeffs calls it diffusion_coeff.
!  - eq_delta_angular_momentum (DJ): equal_grid_to_model calls it
!    delta_angular_momentum; solve_angular_momentum_transport calls it dj.
!  - eq_mass (EM): diffuse_composition_driver calls it equally_spaced_mass.
!  - diffusion_converged (LOKAD): check_angular_momentum calls it
!    already_converged_flag.
!  - wind_loss_explicit/wind_loss_implicit (WIND1/WIND2):
!    magnetic_braking's compute_solid_body_braking calls the same pair
!    domega_start/domega_end.
! The sub-step's interpolated structure is read from rot_scr
! (rotation_scratch_lib, filled by mid_timestep_model each sub-step),
! and the model's arrays and scalars (log_mass, m, dm, j_rot, xa, nz,
! log_L, log_total_mass, star_mass, log_Teff) from star. The arguments
! carry only per-call data: the sub-step, the angular momentum saved at
! its start, the bounds of the transport region, the surface convection
! zone's terms, and the outputs.
module secular_transport
      use rotation_configuration, only: rotation_settings, disk_locking_engaged
      use transport_checks, only: check_angular_momentum, check_composition
      use am_grid, only: am_transport_grid, equal_grid_to_model
      use rotational_mixing, only: diffuse_composition_driver
      use circulation, only: circulation_velocities, diffusion_velocity_scales, zahn_coupling_factor, compute_quadrupole
      use am_transport_coefficients, only: am_diffusion_coeffs, &
           solve_angular_momentum_transport
      use omega_profile, only: solid_body_omega
      use stability_coefficients, only: rotation_stability_setup
      implicit none
      private :: compute_omega_gradients
contains

! ---------------------------------------------------------------
! One sub-step (sub_timestep, s) of angular momentum transport and
! rotational mixing over the interfaces zone_min..zone_max (interface i
! lies between zones i-1 and i; mid_timestep_model's core_boundary_zone
! and envelope_boundary_zone).
! Method: the circulation and instability velocities U at the
! interfaces (diffusion_velocity: Eddington-Sweet circulation, secular
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
! mixing_diffusion_coeff and diffusion_velocity; star%am_diffusion_coeff
! and star%mixing_diffusion_coeff for the profile. A nonzero ierr stops
! the run.
subroutine apply_secular_transport(s, rot_scr, mid_scr, circ_scr, &
     sub_timestep, specific_angular_momentum_saved, &
     zone_min, zone_max, redo_flag, cut_count, &
     cz_moment_of_inertia, cz_mass_bottom, cz_mass_top, omega_surface, &
     surface_cz_active, mixing_diffusion_coeff, diffusion_velocity, &
     diffusion_solve_ok, ierr)
      use rotation_scratch_lib, only: rotation_diffusion_state, &
           midstep_physics_state, circulation_velocity_state
      use magnetic_braking, only: brake_convection_zone, compute_solid_body_braking, wind_conditions
      use convection_lib, only: homogenize_convection_zones
      use star_info_lib, only: star_info, json
      use physical_constants
      use numerics_lib
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s
      type(rotation_diffusion_state), intent(inout) :: rot_scr
      type(midstep_physics_state), intent(inout) :: mid_scr
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
      double precision, intent(inout) :: diffusion_velocity(json)
      logical, intent(inout) :: diffusion_solve_ok
      integer, intent(out) :: ierr

! --- locals (arrays of size json) ---
      double precision :: dynamical_shear_omega_limit(json), eq_mass(json), &
           am_diffusion_coeff(json), eq_mixing_diffusion_coeff(json), &
           eq_moment_of_inertia(json), eq_angular_momentum(json), &
           eq_am_diffusion_coeff(json), radius_unlogged(json), radius_mid(json), &
           dlnomega_dlnr(json), eq_delta_angular_momentum(json), &
           specific_angular_momentum_prev(json), omega_start(json), &
           eq_omega(json)

! The tridiagonal system of the angular momentum diffusion equation:
! am_diffusion_coeffs fills it and solve_angular_momentum_transport
! solves it (Thomas algorithm). The solve's new omega (its solution) is
! not used here, only eq_delta_angular_momentum and
! sum_delta_angular_momentum, so it goes into an unused local.
! surface_wind_loss_term is the wind's loss of angular momentum from the
! outermost grid point, from am_diffusion_coeffs to the solve.
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), unused_tridia_solution(json), surface_wind_loss_term

! --- other locals ---
      logical :: disk_lock_active
      logical :: lcz_first_zone, lcz_last_zone
      double precision :: total_luminosity
      integer :: i, j
      integer :: iteration
      double precision :: omega_surface_start
      logical :: wind_loss_active
      integer :: solid_body_zone_start, solid_body_zone_end
      logical :: in_unstable_region
      integer :: scan_start_zone
      logical :: unstable_zone_found
      integer :: zone_begin, zone_end
      double precision :: grid_spacing
      logical :: single_interface_flag
      logical :: fix_omega_at_surface
      logical :: diffusion_converged
      double precision :: wind_loss_explicit, wind_loss_implicit
      double precision :: sum_delta_angular_momentum
      double precision :: solid_cz_mass_bottom, solid_cz_mass_top
      integer :: solid_start_zone
      integer :: species_begin, species_end

      ierr = 0
! A sub-step that is stable from the start returns before any check
! sets redo_flag, so it starts out false (otherwise a previous redo's
! .true. would survive and rewind the clock again).
      redo_flag = .false.
      call refuse_advection_diffusion()
      if (ierr /= 0) return

! Disk locking (MHP 9/94): while it is engaged, a transport region that
! reaches the surface holds its surface omega fixed, and an uncoupled
! surface convection zone feels no wind.
      disk_lock_active = .false.
      if(disk_locking_engaged(s%disk_gate_age_gyr)) disk_lock_active = .true.
      call save_end_zone_flags()
      call set_up_interfaces()
      if (ierr /= 0) return
      call compute_unperturbed_velocities()
      call save_substep_start()

! The iterations. The first uses the velocities of the unperturbed
! model. Each later one recomputes them from omega averaged between the
! sub-step's start and the previous iteration's result
! (restart_iteration), averaged in turn with the previous velocities,
! for a more accurate run of diffusion coefficients.
      do iteration = 1,rotation_settings%max_diffusion_iters
         omega_surface = s%omega_mid(s%nz)
         if(iteration.gt.1) call restart_iteration()
         call spin_down_uncoupled_surface()
         if (ierr /= 0) return
! Stable against every mechanism (no U > 0, legacy LDO false): nothing
! to transport.
         if(.not.unstable_zone_found) return
         call flag_end_zones_convective()
         call set_diffusion_coefficients()
         call transport_unstable_regions()
         if (ierr /= 0) return
         call restore_end_zone_flags()
         call check_iterate()
         if (ierr /= 0) return
! Converged.
         if(diffusion_solve_ok)exit
! A check halved the sub-step: the caller redoes it.
         if(redo_flag)return
      end do

! Converged, or at the iteration limit, where the last iterate is
! accepted (no redo, no error; check_angular_momentum only prints).
! diffusion_solve_ok = true also selects diffuse_composition_driver's
! final mode, which rescans the interfaces and builds each region's
! grid itself.
      diffusion_solve_ok = .true.
      call flag_end_zones_convective()
      call mix_remaining_species()
      if (ierr /= 0) return
      call restore_end_zone_flags()

contains

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
subroutine refuse_advection_diffusion()
      if (rotation_settings%use_diffusion_advection_transport) then
         write(*,'(a)') 'apply_secular_transport: use_diffusion_advection_transport' &
              //' (LDIFAD) is not available -- the advection-diffusion' &
              //' transport is commented out. Set it .false.'
         ierr = 1
      end if
end subroutine refuse_advection_diffusion

! ---------------------------------------------------------------
! The diffusion's boundary conditions treat the first and last zones as
! convective, a crude central and surface boundary condition
! (flag_end_zones_convective). Their own flags are kept here and
! restored after each pass (restore_end_zone_flags).
subroutine save_end_zone_flags()
      lcz_first_zone = rot_scr%am_transport_convective_flag_mid(1)
      lcz_last_zone = rot_scr%am_transport_convective_flag_mid(s%nz)
end subroutine save_end_zone_flags

! ---------------------------------------------------------------
! The central and surface zones count as convective in the solve.
subroutine flag_end_zones_convective()
      rot_scr%am_transport_convective_flag_mid(1) = .true.
      rot_scr%am_transport_convective_flag_mid(s%nz) = .true.
end subroutine flag_end_zones_convective

! ---------------------------------------------------------------
! The first and last zones get their own convective flags back.
subroutine restore_end_zone_flags()
      rot_scr%am_transport_convective_flag_mid(1) = lcz_first_zone
      rot_scr%am_transport_convective_flag_mid(s%nz) = lcz_last_zone
end subroutine restore_end_zone_flags

! ---------------------------------------------------------------
! Once per sub-step, what every iteration shares: the structure at the
! interfaces and the parts of the velocities that depend only on it
! (rotation_stability_setup), with the unlogged radii radius_unlogged
! and the dynamical-shear limit on |d omega/d ln r|; the surface L,
! total_luminosity (erg/s); d omega/d ln r of the sub-step's omega; and
! the quadrupole moment of the distorted star, for the correction of
! the Eddington-Sweet circulation (compute_quadrupole, which sets ierr
! on a zero pivot).
subroutine set_up_interfaces()
      call rotation_stability_setup(rot_scr, mid_scr, &
           rot_scr%log_density_mid,rot_scr%mean_gravity_mid,rot_scr%luminosity_lsun_mid,rot_scr%log_pressure_mid, &
           rot_scr%log_radius_mid,s%m,rot_scr%log_temperature_mid,zone_min,zone_max, &
           s%nz,radius_unlogged,dynamical_shear_omega_limit)
      total_luminosity = s%solar_luminosity_cgs*rot_scr%luminosity_lsun_mid(s%nz)
      call compute_omega_gradients(zone_min,zone_max,rot_scr%radius_interface, &
           rot_scr%log_radius_mid,s%omega_mid,dlnomega_dlnr)
      call compute_quadrupole(rot_scr, &
           rot_scr%log_density_mid,rot_scr%mean_gravity_mid,radius_unlogged,s%omega_mid,s%nz,ierr)
end subroutine set_up_interfaces

! ---------------------------------------------------------------
! The stability of the unperturbed model, and its velocities (the first
! iteration's): secular and dynamical shear, Eddington-Sweet
! circulation, the GSF instability. The Zahn coupling factor star%vfc
! is zero unless compute_circulation sets it.
subroutine compute_unperturbed_velocities()
      do i = 1,s%nz
         s%vfc(i) = 0.0D0
      end do
      iteration = 1
      call compute_circulation()
end subroutine compute_unperturbed_velocities

! ---------------------------------------------------------------
! The state at the start of the sub-step: j/m (check_angular_momentum's
! first previous iterate), omega, and the surface omega for the wind
! torque. No iterate has been accepted yet.
subroutine save_substep_start()
      do i =1,s%nz
         specific_angular_momentum_prev(i) = specific_angular_momentum_saved(i)
         omega_start(i) = s%omega_mid(i)
      end do
      omega_surface_start = rot_scr%omega_substep_start(s%nz)
      diffusion_solve_ok = .false.
end subroutine save_substep_start

! ---------------------------------------------------------------
! The circulation and instability velocities of this iteration
! (diffusion_velocity) and whether any transport acts
! (unstable_zone_found). Under use_zahn_coupling_factor (legacy LVFC),
! also the Zahn (1992) coupling factor star%vfc, by which
! diffusion_velocity_scales scales D_mix (gated on LVFC as in the F77
! seculr.f, not on LDIFAD).
subroutine compute_circulation()
      call circulation_velocities(rot_scr, mid_scr, circ_scr, &
           rot_scr%log_radius_mid,radius_unlogged,zone_min,zone_max,iteration, &
           rot_scr%am_transport_convective_flag_mid,s%nz,s%omega_mid,unstable_zone_found, &
           dlnomega_dlnr,dynamical_shear_omega_limit,diffusion_velocity, &
           total_luminosity,rot_scr%log_pressure_mid,s%es_circulation_velocity, &
           s%gsf_circulation_velocity,s%secular_shear_velocity)
      if(rotation_settings%use_zahn_coupling_factor) &
           call zahn_coupling_factor(rot_scr, &
                rot_scr%log_density_mid,radius_unlogged,diffusion_velocity, &
           zone_min,zone_max,s%omega_mid,s%vfc)
end subroutine compute_circulation

! ---------------------------------------------------------------
! The start of the second and later iterations. omega is averaged with
! its value at the start of the sub-step,
!   omega_mid = (omega_mid + omega_start)/2,
! and the velocities are recomputed from it. Then j/m goes back to its
! start, composition slots 1-4 to rot_scr%composition_snapshot, and the
! convection zones are remixed to their proper depth (MHP 10/91):
! otherwise a deepening convection zone takes the local abundance at its
! old base.
subroutine restart_iteration()
      do i = 1,s%nz
         s%omega_mid(i) = 0.5D0*(s%omega_mid(i)+omega_start(i))
      end do
      call compute_omega_gradients(zone_min,zone_max,rot_scr%radius_interface, &
           rot_scr%log_radius_mid,s%omega_mid,dlnomega_dlnr)
      call compute_circulation()
      do i = 1,s%nz
         s%j_rot(i) = specific_angular_momentum_saved(i)
         do j = 1,4
            s%xa(j,i) = rot_scr%composition_snapshot(j,i)
         end do
      end do
      call homogenize_convection_zones(s%xa,s%dm,rot_scr%am_transport_convective_flag_mid, &
           s%nz)
end subroutine restart_iteration

! ---------------------------------------------------------------
! The wind's spin-down (brake_convection_zone, which removes the lost
! j/m from every zone it brakes):
! 1. of a surface convection zone the diffusion does not couple to the
!    interior (U = 0 at its base, interface zone_max), when disk
!    locking is not engaged: zones zone_max..nz;
! 2. otherwise, of the whole star when it rotates as a solid body
!    (solid_body_mode_flag = 1 acts like force_solid_body_rotation),
!    whose omega and shape are then recomputed. This branch does not
!    test disk locking.
subroutine spin_down_uncoupled_surface()
      if(surface_cz_active.and.diffusion_velocity(zone_max).eq.0.0D0 &
           .and. .not.disk_lock_active)then
! The moment of inertia of the surface convection zone.
         cz_moment_of_inertia = 0.0D0
         do i = zone_max,s%nz
            cz_moment_of_inertia = cz_moment_of_inertia + rot_scr%moment_of_inertia_mid(i)
         end do
         wind_loss_active = rotation_settings%use_wind_torque
         call brake_convection_zone(s%log_L,sub_timestep,cz_mass_bottom, &
              cz_mass_top,zone_max,s%nz,wind_loss_active,omega_surface, &
              s%star_mass,s%log_Teff,cz_moment_of_inertia, &
              s%j_rot, wind_now(), ierr)
         if (ierr /= 0) return
      else if(surface_cz_active .and. (rotation_settings%force_solid_body_rotation .or. &
           (rotation_settings%solid_body_mode_flag.eq.1)))then
! The torque on the whole star: from mass 0 to the total mass, zones
! 1..nz.
         wind_loss_active = rotation_settings%use_wind_torque
         solid_cz_mass_bottom = 0.0D0
         solid_cz_mass_top = exp(ln10*s%log_total_mass)
         solid_start_zone = 1
         call brake_convection_zone(s%log_L,sub_timestep,solid_cz_mass_bottom, &
              solid_cz_mass_top,solid_start_zone,s%nz,wind_loss_active, &
              omega_surface,s%star_mass,s%log_Teff,cz_moment_of_inertia, &
              s%j_rot, wind_now(), ierr)
         if (ierr /= 0) return
! One omega for the whole star from its new j.
         solid_body_zone_start = 1
         solid_body_zone_end = s%nz
         call solid_body_omega(rot_scr%log_density_mid,s%j_rot,rot_scr%log_radius_mid, &
              s%log_mass,s%dm,solid_body_zone_start,solid_body_zone_end, &
              rot_scr%eta_squared_mid,rot_scr%moment_of_inertia_mid,s%omega_mid,rot_scr%di_domega_mid,rot_scr%mean_radius_mid,s%nz,s%disk_gate_age_gyr)
      endif
end subroutine spin_down_uncoupled_surface

! ---------------------------------------------------------------
! The diffusion coefficients at every interface from the velocities
! (diffusion_velocity_scales): D_omega (am_diffusion_coeff, legacy COD)
! and D_mix (mixing_diffusion_coeff). The model-grid coefficients are
! also profile columns (D_omega, D_mix_rot), which show the step's last
! sub-step.
subroutine set_diffusion_coefficients()
      call diffusion_velocity_scales(radius_unlogged,s%nz,radius_mid, &
           am_diffusion_coeff,mixing_diffusion_coeff,s%es_circulation_velocity, &
           s%gsf_circulation_velocity,s%secular_shear_velocity,s%vfc)
      do i = 1, s%nz
         s%am_diffusion_coeff(i) = am_diffusion_coeff(i)
         s%mixing_diffusion_coeff(i) = mixing_diffusion_coeff(i)
      end do
end subroutine set_diffusion_coefficients

! ---------------------------------------------------------------
! Each unstable region is solved separately, scanning outward from
! zone_min: on its own grid, the angular momentum transport (none under
! solid-body rotation: force_solid_body_rotation, or
! solid_body_mode_flag = 1; MHP 3/09), then the mixing of composition
! slots 1-4 (X, Y, Z, He3), the species that set the gradient of the
! mean molecular weight; the other species mix after the iterations
! (mix_remaining_species).
subroutine transport_unstable_regions()
      in_unstable_region = .false.
      scan_start_zone = zone_min
      region_loop: do
         call find_unstable_region()
! No further unstable region.
         if(.not.unstable_zone_found) exit region_loop
         call grid_unstable_region()
! A region of a single unstable interface (zone_end - zone_begin <= 1;
! legacy LDUM2) is skipped.
         if(single_interface_flag) then
            if(scan_start_zone.le.zone_max) then
               cycle region_loop
            else
               exit region_loop
            endif
         endif
         if(.not.rotation_settings%force_solid_body_rotation .and. &
              (rotation_settings%solid_body_mode_flag.ne.1))then
            call transport_angular_momentum_in_region()
            if (ierr /= 0) return
         endif
         species_begin = 1
         species_end = 4
         call mix_composition()
         if (ierr /= 0) return
! On to the next region, unless the scan has passed zone_max.
         if(scan_start_zone.gt.zone_max) exit region_loop
      end do region_loop
end subroutine transport_unstable_regions

! ---------------------------------------------------------------
! The next unstable region (a run of interfaces with U > 0), scanning
! from scan_start_zone to zone_max: zone_begin is the zone below its
! first unstable interface, zone_end the zone above its last, and
! scan_start_zone moves past it. unstable_zone_found is false when the
! scan meets no U > 0.
subroutine find_unstable_region()
      unstable_zone_found = .false.
      do j = scan_start_zone,zone_max
         if(diffusion_velocity(j).gt.0.0D0) then
            unstable_zone_found = .true.
            if(.not.in_unstable_region) then
               in_unstable_region = .true.
               zone_begin = j - 1
            endif
         else if(in_unstable_region) then
            zone_end = j - 1
            in_unstable_region = .false.
            scan_start_zone = j + 1
            exit
         endif
      end do
      if (j .gt. zone_max) then
! The scan reached zone_max: a region still open ends there.
      if(in_unstable_region) zone_end = zone_max
      scan_start_zone = zone_max + 1
      end if
end subroutine find_unstable_region

! ---------------------------------------------------------------
! The region zone_begin..zone_end on a grid equally spaced in chi, a
! normalized combination of the variables the Henyey mesh is placed by,
!   chi = log M/dm + L/(L_tot dl) - log P/dp
! (rotation_scratch_lib), so that the angular momentum grid stays close
! to the structure's: D_omega, D_mix, I, J, m and omega on that grid
! (am_transport_grid), with its spacing grid_spacing.
! single_interface_flag marks a region of a single unstable interface.
subroutine grid_unstable_region()
      call am_transport_grid(rot_scr, mid_scr, &
           am_diffusion_coeff,mixing_diffusion_coeff,rot_scr%log_density_mid, &
           rot_scr%moment_of_inertia_mid,specific_angular_momentum_saved,rot_scr%luminosity_lsun_mid, &
           rot_scr%log_pressure_mid,rot_scr%log_radius_mid,radius_unlogged,s%log_mass,s%m, &
           s%dm,s%log_total_mass,zone_begin,zone_end, &
           rot_scr%am_transport_convective_flag_mid,s%nz,omega_start,grid_spacing, &
           eq_am_diffusion_coeff,eq_mixing_diffusion_coeff, &
           eq_moment_of_inertia,eq_angular_momentum,eq_mass,eq_omega, &
           single_interface_flag)
end subroutine grid_unstable_region

! ---------------------------------------------------------------
! Angular momentum transport across one unstable region, on its equally
! spaced grid of n = rot_scr%num_eq_points points: the implicit
! diffusion equation for omega at the end of the sub-step (primed),
!   I_k (omega'_k - omega_k) = (dt/dchi) [D_{k+1} (omega'_{k+1} - omega'_k)
!                                        - D_k (omega'_k - omega'_{k-1})],
! with D_k = eq_am_diffusion_coeff(k) at the edge between points k-1
! and k and no flux through the region's ends (am_diffusion_coeffs).
! When the region reaches the surface convection zone, its outermost
! point also loses the wind's I_n (wind_loss_explicit +
! wind_loss_implicit)/2, or, under disk locking, keeps its omega
! (set_surface_wind_terms). The change of angular momentum is then
! mapped back to the model. (The commented-out advection-diffusion
! transport solved a band system for omega and its derivatives here
! instead.)
subroutine transport_angular_momentum_in_region()
      call set_surface_wind_terms()
      if (ierr /= 0) return
      diffusion_converged = .false.
! The diffusion equation's tridiagonal system.
      call am_diffusion_coeffs(eq_am_diffusion_coeff,grid_spacing,sub_timestep, &
           eq_moment_of_inertia,eq_angular_momentum,eq_omega,rot_scr%num_eq_points, &
           wind_loss_explicit,wind_loss_implicit,fix_omega_at_surface, &
           sub_diag,diag,super_diag,rhs,surface_wind_loss_term)
! Its solve: eq_delta_angular_momentum = I (omega' - omega) at each
! grid point, sum_delta_angular_momentum their sum.
      call solve_angular_momentum_transport(rot_scr%num_eq_points,eq_moment_of_inertia, &
           eq_delta_angular_momentum,sum_delta_angular_momentum, &
           sub_diag,diag,super_diag,rhs,unused_tridia_solution, &
           surface_wind_loss_term, ierr)
      if (ierr /= 0) return
! (commented out: the advection-diffusion branch)
!      else
!         call am_advection_diffusion_coeffs(grid_spacing,sub_timestep,eq_moment_of_inertia, &
!              eq_omega,rot_scr%num_eq_points,wind_loss_explicit,wind_loss_implicit, &
!              eq_delta_angular_momentum,eq_mixing_diffusion_coeff, &
!              sum_delta_angular_momentum,fix_omega_at_surface, &
!              diffusion_converged, ierr)
!         if (ierr /= 0) return
!      endif
! The fractional change dJ/J interpolated back to the region's model
! zones (and any convection zone adjoining it), rescaled if needed so
! that the total change is sum_delta_angular_momentum
! (equal_grid_to_model).
      call equal_grid_to_model(rot_scr, eq_delta_angular_momentum,eq_angular_momentum, &
           s%dm,zone_begin,zone_end,rot_scr%am_transport_convective_flag_mid, &
           s%nz,sum_delta_angular_momentum,s%j_rot)
end subroutine transport_angular_momentum_in_region

! ---------------------------------------------------------------
! The wind's terms for the region's outermost grid point. They act only
! when the region reaches the surface convection zone (zone_end =
! zone_max): under disk locking the point keeps its omega
! (fix_omega_at_surface) and loses nothing; otherwise, with
! use_wind_torque, it loses I_n (wind_loss_explicit +
! wind_loss_implicit)/2. Elsewhere both terms are zero.
subroutine set_surface_wind_terms()
      if(surface_cz_active.and.zone_end.eq.zone_max)then
         if(disk_lock_active)then
            wind_loss_explicit = 0.0D0
            wind_loss_implicit = 0.0D0
            fix_omega_at_surface = .true.
         else
            fix_omega_at_surface = .false.
            if(rotation_settings%use_wind_torque)then
! The drop of the surface omega over the sub-step from the torque at
! its start (wind_loss_explicit) and at its current end
! (wind_loss_implicit, averaged with its previous value from the second
! iteration on), for the outermost grid point's moment of inertia.
               cz_moment_of_inertia = eq_moment_of_inertia(rot_scr%num_eq_points)
               omega_surface = s%omega_mid(s%nz)
               call compute_solid_body_braking(s%log_L,sub_timestep, &
                    cz_moment_of_inertia,iteration,omega_surface, &
                    s%star_mass,s%log_Teff,omega_surface_start, &
                    wind_loss_explicit,wind_loss_implicit, wind_now(), ierr)
               if (ierr /= 0) return
            else
               wind_loss_explicit = 0.0D0
               wind_loss_implicit = 0.0D0
            endif
         endif
      else
         fix_omega_at_surface = .false.
         wind_loss_explicit = 0.0D0
         wind_loss_implicit = 0.0D0
      endif
end subroutine set_surface_wind_terms

! ---------------------------------------------------------------
! The end of an iteration. check_angular_momentum stops the run on a
! negative j/m, tests convergence against the previous iterate, finds
! omega and the shape from the new j, and smooths any new reversal of
! the omega gradient. Unless it asked for a redo, check_composition
! checks the mixed composition (an anomalous abundance halves
! sub_timestep and asks for a redo) and recomputes the mean molecular
! weight. Mixing without angular momentum transport
! (no_am_transport_in_core; diffusion_velocity_scales zeros D_omega)
! accepts the first iteration (MHP 9/93).
subroutine check_iterate()
      call check_angular_momentum(s, rot_scr%log_density_mid,specific_angular_momentum_prev, &
           specific_angular_momentum_saved,rot_scr%log_radius_mid,s%log_mass,s%dm, &
           iteration, &
           rot_scr%am_transport_convective_flag_mid,s%nz,sub_timestep,rot_scr%eta_squared_mid, &
           rot_scr%moment_of_inertia_mid,s%j_rot,cut_count, &
           diffusion_solve_ok,redo_flag,s%omega_mid,rot_scr%di_domega_mid,rot_scr%mean_radius_mid, &
           diffusion_converged, ierr)
      if (ierr /= 0) return
      if(.not.redo_flag)call check_composition(s, mid_scr, s%xa,iteration, &
           s%nz,sub_timestep,cut_count, &
           diffusion_solve_ok,redo_flag, ierr)
      if (ierr /= 0) return
      if(rotation_settings%no_am_transport_in_core)diffusion_solve_ok = .true.
end subroutine check_iterate

! ---------------------------------------------------------------
! The species the iterations leave alone (slots 5 to 11, or to 15 with
! the extended composition), mixed once across every unstable region.
subroutine mix_remaining_species()
      species_begin = 5
      if(rotation_settings%use_extended_composition)then
         species_end = 15
      else
         species_end = 11
      endif
      call mix_composition()
end subroutine mix_remaining_species

! ---------------------------------------------------------------
! The rotational mixing of species species_begin..species_end over the
! sub-step (diffuse_composition_driver): across the current region, on
! its grid, during the iterations; across every region once
! diffusion_solve_ok is set.
subroutine mix_composition()
      call diffuse_composition_driver(rot_scr, mid_scr, sub_timestep,mixing_diffusion_coeff, &
           eq_mixing_diffusion_coeff,eq_mass,rot_scr%log_density_mid,rot_scr%luminosity_lsun_mid, &
           rot_scr%log_pressure_mid,rot_scr%log_radius_mid,s%log_mass,s%m,s%dm, &
           s%log_total_mass,diffusion_velocity,zone_begin,zone_end,zone_max, &
           zone_min,rot_scr%am_transport_convective_flag_mid,diffusion_solve_ok, &
           s%nz,s%xa,species_begin,species_end, ierr)
end subroutine mix_composition


! ---------------------------------------------------------------
! The step's state the wind laws read: tau_cz and log P_phot at the
! start and the end of the step, and the fraction of it elapsed.
function wind_now() result(conditions)
      type(wind_conditions) :: conditions
      conditions = wind_conditions( &
           turnover_timescale_start=s%convective_turnover_timescale_old, &
           turnover_timescale_end=s%convective_turnover_timescale, &
           log_pphot_start=s%pphot0, log_pphot_end=s%pphot, &
           step_fraction=s%fracstep)
end function wind_now

end subroutine apply_secular_transport

! ---------------------------------------------------------------
! d omega/d ln r (not d ln omega/d ln r) at each interface i =
! zone_min..zone_max (between zones i-1 and i), centred on the interface
! radius r_c = radius_interface(i) (cm):
!   dlnomega_dlnr(i) = (omega_i - omega_{i-1})/4
!                      * (1/(ln r_i - ln r_c) + 1/(ln r_c - ln r_{i-1})),
! with log_radius = log10 r and omega (rad/s) at the zone centres. This
! is (omega_i - omega_{i-1})/(ln r_i - ln r_{i-1}) when r_c lies midway
! in ln r. Entries outside zone_min..zone_max are left as they are.
pure subroutine compute_omega_gradients(zone_min, zone_max, radius_interface, log_radius, &
     omega, dlnomega_dlnr)
      use physical_constants, only: ln10
      use math_lib
      integer, intent(in) :: zone_min, zone_max
      double precision, intent(in) :: radius_interface(:), log_radius(:), omega(:)
      double precision, intent(inout) :: dlnomega_dlnr(:)
      integer :: i
      double precision :: log_radius_center
! dlnr_weight (legacy DR1): the sum of the two one-sided inverse ln r
! spacings.
      double precision :: dlnr_weight
      do i = zone_min,zone_max
         log_radius_center = log10(radius_interface(i))
         dlnr_weight = 1.0D0/ln10/(log_radius(i)-log_radius_center)+ &
              1.0D0/ln10/(log_radius_center-log_radius(i-1))
         dlnomega_dlnr(i) = 0.25D0*(omega(i)-omega(i-1))*dlnr_weight
      end do
end subroutine compute_omega_gradients

end module secular_transport

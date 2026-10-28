!----------------------------------------------------------------------
! evolve_rotation
!----------------------------------------------------------------------
! The rotation of the star over a timestep: evolve_angular_momentum
! sub-steps the transport of angular momentum and the mixing it drives,
! with the wind torque, and updates omega.
module evolve_rotation
      use shell_physics, only: compute_shell_physics
      use rotation_configuration, only: rotation_settings, disk_locking_engaged
      implicit none
      private
      public :: evolve_angular_momentum

contains

! ---------------------------------------------------------------
! 11/91 MHP revised subroutine (replaces old FINDW).
! GETW evolves the angular momentum distribution in time.
! It enforces local conservation of angular momentum in radiative regions,
! and enforces the user-specified rotation law in convective regions.
! Both meridional circulation and differential rotation with depth can
! cause the radial transport of angular momentum and material.
! GETW solves for the transport of angular momentum and the associated mixing
! if desired (LINSTB=T), and accounts for angular momentum loss caused by a
! magnetic stellar wind if desired (LJDOT0=T).
! If instabilities are being treated, the burning of the light elements
! lithium and beryllium are considered here if desired (LEXCOM=T).
! For a discussion of the overall technique used see Pinsonneault, Kawaler,
! Sofia, and Demarque (1989), Ap.J. vol. 338, p.424.
! 11/91 JENV0 added to call.
!
! NOTE: this routine is only reached inside an IF(LROT) gate; within
! it, whether secular_transport (the rotation-diffusion instability
! solver) is called is gated by instability_transport_active (LINSTB).
subroutine evolve_angular_momentum(s, full_timestep, max_domega_step, wind_loss_active, &
     envelope_boundary_zone_prev, ierr)
      use am_regions, only: am_convective_regions
      use omega_profile, only: enforce_rotation_profile, omega_from_j
      use rotation_scratch_lib, only: rotation_diffusion_state, &
           midstep_physics_state, circulation_velocity_state
      use mid_timestep_model
      use secular_transport
      use magnetic_braking, only: brake_convection_zone, wind_conditions
      use convection_lib, only: find_convection_zones, homogenize_convection_zones
      use star_info_lib, only: star_info, json
      use physical_constants
      use burn_step, only: store_end_li_be_rates, burn_li_be
      use convection_geometry, only: pressure_scale_height
      use math_lib
      implicit none
! di_domega of omega_from_j and enforce_rotation_profile, which the
! step does not use
      double precision :: qiw(json)
      type(star_info), intent(inout) :: s
! The step's rotation scratch (rotation_scratch_lib): the solver's
! working arrays and the structure interpolated to the sub-step
! (rot_scr), the physics interpolated to it (mid_scr) and the
! circulation's previous iterates (circ_scr). On the heap for the step,
! passed down as arguments.
      type(rotation_diffusion_state), allocatable :: rot_scr
      type(midstep_physics_state), allocatable :: mid_scr
      type(circulation_velocity_state), allocatable :: circ_scr

      double precision, intent(in) :: full_timestep
      double precision, intent(inout) :: max_domega_step
      logical, intent(inout) :: wind_loss_active
      integer, intent(inout) :: envelope_boundary_zone_prev
      integer, intent(out) :: ierr

! am_transport_convective_flag (originally LCZ): convective-for-AM-
! transport-purposes flag, set by am_convective_regions below.
      logical :: am_transport_convective_flag(json)
      double precision :: specific_angular_momentum_saved(json)
! The "_mid" midpoint-in-time structure arrays (computed by
! mid_timestep_model each diffusion sub-step) live in rot_scr
! (rotation_scratch_lib) -- see rot_scr%*_mid uses below.
! cod2 and diffusion_velocity (HV): produced by secular_transport for
! the commented-out burn_settle_mix, their only consumer.
      double precision :: cod2(json), diffusion_velocity(json)
      integer :: convective_zone_bounds(12,2), &
           convective_zone_bounds_burn(12,2), radiative_zone_bounds(13,2)
      double precision :: cz_mass_bottom, cz_mass_top
      integer :: num_species_tracked
      logical :: disk_lock_engaged
      double precision :: omega_surface
      double precision :: moment_of_inertia_cz
! envelope_boundary_zone/core_boundary_zone (originally IMAX/IMIN) are
! reused across both the early (LINSTB=.false.) return branch, where
! envelope_boundary_zone is set from convective_zone_bounds directly,
! and the main diffusion loop, where MIDMOD sets both from the
! midpoint-in-time structure -- the two branches never execute in the
! same call (the early branch exits via GOTO 9999), so this matches
! mid_timestep_model.f90's own naming for the same physical quantities throughout.
      integer :: envelope_boundary_zone, core_boundary_zone
      integer :: zone_index, species_index
      integer :: redo_count
      integer :: num_diffusion_steps, num_wind_diffusion_steps
      double precision :: sub_timestep
      double precision :: elapsed_substep_time
      logical :: first_call
      double precision :: delts_test
      double precision :: fx
      logical :: surface_cz_active, fully_convective_flag
      integer :: num_radiative_zones, num_convective_zones
      logical :: redo_needed_flag
      logical :: diffusion_solve_ok
      integer :: core_boundary_zone_cur, envelope_boundary_zone_cur
      integer :: num_convective_zones_burn
      integer :: iend
! (commented out with the advection-diffusion transport's theta term)
!      double precision :: omega_avg, domega_dr, delta_radius_step

      ierr = 0
      num_species_tracked = 11
      if(rotation_settings%use_extended_composition) num_species_tracked = 15
! The regions for angular momentum transport: the convection zones with
! their overshoot (am_convective_regions). With instability transport
! on, the overshoot region needs at least one shell, to couple the
! radiative and convective regions in the diffusion equations.
      call am_convective_regions(s%xa,s%logRho,s%logP,s%logR,s%log_mass, &
           s%logT,s%convective_flag,s%nz, &
           am_transport_convective_flag,radiative_zone_bounds, &
           convective_zone_bounds,num_radiative_zones,num_convective_zones, s%iov1, s%iov2, ierr)
      if (ierr /= 0) return
      disk_lock_engaged = .false.
      if(disk_locking_engaged(s%disk_gate_age_gyr)) disk_lock_engaged = .true.
      if(.not.rotation_settings%instability_transport_active)then
         call spin_down_without_transport()
         return
      end if
! The rotation curve of the angular momentum distribution, convection
! zones as solid bodies.
      call omega_from_j(s%logRho,s%j_rot,s%logR,s%log_mass, &
           s%dm,am_transport_convective_flag,s%nz,s%eta_squared, &
           s%i_rot,s%omega,qiw,s%mean_radius,s%disk_gate_age_gyr)
      if (.not. full_timestep.gt.0.0D0) return
      call count_substeps()
! The quantities both kinds of instability need: the temperature
! gradients, opacities, viscosities, thermal diffusivity and heat
! capacity.
      call compute_shell_physics(s, s%fp_rot,s%ft_rot,s%xa,s%logRho,s%luminosity_lsun,s%logP, &
           s%logR,s%log_mass,s%logT,s%convective_flag,s%nz,s%log_Teff, ierr)
      if (ierr /= 0) return
      allocate(rot_scr, mid_scr, circ_scr)
! The sub-steps: the structure is interpolated between the start and the
! end of the step; a sub-step whose transport asks for a redo is
! repeated from its start, halved.
      elapsed_substep_time = 0.0D0
      first_call = .true.
      do
      redo_count = 0
      retry: do
      call advance_substep_clock()
      call record_substep_start_omega()
      fx = elapsed_substep_time/full_timestep
      s%fracstep = fx
      call advance_rotation_substep(s, rot_scr, mid_scr, &
           full_timestep, sub_timestep, fx, first_call, &
           moment_of_inertia_cz, cz_mass_bottom, cz_mass_top, &
           core_boundary_zone, envelope_boundary_zone, &
           fully_convective_flag, surface_cz_active, &
           radiative_zone_bounds, convective_zone_bounds, &
           num_radiative_zones, num_convective_zones, ierr)
      if (ierr /= 0) return
      if(.not.disk_lock_engaged .and. wind_loss_active .and. surface_cz_active) then
         call brake_surface_convection_zone()
         if (ierr /= 0) return
      endif
      if(rotation_settings%instability_transport_active.and..not.fully_convective_flag) then
         call apply_instability_transport()
         if (ierr /= 0) return
         if (redo_needed_flag) cycle retry
      endif
! the start of the next sub-step (burning is accounted for in fully
! convective models and those without instabilities too)
      do zone_index = 1,s%nz
         do species_index = 1,11
            s%xa_start(species_index,zone_index) = s%xa(species_index,zone_index)
         end do
      end do
      first_call = .false.
      if(rotation_settings%use_extended_composition)then
         call burn_light_elements_over_substep()
         if (ierr /= 0) return
      endif
! The commented-out Bulirsch-Stoer path (MHP 6/00) extrapolated mixing
! plus burning to zero timestep here, when instability transport was on
! and the model not fully convective; while it was on, the xa_start
! updates above were skipped:
!     if(rotation_settings%instability_transport_active .and. .not.fully_convective_flag)then
!       call burn_settle_mix(cod2,sub_timestep,star%xa,rot_scr%log_density_mid,rot_scr%luminosity_lsun_mid, &
!            rot_scr%log_pressure_mid,rot_scr%log_radius_mid,star%log_mass,star%m,star%dm,star%log_total_mass, &
!            rot_scr%log_temperature_mid,diffusion_velocity,envelope_boundary_zone,core_boundary_zone, &
!            envelope_boundary_zone_prev,envelope_boundary_zone_cur,diffusion_solve_ok, &
!            rot_scr%am_transport_convective_flag_mid,star%nz,radiative_zone_bounds, &
!            convective_zone_bounds,num_radiative_zones,num_convective_zones, ierr)
!       if (ierr /= 0) return
!     endif
      exit retry
      end do retry
      if (elapsed_substep_time.ge.full_timestep) exit
      end do
! omega from the new angular momentum distribution, self-consistently
      do zone_index = 1,s%nz
         s%omega(zone_index) = s%omega_mid(zone_index)
      end do
      call omega_from_j(s%logRho,s%j_rot,s%logR,s%log_mass, &
           s%dm,am_transport_convective_flag,s%nz,s%eta_squared, &
           s%i_rot,s%omega,qiw,s%mean_radius,s%disk_gate_age_gyr)

contains

! ---------------------------------------------------------------
! Without instability transport: the rotation curve (solid body
! everywhere if forced), and the wind's spin-down of the surface
! convection zone over the whole step, if the step has a duration (a
! model solved but not aged gets a zero step and loses nothing; until
! 2026-10 it got its step negated, which the brake applied as a spin-up,
! docs/bug_sweep.md O24).
subroutine spin_down_without_transport()
      omega_surface = s%omega(s%nz)
! solid_body_mode_flag = 1 acts like force_solid_body_rotation
      if(rotation_settings%force_solid_body_rotation .or. (rotation_settings%solid_body_mode_flag.eq.1))then
         do zone_index = 1,s%nz
            am_transport_convective_flag(zone_index) = .true.
         end do
      endif
      call omega_from_j(s%logRho,s%j_rot,s%logR, &
           s%log_mass,s%dm,am_transport_convective_flag,s%nz, &
           s%eta_squared,s%i_rot,s%omega,qiw,s%mean_radius,s%disk_gate_age_gyr)
      if(full_timestep.gt.0.0D0 .and. .not.disk_lock_engaged .and. &
           rotation_settings%use_wind_torque .and. s%convective_flag(s%nz)) then
! the moment of inertia of the surface convection zone (the whole star
! under solid-body rotation), and the masses at its bottom and top
         moment_of_inertia_cz = 0.0D0
         if(rotation_settings%force_solid_body_rotation.or. (rotation_settings%solid_body_mode_flag.eq.1))then
           envelope_boundary_zone = 1
         else
           envelope_boundary_zone = convective_zone_bounds(num_convective_zones,1)
         endif
         do zone_index = envelope_boundary_zone,s%nz
            moment_of_inertia_cz = moment_of_inertia_cz + s%i_rot(zone_index)
         end do
         cz_mass_top = exp(ln10*s%log_total_mass)
         if(envelope_boundary_zone.gt.1)then
            cz_mass_bottom = 0.5D0*(s%m(envelope_boundary_zone)+ &
                 s%m(envelope_boundary_zone-1))
         else
            cz_mass_bottom = 0.0D0
         endif
         wind_loss_active = rotation_settings%use_wind_torque
! the loss rate averaged over the step
         s%fracstep = 0.5
         call brake_convection_zone(s%log_L,full_timestep,cz_mass_bottom, &
              cz_mass_top,envelope_boundary_zone,s%nz,wind_loss_active, &
              omega_surface, &
              s%star_mass,s%log_Teff,moment_of_inertia_cz, &
              s%j_rot, wind_now(), ierr)
         if (ierr /= 0) return
! the convection zone's new rotation after the wind
         iend = s%nz
         call enforce_rotation_profile(s%logRho,s%j_rot,s%logR, &
              s%log_mass,s%dm,envelope_boundary_zone,iend,s%eta_squared, &
              s%i_rot,s%omega,qiw,s%mean_radius,s%nz,s%disk_gate_age_gyr)
      endif
end subroutine spin_down_without_transport

! ---------------------------------------------------------------
! The number of sub-steps: enough that omega changes by at most
! max_domega_step (max_domega_global if unset) per diffusion_timestep_factor.
subroutine count_substeps()
      if(max_domega_step.eq.0.0D0) max_domega_step = rotation_settings%max_domega_global
      num_diffusion_steps = int(max_domega_step/rotation_settings%diffusion_timestep_factor)
      if (mod(max_domega_step,rotation_settings%diffusion_timestep_factor).ne.0.0D0) num_diffusion_steps = num_diffusion_steps + 1
      num_wind_diffusion_steps = int(rotation_settings%max_domega_global/rotation_settings%diffusion_timestep_factor)
      if (mod(rotation_settings%max_domega_global,rotation_settings%diffusion_timestep_factor).ne.0.0D0) num_wind_diffusion_steps = num_wind_diffusion_steps + 1
      num_diffusion_steps = min(num_diffusion_steps,num_wind_diffusion_steps)
      sub_timestep = full_timestep/dfloat(num_diffusion_steps)
end subroutine count_substeps

! ---------------------------------------------------------------
! Advance the sub-step clock; the last sub-step ends exactly at the end
! of the step (within 1e-6 of it counts, against round-off).
subroutine advance_substep_clock()
      delts_test = full_timestep - full_timestep*1.0D-6
      if(elapsed_substep_time+sub_timestep.ge.delts_test)then
         sub_timestep = full_timestep - elapsed_substep_time
         elapsed_substep_time = full_timestep
      else
         elapsed_substep_time = elapsed_substep_time + sub_timestep
      endif
end subroutine advance_substep_clock

! ---------------------------------------------------------------
! omega at the start of the sub-step (the step's starting omega on the
! first, the current omega afterwards), for the wind torque. The
! commented-out advection-diffusion transport (Maeder & Zahn 1998) also
! carried its theta term forward here.
subroutine record_substep_start_omega()
      if(first_call)then
!         rot_scr%theta_prev(1) = rot_scr%theta_rezoned(1)
         rot_scr%omega_substep_start(1) = s%old_omega(1)
         do zone_index = 2,s%nz
!            rot_scr%theta_prev(zone_index) = rot_scr%theta_rezoned(zone_index)
            rot_scr%omega_substep_start(zone_index) = s%old_omega(zone_index)
         end do
      else
         rot_scr%omega_substep_start(1) = s%omega(1)
         do zone_index = 2,s%nz
!            omega_avg = 0.5D0*(star%omega(zone_index)+star%omega(zone_index-1))
!            delta_radius_step = exp10(rot_scr%log_radius_mid(zone_index))- &
!                 exp10(rot_scr%log_radius_mid(zone_index-1))
!            domega_dr = (star%omega(zone_index)-star%omega(zone_index-1))/delta_radius_step
!            rot_scr%theta_prev(zone_index) = rot_scr%theta_mean(zone_index)*omega_avg*domega_dr
            rot_scr%omega_substep_start(zone_index) = s%omega(zone_index)
         end do
      endif
end subroutine record_substep_start_omega

! ---------------------------------------------------------------
! The wind's spin-down of the surface convection zone over the sub-step
! (Weber-Davis magnetic wind), when transport is not doing it.
subroutine brake_surface_convection_zone()
      omega_surface = s%omega_mid(s%nz)
      if(fully_convective_flag.or..not.rotation_settings%instability_transport_active) then
         call brake_convection_zone(s%log_L,sub_timestep,cz_mass_bottom,cz_mass_top, &
              envelope_boundary_zone,s%nz,wind_loss_active,omega_surface, &
              s%star_mass,s%log_Teff,moment_of_inertia_cz, &
              s%j_rot, wind_now(), ierr)
         if (ierr /= 0) return
         iend = s%nz
         call enforce_rotation_profile(rot_scr%log_density_mid,s%j_rot,rot_scr%log_radius_mid, &
              s%log_mass,s%dm,envelope_boundary_zone,iend,rot_scr%eta_squared_mid, &
              rot_scr%moment_of_inertia_mid,s%omega_mid,rot_scr%di_domega_mid,rot_scr%mean_radius_mid,s%nz,s%disk_gate_age_gyr)
      endif
end subroutine brake_surface_convection_zone

! ---------------------------------------------------------------
! The secular instabilities in the radiative regions, as diffusion
! (apply_secular_transport), after the convection zones are mixed. When
! the solve asks for a redo (redo_needed_flag), the angular momentum,
! mean molecular weight and composition go back to the sub-step's start
! and the clock back two sub-steps.
subroutine apply_instability_transport()
      call homogenize_convection_zones(s%xa,s%dm,rot_scr%am_transport_convective_flag_mid,s%nz)
      do zone_index = 1,s%nz
         specific_angular_momentum_saved(zone_index) = s%j_rot(zone_index)
      end do
      call apply_secular_transport(s, rot_scr, mid_scr, circ_scr, &
           sub_timestep, specific_angular_momentum_saved, &
           core_boundary_zone, envelope_boundary_zone, &
           redo_needed_flag, redo_count, &
           moment_of_inertia_cz, cz_mass_bottom, cz_mass_top, &
           omega_surface, surface_cz_active, &
           cod2, diffusion_velocity, diffusion_solve_ok, ierr)
      if (ierr /= 0) return
      if(redo_needed_flag)then
         fx = 2.0D0*sub_timestep/full_timestep
         elapsed_substep_time = elapsed_substep_time - 2.0D0*sub_timestep
         do zone_index = 1,s%nz
            s%j_rot(zone_index) = specific_angular_momentum_saved(zone_index)
            mid_scr%mu(zone_index) = mid_scr%mu(zone_index) - fx*(s%mu(zone_index)-s%old_mu(zone_index))
            do species_index = 1,num_species_tracked
               s%xa(species_index,zone_index) = s%xa_start(species_index,zone_index)
            end do
         end do
      endif
end subroutine apply_instability_transport

! ---------------------------------------------------------------
! Light-element burning over the sub-step: the surface convection zone
! at its end (with envelope overshoot under lovstm), the end rates from
! the mid-sub-step model, the burn, then the convection zones mixed, the
! tiny abundances zeroed and the start of the next sub-step stored.
subroutine burn_light_elements_over_substep()
      call find_convection_zones(s%xa,rot_scr%log_density_mid,rot_scr%log_pressure_mid,rot_scr%log_radius_mid, &
           s%log_mass,rot_scr%log_temperature_mid,rot_scr%convective_flag_mid,s%nz, &
           radiative_zone_bounds,convective_zone_bounds, &
           convective_zone_bounds_burn,core_boundary_zone_cur, &
           envelope_boundary_zone_cur,num_radiative_zones,num_convective_zones, &
           num_convective_zones_burn, s%iov1, s%iov2, ierr)
      if (ierr /= 0) return
      if(rotation_settings%intermediate_overshoot .and. rot_scr%convective_flag_mid(s%nz))then
         s%envelope_overshoot_depth_end = rotation_settings%overshoot_alpha_envelope* &
              pressure_scale_height(rot_scr%log_pressure_mid(envelope_boundary_zone_cur), &
              rot_scr%log_radius_mid(envelope_boundary_zone_cur), &
              rot_scr%log_density_mid(envelope_boundary_zone_cur), &
              s%log_mass(envelope_boundary_zone_cur))
      else
         s%envelope_overshoot_depth_end = 0.0D0
      endif
      call store_end_li_be_rates(s, s%xa,rot_scr%log_density_mid,rot_scr%log_temperature_mid,s%nz)
! The base of the surface mixed region is located on the sub-step's
! interpolated model's gradients with rotational mixing, else on the
! model's own.
      if (rotation_settings%instability_transport_active) then
         call burn_li_be(s, sub_timestep,s%xa,rot_scr%log_radius_mid,s%m, &
              s%dm,rot_scr%log_temperature_mid,mid_scr%grada,mid_scr%gradr, &
              envelope_boundary_zone_cur,envelope_boundary_zone_prev,s%nz)
      else
         call burn_li_be(s, sub_timestep,s%xa,rot_scr%log_radius_mid,s%m, &
              s%dm,rot_scr%log_temperature_mid,s%grada,s%gradr, &
              envelope_boundary_zone_cur,envelope_boundary_zone_prev,s%nz)
      end if
      envelope_boundary_zone_prev = envelope_boundary_zone_cur
      s%envelope_overshoot_depth_start = s%envelope_overshoot_depth_end
      do zone_index = 1,s%nz
         s%rate_li6_start(zone_index) = s%rate_li6(zone_index)
         s%rate_li7_start(zone_index) = s%rate_li7(zone_index)
         s%rate_be9_start(zone_index) = s%rate_be9(zone_index)
      end do
! The commented-out Bulirsch-Stoer path (rotation/burn_settle_mix.f90)
! burned here instead, when instability transport was on and the model
! not fully convective, and left the start-of-step rates alone:
!           call burn_li_be_midpoint(sub_timestep,star%xa,rot_scr%log_radius_mid,star%m, &
!                star%dm,rot_scr%log_temperature_mid,envelope_boundary_zone_cur, &
!                envelope_boundary_zone_prev,star%nz)
      call homogenize_convection_zones(s%xa,s%dm,rot_scr%am_transport_convective_flag_mid,s%nz)
      do zone_index = 1,s%nz
         do species_index = 12,15
            if(s%xa(species_index,zone_index).lt.1.0D-24)s%xa(species_index,zone_index)=0.0D0
         end do
      end do
      do zone_index = 1,s%nz
         do species_index = 12,15
            s%xa_start(species_index,zone_index) = s%xa(species_index,zone_index)
         end do
      end do
end subroutine burn_light_elements_over_substep


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

end subroutine evolve_angular_momentum

end module evolve_rotation

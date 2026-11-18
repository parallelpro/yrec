!----------------------------------------------------------------------
! step_driver
!----------------------------------------------------------------------
! One model: evolve_step advances the star by one timestep (the mass
! change, the composition step, the rezone, the structure solve and its
! coupling to rotation, the updates after convergence, the output and
! the next timestep). A solve that diverges stops the run; nothing cuts
! the timestep.
!
! The stages are the module procedures below evolve_step, in the order
! it calls them. Each takes the star, the step's record (step_state:
! what one stage hands another) where it needs it, and ierr; everything
! else a stage uses is its own local. A new stage is one more such
! procedure and its call (core/README.md, "Add a stage of the step").
module step_driver
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star_info
      implicit none
      private
      public :: evolve_step

! What the stages of one step hand each other. evolve_step owns the one
! instance, new at every call (the defaults below): nothing in it
! outlives the step, and what does is in the star.
      type :: step_state
! The model is aged this step (advance_composition_and_age). False for a
! model that is solved but not aged (a rescale card): its stages get a
! zero step.
         logical :: evolve_model = .false.
! change_mass accreted this step, so that after convergence the envelope
! is rebuilt down to target_envelope_mass, log10 of m(nz)/M before the
! change (rebuild_envelope_after_accretion).
         logical :: new_atmosphere_fit_needed = .false.
         real(dp) :: target_envelope_mass = 0.0_dp
! The base zone of the surface convection zone at the start of the step
! (rezone_or_snapshot, with the extended composition; 0 otherwise), for
! the angular momentum step and the Li/Be burn.
         integer :: envelope_cz_zone_start = 0
! The structure solve so far (solve_level): the Newton iterations done
! over the levels, and whether the last level that ran converged.
         integer :: iterations_done = 0
         logical :: converged = .false.
      end type step_state

! A rescale-only card ages its model only while the centre is cooler
! than this log10 Tc (K): pre-main-sequence models are rescaled and aged
! at the same time (MHP 7/98).
      real(dp), parameter :: rescale_card_aging_max_log_tc = 6.6_dp
! The timestep written to a model that reached its card's end age is at
! least this fraction of the model's age.
      real(dp), parameter :: end_age_min_step_fraction = 1.0e-3_dp
! Gyr and years.
      real(dp), parameter :: years_per_gyr = 1.0e9_dp, gyr_per_year = 1.0e-9_dp

contains

!----------------------------------------------------------------------
! evolve_step (formerly the body of run_yrec's model loop)
!----------------------------------------------------------------------
! Advances the model by one timestep, in named stages (the module
! procedures below):
!   reload_model_if_diverged     the stop after a diverged solve
!   advance_composition_and_age  the mass change, the composition step
!                                (predictor pass) and the new age
!   rezone_or_snapshot           the new mesh and the start-of-step state
!   solve_structure              Henyey corrector levels 1-3 (solve_level)
!   converge_with_rotation       level 4 inside the structure-rotation
!                                iteration (evolve_angular_momentum)
! and then, on the converged model: the shell boundaries
! (locate_shell_boundaries), the Li/Be burn (burn_light_elements), the
! envelope rebuild (rebuild_envelope_after_accretion), the next timestep
! (set_next_timestep), the stitched model, the observables, the turnover
! time (refresh_turnover_timescale), the output and the stop conditions.
! output is the job's output state (run_yrec's). model_iteration is the
! model's count within the run card (the star calibration's check
! starts at the second). The run-card index is star%current_card. Exits:
!   step_status = step_continue        advance accepted, continue
!   step_status = step_kind_card_done  an age, abundance or structure-limit
!                                      stop reached
!   step_status = step_leave_run_loop  target radius crossed
!   (named constants from core/run/stop_conditions.f90)
!   ierr /= 0                          error; run_yrec returns it to the
!                                      caller
!
! Nothing is carried from one step to the next: the step's record is new
! at every call, and what outlives the step is in the star.
!   star reads: control%rewind_short_file, model_diverged_flag
!   star writes: none
subroutine evolve_step(s, output, model_iteration, step_status, ierr)
      use core_messages, only: rewind_core_log
      use shell_boundaries, only: locate_shell_boundaries
      use yrec_output, only: output_state, output_write_model
      use observables, only: compute_observables
      use stitched_model, only: build_stitched_model
      use stop_conditions, only: check_stop_conditions, &
           converged_model_is_nan, step_continue
      type(star_info), intent(inout) :: s
! the job's output state (yrec_output), which run_yrec owns
      type(output_state), intent(inout) :: output
      integer, intent(in) :: model_iteration
      integer, intent(out) :: step_status, ierr
      type(step_state) :: step

      step_status = step_continue
      ierr = 0

! With control%rewind_short_file (legacy LRWSH), rewind the run log so
! that it holds only the current model's solver diagnostics.
      if (s%control%rewind_short_file) call rewind_core_log()

! A diverged solve sets star%model_diverged_flag and cycles retry_step,
! but the step is never retried: nothing cuts the timestep, and on the
! next pass reload_model_if_diverged hands the divergence flag to
! read_starting_model, which stops the run ("ERROR: starting model: the
! previous model failed to converge; run stopped", ierr = 1), as the
! legacy code did. (The loop replaces the legacy backward goto 15.)
      retry_step: do
         call reload_model_if_diverged(s, ierr)
         if (ierr /= 0) return
         call advance_composition_and_age(s, step, ierr)
         if (ierr /= 0) return
         call rezone_or_snapshot(s, step, ierr)
         if (ierr /= 0) return
         call solve_structure(s, step, ierr)
         if (ierr /= 0) return
         if (s%model_diverged_flag) cycle retry_step
         call converge_with_rotation(s, step, ierr)
         if (ierr /= 0) return
         if (s%model_diverged_flag) cycle retry_step
         exit retry_step
      end do retry_step
! A non-finite converged structure is an error, not a model
! (converged_model_is_nan in stop_conditions).
      if (converged_model_is_nan(s)) then
         ierr = 1
         return
      end if
! Locate the hydrogen-burning shell and the edges of the convective core
! and of the surface convection zone, where they exist.
      call locate_shell_boundaries(s)
! Burn Li6, Li7 and Be9 on the converged model.
      call burn_light_elements(s, step, ierr)
      if (ierr /= 0) return
! When change_mass accreted, move the fitting point.
      call rebuild_envelope_after_accretion(s, step, ierr)
      if (ierr /= 0) return
! The timestep of the next model, and the model's count.
      call set_next_timestep(s)
! Assemble the whole converged star (interior, envelope and atmosphere)
! into the stitched arrays before the observables: the observables and
! the turnover time walk them, star%photosphere_logP comes from this
! build, and the output writers below only copy; nothing downstream
! integrates anything.
      call build_stitched_model(s, ierr)
      if (ierr /= 0) return
! The per-model observables, computed in the star layer (one procedure
! per family of quantities, core/diagnostics/observables.f90), including
! the renormalized luminosity_breakdown_lsun; the output below only reads
! them.
      call compute_observables(s, ierr)
      if (ierr /= 0) return
! The turnover time, physics state that the next step's deuterium cap
! and wind read, and its step-start copies.
      call refresh_turnover_timescale(s, ierr)
      if (ierr /= 0) return
! Write the model's output (formerly wrtout).
      call output_write_model(s, output, ierr)
      if (ierr /= 0) return

! The stop conditions (MHP 10/24), in order: the end age, the table of
! central D/X/Y abundances, the structure limits (log L, Teff, log g,
! nu_max) and the star-calibration target radius.
      call check_stop_conditions(s, model_iteration, step_status)
end subroutine evolve_step

! ---------------------------------------------------------------
! If the previous attempt diverged, call read_starting_model (formerly
! starin) with the divergence flag set, which ends the run with ierr = 1
! (no timestep is cut; see the note at retry_step). Otherwise nothing
! happens here, except that star%punch_pending_flag is always set; with
! the helium-flash option, the next run card stores its starting model
! on last_model_unit when the flag is set.
!   star reads: model_diverged_flag
!   star writes: punch_pending_flag
subroutine reload_model_if_diverged(s, ierr)
      use starting_model, only: read_starting_model
      use stop_conditions, only: check_rotation_initialised
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr

      ierr = 0
      if (s%model_diverged_flag) then
         call read_starting_model(s, ierr)
         if (ierr /= 0) return
         call check_rotation_initialised(s, ierr)
         if (ierr /= 0) return
      end if
      s%punch_pending_flag = .true.
end subroutine reload_model_if_diverged

! ---------------------------------------------------------------
! The composition half of the step, when the model is to be aged
! (step%evolve_model, set here): the mass change at the surface
! (change_mass, which sets step%target_envelope_mass and, when it
! accretes, step%new_atmosphere_fit_needed), the start-of-step
! composition, the composition step (burning and mixing; the predictor
! pass of composition_step), and the new age. ierr: change_mass' or the
! composition step's code.
!   star reads: card_kind, current_card, dt, logT, model_number, nz,
!     physics%use_extended_composition, xa
!   star writes: age_gyr, core_cz_top_index, envelope_cz_bottom_index,
!     timestep_yr, xa_start
subroutine advance_composition_and_age(s, step, ierr)
      use star_info_lib, only: max_convective_zones, n_species_basic, &
           n_species_extended
      use run_cards, only: card_rescale
      use physical_constants, only: seconds_per_year
      use mass_change, only: change_mass
      use composition_step, only: advance_composition, predictor_pass
      type(star_info), intent(inout) :: s
      type(step_state), intent(inout) :: step
      integer, intent(out) :: ierr
! an output of the composition step that the step does not use: the
! convection zones without overshoot
      integer :: mixed_zone_bounds_no_overshoot(max_convective_zones,2)
      integer :: num_species, i, j

      ierr = 0
! The model is aged when model_number >= 0, except on a rescale-only
! card, where it is aged only while the centre is cooler than
! rescale_card_aging_max_log_tc.
      step%evolve_model = s%model_number >= 0 .and. &
           (s%card_kind(s%current_card) /= card_rescale .or. &
           s%logT(1) < rescale_card_aging_max_log_tc)
      step%new_atmosphere_fit_needed = .false.
      if (.not. step%evolve_model) return
! Mass loss or accretion at the surface.
      call change_mass(s, step%target_envelope_mass, &
           step%new_atmosphere_fit_needed, ierr)
      if (ierr /= 0) return
! Store the composition at the start of the timestep (xa_start).
      num_species = n_species_basic
      if (s%physics%use_extended_composition) num_species=n_species_extended
      do i = 1,s%nz
         do j = 1,num_species
            s%xa_start(j,i) = s%xa(j,i)
         end do
      end do
      s%timestep_yr = s%dt/seconds_per_year
! The composition step on the start-of-step structure. Its zone outputs
! are arguments, not star members set inside it, because the corrector
! passes at levels >= 3 (solve_level) pass their own local storage,
! deliberately separate from these star copies.
      call advance_composition(s, s%dt, predictor_pass, &
           s%core_cz_top_index, s%envelope_cz_bottom_index, &
           mixed_zone_bounds_no_overshoot, ierr)
      if (ierr /= 0) return
! The new age (Gyr).
      s%age_gyr = s%age_gyr + gyr_per_year*s%timestep_yr
end subroutine advance_composition_and_age

! ---------------------------------------------------------------
! Rezone the new model (rezone, formerly hpoint), which also keeps the
! start-of-step structure in star%*_start; with the helium-flash option
! there is no rezoning, and the start-of-step structure is copied here
! instead. Then store the start-of-step rotational kinetic energy and,
! with the light elements, their start-of-step state (the surface
! convection zone's base, kept in step%envelope_cz_zone_start, and the
! Li/Be burning rates). ierr: the rezone's or the overshoot depth's code.
!   star reads: envelope_cz_bottom_index, kinetic_energy_rot, log_mass,
!     log_Teff, logP, logR, logRho, logT, luminosity_lsun, nz,
!     physics%helium_flash_active, physics%rotation_active,
!     physics%use_extended_composition, xa
!   star writes: cz_base_radius_last_call,
!     envelope_cz_bottom_index_start, envelope_overshoot_depth_start,
!     kinetic_energy_rot_start, log_Teff_start, logP_start, logR_start,
!     logRho_start, logT_start, luminosity_lsun_start,
!     surface_luminosity_lsun_start
subroutine rezone_or_snapshot(s, step, ierr)
      use rezoning, only: rezone
      use convection_lib, only: envelope_overshoot_depth
      use burn_step, only: store_start_li_be_rates
      type(star_info), intent(inout) :: s
      type(step_state), intent(inout) :: step
      integer, intent(out) :: ierr
      integer :: i

      ierr = 0
! Rezone the new model, except with the helium-flash option.
      if (.not.s%physics%helium_flash_active) then
         call rezone(s, ierr)
         if (ierr /= 0) return
! The base zone of the surface convection zone, where deuterium burning
! is capped (compute_zone_energy, compute_burn_rates).
         s%envelope_cz_bottom_index_start = s%envelope_cz_bottom_index
      else
! Keep the start-of-step structure, which the structure-change timestep
! limit (structure_change_limit, formerly ptime) compares against.
         do i=1, s%nz
            s%logP_start(i) = s%logP(i)
            s%logT_start(i) = s%logT(i)
            s%logR_start(i) = s%logR(i)
            s%luminosity_lsun_start(i) = s%luminosity_lsun(i)
            s%logRho_start(i) = s%logRho(i)
         end do
! Teff and the surface L are kept as well.
         s%log_Teff_start = s%log_Teff
         s%surface_luminosity_lsun_start = s%luminosity_lsun(s%nz)
      end if
! The start-of-step rotational kinetic energy of each zone.
      if (s%physics%rotation_active) then
         do i = 1,s%nz
            s%kinetic_energy_rot_start(i) = s%kinetic_energy_rot(i)
         end do
      end if
! The start-of-step depth of the surface convection zone, for the Li/Be
! burn: its base zone, and the envelope overshoot depth alpha_ov H_P
! below it (envelope_overshoot_depth_start).
! cz_base_radius_last_call = 0 tells locate_surface_mixing_base that no base
! radius has been stored yet this step.
      if (s%physics%use_extended_composition) then
         s%cz_base_radius_last_call = 0.0D0
         step%envelope_cz_zone_start = s%envelope_cz_bottom_index
         call envelope_overshoot_depth(s%xa, s%logRho, s%logP, s%logR, s%log_mass, &
              s%logT, s%envelope_cz_bottom_index, s%envelope_overshoot_depth_start, ierr)
         if (ierr /= 0) return
! The Li/Be burning rates of the start-of-step model.
         call store_start_li_be_rates(s)
      end if
end subroutine rezone_or_snapshot

! ---------------------------------------------------------------
! The structure solve, Henyey corrector levels 1-3 (solve_level): the
! first guess (the previous model or, with
! physics%improved_first_guess_flag, the previous model extrapolated by its
! rates of change); level 1 on the previous model's envelope triangle;
! level 2 after re-checking the surface boundary condition; the rates of
! change for the next model's first guess; level 3. Returns early with
! star%model_diverged_flag set when a level diverges (evolve_step then
! cycles retry_step); ierr on error. step carries the aging flag in and
! the iteration count and the convergence verdict out.
!   star reads: dt, logR_start, luminosity_lsun_start, nz,
!     physics%improved_first_guess_flag, physics%max_iter_level1,
!     physics%max_iter_level2, physics%max_iter_level3,
!     physics%new_triangle_init_flag,
!     physics%solver_underrelax_increment,
!     physics%solver_underrelax_start
!   star writes: dln_luminosity_dt, dlog_pressure_dt, dlog_radius_dt,
!     dlog_temperature_dt, eps_grav, log_pressure_delta,
!     log_temperature_delta, logP, logR, logT, luminosity_lsun,
!     model_diverged_flag, newton_correction_factor,
!     recompute_envelope_triangle
subroutine solve_structure(s, step, ierr)
      type(star_info), intent(inout) :: s
      type(step_state), intent(inout) :: step
      integer, intent(out) :: ierr
      real(dp) :: delta_lum_step, delta_pressure_step, &
           delta_radius_step, delta_temp_step
      integer :: i

      ierr = 0
! With physics%new_triangle_init_flag, a new envelope triangle is
! computed on the first iteration (the previous model's is ignored).
      if (s%physics%new_triangle_init_flag) s%recompute_envelope_triangle = .true.
! A model that is solved but not aged (step%evolve_model false) keeps
! its s%dt; the stages that age read the flag: the first guess below,
! the rates of change after level 2, the composition rerun and the
! Henyey pass (solve_level, which hands the pass a zero step), the
! light-element burn and the rotation step. (Until 2026-10 the flag was
! encoded as a negative s%dt, made positive again before
! compute_timestep.)
! The under-relaxation factor starts at |fcorr0| - fcorri. The sign of
! physics%solver_underrelax_start selects whether henyey_iterate raises
! it by fcorri per iteration, up to 1 (fcorr0 > 0), or holds it fixed.
      s%newton_correction_factor = abs(s%physics%solver_underrelax_start) - &
           s%physics%solver_underrelax_increment
      step%iterations_done = 0
      s%model_diverged_flag = .false.
      step%converged = .false.
! Without the improved first guess, or for a model not aged, start from
! the previous model: zero the changes of log T and log P over the step
! (from which the energy equation builds the entropy term T dS/dt) and
! the rates of change of the previous step.
      if (.not.s%physics%improved_first_guess_flag .or. .not.step%evolve_model &
           .or. s%dt <= 0.0D0) then
         do i = 1,s%nz
            s%log_temperature_delta(i) = 0.0D0
            s%log_pressure_delta(i) = 0.0D0
            s%dlog_temperature_dt(i) = 0.0D0
            s%dlog_pressure_dt(i) = 0.0D0
            s%dln_luminosity_dt(i) = 0.0D0
            s%dlog_radius_dt(i) = 0.0D0
! Zero the gravitational luminosity.
            s%eps_grav(i) = 0.0D0
         end do
      else
! Improved first guess: extrapolate each variable by its rate of change
! over the previous step (dlog_temperature_dt, dlog_pressure_dt,
! dln_luminosity_dt and dlog_radius_dt):
!   log T += (d log T/dt) dt,  log P += (d log P/dt) dt,
!   L += L (d ln L/dt) dt,     log R += (d log R/dt) dt,
! and start the changes of log T and log P over the step at those values.
         do i = 1,s%nz
            delta_temp_step = s%dlog_temperature_dt(i)*s%dt
            delta_pressure_step = s%dlog_pressure_dt(i)*s%dt
            delta_lum_step = s%luminosity_lsun(i)*s%dln_luminosity_dt(i)*s%dt
            delta_radius_step = s%dlog_radius_dt(i)*s%dt
            s%log_temperature_delta(i) = delta_temp_step
            s%log_pressure_delta(i) = delta_pressure_step
            s%logT(i) = s%logT(i) + delta_temp_step
            s%logP(i) = s%logP(i) + delta_pressure_step
            s%luminosity_lsun(i) = s%luminosity_lsun(i) + delta_lum_step
            s%logR(i) = s%logR(i) + delta_radius_step
! Zero the gravitational luminosity.
            s%eps_grav(i) = 0.0D0
         end do
      end if

! Level 1 uses the previous model's envelope triangle; for the first
! model of a run the triangle is generated here. The four levels differ
! only in their arguments to solve_level: the level, its iteration limit
! and whether the surface boundary condition is re-checked.
      call solve_level(s, step, 1, s%physics%max_iter_level1, .false., ierr)
      if (ierr /= 0) return
      if (s%model_diverged_flag) return
! Level 2 checks the envelope triangle (the surface boundary condition)
! before iterating.
      call solve_level(s, step, 2, s%physics%max_iter_level2, .true., ierr)
      if (ierr /= 0) return
      if (s%model_diverged_flag) return
! The rates of change of the structure over this step, from the level-2
! solution, for the improved first guess of the next model
! (physics%improved_first_guess_flag, legacy LNEWS; 7/91):
!   d log T/dt = log_temperature_delta/dt,  d log P/dt = log_pressure_delta/dt,
!   d ln L/dt = 2 (L - L_start)/(L + L_start)/dt,
!   d log R/dt = (log R - log R_start)/dt.
      if (step%evolve_model .and. s%dt > 0.0D0) then
         do i = 1,s%nz
            s%dlog_temperature_dt(i)=s%log_temperature_delta(i)/s%dt
            s%dlog_pressure_dt(i)=s%log_pressure_delta(i)/s%dt
            s%dln_luminosity_dt(i)=2.0D0*(s%luminosity_lsun(i)-s%luminosity_lsun_start(i))/ &
                 (s%luminosity_lsun(i)+s%luminosity_lsun_start(i))/s%dt
            s%dlog_radius_dt(i)=(s%logR(i)-s%logR_start(i))/s%dt
         end do
      end if
! Level 3: from here on, solve_level redoes the composition step on the
! updated structure before iterating.
      call solve_level(s, step, 3, s%physics%max_iter_level3, .false., ierr)
end subroutine solve_structure

! ---------------------------------------------------------------
! The structure-rotation iteration (num_rotation_structure_passes
! passes; one without rotation). Each pass: corrector level 4, the
! convergence check, the convection zones mixed fully
! (homogenize_convection_zones, formerly mixcz), then the angular
! momentum update (evolve_angular_momentum, formerly getw) and the new
! shape factors (rotation_shape_factors, formerly fpft). The
! start-of-step composition and angular momentum are restored on repeat
! passes, so that nothing is counted twice. A pass whose level 4 has not
! converged (step%converged) sets star%model_diverged_flag and returns.
!   star reads: age_gyr, dm, dt, log_mass, logR, logRho, nz, omega,
!     physics%max_iter_level4, physics%rotation_active,
!     physics%use_wind_torque
!   star writes: convective_flag, disk_gate_age_gyr, eta_squared,
!     fp_rot, ft_rot, j_rot, j_rot_iteration_start, max_domega_frac,
!     mean_gravity, mean_radius, model_diverged_flag,
!     num_rotation_structure_passes, xa, xa_iteration_start, xa_start
subroutine converge_with_rotation(s, step, ierr)
      use star_info_lib, only: n_species_extended
      use evolve_rotation, only: evolve_angular_momentum
      use rotating_shape, only: rotation_shape_factors
      use convection_lib, only: homogenize_convection_zones
      type(star_info), intent(inout) :: s
      type(step_state), intent(inout) :: step
      integer, intent(out) :: ierr
      integer :: i, j, itrot
      logical :: wind_loss_active

      ierr = 0
      if (.not.s%physics%rotation_active) then
         s%num_rotation_structure_passes = 1
      end if
! When the structure and rotation are iterated (MHP 05/02), keep copies
! of the start-of-step specific angular momentum (j_rot) and composition
! (xa_start, legacy HCOMPP), so that every pass starts from them.
      if (s%num_rotation_structure_passes > 1) then
         do i = 1,s%nz
            s%j_rot_iteration_start(i) = s%j_rot(i)
            do j = 1,n_species_extended
               s%xa_iteration_start(j,i) = s%xa_start(j,i)
            end do
         end do
      end if
      do itrot = 1, s%num_rotation_structure_passes
! On repeat passes, restore the start-of-step composition (MHP 05/02).
         if (itrot > 1) then
            do i = 1,s%nz
               do j = 1,n_species_extended
                  s%xa_start(j,i) = s%xa_iteration_start(j,i)
               end do
            end do
         end if
! Level 4 repeats the iteration between the mixing (the composition
! step, redone in solve_level) and the structure (7/91). Use it for
! high-precision work such as solar models, but not where semiconvection
! matters: iterating between the burning and the structure then
! oscillates. Unlike levels 1-3 it runs inside the structure-rotation
! iteration. The surface boundary condition is not re-checked here
! (check_surface_bc = .false.), although level 3 changed the composition
! and hence the structure; the legacy comment left that check for later.
         call solve_level(s, step, 4, s%physics%max_iter_level4, .false., ierr)
         if (ierr /= 0) return
         if (.not.step%converged) then
! Not converged within the iterations of levels 1-4 (max_iter_level1 to
! max_iter_level4): flag the divergence.
            s%model_diverged_flag = .true.
            return
         end if

! The model has converged. Mix the convection zones fully.
         call homogenize_convection_zones(s%xa,s%dm,s%convective_flag,s%nz)

! The age (Gyr) that the disk-locking test reads.
         s%disk_gate_age_gyr = s%age_gyr
         if (s%physics%rotation_active) then
! On repeat passes, restore the start-of-step angular momentum.
            if (itrot > 1) then
               do i = 1,s%nz
                  s%j_rot(i) = s%j_rot_iteration_start(i)
               end do
            end if
! Evolve the angular momentum over the step: the new run of omega (with
! the wind torque when physics%use_wind_torque is set). A model not aged
! gets a zero step, as the Henyey pass does: it spans no time, so no
! angular momentum moves or leaves. (Until 2026-10 it got its step
! negated, and without instability transport the wind brake applied it,
! spinning the surface convection zone up; docs/bug_sweep.md O24.)
            wind_loss_active = s%physics%use_wind_torque
            call evolve_angular_momentum(s, merge(s%dt, 0.0d0, step%evolve_model), &
                 s%max_domega_frac, wind_loss_active, &
                 step%envelope_cz_zone_start, ierr)
            if (ierr /= 0) return
! The rotational shape factors f_P and f_T for the new omega.
            call rotation_shape_factors(s%logRho, s%logR, s%log_mass, s%nz, s%omega, &
                 s%eta_squared, s%fp_rot, s%ft_rot, s%mean_gravity, s%mean_radius, ierr)
            if (ierr /= 0) return
         end if
      end do
end subroutine converge_with_rotation

! ---------------------------------------------------------------
! One level of the corrector ladder. Every level runs the same three
! stages; only the level number, its iteration limit and whether the
! surface boundary condition is re-checked differ:
!   1. refresh_surface_bc: the surface boundary condition (before the
!      composition step, because it reads the surface composition);
!   2. at levels >= 3, for a model being aged (step%evolve_model): the
!      composition step again (corrector pass), from the start-of-step
!      composition, on the updated structure;
!   3. henyey_iterate: the Newton iterations, which add to
!      step%iterations_done and set step%converged.
! The caller checks ierr and the divergence flag after each level. A
! level with an iteration limit <= 0 is skipped, and leaves
! step%converged as the level before set it.
!   star reads: dt, envelope_cz_bottom_index, nz,
!     physics%use_extended_composition, xa_start
!   star writes: xa
subroutine solve_level(s, step, level, level_max_iterations, check_surface_bc, ierr)
      use star_info_lib, only: max_convective_zones, n_species_basic, &
           n_species_extended
      use composition_step, only: advance_composition, corrector_pass
      use henyey_iteration, only: refresh_surface_bc, henyey_iterate
      type(star_info), intent(inout) :: s
      type(step_state), intent(inout) :: step
      integer, intent(in) :: level, level_max_iterations
      logical, intent(in) :: check_surface_bc
      integer, intent(out) :: ierr
      integer :: saha_iterations, envelope_zone_index, num_species, i, j
      logical :: envelope_recomputed_flag, recompute_surface_bc
! The zone outputs of the redone composition step. Only the envelope
! edge is used (by the Newton iterations, below); the others are
! discarded, so star's copies keep their values from the predictor
! pass.
      integer :: core_cz_edge_rerun
      integer :: mixed_zone_bounds_rerun(max_convective_zones,2)

      ierr = 0
      recompute_surface_bc = check_surface_bc
      if (level_max_iterations <= 0) return

      call refresh_surface_bc(s, level, recompute_surface_bc, &
           saha_iterations, envelope_recomputed_flag, ierr)
      if (ierr /= 0) return

! envelope_zone_index gates the accretion correction to the entropy term
! of the energy equation (compute_henyey_coefficients adds it from that
! zone out, the surface convection zone the accreted matter mixes into).
! It is the predictor's envelope edge (star%envelope_cz_bottom_index) at
! levels 1-2 and the redone composition step's at levels >= 3.
      envelope_zone_index = s%envelope_cz_bottom_index
      core_cz_edge_rerun = 0
      mixed_zone_bounds_rerun = 0
      if (level > 2 .and. step%evolve_model .and. s%dt > 0.0d0) then
         num_species = n_species_basic
         if (s%physics%use_extended_composition) num_species = n_species_extended
         do i = 1,s%nz
            do j = 1,num_species
               s%xa(j,i) = s%xa_start(j,i)
            end do
         end do
         call advance_composition(s, s%dt, corrector_pass, &
              core_cz_edge_rerun, envelope_zone_index, &
              mixed_zone_bounds_rerun, ierr)
         if (ierr /= 0) return
      end if

! The pass gets the step only for a model being aged (a zero step: no
! time derivatives).
      call henyey_iterate(s, merge(s%dt, 0.0d0, step%evolve_model), &
           level_max_iterations, step%converged, recompute_surface_bc, &
           envelope_recomputed_flag, saha_iterations, envelope_zone_index, &
           step%iterations_done, ierr)
end subroutine solve_level

! ---------------------------------------------------------------
! The Li6, Li7 and Be9 burn over the step, on the converged model, for
! models without rotation or without rotational mixing (with it,
! evolve_angular_momentum burns them in its sub-steps): the end-of-step
! convection zones (find_convection_zones, formerly convec) and overshoot
! depth, the end-of-step burning rates (store_end_li_be_rates), then the
! burn (burn_li_be), from the start-of-step base zone
! step%envelope_cz_zone_start. Only with the extended composition, for a
! model being aged (step%evolve_model). Deuterium burns in the
! composition step instead. ierr: the zone search's or the overshoot
! depth's code.
!   star reads: dm, dt, grada, gradr, log_mass, logP, logR, logRho,
!     logT, m, model_number, nz, physics%instability_transport_active,
!     physics%rotation_active, physics%use_extended_composition
!   star writes: convective_flag, core_cz_top_index,
!     envelope_cz_bottom_index, envelope_overshoot_bottom_index,
!     envelope_overshoot_depth_end, envelope_overshoot_top_index,
!     mixed_zone_bounds, xa
subroutine burn_light_elements(s, step, ierr)
      use star_info_lib, only: max_convective_zones, max_radiative_zones
      use burn_step, only: store_end_li_be_rates, burn_li_be
      use convection_lib, only: find_convection_zones, envelope_overshoot_depth
      type(star_info), intent(inout) :: s
      type(step_state), intent(in) :: step
      integer, intent(out) :: ierr
! outputs of the zone search that the step does not use: the convection
! zones without overshoot, the radiative zones and the counts
      integer :: mixed_zone_bounds_no_overshoot(max_convective_zones,2), &
           radiative_zone_bounds(max_radiative_zones,2)
      integer :: num_mixed_zones, num_mixed_zones_no_overshoot, &
           num_radiative_zones, envelope_cz_zone_end

      ierr = 0
      if (.not. (s%physics%use_extended_composition .and. s%model_number >= 0 .and. &
           step%evolve_model .and. s%dt > 0.0D0)) return
! Only without rotation, or without rotational mixing.
      if (s%physics%rotation_active .and. s%physics%instability_transport_active) return
! The convection zones at the end of the step.
      call find_convection_zones(s%xa, s%logRho, s%logP, s%logR, s%log_mass, s%logT, &
           s%convective_flag, s%nz, radiative_zone_bounds, s%mixed_zone_bounds, &
           mixed_zone_bounds_no_overshoot, s%core_cz_top_index, s%envelope_cz_bottom_index, &
           num_radiative_zones, num_mixed_zones, num_mixed_zones_no_overshoot, &
           s%envelope_overshoot_bottom_index, s%envelope_overshoot_top_index, ierr)
      if (ierr /= 0) return
! The end-of-step base zone of the surface convection zone and the
! envelope overshoot depth below it (alpha_ov H_P).
      envelope_cz_zone_end = s%envelope_cz_bottom_index
      call envelope_overshoot_depth(s%xa, s%logRho, s%logP, s%logR, s%log_mass, &
           s%logT, s%envelope_cz_bottom_index, s%envelope_overshoot_depth_end, ierr)
      if (ierr /= 0) return
! The burning rates at the end of the step, from the converged model;
! then the burn.
      call store_end_li_be_rates(s, s%xa,s%logRho,s%logT,s%nz)
      call burn_li_be(s, s%dt,s%xa,s%logR,s%m,s%dm,s%logT,s%grada,s%gradr, &
           envelope_cz_zone_end,step%envelope_cz_zone_start,s%nz)
end subroutine burn_light_elements

! ---------------------------------------------------------------
! When change_mass accreted (step%new_atmosphere_fit_needed), move the
! fitting point: rebuild the envelope down to step%target_envelope_mass
! (rebuild_envelope), then the rotational shape factors f_P and f_T
! (fp_rot, ft_rot) for the current omega on the new points.
!   star reads: log_mass, logR, logRho, nz, omega
!   star writes: eta_squared, fp_rot, ft_rot, mean_gravity, mean_radius
subroutine rebuild_envelope_after_accretion(s, step, ierr)
      use envelope_rebuild, only: rebuild_envelope
      use rotating_shape, only: rotation_shape_factors
      type(star_info), intent(inout) :: s
      type(step_state), intent(inout) :: step
      integer, intent(out) :: ierr

      ierr = 0
      if (.not. step%new_atmosphere_fit_needed) return
      call rebuild_envelope(s, step%target_envelope_mass, ierr)
      if (ierr /= 0) return
      call rotation_shape_factors(s%logRho, s%logR, s%log_mass, s%nz, s%omega, &
           s%eta_squared, s%fp_rot, s%ft_rot, s%mean_gravity, s%mean_radius, ierr)
      if (ierr /= 0) return
      step%new_atmosphere_fit_needed = .false.
end subroutine rebuild_envelope_after_accretion

! ---------------------------------------------------------------
! The timestep of the next model (compute_timestep, from the magnitude
! of the step just taken, as a stored model's timestep is taken), and
! the model's count.
! With a fixed end age (card_max_age_yr > 0) that this model reached,
! compute_timestep has clipped the next timestep to (nearly) zero. The
! timestep just taken (dt_saved, floored at end_age_min_step_fraction of
! the age) is written to the model instead, so that the sequence can be
! continued. Other stops keep the new timestep. (The floor read
! 1e-3*age_gyr years, the age in Gyr taken as years, until 2026-10, as
! in main.f.)
!   star reads: age_gyr, card_kind, card_max_age_yr, card_stops_active,
!     current_card
!   star writes: dt, dt_saved, model_number, timestep_yr
subroutine set_next_timestep(s)
      use run_cards, only: card_rescale
      use physical_constants, only: seconds_per_year
      use timestep_control, only: compute_timestep
      use stop_conditions, only: reached_end_age
      type(star_info), intent(inout) :: s

      s%dt = abs(s%dt)
      s%dt_saved = s%dt
      call compute_timestep(s)
      if (s%card_stops_active(s%current_card) .and. s%card_max_age_yr(s%current_card) > 0.0D0) then
         if (reached_end_age(s)) then
            s%dt = max(s%dt_saved, &
                 end_age_min_step_fraction*s%age_gyr*years_per_gyr*seconds_per_year)
            s%timestep_yr = s%dt/seconds_per_year
         else
            s%dt_saved = s%dt
         end if
      else
         s%dt_saved = s%dt
      end if
! A rescale-only card does not count models.
      if (s%card_kind(s%current_card) /= card_rescale) s%model_number = s%model_number+1
end subroutine set_next_timestep

! ---------------------------------------------------------------
! The convective turnover time of the converged model, from the stitched
! model (compute_turnover_timescale), and the copies of it and of
! photosphere_logP in their _prev members: the next step's start values,
! from which the wind interpolates to the step's end values
! (evolve_rotation, secular_transport). compute_turnover_timescale also
! writes star%envelope_radius_rsun, the base radius it finds, over the
! value compute_observables gave it. read_starting_model does the same
! for a model it reads.
!   star reads: convective_turnover_timescale, photosphere_logP
!   star writes: convective_turnover_timescale_prev,
!     photosphere_logP_prev
subroutine refresh_turnover_timescale(s, ierr)
      use turnover_timescale, only: compute_turnover_timescale
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      call compute_turnover_timescale(s, ierr)
      s%convective_turnover_timescale_prev = &
           s%convective_turnover_timescale
      s%photosphere_logP_prev = s%photosphere_logP
end subroutine refresh_turnover_timescale

end module step_driver

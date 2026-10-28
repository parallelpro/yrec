!----------------------------------------------------------------------
! step_driver
!----------------------------------------------------------------------
! One model: evolve_step advances the star by one timestep (the mass
! change, the composition step, the rezone, the structure solve and its
! coupling to rotation, the updates after convergence, the output and
! the next timestep). A solve that diverges stops the run; nothing cuts
! the timestep.
module step_driver
      use shell_boundaries, only: locate_shell_boundaries
      use rezoning, only: rezone
      use starting_model, only: read_starting_model
      use envelope_rebuild, only: rebuild_envelope
      use star_info_lib, only: star_info, json
      implicit none
      private
      public :: evolve_step

contains

!----------------------------------------------------------------------
! evolve_step (formerly the body of run_yrec's model loop)
!----------------------------------------------------------------------
! Advances the model by one timestep, in named stages (the internal
! procedures below):
!   advance_composition_and_age  the mass change, the composition step
!                                (predictor pass) and the new age
!   rezone_or_snapshot           the new mesh and the start-of-step state
!   solve_structure              Henyey corrector levels 1-3 (solve_level)
!   converge_with_rotation       level 4 inside the structure-rotation
!                                iteration (evolve_angular_momentum)
! and then, on the converged model: the shell boundaries, the Li/Be burn,
! the envelope rebuild, the next timestep (compute_timestep), the
! stitched model, the observables, the output and the stop conditions.
! model_iteration is the model's count within the run card (the star
! calibration's check starts at the second). Exits:
!   step_status = step_continue        advance accepted, continue
!   step_status = step_kind_card_done  an age, abundance or structure-limit
!                                      stop reached
!   step_status = step_leave_run_loop  target radius crossed
!   (named constants from core/run/stop_conditions.f90)
!   ierr /= 0                          error; run_yrec returns it to the
!                                      caller
!
! The blanket save is deliberate: these locals were statics of the legacy
! main program, and several carry state from one step to the next (the
! convergence bookkeeping, the surface convection zone's base zone of the
! previous step). The reset block below (evolve_step_reset_pending)
! clears them for a second run in the same process, alongside
! evolve_state.
subroutine evolve_step(s, model_iteration, step_status, ierr)
      use timestep_control, only: compute_timestep
      use evolve_rotation, only: evolve_angular_momentum
      use rotating_shape, only: rotation_shape_factors

      use star_info_lib, only: evolve_step_reset_pending, i_h1, i_h2, &
           max_convective_zones, max_radiative_zones, &
           i_he4, n_species_basic, n_species_extended
      use io_units
      use physical_constants
      use burn_step, only: store_start_li_be_rates, store_end_li_be_rates, &
           burn_li_be
      use composition_step, only: advance_composition, predictor_pass, corrector_pass
      use henyey_iteration, only: refresh_surface_bc, henyey_iterate
      use convection_geometry, only: pressure_scale_height
      use yrec_output, only: output_write_model
      use observables, only: compute_observables
      use stitched_model, only: build_stitched_model
      use mass_change, only: change_mass
      use convection_lib, only: find_convection_zones, homogenize_convection_zones
      use stop_conditions
      implicit none
! outputs of the calls below that the step does not use: the timestep
! compute_timestep chose, the model read's ikut flag and convective
! velocity, the rotation's di_domega, the convection zones without
! overshoot and the radiative zones
      double precision :: chosen_dt, convective_velocity
      integer :: ikut_flag
      double precision :: qiw(json)
      integer :: mixed_zone_bounds_no_overshoot(max_convective_zones,2), &
           radiative_zone_bounds(max_radiative_zones,2)
      type(star_info), intent(inout) :: s

! The run-card index is star%job%nk (legacy NK in /ZRAMP/), not an
! argument.
      integer, intent(in) :: model_iteration
      integer, intent(out) :: step_status, ierr

      double precision :: delta_lum_step, delta_pressure_step, &
           delta_radius_step, delta_temp_step, &
           target_envelope_mass
      integer :: envelope_cz_zone_end, envelope_cz_zone_prev, &
           iteration_level, iterations_done, max_iterations, &
           num_mixed_zones, num_mixed_zones_no_overshoot, &
           num_radiative_zones, num_species, i, j, ii, itrot, jerr
      logical :: converged, evolve_model_flag, &
           new_atmosphere_fit_needed, recompute_surface_bc, &
           wind_loss_active
! Set by compute_timestep when it switches the structure-based timestep
! limits off (log Tc > 7.1 with L_grav < 0); the driver then clears
! star%job%use_structure_dt_limits.
      logical :: disable_structure_dt_limits
      save
   ! INTENTIONAL: cross-step driver state (see the header); reset via
   ! evolve_step_reset_pending
      step_status = step_continue
      ierr = 0

! On a second run_yrec call in the same process, put the saved locals
! back to their process-start value (zero). yrec_reset sets the flag; a
! fresh process never sees it, so the first run is untouched.
      if (evolve_step_reset_pending) then
         delta_lum_step = 0.0d0
         delta_pressure_step = 0.0d0
         delta_radius_step = 0.0d0
         delta_temp_step = 0.0d0
         target_envelope_mass = 0.0d0
         envelope_cz_zone_end = 0
         envelope_cz_zone_prev = 0
         iteration_level = 0
         iterations_done = 0
         max_iterations = 0
         num_mixed_zones = 0
         num_mixed_zones_no_overshoot = 0
         num_radiative_zones = 0
         num_species = 0
         converged = .false.
         evolve_model_flag = .false.
         new_atmosphere_fit_needed = .false.
         recompute_surface_bc = .false.
         wind_loss_active = .false.
         evolve_step_reset_pending = .false.
      end if

! The stages are the internal procedures below. Each sets ierr on an
! error, checked after each call; a diverged solve sets
! star%model_diverged_flag and cycles retry_step.
      call update_output_flags_for_step

! A diverged solve cycles retry_step, but the step is never retried:
! nothing cuts the timestep, and on the next pass reload_model_if_diverged
! hands the divergence flag to read_starting_model, which stops the run
! ("MODEL FAILED TO CONVERGE", ierr = 1), as the legacy code did. (The
! loop replaces the legacy backward goto 15.)
      retry_step: do
         call reload_model_if_diverged
         if (ierr /= 0) return
         call advance_composition_and_age
         if (ierr /= 0) return
         call rezone_or_snapshot
         if (ierr /= 0) return
         call solve_structure
         if (ierr /= 0) return
         if (s%model_diverged_flag) cycle retry_step
         call converge_with_rotation
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
       call locate_shell_boundaries(s, s%xa,s%luminosity_lsun,s%convective_flag,s%nz, &
              s%core_cz_top_index,s%envelope_cz_bottom_index,s%h_shell_zone_begin,s%h_shell_end_index,s%h_shell_midpoint_zone, &
              s%has_h_shell)
! Burn Li6, Li7 and Be9 on the converged model (burn_light_elements,
! below).
         call burn_light_elements
! When change_mass accreted (new_atmosphere_fit_needed), move the
! fitting point: rebuild the envelope down to target_envelope_mass.
         if (new_atmosphere_fit_needed) then
            call rebuild_envelope(s, target_envelope_mass,s%xa,s%logRho,s%luminosity_lsun,s%logP,s%logR,s%log_mass,s%m,s%dm, &
                            s%log_total_mass,s%logT,s%convective_flag,s%eta_squared,s%i_rot,s%j_rot,qiw,s%mean_radius, &
                            s%kinetic_energy_rot,s%log_L,s%total_angular_momentum,s%total_rotational_ke,s%log_Teff,s%nz,s%recompute_envelope_triangle,ierr)
            if (ierr /= 0) return
! The rotational shape factors f_P and f_T (fp_rot, ft_rot) for the
! current omega on the new points.
            call rotation_shape_factors(s%logRho,s%logR,s%log_mass,s%nz,s%omega,s%eta_squared,s%fp_rot,s%ft_rot,s%mean_gravity,s%mean_radius,ierr)
            if (ierr /= 0) return
            new_atmosphere_fit_needed = .false.
         endif
! The timestep of the next model (by its magnitude, as a stored model's
! timestep is taken).
       s%dt = dabs(s%dt)
       s%dt_saved = s%dt
       disable_structure_dt_limits = .false.
       call compute_timestep(s, s%dt,chosen_dt,s%nz,s%logRho,s%luminosity_lsun,s%m,s%logT,s%xa,s%core_cz_top_index, &
              s%h_shell_midpoint_zone,s%luminosity_breakdown,s%age_gyr,s%timestep_yr,s%job%nk,s%logP,s%logR,s%omega, &
              s%max_domega_frac,s%h_shell_zone_begin,s%log_Teff, disable_structure_dt_limits)
       if (disable_structure_dt_limits) s%job%use_structure_dt_limits = .false.
! With a fixed end age (target_end_age > 0) that this model reached,
! compute_timestep has clipped the next timestep to (nearly) zero. Write
! the timestep just taken (dt_saved, floored at 1e-3 of the age, dage
! in Gyr) to the model instead, so that the sequence can be continued.
! Other stops keep the new timestep. (The floor read 1e-3*dage years,
! dage taken as years, until 2026-10, as in main.f.)
       if (s%job%end_age_stop_active(s%job%nk) .and. s%job%target_end_age(s%job%nk).gt.0.0D0) then
          if (reached_end_age(s, s%job%nk)) then
             s%dt = max(s%dt_saved,1.0D-3*s%age_gyr*1.0D9*seconds_per_year)
             s%timestep_yr = s%dt/seconds_per_year
          else
             s%dt_saved = s%dt
          endif
       else
          s%dt_saved = s%dt
       endif
! A rescale-only card (rescale_kind = 2) does not count models.
       if (s%job%rescale_kind(s%job%nk).ne.2) s%model_number = s%model_number+1
! Assemble the whole converged star (interior, envelope and atmosphere)
! into the stitched arrays before the observables: the turnover-time
! observable walks them, star%pphot comes from this build, and the output
! writers below only copy; nothing downstream integrates anything.
       call build_stitched_model(s)
! The per-model observables, computed in the star layer (one procedure
! per family of quantities, core/diagnostics/observables.f90), including
! the renormalized luminosity_breakdown; the output below only reads
! them.
       call compute_observables(s, ierr)
       if (ierr /= 0) return
! Write the model's output (formerly wrtout).
       call output_write_model(s, ierr)
       if (ierr /= 0) return

! The stop conditions (MHP 10/24), in order: the end age, the table of
! central D/X/Y abundances, the structure limits (log L, Teff, log g,
! nu_max) and the star-calibration target radius.
      call check_stop_conditions(s, model_iteration, step_status)
      return
contains

! ---------------------------------------------------------------
! Per-step run-log bookkeeping: with ctrl%rewind_short_file (legacy
! LRWSH), rewind the run log (run_log_unit, log_output_file / run.log)
! so that it holds only the current model's solver diagnostics.
subroutine update_output_flags_for_step

          if (s%ctrl%rewind_short_file) then
             rewind(run_log_unit)
          endif
end subroutine update_output_flags_for_step

! ---------------------------------------------------------------
! If the previous attempt diverged, call read_starting_model (formerly
! starin) with the divergence flag set, which ends the run with ierr = 1
! (no timestep is cut; see the note at retry_step). Otherwise nothing
! happens here, except that star%punch_pending_flag is always set; with
! the helium-flash option, the next run card stores its starting model
! on last_model_unit when the flag is set.
subroutine reload_model_if_diverged
            if (s%model_diverged_flag) then
             call read_starting_model(s, s%timestep_yr, s%dt, chosen_dt, &
                  s%triangle_orientation, ikut_flag, s%stored_vertex_index, &
                  s%model_diverged_flag, s%recompute_envelope_triangle, s%job%nk, &
                  s%total_angular_momentum, &
                  s%total_rotational_ke, convective_velocity, &
                  s%job%mixture_weights, ierr)
             if (ierr /= 0) return
             call check_rotation_initialised(s, ierr)
             if (ierr /= 0) return
          endif
          s%punch_pending_flag = .true.
end subroutine reload_model_if_diverged

! ---------------------------------------------------------------
! The composition half of the step, when the model is to be aged: the
! mass change at the surface (change_mass), the start-of-step
! composition, the composition step (burning and mixing; the predictor
! pass of composition_step), and the new age. Sets ierr on failure.
subroutine advance_composition_and_age
! The model is aged when model_number >= 0, except on a rescale-only
! card (rescale_kind = 2), where it is aged only while the centre is
! cooler than log Tc = 6.6: pre-main-sequence models are rescaled and
! aged at the same time (MHP 7/98).
            evolve_model_flag = s%model_number.ge.0 .and. &
                 (s%job%rescale_kind(s%job%nk).ne.2 .or. s%logT(1).lt.6.6D0)
            new_atmosphere_fit_needed = .false.
            if (evolve_model_flag) then
! Mass loss or accretion at the surface. change_mass sets
! target_envelope_mass (log10 of m(nz)/M before the change) and, when it
! accretes, new_atmosphere_fit_needed, for the envelope rebuild after
! convergence.
               call change_mass(s, target_envelope_mass, new_atmosphere_fit_needed, ierr)
               if (ierr /= 0) return
! Store the composition at the start of the timestep (xa_start).
               num_species = n_species_basic
               if (s%job%use_extended_composition) num_species=n_species_extended
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
                    mixed_zone_bounds_no_overshoot, jerr)
               if (jerr /= 0) then
               ! Return the error to the caller.
                  ierr = jerr
                  return
               end if
! The new age, dage in Gyr.
             s%age_gyr = s%age_gyr + 1.0D-9*s%timestep_yr
            endif
end subroutine advance_composition_and_age

! ---------------------------------------------------------------
! Rezone the new model (rezone, formerly hpoint), which also keeps the
! start-of-step structure in star%*_start; with the helium-flash option
! there is no rezoning, and the start-of-step structure is copied here
! instead. Then store the start-of-step rotational kinetic energy and,
! with the light elements, their start-of-step state (the surface
! convection zone's base and the Li/Be burning rates). Sets ierr on
! failure.
subroutine rezone_or_snapshot
      use math_lib
! Rezone the new model, except with the helium-flash option.
          if (.not.s%ctrl%helium_flash_active) then
             call rezone(s, s%stored_vertex_index, s%reset_triangle, s%h_shell_zone_begin, &
                  s%has_h_shell, s%total_angular_momentum, &
                  s%total_rotational_ke, ierr)
             if (ierr /= 0) return
! The base zone of the surface convection zone, where deuterium burning
! is capped (compute_zone_energy, compute_burn_rates).
               s%jcz = s%envelope_cz_bottom_index
            else
! Keep the start-of-step structure, which the structure-change timestep
! limit (timestep_limit_structure, formerly ptime) compares against.
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
          endif
! The start-of-step rotational kinetic energy of each zone.
            if (s%job%rotation_active) then
               do i = 1,s%nz
                  s%kinetic_energy_rot_old(i) = s%kinetic_energy_rot(i)
               end do
            endif
! The start-of-step depth of the surface convection zone, for the Li/Be
! burn: its base zone, and the envelope overshoot depth alpha_ov H_P
! below it (stored as pressure_scale_height_start, though it is not H_P).
! cz_base_radius_prev = 0 tells locate_surface_mixing_base that no base
! radius has been stored yet this step.
            if (s%job%use_extended_composition) then
               s%cz_base_radius_prev = 0.0D0
               envelope_cz_zone_prev = s%envelope_cz_bottom_index
               s%envelope_overshoot_depth_start = envelope_overshoot_depth()
! The Li/Be burning rates of the start-of-step model.
               call store_start_li_be_rates(s)
            endif
end subroutine rezone_or_snapshot

! ---------------------------------------------------------------
! The structure solve, Henyey corrector levels 1-3 (solve_level): the
! first guess (the previous model or, with
! ctrl%improved_first_guess_flag, the previous model extrapolated by its
! rates of change); level 1 on the previous model's envelope triangle;
! level 2 after re-checking the surface boundary condition; the rates of
! change for the next model's first guess; level 3. Returns early with
! star%model_diverged_flag set when a level diverges (evolve_step then
! cycles retry_step); ierr on error.
subroutine solve_structure
! Flags for henyey_iterate: with ctrl%lnew0, a new envelope triangle is
! computed on the first iteration (the previous model's is ignored);
! converged and star%model_diverged_flag report the outcome.
          if (s%ctrl%lnew0) s%recompute_envelope_triangle = .true.
! A model that is solved but not aged (evolve_model_flag false) keeps
! its s%dt; the stages that age read the flag: the first guess below,
! the rates of change after level 2, the composition rerun and the
! Henyey pass (solve_level, which hands the pass a zero step), the
! light-element burn and the rotation step. (Until 2026-10 the flag was
! encoded as a negative s%dt, made positive again before
! compute_timestep.)
! The under-relaxation factor starts at |fcorr0| - fcorri. The sign of
! ctrl%fcorr0 selects whether henyey_iterate raises it by fcorri per
! iteration, up to 1 (fcorr0 > 0), or holds it fixed.
            s%job%fcorr = dabs(s%ctrl%fcorr0) - s%ctrl%fcorri
            iterations_done = 0
            s%model_diverged_flag = .false.
            converged = .false.
! Without the improved first guess, or for a model not aged, start from
! the previous model: zero the changes of log T and log P over the step
! (from which the energy equation builds the entropy term T dS/dt) and
! the rates of change of the previous step.
            if (.not.s%ctrl%improved_first_guess_flag .or. .not.evolve_model_flag &
                 .or. s%dt.le.0.0D0) then
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
! over the previous step (the *_entropy_term arrays, which are rates, not
! entropies):
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
            endif

! Level 1 uses the previous model's envelope triangle; for the first
! model of a run the triangle is generated here. The four levels differ
! only in their arguments to solve_level: the level, its iteration limit
! and whether the surface boundary condition is re-checked.
            call solve_level(1, s%ctrl%max_iter_level1, .false.)
            if (ierr /= 0) return
! Level 2 checks the envelope triangle (the surface boundary condition)
! before iterating.
            if (s%model_diverged_flag) return   ! (host cycles retry_step on the flag)
            call solve_level(2, s%ctrl%max_iter_level2, .true.)
            if (ierr /= 0) return
            if (s%model_diverged_flag) return   ! (host cycles retry_step on the flag)
! The rates of change of the structure over this step, from the level-2
! solution, for the improved first guess of the next model
! (ctrl%improved_first_guess_flag, legacy LNEWS; 7/91):
!   d log T/dt = log_temperature_delta/dt,  d log P/dt = log_pressure_delta/dt,
!   d ln L/dt = 2 (L - L_start)/(L + L_start)/dt,
!   d log R/dt = (log R - log R_start)/dt.
            if (evolve_model_flag .and. s%dt.gt.0.0D0) then
               do ii = 1,s%nz
                  s%dlog_temperature_dt(ii)=s%log_temperature_delta(ii)/s%dt
                  s%dlog_pressure_dt(ii)=s%log_pressure_delta(ii)/s%dt
                  s%dln_luminosity_dt(ii)=2.0D0*(s%luminosity_lsun(ii)-s%luminosity_lsun_start(ii))/(s%luminosity_lsun(ii)+s%luminosity_lsun_start(ii))/s%dt
                  s%dlog_radius_dt(ii)=(s%logR(ii)-s%logR_start(ii))/s%dt
               end do
            endif
! Level 3: from here on, solve_level redoes the composition step on the
! updated structure before iterating.
            call solve_level(3, s%ctrl%max_iter_level3, .false.)
            if (ierr /= 0) return
            if (s%model_diverged_flag) return   ! (host cycles retry_step on the flag)
end subroutine solve_structure

! ---------------------------------------------------------------
! The structure-rotation iteration (num_rotation_structure_iters
! passes; one without rotation). Each pass: corrector level 4, the
! convergence check, the convection zones mixed fully
! (homogenize_convection_zones, formerly mixcz), then the angular
! momentum update (evolve_angular_momentum, formerly getw) and the new
! shape factors (rotation_shape_factors, formerly fpft). The
! start-of-step composition and angular momentum are restored on repeat
! passes, so that nothing is counted twice.
subroutine converge_with_rotation
            if (.not.s%job%rotation_active) then
               s%job%num_rotation_structure_iters = 1
            endif
! When the structure and rotation are iterated (MHP 05/02), keep copies
! of the start-of-step specific angular momentum (j_rot) and composition
! (xa_start, legacy HCOMPP), so that every pass starts from them.
            if (s%job%num_rotation_structure_iters.gt.1) then
               do i = 1,s%nz
                  s%orig_specific_angular_momentum(i) = s%j_rot(i)
                  do j = 1,n_species_extended
                     s%orig_composition(j,i) = s%xa_start(j,i)
                  end do
               end do
            endif
            do itrot = 1, s%job%num_rotation_structure_iters
! On repeat passes, restore the start-of-step composition (MHP 05/02).
               if (itrot.gt.1) then
                  do i = 1,s%nz
                     do j = 1,n_species_extended
                        s%xa_start(j,i) = s%orig_composition(j,i)
                     end do
                  end do
               endif
! Level 4 repeats the iteration between the mixing (the composition
! step, redone in solve_level) and the structure (7/91). Use it for
! high-precision work such as solar models, but not where semiconvection
! matters: iterating between the burning and the structure then
! oscillates. Unlike levels 1-3 it runs inside the structure-rotation
! iteration. The surface boundary condition is not re-checked here
! (check_surface_bc = .false.), although level 3 changed the composition
! and hence the structure; the legacy comment left that check for later.
            call solve_level(4, s%ctrl%max_iter_level4, .false.)
            if (ierr /= 0) return
            if (.not.converged) then
! Not converged within the iterations of levels 1-4 (max_iter_level1 to
! max_iter_level4): flag the divergence.
               s%model_diverged_flag = .true.
               return   ! (host cycles retry_step on the flag)
            endif

! The model has converged. Mix the convection zones fully.
      call homogenize_convection_zones(s%xa,s%dm,s%convective_flag,s%nz)

! The age (Gyr) that the disk-locking test reads.
            s%disk_gate_age_gyr = s%age_gyr
            if (s%job%rotation_active) then
! On repeat passes, restore the start-of-step angular momentum.
               if (itrot.gt.1) then
                  do i = 1,s%nz
                     s%j_rot(i) = s%orig_specific_angular_momentum(i)
                  end do
               endif
! Evolve the angular momentum over the step: the new run of omega (with
! the wind torque when job%use_wind_torque is set). A model not aged
! passes its step negated, as the negative-dt encoding did: the rotation
! step skips its sub-steps for a non-positive step, but without
! instability transport the wind brake still applies it, which spins the
! surface convection zone up (inherited from GETW; docs/bug_sweep.md O24,
! open).
               wind_loss_active = s%job%use_wind_torque
               call evolve_angular_momentum(s, merge(s%dt, -dabs(s%dt), evolve_model_flag), &
                    s%max_domega_frac, wind_loss_active, &
                    envelope_cz_zone_prev, jerr)
               if (jerr /= 0) then
               ! Return the error to the caller.
                  ierr = jerr
                  return
               end if
! The rotational shape factors f_P and f_T for the new omega.
               call rotation_shape_factors(s%logRho,s%logR,s%log_mass,s%nz,s%omega,s%eta_squared,s%fp_rot,s%ft_rot,s%mean_gravity,s%mean_radius,ierr)
               if (ierr /= 0) return
            endif
            end do
end subroutine converge_with_rotation

! ---------------------------------------------------------------
! The Li6, Li7 and Be9 burn over the step, on the converged model, for
! models without rotation or without rotational mixing (with it,
! evolve_angular_momentum burns them in its sub-steps): the end-of-step
! convection zones (find_convection_zones, formerly convec) and overshoot
! depth, the end-of-step burning rates (store_end_li_be_rates), then the
! burn (burn_li_be). Deuterium burns in the composition step instead.
subroutine burn_light_elements
      use math_lib
! Only with the extended composition, for a model being aged.
         if (s%job%use_extended_composition .and. s%model_number.ge.0 .and. &
              evolve_model_flag .and. s%dt.gt.0.0D0) then
! Only without rotation, or without rotational mixing.
            if (.not.s%job%rotation_active .or. .not.s%job%instability_transport_active) then
! The convection zones at the end of the step.
               call find_convection_zones(s%xa,s%logRho,s%logP,s%logR,s%log_mass,s%logT,s%convective_flag,s%nz,radiative_zone_bounds,s%mixed_zone_bounds, &
                            mixed_zone_bounds_no_overshoot,s%core_cz_top_index,s%envelope_cz_bottom_index,num_radiative_zones,num_mixed_zones,num_mixed_zones_no_overshoot, s%iov1, s%iov2, ierr)
               if (ierr /= 0) return
! The end-of-step base zone of the surface convection zone and the
! envelope overshoot depth below it (alpha_ov H_P).
               envelope_cz_zone_end = s%envelope_cz_bottom_index
               s%envelope_overshoot_depth_end = envelope_overshoot_depth()
! The burning rates at the end of the step, from the converged model;
! then the burn.
               call store_end_li_be_rates(s, s%xa,s%logRho,s%logT,s%nz)
               call burn_li_be(s, s%dt,s%xa,s%logR,s%m,s%dm,s%logT,s%grada,s%gradr, &
                    envelope_cz_zone_end,envelope_cz_zone_prev,s%nz)
            endif
         endif
end subroutine burn_light_elements


! ---------------------------------------------------------------
! One level of the corrector ladder. Every level runs the same three
! stages; only the level number, its iteration limit and whether the
! surface boundary condition is re-checked differ:
!   1. refresh_surface_bc: the surface boundary condition (before the
!      composition step, because it reads the surface composition);
!   2. at levels >= 3, for a model being aged: the composition step
!      again (corrector pass), from the start-of-step composition, on
!      the updated structure;
!   3. henyey_iterate: the Newton iterations.
! Everything else comes from the host's locals and star, by host
! association. The caller checks ierr and the divergence flag after each
! level. A level with an iteration limit <= 0 is skipped.
subroutine solve_level(level, level_max_iterations, check_surface_bc)
      use star_info_lib, only: max_convective_zones
      integer, intent(in) :: level, level_max_iterations
      logical, intent(in) :: check_surface_bc
      integer :: saha_state, envelope_zone_index, num_species, i, j
      logical :: envelope_recomputed_flag
! The zone outputs of the redone composition step. Only the envelope
! edge is used (by the Newton iterations, below); the others are
! discarded, so star's copies keep their values from the predictor
! pass.
      integer :: core_cz_edge_rerun
      integer :: mixed_zone_bounds_rerun(max_convective_zones,2)

      ierr = 0
      iteration_level = level
      max_iterations = level_max_iterations
      recompute_surface_bc = check_surface_bc
      if (max_iterations.le.0) return

      call refresh_surface_bc(s, iteration_level, s%recompute_envelope_triangle, &
           s%reset_triangle, recompute_surface_bc, s%triangle_orientation, &
           s%stored_vertex_index, saha_state, envelope_recomputed_flag, ierr)
      if (ierr /= 0) return

! envelope_zone_index gates the accretion correction to the entropy term
! of the energy equation (compute_henyey_coefficients adds it from that
! zone out). It is the redone composition step's envelope edge at levels
! >= 3, and 0 (every zone) at levels 1-2. That reproduces the legacy
! code, where it was an unset local zeroed by -finit-local-zero
! (docs/ROADMAP.md section 13).
      envelope_zone_index = 0
      core_cz_edge_rerun = 0
      mixed_zone_bounds_rerun = 0
      if (iteration_level.gt.2 .and. evolve_model_flag .and. s%dt.gt.0.0d0) then
         num_species = n_species_basic
         if (s%job%use_extended_composition) num_species = n_species_extended
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
      call henyey_iterate(s, merge(s%dt, 0.0d0, evolve_model_flag), max_iterations, converged, &
           s%model_diverged_flag, recompute_surface_bc, &
           envelope_recomputed_flag, saha_state, envelope_zone_index, &
           iterations_done, ierr)
end subroutine solve_level

! ---------------------------------------------------------------
! The envelope-overshoot depth below the base of the surface convection
! zone, in cm: d_ov = alpha_ov H_P, with alpha_ov =
! ctrl%overshoot_alpha_envelope and H_P the pressure scale height at the
! base zone (envelope_cz_bottom_index); zero when envelope overshoot is
! off. The same formula serves the start and the end of the step
! (rezone_or_snapshot, burn_light_elements).
double precision function envelope_overshoot_depth()
      use math_lib
      if (s%job%envelope_overshoot_active) then
         envelope_overshoot_depth = s%ctrl%overshoot_alpha_envelope* &
              pressure_scale_height(s%logP(s%envelope_cz_bottom_index), &
              s%logR(s%envelope_cz_bottom_index), &
              s%logRho(s%envelope_cz_bottom_index), &
              s%log_mass(s%envelope_cz_bottom_index))
      else
         envelope_overshoot_depth = 0.0D0
      endif
end function envelope_overshoot_depth

end subroutine evolve_step

end module step_driver

!----------------------------------------------------------------------
! envelope_integration
!----------------------------------------------------------------------
! The star-layer driver of the outer layers' integration. atm_get
! integrates the atmosphere and the envelope from the surface down to the
! fitting point: it fills the kernel's settings (envint_config) from the
! star's controls and atm_choice_in_use, calls the star-blind kernel
! (envelope_integration_kernel.f90, whose integrands in
! envelope_integrands.f90 call the eos and kap facades, through numerics'
! bsstep), applies the integration's two effects on the star
! (star%photosphere_logP, and the gray atmosphere after a failed Allard
! lookup), and reports the kernel's events (report_envint_events). The
! points it integrates go to the star's env_points, atm_points and
! atm_point. The tables and their surface lookups are atm's. The helpers
! fixed_envint_step, surface_log10_radius_cm and ambient_log10_teff serve
! atm_get's callers (compute_surface_boundary, rebuild_envelope,
! read_starting_model, build_stitched_model).
module envelope_integration
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star_info
      implicit none
      private
      public :: atm_get, fixed_envint_step, surface_log10_radius_cm, &
           ambient_log10_teff

! A Bulirsch-Stoer step-size triple (begin, min, max) a driver can hand
! atm_get in place of the controls' {atm,env}_step_* values, for that
! call only (rebuild_envelope, read_starting_model, build_stitched_model).
      type, public :: envint_step_config
            real(dp) :: step_begin, step_min, step_max
      end type envint_step_config

contains

! all three step sizes equal to h (the only pattern the drivers use)
pure function fixed_envint_step(h) result(steps)
      real(dp), intent(in) :: h
      type(envint_step_config) :: steps
      steps%step_begin = h
      steps%step_min = h
      steps%step_max = h
end function fixed_envint_step

! log10 of the surface radius (cm) from log10 L/Lsun and log10 Teff,
! L = 4 pi R**2 sigma Teff**4, as rebuild_envelope, compute_surface_boundary
! and build_stitched_model derive it before integrating. The operand
! order is theirs; physical_constants' log10_radius_from_l_teff adds the
! - 4 log Teff term last, which can differ in the last bit, so the two
! forms are kept apart (merging them moves results).
!   star reads: log10_solar_luminosity
!   star writes: none
real(dp) function surface_log10_radius_cm(s, log_luminosity_lsun, log_teff)
      use physical_constants, only: c4pil, csigl
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: log_luminosity_lsun, log_teff
      surface_log10_radius_cm = 0.5d0*(log_luminosity_lsun + s%log10_solar_luminosity &
           - 4.0d0*log_teff - c4pil - csigl)
end function surface_log10_radius_cm

! log10 Teff of the unspotted ("ambient") surface for a spotted run
! (G Somers 10/14): the caller has checked that spots are on
! (spot_filling_factor /= 0, spot_temp_contrast /= 1) and that the
! surface shell is convective. The single-precision literals
! (0.25, 4.0, 1.0) are the original ones; build_stitched_model keeps
! its own copy with d0 literals, which is not token-identical.
!   star reads: physics%spot_filling_factor, physics%spot_temp_contrast
!   star writes: none
real(dp) function ambient_log10_teff(s, log_teff)
      use math_lib
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: log_teff
      ambient_log10_teff = log_teff - 0.25*log10(s%physics%spot_filling_factor * &
           pow(s%physics%spot_temp_contrast, 4.0) + 1.0 - s%physics%spot_filling_factor)
end function ambient_log10_teff

!   star reads: log10_solar_luminosity, log_fit_point_mass_fraction,
!     log_total_mass, physics%atm_error_tol, physics%atm_step_begin,
!     physics%atm_step_initial, physics%atm_step_max,
!     physics%atm_step_min, physics%calc_envelope_flag,
!     physics%env_error_tol, physics%env_step_begin,
!     physics%env_step_max, physics%env_step_min,
!     physics%envelope_store_log_temperature, physics%extrap_order,
!     physics%max_stage_index, physics%tolerance_fraction
!   star writes: atm_choice_in_use, atm_point, atm_points, env_points,
!     photosphere_logP
subroutine atm_get(s, luminosity_linear, pressure_rotation_factor, &
     temperature_rotation_factor, log10_gravity, log10_star_mass, &
     print_flag, save_boundary_flag, log10_pressure_limit, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     surface_composition, saha_iterations, ierr, &
     vertex_index, stored_envelope_state, stored_vertex_index, &
     vtx_logp, vtx_logr, vtx_logt, atm_steps, env_steps)

! atm_get: one integration of the atmosphere and the envelope for the
! star, at L (Lsun, luminosity_linear), f_P and f_T, log g, log M, log R,
! log Teff and the surface X, Z and composition. The fitting-vertex
! arguments (vertex_index, stored_envelope_state, stored_vertex_index,
! vtx_logp/logr/logt) are optional and go together: only
! compute_surface_boundary passes them, with save_boundary_flag; the
! other drivers call without them and save_boundary_flag false, under
! which the kernel never writes them, so they get zeroed stand-ins here.
! atm_steps and env_steps, when present, replace the controls' step
! sizes for this call. print_flag adds the fitted vertex to the solver
! diagnostics.
      use envelope_integration_kernel, only: envint_config, envint_events, &
           integrate_envelope_atmosphere
      use physical_constants, only: cgl
      use atm_lib, only: hsra_t_tau_offset, atm_eddington_gray
      implicit none
      type(star_info), intent(inout) :: s

      real(dp), intent(inout) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, &
           log10_gravity
      real(dp), intent(in) :: log10_star_mass
      logical, intent(in) :: print_flag
      logical, intent(in) :: save_boundary_flag
      real(dp), intent(in) :: log10_pressure_limit
      real(dp), intent(inout) :: log10_radius
      real(dp), intent(inout) :: log10_teff
      real(dp), intent(inout) :: hydrogen_fraction, metal_fraction
      ! surface_composition: the outer shell's composition (star%xa layout),
      ! for the EOS's Debye-Hueckel correction.
      real(dp), intent(in) :: surface_composition(:)
      integer, intent(inout) :: saha_iterations
      integer, intent(out) :: ierr
      integer, intent(in), optional :: vertex_index
      real(dp), intent(inout), optional :: stored_envelope_state(:)
      integer, intent(inout), optional :: stored_vertex_index
      real(dp), intent(inout), optional :: vtx_logp(:), vtx_logr(:), &
           vtx_logt(:)
      type(envint_step_config), intent(in), optional :: atm_steps, env_steps

      type(envint_config) :: cfg
      type(envint_events) :: events
      logical :: switched_to_gray
      real(dp) :: log10_photo_pressure
! scratch stand-ins for the absent vertex group (never read back)
      real(dp) :: scratch_state(4), scratch_logp(3), scratch_logr(3), &
           scratch_logt(3)
      integer :: scratch_vertex_index, scratch_stored_vertex_index

      cfg%atm_choice = s%atm_choice_in_use
      if (present(atm_steps)) then
         cfg%atm_step_begin = atm_steps%step_begin
         cfg%atm_step_min = atm_steps%step_min
         cfg%atm_step_max = atm_steps%step_max
      else
         cfg%atm_step_begin = s%physics%atm_step_begin
         cfg%atm_step_min = s%physics%atm_step_min
         cfg%atm_step_max = s%physics%atm_step_max
      end if
      if (present(env_steps)) then
         cfg%env_step_begin = env_steps%step_begin
         cfg%env_step_min = env_steps%step_min
         cfg%env_step_max = env_steps%step_max
      else
         cfg%env_step_begin = s%physics%env_step_begin
         cfg%env_step_min = s%physics%env_step_min
         cfg%env_step_max = s%physics%env_step_max
      end if
      cfg%atm_error_tol = s%physics%atm_error_tol
      cfg%env_error_tol = s%physics%env_error_tol
      cfg%atm_step_initial = s%physics%atm_step_initial
! the HSRA T(tau) fit at tau = 2/3, for atm_choice = 2
      cfg%hsra_log_t_offset = hsra_t_tau_offset()
      cfg%log10_total_mass = s%log_total_mass
      cfg%surface_composition = surface_composition
      cfg%log_fit_point_mass_fraction = s%log_fit_point_mass_fraction
      cfg%envelope_store_log_temperature = s%physics%envelope_store_log_temperature
      cfg%log10_solar_luminosity = s%log10_solar_luminosity
      cfg%log10_gravitational_constant = cgl
      cfg%calc_envelope = s%physics%calc_envelope_flag
      cfg%tolerance_fraction = s%physics%tolerance_fraction
      cfg%max_stage_index = s%physics%max_stage_index
      cfg%extrap_order = s%physics%extrap_order

      if (present(vertex_index)) then
         call run_kernel(vertex_index, stored_envelope_state, &
              stored_vertex_index, vtx_logp, vtx_logr, vtx_logt)
      else
         scratch_vertex_index = 0
         scratch_stored_vertex_index = 0
         scratch_state = 0.0d0
         scratch_logp = 0.0d0
         scratch_logr = 0.0d0
         scratch_logt = 0.0d0
         call run_kernel(scratch_vertex_index, scratch_state, &
              scratch_stored_vertex_index, scratch_logp, scratch_logr, &
              scratch_logt)
      end if
! The kernel's notes, diagnostics and failure, as it used to print them.
      call report_envint_events(events, switched_to_gray, print_flag)

! star-state consequences of the integration, applied HERE and only
! here (the kernel never sees star_info):
      if (log10_photo_pressure > -1.0d98) s%photosphere_logP = log10_photo_pressure
      if (switched_to_gray) then
! Allard table lookup out of range: the kernel fell back to a gray
! atmosphere for this and all future calls (historical behavior)
         s%atm_choice_in_use = atm_eddington_gray
      end if

contains

      subroutine run_kernel(vidx, state, sidx, logp, logr, logt)
            integer, intent(in) :: vidx
            real(dp), intent(inout) :: state(:)
            integer, intent(inout) :: sidx
            real(dp), intent(inout) :: logp(:), logr(:), logt(:)
            call integrate_envelope_atmosphere(cfg, switched_to_gray, &
                 log10_photo_pressure, &
                 luminosity_linear, pressure_rotation_factor, &
                 temperature_rotation_factor, log10_gravity, log10_star_mass, &
                 vidx, save_boundary_flag, log10_pressure_limit, &
                 log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
                 state, sidx, saha_iterations, logp, logr, logt, &
                 s%env_points, s%atm_points, s%atm_point, events, ierr)
      end subroutine run_kernel

end subroutine atm_get

! ---------------------------------------------------------------
! The messages of one envelope integration (integrate_envelope_atmosphere
! prints nothing): the switch to a gray atmosphere after a failed Allard
! lookup, to standard output and the run log; then the kernel's events
! in order. The completed integrations' step counts and errors go to the
! run log with the solver diagnostics on, and the fitted vertex as well
! when print_flag asks for it; a starting optical depth past the
! photosphere goes to the run log; the integrations' failures go to the
! terminal and the run log, except a failed Bulirsch-Stoer step, whose
! status goes to standard output.
subroutine report_envint_events(events, switched_to_gray, print_flag)
      use envelope_integration_kernel, only: envint_events, &
           envint_start_retry, envint_atmosphere_done, &
           envint_atmosphere_maxstp, envint_bsstep_limits, &
           envint_bsstep_failed, envint_vertex, envint_envelope_done, &
           envint_envelope_maxstp
      use core_messages, only: core_report, core_log, core_stdout, &
           core_stdout_and_log, log_envelope_start_retry, message_length, &
           unwritten_line
      use run_messages, only: solver_diagnostics
      use numerics_lib, only: numerics_status_text, bsstep_max_extrap_order, &
           bsstep_max_stages
      type(envint_events), intent(in) :: events
      logical, intent(in) :: switched_to_gray, print_flag
      character(len=message_length) :: lines(3)
      integer :: k

      if (switched_to_gray) then
         lines = unwritten_line
         lines(1) = ' note: Allard table lookup out of range;' &
              // ' switching to a gray atmosphere'
         call core_stdout_and_log(lines)
      end if
      do k = 1, events%num_events
         lines = unwritten_line
         associate (e => events%event(k))
         select case (e%kind)
         case (envint_start_retry)
            call log_envelope_start_retry(e%values(1), e%values(2))
         case (envint_atmosphere_done)
            if (solver_diagnostics()) then
               write(lines,35) e%counts(1), e%counts(2), e%values(1)
   35          format(1X,'ATMOSPHERE INTEGRATION COMPLETE',1X, &
                    'NUMBER OF STEPS ACCEPTED',I5,' REJECTED', &
                    I5/5X,'MAXIMUM RELATIVE ERROR IN P ',1PE22.13)
               call core_log(lines)
            end if
         case (envint_atmosphere_maxstp)
            write(lines,50)
   50       format(5X,'ATMOSPHERE INTEGRATION FAILED AFTER MAXSTP',1X, &
                 'INTEGRATIONS.I QUIT.')
            call core_report(lines)
         case (envint_bsstep_limits)
            write(lines,'(1x,a,i0,a,i0,a,i0,a,i0)') &
                 'ENVINT: Bulirsch-Stoer limits out of range: imax = ', &
                 e%counts(1), ', nuse = ', e%counts(2), '; need 2 <= nuse <= ', &
                 bsstep_max_extrap_order, ' and nuse <= imax <= ', bsstep_max_stages
            call core_report(lines)
         case (envint_bsstep_failed)
            lines(1) = ' bsstep: '//numerics_status_text(e%counts(1))
            call core_stdout(lines)
         case (envint_vertex)
            if (print_flag .and. solver_diagnostics()) then
               write(lines,230) e%values(1:4)
 230           format(4X,3F16.12,8X,F16.12)
               call core_log(lines)
            end if
         case (envint_envelope_done)
            if (solver_diagnostics()) then
               write(lines,215) e%counts(1), e%counts(2), e%values(1:5)
 215           format(1X,'ENVELOPE INTEGRATION COMPLETE',1X, &
                    'NUMBER OF STEPS ACCEPTED',I5,' REJECTED', &
                    I5/5X,'SENV-LAST M=',1PE22.13,'  LAST M=',E22.13/ &
                    5X,'MAX RELATIVE ERRORS:M ',1PE14.5,'  T ',E14.5, &
                    '  R ',E14.5)
               call core_log(lines)
            end if
         case (envint_envelope_maxstp)
            write(lines,911)
 911        format(5X,'ENVELOPE INTEGRATION FAILED AFTER MAXSTP TRIES.',1X, &
                 'I QUIT')
            call core_report(lines)
         end select
         end associate
      end do
end subroutine report_envint_events

end module envelope_integration

! Sources
! - The integration was atm_lib's until the physics-purity pass of 2026
!   (docs/ROADMAP.md, the atm split), which moved it into core as
!   star-solver code (MESA's create_atm and envelope integration), with
!   the kernel split from this driver (envint purity split); atm keeps
!   the tables and their surface lookups.
! - The step-size triple replaced a save, override and restore of the
!   controls around the call in rebuild_envelope, read_starting_model
!   and build_stitched_model (2026, readability wave 2).
! - The kernel prints nothing since round 2 (core_messages; this driver
!   reports its events).

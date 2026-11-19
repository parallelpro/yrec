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
! atm_point. gray_photosphere integrates the Eddington gray atmosphere
! alone, at a given Teff and g, into its own points (the accretion
! surface of the retired accretion-entropy path, after a failed Allard
! lookup: core/step/retired/accretion_entropy.f90.txt; no caller in the
! run since, bug_sweep M323). The tables and their
! surface lookups are atm's. The helpers fixed_envint_step,
! surface_log10_radius_cm and ambient_log10_teff serve atm_get's callers
! (compute_surface_boundary, rebuild_envelope, read_starting_model,
! build_stitched_model).
module envelope_integration
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star_info
      use envelope_integration_kernel, only: fitting_vertex
      implicit none
      private
      public :: atm_get, gray_photosphere, fixed_envint_step, &
           surface_log10_radius_cm, ambient_log10_teff
! the vertex record atm_get takes (the kernel's), for its caller
      public :: fitting_vertex

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
! and build_stitched_model derive it before integrating (and the stitched
! model's atmosphere takes it as its base when no envelope point is
! stored). The operand
! order is theirs; physical_constants' log10_radius_from_l_teff adds the
! - 4 log Teff term last, which can differ in the last bit, so the two
! forms are kept apart (merging them moves results).
!   star reads: log10_solar_luminosity
!   star writes: none
real(dp) function surface_log10_radius_cm(s, log_luminosity_lsun, log_teff)
      use physical_constants, only: c4pil, csigl
      type(star_info), intent(in) :: s
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
      type(star_info), intent(in) :: s
      real(dp), intent(in) :: log_teff
      ambient_log10_teff = log_teff - 0.25*log10(s%physics%spot_filling_factor * &
           pow(s%physics%spot_temp_contrast, 4.0) + 1.0 - s%physics%spot_filling_factor)
end function ambient_log10_teff

! ---------------------------------------------------------------
! atm_get
! ---------------------------------------------------------------
! One integration of the atmosphere and the envelope for the star, at
! the surface the caller gives: L (Lsun, luminosity_linear), the
! rotational distortion factors f_P and f_T, log g, log R, log Teff, X,
! Z and the outer shell's composition (star%xa layout, for the EOS's
! Debye-Hueckel correction). log10_pressure_limit is the guess at the
! fitting point's pressure (the kernel's header); saha_iterations is the
! EOS's iteration counter. The mass is the star's (log_total_mass).
! vertex, when present, is the vertex of the surface boundary's triangle
! this integration fits (the kernel's fitting_vertex, which
! compute_surface_boundary fills from the star's triangle and stores
! back): the integration writes its fitting point and stored envelope
! state. atm_steps and env_steps, when present, replace the controls'
! step sizes for this call. ierr is the kernel's.
!   star reads: none
!   star writes: atm_choice_in_use, atm_point, atm_points, env_points,
!     photosphere_logP
subroutine atm_get(s, luminosity_linear, pressure_rotation_factor, &
     temperature_rotation_factor, log10_gravity, log10_radius, log10_teff, &
     hydrogen_fraction, metal_fraction, surface_composition, &
     log10_pressure_limit, saha_iterations, ierr, vertex, atm_steps, env_steps)
      use envelope_integration_kernel, only: envint_config, envint_events, &
           integrate_envelope_atmosphere, no_photosphere_log10_pressure
      use atm_lib, only: atm_eddington_gray
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: luminosity_linear, pressure_rotation_factor, &
           temperature_rotation_factor, log10_gravity, log10_radius, log10_teff
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(in) :: surface_composition(:)
      real(dp), intent(in) :: log10_pressure_limit
      integer, intent(inout) :: saha_iterations
      integer, intent(out) :: ierr
      type(fitting_vertex), intent(inout), optional :: vertex
      type(envint_step_config), intent(in), optional :: atm_steps, env_steps

      type(envint_config) :: cfg
      type(envint_events) :: events
      logical :: switched_to_gray
      real(dp) :: log10_photo_pressure
      integer :: vertex_index

      cfg = envint_config_from_star(s, surface_composition)
      if (present(atm_steps)) then
         cfg%atm_step_begin = atm_steps%step_begin
         cfg%atm_step_min = atm_steps%step_min
         cfg%atm_step_max = atm_steps%step_max
      end if
      if (present(env_steps)) then
         cfg%env_step_begin = env_steps%step_begin
         cfg%env_step_min = env_steps%step_min
         cfg%env_step_max = env_steps%step_max
      end if

      call integrate_envelope_atmosphere(cfg, luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, &
           log10_gravity, log10_radius, log10_teff, hydrogen_fraction, &
           metal_fraction, log10_pressure_limit, saha_iterations, &
           s%env_points, s%atm_points, s%atm_point, switched_to_gray, &
           log10_photo_pressure, events, ierr, vertex)
! The kernel's notes, diagnostics and failure, with the vertex being
! fitted (0: none).
      vertex_index = 0
      if (present(vertex)) vertex_index = vertex%index
      call report_envint_events(events, switched_to_gray, vertex_index)

! The integration's two effects on the star, applied here and only here
! (the kernel never sees star_info): the photospheric pressure, when the
! atmosphere finished, and the gray atmosphere after an Allard lookup
! out of the tables' range, which stays the star's working choice until
! switch_atmosphere_choice (surface_boundary) returns to the tables.
      if (log10_photo_pressure > no_photosphere_log10_pressure) &
           s%photosphere_logP = log10_photo_pressure
      if (switched_to_gray) s%atm_choice_in_use = atm_eddington_gray
end subroutine atm_get

! ---------------------------------------------------------------
! The kernel's settings as the star's controls give them (envint_config):
! the atmosphere choice in use, the controls' step sizes and
! tolerances, the fitting mass, and calc_envelope_flag. atm_get and
! gray_photosphere change what their call needs.
!   star reads: atm_choice_in_use, log10_solar_luminosity,
!     log_fit_point_mass_fraction, log_total_mass,
!     physics%atm_error_tol, physics%atm_step_begin,
!     physics%atm_step_initial, physics%atm_step_max,
!     physics%atm_step_min, physics%calc_envelope_flag,
!     physics%env_error_tol, physics%env_step_begin,
!     physics%env_step_max, physics%env_step_min,
!     physics%envelope_store_log_temperature, physics%extrap_order,
!     physics%max_stage_index, physics%tolerance_fraction
!   star writes: none
function envint_config_from_star(s, surface_composition) result(cfg)
      use envelope_integration_kernel, only: envint_config
      use physical_constants, only: cgl
      use atm_lib, only: hsra_t_tau_offset
      type(star_info), intent(in) :: s
      real(dp), intent(in) :: surface_composition(:)
      type(envint_config) :: cfg

      cfg%atm_choice = s%atm_choice_in_use
      cfg%atm_step_begin = s%physics%atm_step_begin
      cfg%atm_step_min = s%physics%atm_step_min
      cfg%atm_step_max = s%physics%atm_step_max
      cfg%env_step_begin = s%physics%env_step_begin
      cfg%env_step_min = s%physics%env_step_min
      cfg%env_step_max = s%physics%env_step_max
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
end function envint_config_from_star

!----------------------------------------------------------------------
! gray_photosphere
!----------------------------------------------------------------------
! log P and log T at the photosphere (tau = 2/3, where T = Teff) of the
! Eddington gray atmosphere (atm_choice 0) at log Teff and log g, for
! the surface X, Z and composition: one integration of the atmosphere
! alone (calc_envelope off), with the controls' atmosphere steps and
! tolerance and no rotation (f_P = f_T = 1). Its points go to locals:
! the star's atmosphere and envelope points, atm_point,
! photosphere_logP and atm_choice_in_use are not written. The
! integration's events are reported as atm_get's are; ierr is the
! kernel's.
!   star reads: log_L_lsun
!   star writes: none
subroutine gray_photosphere(s, log10_teff, log10_gravity, &
     hydrogen_fraction, metal_fraction, surface_composition, &
     log10_pressure, log10_temperature, ierr)
      use envelope_integration_kernel, only: envint_config, envint_events, &
           integrate_envelope_atmosphere
      use envelope_points, only: envelope_structure_state, &
           atmosphere_structure_state, atmosphere_point_state
      use atm_lib, only: atm_eddington_gray
      use math_lib
      type(star_info), intent(in) :: s
      real(dp), intent(in) :: log10_teff, log10_gravity
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(in) :: surface_composition(:)
      real(dp), intent(out) :: log10_pressure, log10_temperature
      integer, intent(out) :: ierr
! No rotational distortion, and no guess at the fitting point's
! pressure (no envelope is integrated).
      real(dp), parameter :: no_distortion = 1.0_dp, unused_pressure_limit = 0.0_dp
      type(envint_config) :: cfg
      type(envint_events) :: events
! the integration's points, on the heap for their size
      type(envelope_structure_state), allocatable :: env_points
      type(atmosphere_structure_state), allocatable :: atm_points
      type(atmosphere_point_state) :: atm_point
      logical :: switched_to_gray
      real(dp) :: log10_photo_pressure, luminosity, radius
      integer :: saha_iterations

      cfg = envint_config_from_star(s, surface_composition)
      cfg%atm_choice = atm_eddington_gray
      cfg%calc_envelope = .false.
      allocate(env_points, atm_points)
      env_points%num_env_points = 0
      atm_points%num_atm_points = 0
      luminosity = exp10(s%log_L_lsun)
      radius = surface_log10_radius_cm(s, s%log_L_lsun, log10_teff)
      saha_iterations = 0
      log10_pressure = 0.0d0
      log10_temperature = 0.0d0
      call integrate_envelope_atmosphere(cfg, luminosity, no_distortion, &
           no_distortion, log10_gravity, radius, log10_teff, &
           hydrogen_fraction, metal_fraction, unused_pressure_limit, &
           saha_iterations, env_points, atm_points, atm_point, &
           switched_to_gray, log10_photo_pressure, events, ierr)
      call report_envint_events(events, switched_to_gray, 0)
      if (ierr /= 0) return
      log10_pressure = atm_point%log10_pressure
      log10_temperature = atm_point%log10_temperature
end subroutine gray_photosphere

! ---------------------------------------------------------------
! The messages of one envelope integration (integrate_envelope_atmosphere
! prints nothing): the switch to a gray atmosphere after a failed Allard
! lookup, to the terminal and the run log; then the kernel's events in
! order. The completed integrations' step counts and errors go to the
! run log with the solver diagnostics on; a starting optical depth past
! the photosphere goes to the run log; the integrations' failures (an
! atmosphere or an envelope out of steps, an atmosphere out of starting
! densities, a failed Bulirsch-Stoer step) go to the terminal and the
! run log, with the vertex being fitted (vertex, 0 when none). The fitted vertex (envint_vertex) has no line: the
! triangle's diagnostics (compute_surface_boundary's) list every vertex.
subroutine report_envint_events(events, switched_to_gray, vertex)
      use envelope_integration_kernel, only: envint_events, &
           envint_start_retry, envint_atmosphere_done, &
           envint_atmosphere_maxstp, envint_bsstep_limits, &
           envint_bsstep_failed, envint_vertex, envint_envelope_done, &
           envint_envelope_maxstp, envint_start_exhausted
      use core_messages, only: core_report, core_log, message_length, &
           unwritten_line
      use run_messages, only: solver_diagnostics
      use numerics_lib, only: numerics_status_text, bsstep_max_extrap_order, &
           bsstep_max_stages
      type(envint_events), intent(in) :: events
      logical, intent(in) :: switched_to_gray
      integer, intent(in) :: vertex
      character(len=message_length) :: lines(3)
      character(len=40) :: where
      integer :: k

      if (switched_to_gray) then
         lines = unwritten_line
         lines(1) = ' note: Allard table lookup out of range;' &
              // ' switching to a gray atmosphere'
         call core_report(lines)
      end if
      where = ''
      if (vertex > 0) write(where,'(a,i0,a)') ' (vertex ', vertex, ')'
      do k = 1, events%num_events
         lines = unwritten_line
         associate (e => events%event(k))
         select case (e%kind)
         case (envint_start_retry)
            write(lines,'(1x,a,f9.4,a,f9.4,a)') 'envelope integration: the atmosphere''s' &
                 //' starting log tau, ', e%values(1), ', lies past the photosphere''s, ', &
                 e%values(2), '; it starts again from a tenth of the density'
            call core_log(lines)
         case (envint_atmosphere_done)
            if (solver_diagnostics()) then
               write(lines,35) e%counts(1), e%counts(2), e%values(1)
   35          format(1X,'ATMOSPHERE INTEGRATION COMPLETE',1X, &
                    'NUMBER OF STEPS ACCEPTED',I5,' REJECTED', &
                    I5/5X,'MAXIMUM RELATIVE ERROR IN P ',1PE22.13)
               call core_log(lines)
            end if
         case (envint_atmosphere_maxstp)
            write(lines,'(1x,a,a,a,i0,a)') 'ERROR: envelope integration', trim(where), &
                 ': the atmosphere did not reach the photosphere in ', e%counts(1), ' steps'
            call core_report(lines)
         case (envint_bsstep_limits)
            write(lines,'(1x,a,i0,a,i0,a,i0,a,i0)') &
                 'ERROR: envelope integration: the Bulirsch-Stoer limits are out of range: imax = ', &
                 e%counts(1), ', nuse = ', e%counts(2), '; need 2 <= nuse <= ', &
                 bsstep_max_extrap_order, ' and nuse <= imax <= ', bsstep_max_stages
            call core_report(lines)
         case (envint_bsstep_failed)
            if (nint(e%values(2)) == 1) then
               write(lines,'(1x,a,a,a,f9.4,a,a)') 'envelope integration', trim(where), &
                    ': a Bulirsch-Stoer step of the atmosphere from log tau = ', e%values(1), &
                    ' failed: ', trim(numerics_status_text(e%counts(1)))
            else
               write(lines,'(1x,a,a,a,f9.4,a,a)') 'envelope integration', trim(where), &
                    ': a Bulirsch-Stoer step of the envelope from log P = ', e%values(1), &
                    ' failed: ', trim(numerics_status_text(e%counts(1)))
            end if
            call core_report(lines)
         case (envint_vertex)
            continue
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
            write(lines,'(1x,a,a,a,i0,a)') 'ERROR: envelope integration', trim(where), &
                 ': the envelope did not reach the fitting mass in ', e%counts(1), ' steps'
            call core_report(lines)
         case (envint_start_exhausted)
            write(lines,'(1x,a,a,a,f9.4,a,f9.4,a,i0,a,1pe9.2,a)') &
                 'ERROR: envelope integration', trim(where), &
                 ': the atmosphere''s starting log tau, ', e%values(1), &
                 ', still lies past the photosphere''s, ', e%values(2), ', after ', &
                 e%counts(1), ' starting densities (the last ', e%values(3), &
                 ' g/cm^3); run stopped'
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
! - Round 3 (2026-11): atm_get takes its inputs intent(in) and the
!   triangle's vertex as one record (it had 24 arguments, the vertex in
!   six of them); the line of the fitted vertex, which no caller could
!   reach, went with its print_flag.

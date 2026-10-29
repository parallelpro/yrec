!----------------------------------------------------------------------
! controls_check
!----------------------------------------------------------------------
! The controls' validation, errors and warnings, and the one way a
! controls message reaches the terminal and the run log:
!   check_controls              the errors: values and combinations the
!                               run cannot use (ierr = 1); the controls
!                               read calls it (read_input), before
!                               derive_controls
!   warn_inert_controls         legacy controls a deck sets that nothing
!                               uses; the controls read calls it
!   warn_inconsistent_controls  combinations that are silently inert or
!                               cancel each other; run_yrec calls it
!                               after prepare_run
!   report_control              one line on the terminal and in the run
!                               log (derive_controls' messages use it too)
! An error goes to the terminal and the run log. A warning goes to the
! terminal only, so that the run log stays byte-pinnable, and never
! stops the run: an inert control is a probable mistake, not a certain
! one. Each rule names the code it is checked against.
module controls_check
      use io_units, only: terminal_unit, run_log_unit
      use numerics_lib, only: bsstep_max_stages, bsstep_max_extrap_order
      use star_info_lib, only: star_info
      use controls_buffer, only: max_runs
      implicit none
      private
      public :: check_controls, warn_inert_controls, warn_inconsistent_controls, &
           report_control

! The rotation transport's iteration history holds 16 iterations
! (rotation/seculr/transport_checks' max_recorded_iterations;
! docs/bug_sweep.md O14).
      integer, parameter :: max_transport_iterations = 16
! The Bulirsch-Stoer integrator's limits are numerics_lib's
! bsstep_max_stages (11, its stage sequence) and bsstep_max_extrap_order
! (7, its extrapolation tableau's columns); order 1 would index the stage
! sequence at 0, and an order above the stage count leaves its step-size
! retry loop without an end.
      integer, parameter :: message_length = 512

contains

! ---------------------------------------------------------------
! line, trimmed, on the terminal and in the run log
subroutine report_control(line)
      character(len=*), intent(in) :: line
      write(terminal_unit,'(a)') trim(line)
      write(run_log_unit,'(a)') trim(line)
end subroutine report_control

! ---------------------------------------------------------------
! a warning, on the terminal only
subroutine warn(text)
      character(len=*), intent(in) :: text
      write(terminal_unit,'(2a)') ' WARNING (controls): ', text
end subroutine warn

!----------------------------------------------------------------------
! check_controls: the errors, from the controls as read (the read stores
! them into s before it opens the run's files). ierr = 1 when
!   - semiconvection is on with any overshoot (rejected since the F77
!     PARMIN);
!   - the envelope triangle's size is not positive (with zero area its
!     walk never contains the model, surface_boundary);
!   - a rotating run's max_diffusion_iters is outside 1 to 16;
!   - the Bulirsch-Stoer controls are outside the integrator's limits
!     (1 <= max_stage_index <= 11, 2 <= extrap_order <= 7, extrap_order
!     <= max_stage_index) or the tolerance is not positive;
!   - num_runs is outside 0 to max_runs (the run cards' arrays).
!   star reads: ctrl%core_overshoot_active,
!     ctrl%envelope_overshoot_active, ctrl%extrap_order, ctrl%lovstm,
!     ctrl%max_diffusion_iters, ctrl%max_stage_index,
!     ctrl%rotation_active, ctrl%tolerance_fraction,
!     ctrl%tri_delta_logl, ctrl%tri_delta_teffl,
!     ctrl%use_semiconvection, job%num_runs
!   star writes: none
subroutine check_controls(s, ierr)
      type(star_info), intent(in) :: s
      integer, intent(inout) :: ierr
      character(len=message_length) :: line

      if (s%ctrl%use_semiconvection .and. (s%ctrl%core_overshoot_active .or. &
           s%ctrl%envelope_overshoot_active .or. s%ctrl%lovstm)) then
         write(line,'(1x,a,3l2)') 'ERROR: SEMI-CONVECTION and overshoot are both on:' &
              //' use_semiconvection (LSEMIC) with core, envelope, intermediate overshoot' &
              //' (LOVSTC, LOVSTE, LOVSTM) =', s%ctrl%core_overshoot_active, &
              s%ctrl%envelope_overshoot_active, s%ctrl%lovstm
         call report_control(line)
         ierr = 1
         return
      end if
      if (s%ctrl%tri_delta_teffl <= 0.0d0 .or. s%ctrl%tri_delta_logl <= 0.0d0) then
         write(line,'(1x,a,1pe11.3,a,1pe11.3)') 'ERROR: the envelope triangle''s size must be' &
              //' positive: tri_delta_teffl (TRIDT) =', s%ctrl%tri_delta_teffl, &
              ', tri_delta_logl (TRIDL) =', s%ctrl%tri_delta_logl
         call report_control(line)
         ierr = 1
         return
      end if
      if (s%ctrl%rotation_active .and. (s%ctrl%max_diffusion_iters < 1 .or. &
           s%ctrl%max_diffusion_iters > max_transport_iterations)) then
         write(line,'(1x,a,i0)') 'ERROR: max_diffusion_iters (ITDIF2) must be 1 to 16, not ', &
              s%ctrl%max_diffusion_iters
         call report_control(line)
         ierr = 1
         return
      end if
      if (s%ctrl%tolerance_fraction <= 0.0d0 .or. s%ctrl%max_stage_index < 1 .or. &
           s%ctrl%max_stage_index > bsstep_max_stages .or. s%ctrl%extrap_order < 2 &
           .or. s%ctrl%extrap_order > bsstep_max_extrap_order &
           .or. s%ctrl%extrap_order > s%ctrl%max_stage_index) then
         write(line,'(1x,a,1pe11.3,a,i0,a,i0)') 'ERROR: the Bulirsch-Stoer controls need' &
              //' tolerance_fraction (STOLR0) > 0, max_stage_index (IMAX) 1 to 11 and' &
              //' extrap_order (NUSE) 2 to 7 and at most IMAX, not', &
              s%ctrl%tolerance_fraction, ', ', s%ctrl%max_stage_index, ', ', s%ctrl%extrap_order
         call report_control(line)
         ierr = 1
         return
      end if
      if (s%job%num_runs < 0 .or. s%job%num_runs > max_runs) then
         write(line,'(1x,a,i0,a,i0)') 'ERROR: num_runs (NUMRUN) must be 0 to ', max_runs, &
              ', not ', s%job%num_runs
         call report_control(line)
         ierr = 1
         return
      end if
end subroutine check_controls

!----------------------------------------------------------------------
! warn_inert_controls: legacy controls that a deck may set and that
! nothing uses (each warns only when set away from its default).
subroutine warn_inert_controls
      use legacy_controls, only: lsolwind, lpmm, pmmwmax, zopal952, zalex2, fscomp
! core/step/mass_change.f90: the rotation-scaled solar wind is commented out
      if (lsolwind) call warn('use_rotation_scaled_solar_wind (LSOLWIND) is on, but the' &
           //' rotation-scaled solar wind is commented out: it does nothing')
! the Matt law is use_pmm_wind_law (LMWIND); lpmm is read and used nowhere
      if (lpmm) call warn('use_pmm_wind (LPMM) does nothing; the Matt wind law is' &
           //' use_pmm_wind_law (LMWIND)')
! the saturation is wind_saturation_omega or solar_omega_sat
      if (pmmwmax /= 0.0d0) call warn('pmm_omega_sat (PMMWMAX) does nothing; the wind''s' &
           //' saturation is wind_saturation_omega or solar_omega_sat')
! the OPAL 1995 and the Alexander tables take one Z (opal95_single_table_z,
! alex_table_z1)
      if (zopal952 /= 0.0d0) call warn('opal95_z2 (ZOPAL952) does nothing; the OPAL 1995' &
           //' table takes one Z, opal95_single_table_z')
      if (zalex2 /= 0.0d0) call warn('alex_z2 (ZALEX2) does nothing; the Alexander tables' &
           //' take one Z, alex_table_z1')
! the path is expanded (run_files) and opened nowhere
      if (len_trim(fscomp) > 0 .and. fscomp(1:1) /= achar(0)) &
           call warn('surface_composition_file (FSCOMP) is set, but no surface' &
           //' composition file is written')
end subroutine warn_inert_controls

!----------------------------------------------------------------------
! warn_inconsistent_controls: combinations that are silently inert or
! self-cancelling (the class of mistake that shipped a template with
! overshoot_alpha_envelope = 0.5 and envelope_overshoot_active unset for
! years). Each rule is checked against the gating code named above it.
!   star reads: ctrl%betac, ctrl%calibrate_solar_model,
!     ctrl%calibrate_star_flag, ctrl%core_overshoot_active,
!     ctrl%diffuse_lithium, ctrl%envelope_overshoot_active, ctrl%lovmax,
!     ctrl%overshoot_alpha_core, ctrl%overshoot_alpha_envelope,
!     ctrl%use_new_diffusion_routines
!   star writes: none
subroutine warn_inconsistent_controls(s)
      type(star_info), intent(in) :: s

! overshoot extents are consumed only inside the active-flag gates
! (convection/convection_zones.f90's extend_by_overshoot)
      if (s%ctrl%overshoot_alpha_core > 0.0d0 .and. &
           .not. s%ctrl%core_overshoot_active) &
           call warn('overshoot_alpha_core is set but core_overshoot_active' &
           //' is off -- core overshoot is NOT applied')
      if (s%ctrl%overshoot_alpha_envelope > 0.0d0 .and. &
           .not. s%ctrl%envelope_overshoot_active) &
           call warn('overshoot_alpha_envelope is set but ' &
           //'envelope_overshoot_active is off -- envelope overshoot is' &
           //' NOT applied')
      if (s%ctrl%core_overshoot_active .and. &
           s%ctrl%overshoot_alpha_core <= 0.0d0) &
           call warn('core_overshoot_active is on but ' &
           //'overshoot_alpha_core is zero -- zero-extent overshoot')
      if (s%ctrl%envelope_overshoot_active .and. &
           s%ctrl%overshoot_alpha_envelope <= 0.0d0) &
           call warn('envelope_overshoot_active is on but ' &
           //'overshoot_alpha_envelope is zero -- zero-extent overshoot')

! lovmax caps the core overshoot extent at betac times the core radius
! (convection/convection_zones.f90's extend_by_overshoot) -- with
! overshoot_beta_core = 0 the cap is zero
      if (s%ctrl%core_overshoot_active .and. s%ctrl%lovmax .and. &
           s%ctrl%betac <= 0.0d0) &
           call warn('overshoot_max_extent_flag is on with ' &
           //'overshoot_beta_core = 0 -- core overshoot is capped to zero')

! light-element settling lives only in the microdiff pipeline
! (element_diffusion/; core/step/composition_step.f90's settle hands it
! to diffuse_elements, which takes that path only with
! use_new_diffusion_routines)
      if (s%ctrl%diffuse_lithium .and. &
           .not. s%ctrl%use_new_diffusion_routines) &
           call warn('diffuse_lithium is on but use_new_diffusion_routines' &
           //' is off -- light-element settling is NOT computed')

! the two calibration modes expand the run list with incompatible card
! protocols (core/run/calibration.f90: solar triples against star pairs)
      if (s%ctrl%calibrate_solar_model .and. &
           s%ctrl%calibrate_star_flag) &
           call warn('calibrate_solar_model and calibrate_star_flag are' &
           //' both on -- their run-list protocols conflict')
end subroutine warn_inconsistent_controls

end module controls_check

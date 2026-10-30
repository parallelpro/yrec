!----------------------------------------------------------------------
! controls_check
!----------------------------------------------------------------------
! The controls' validation, errors and warnings:
!   check_controls              the errors: values and combinations the
!                               run cannot use (ierr = 1); the controls
!                               read calls it (read_input), before
!                               derive_controls
!   warn_inert_controls         legacy controls a deck sets that nothing
!                               uses; the controls read calls it
!   warn_inconsistent_controls  combinations that are silently inert or
!                               cancel each other; run_yrec calls it
!                               after prepare_run
! Their messages go through io_units (report_error, report_warning): on
! the terminal and in the run log. A warning never stops the run: an
! inert control is a probable mistake, not a certain one. Each rule
! names the code it is checked against.
module controls_check
      use io_units, only: report_error, report_warning
      use numerics_lib, only: bsstep_max_stages, bsstep_max_extrap_order
      use star_info_lib, only: star_info
      use controls_buffer, only: max_runs, use_semiconvection, core_overshoot_active, &
           envelope_overshoot_active, intermediate_overshoot_active, tri_delta_teffl, tri_delta_logl, &
           rotation_active, max_diffusion_iters, tolerance_fraction, max_stage_index, &
           extrap_order, num_runs
      use rotation_configuration, only: max_transport_iterations
      implicit none
      private
      public :: check_controls, warn_inert_controls, warn_inconsistent_controls

! the longest message line
      integer, parameter :: message_length = 512

contains

! ---------------------------------------------------------------
! a warning about the controls
subroutine warn(text)
      character(len=*), intent(in) :: text
      call report_warning(text)
end subroutine warn

!----------------------------------------------------------------------
! check_controls: the errors, from the controls as read. ierr = 1 when
!   - semiconvection is on with any overshoot (rejected since the F77
!     PARMIN);
!   - the envelope triangle's size is not positive (with zero area its
!     walk never contains the model, surface_boundary);
!   - a rotating run's max_diffusion_iters is outside 1 to
!     max_transport_iterations (rotation_configuration, 16: the transport
!     check's iteration history; docs/bug_sweep.md O14);
!   - the Bulirsch-Stoer controls are outside the integrator's limits,
!     numerics_lib's bsstep_max_stages (11, its stage sequence) and
!     bsstep_max_extrap_order (7, its extrapolation tableau's columns):
!     1 <= max_stage_index <= 11, 2 <= extrap_order <= 7 (order 1 would
!     index the stage sequence at 0) and extrap_order <= max_stage_index
!     (above, its step-size retry loop has no end); or the tolerance is
!     not positive;
!   - num_runs is outside 0 to max_runs (the run cards' arrays).
! It reads the controls in the read's buffer, after adopt_canonical_names
! and map_user_inputs.
subroutine check_controls(ierr)
      integer, intent(inout) :: ierr
      character(len=message_length) :: line

      if (use_semiconvection .and. (core_overshoot_active .or. &
           envelope_overshoot_active .or. intermediate_overshoot_active)) then
         write(line,'(a,3l2)') 'SEMI-CONVECTION and overshoot are both on:' &
              //' use_semiconvection (LSEMIC) with core, envelope, intermediate overshoot' &
              //' (LOVSTC, LOVSTE, LOVSTM) =', core_overshoot_active, &
              envelope_overshoot_active, intermediate_overshoot_active
         call report_error(line)
         ierr = 1
         return
      end if
      if (tri_delta_teffl <= 0.0d0 .or. tri_delta_logl <= 0.0d0) then
         write(line,'(a,1pe11.3,a,1pe11.3)') 'the envelope triangle''s size must be' &
              //' positive: tri_delta_teffl (TRIDT) =', tri_delta_teffl, &
              ', tri_delta_logl (TRIDL) =', tri_delta_logl
         call report_error(line)
         ierr = 1
         return
      end if
      if (rotation_active .and. (max_diffusion_iters < 1 .or. &
           max_diffusion_iters > max_transport_iterations)) then
         write(line,'(a,i0,a,i0)') 'max_diffusion_iters (ITDIF2) must be 1 to ', &
              max_transport_iterations, ', not ', max_diffusion_iters
         call report_error(line)
         ierr = 1
         return
      end if
      if (tolerance_fraction <= 0.0d0 .or. max_stage_index < 1 .or. &
           max_stage_index > bsstep_max_stages .or. extrap_order < 2 &
           .or. extrap_order > bsstep_max_extrap_order &
           .or. extrap_order > max_stage_index) then
         write(line,'(a,i0,a,i0,a,1pe11.3,a,i0,a,i0)') 'the Bulirsch-Stoer controls need' &
              //' tolerance_fraction (STOLR0) > 0, max_stage_index (IMAX) 1 to ', &
              bsstep_max_stages, ' and extrap_order (NUSE) 2 to ', bsstep_max_extrap_order, &
              ' and at most IMAX, not', tolerance_fraction, ', ', &
              max_stage_index, ', ', extrap_order
         call report_error(line)
         ierr = 1
         return
      end if
      if (num_runs < 0 .or. num_runs > max_runs) then
         write(line,'(a,i0,a,i0)') 'num_runs (NUMRUN) must be 0 to ', max_runs, &
              ', not ', num_runs
         call report_error(line)
         ierr = 1
         return
      end if
end subroutine check_controls

!----------------------------------------------------------------------
! warn_inert_controls: legacy controls that a deck may set and that
! nothing uses (each warns only when set away from its default).
subroutine warn_inert_controls
      use legacy_controls, only: lsolwind, lpmm, pmmwmax, zopal952, zalex2
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
end subroutine warn_inert_controls

!----------------------------------------------------------------------
! warn_inconsistent_controls: combinations that are silently inert or
! self-cancelling (the class of mistake that shipped a template with
! overshoot_alpha_envelope = 0.5 and envelope_overshoot_active unset for
! years). Each rule is checked against the gating code named above it.
!   star reads: control%calibrate_solar_model, control%calibrate_star_flag,
!     physics%core_overshoot_active, physics%envelope_overshoot_active,
!     physics%overshoot_alpha_core, physics%overshoot_alpha_envelope,
!     physics%overshoot_beta_core, physics%overshoot_max_extent_flag
!   star writes: none
subroutine warn_inconsistent_controls(s)
      type(star_info), intent(in) :: s

! overshoot extents are consumed only inside the active-flag gates
! (convection/convection_zones.f90's extend_by_overshoot)
      if (s%physics%overshoot_alpha_core > 0.0d0 .and. &
           .not. s%physics%core_overshoot_active) &
           call warn('overshoot_alpha_core is set but core_overshoot_active' &
           //' is off -- core overshoot is NOT applied')
      if (s%physics%overshoot_alpha_envelope > 0.0d0 .and. &
           .not. s%physics%envelope_overshoot_active) &
           call warn('overshoot_alpha_envelope is set but ' &
           //'envelope_overshoot_active is off -- envelope overshoot is' &
           //' NOT applied')
      if (s%physics%core_overshoot_active .and. &
           s%physics%overshoot_alpha_core <= 0.0d0) &
           call warn('core_overshoot_active is on but ' &
           //'overshoot_alpha_core is zero -- zero-extent overshoot')
      if (s%physics%envelope_overshoot_active .and. &
           s%physics%overshoot_alpha_envelope <= 0.0d0) &
           call warn('envelope_overshoot_active is on but ' &
           //'overshoot_alpha_envelope is zero -- zero-extent overshoot')

! overshoot_max_extent_flag caps the core overshoot extent at overshoot_beta_core times the core radius
! (convection/convection_zones.f90's extend_by_overshoot) -- with
! overshoot_beta_core = 0 the cap is zero
      if (s%physics%core_overshoot_active .and. s%physics%overshoot_max_extent_flag .and. &
           s%physics%overshoot_beta_core <= 0.0d0) &
           call warn('overshoot_max_extent_flag is on with ' &
           //'overshoot_beta_core = 0 -- core overshoot is capped to zero')

! the two calibration modes expand the run list with incompatible card
! protocols (core/run/calibration.f90: solar triples against star pairs)
      if (s%control%calibrate_solar_model .and. &
           s%control%calibrate_star_flag) &
           call warn('calibrate_solar_model and calibrate_star_flag are' &
           //' both on -- their run-list protocols conflict')
end subroutine warn_inconsistent_controls

end module controls_check

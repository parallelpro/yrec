!----------------------------------------------------------------------
! controls_check
!----------------------------------------------------------------------
! A startup sanity pass over the adopted controls, warning about
! combinations that are silently inert or
! self-cancelling -- the class of mistake that shipped a template
! with overshoot_alpha_envelope = 0.5 and envelope_overshoot_active
! unset for years. Every rule below is verified against the actual
! gating code (file:line in the comment), not inferred physics.
!
! Warnings go to the TERMINAL ONLY (the run log stays byte-pinnable),
! and never abort the run: an inert control is a probable mistake,
! not a certain one. Genuine conflicts that must abort stay where
! they are (e.g. semiconvection with overshoot, in derive_controls).
!
! Called from run_yrec once, after star_setup -- all controls are
! adopted and final by then.
module controls_check
      implicit none
      private
      public :: warn_inconsistent_controls

contains

subroutine warn_inconsistent_controls(s)
      use star_info_lib, only: star_info
      use io_units
      type(star_info), intent(inout) :: s

! overshoot extents are consumed only inside the active-flag gates
! (convection/convection_zones.f90's extend_by_overshoot: core cycle at "SKIP IF NO CORE
! OVERSHOOT", envelope likewise)
      if (s%ctrl%overshoot_alpha_core > 0.0d0 .and. &
           .not. s%job%core_overshoot_active) &
           call warn('overshoot_alpha_core is set but core_overshoot_active' &
           //' is off -- core overshoot is NOT applied')
      if (s%ctrl%overshoot_alpha_envelope > 0.0d0 .and. &
           .not. s%job%envelope_overshoot_active) &
           call warn('overshoot_alpha_envelope is set but ' &
           //'envelope_overshoot_active is off -- envelope overshoot is' &
           //' NOT applied')
      if (s%job%core_overshoot_active .and. &
           s%ctrl%overshoot_alpha_core <= 0.0d0) &
           call warn('core_overshoot_active is on but ' &
           //'overshoot_alpha_core is zero -- zero-extent overshoot')
      if (s%job%envelope_overshoot_active .and. &
           s%ctrl%overshoot_alpha_envelope <= 0.0d0) &
           call warn('envelope_overshoot_active is on but ' &
           //'overshoot_alpha_envelope is zero -- zero-extent overshoot')

! lovmax caps the core overshoot extent at betac*R_edge
! (overshoot_boundaries.f90: pscale_up = min(alpha*H, betac*r)) --
! with overshoot_beta_core = 0 the cap IS zero
      if (s%job%core_overshoot_active .and. s%ctrl%lovmax .and. &
           s%ctrl%betac <= 0.0d0) &
           call warn('overshoot_max_extent_flag is on with ' &
           //'overshoot_beta_core = 0 -- core overshoot is capped to zero')

! light-element settling lives only in the microdiff pipeline
! (element_diffusion/; core/step/composition_step.f90 dispatches there only when
! use_new_diffusion_routines is set)
      if (s%ctrl%diffuse_lithium .and. &
           .not. s%ctrl%use_new_diffusion_routines) &
           call warn('diffuse_lithium is on but use_new_diffusion_routines' &
           //' is off -- light-element settling is NOT computed')

! the two calibration modes expand the run list with incompatible
! card protocols (run_yrec: setcal triples vs setscal pairs)
      if (s%ctrl%calibrate_solar_model .and. &
           s%ctrl%calibrate_star_flag) &
           call warn('calibrate_solar_model and calibrate_star_flag are' &
           //' both on -- their run-list protocols conflict')
end subroutine warn_inconsistent_controls

subroutine warn(text)
      use io_units
      character(len=*), intent(in) :: text
      write(terminal_unit,'(2a)') ' WARNING (controls): ', text
end subroutine warn

end module controls_check

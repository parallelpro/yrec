!----------------------------------------------------------------------
! setcal
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setcal.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Sets up the run-list for an automatic solar calibration: expands the
! two "seed" runs (kind cards 1 and 2) into 48 runs by copying their
! rescaling/mixing-length parameters forward, and hardwires run 2 to
! stop at 1e8 yr and run 3 to stop at the target solar age (both
! scaled by age_scale_factor). Every third run (starting at run 4) is
! reset to rerun the seed rescaling (kind card 1); the runs after that
! alternately reuse run 2's or run 3's age/timestep stop criteria so
! that check_solar_calibration.f90's Newton correction (applied on the even NK runs
! following the odd-NK rescaling run) has fixed X/alpha evaluation
! points to iterate against.
subroutine setup_solar_calibration(age_scale_factor)

      use star_info_lib, only: star
      use controls_lib, only: max_runs
      implicit none

      double precision, intent(in) :: age_scale_factor
! --- locals ---
      integer :: i, j
! 16 three-run calibration cycles (the run arrays in controls_lib are
! dimensioned max_runs = 50).
      integer, parameter :: num_calibration_runs = max_runs-2

! SET UP RUN TO CALIBRATE A SOLAR MODEL.
! THIS CONSISTS OF SETTING THE NUMBER OF RUNS TO 48 (16 THREE-RUN
! CALIBRATION CYCLES), AND COPYING THE RELEVANT PARAMETERS FROM THE
! FIRST RUN TO THE CALIBRATING RUNS.
! mhp 5/96 changed to do solar models in 3 runs rather than 2.
      star%job%num_runs = num_calibration_runs
      do i = 2,num_calibration_runs
         star%job%initial_x_array(i) = star%job%initial_x_array(1)
         star%job%initial_z_array(i) = star%job%initial_z_array(1)
         star%job%mixing_length_array(i) = star%job%mixing_length_array(1)
         star%job%has_senv0_array(i) = star%job%has_senv0_array(1)
         star%job%senv0_array(i) = star%job%senv0_array(1)
      end do
      do i = 4,46,3
         star%job%rescale_kind(i) = star%job%rescale_kind(1)
         star%job%first_call_flag(i) = .true.
         star%job%num_models(i) = star%job%num_models(1)
         star%job%rsclzc(i) = star%job%rsclzc(1)
         star%job%rsclzm1(i) = star%job%rsclzm1(1)
         star%job%rsclzm2(i) = star%job%rsclzm2(1)
         do j = 1,4
            star%job%rescale_params(j,i) = star%job%rescale_params(j,1)
         end do
      end do
! MHP 06/13 HARDWIRE RUN 2 TO 1D8 YEARS AND RUN3 TO CALSOLAGE YEARS
      star%job%target_end_age(2) = 1.0d8*age_scale_factor
      star%job%timestep_override(2) = star%job%timestep_override(2)*age_scale_factor
      star%job%target_end_age(3) = star%ctrl%target_solar_age*age_scale_factor
      star%job%timestep_override(3) = star%job%timestep_override(3)*age_scale_factor
      do i = 5,47,3
         star%job%rescale_kind(i) = 1
         star%job%first_call_flag(i) = .false.
         star%job%num_models(i) = star%job%num_models(2)
         star%job%target_end_age(i) = star%job%target_end_age(2)
         star%job%end_age_stop_active(i) = star%job%end_age_stop_active(2)
         star%job%timestep_override(i) = star%job%timestep_override(2)
         star%job%timestep_override_active(i) = star%job%timestep_override_active(2)
      end do
      do i = 6,num_calibration_runs,3
         star%job%rescale_kind(i) = 1
         star%job%first_call_flag(i) = .false.
         star%job%num_models(i) = star%job%num_models(3)
         star%job%target_end_age(i) = star%job%target_end_age(3)
         star%job%end_age_stop_active(i) = star%job%end_age_stop_active(3)
         star%job%timestep_override(i) = star%job%timestep_override(3)
         star%job%timestep_override_active(i) = star%job%timestep_override_active(3)
      end do
      star%solar_calibration_active = .false.
      return
end subroutine setup_solar_calibration

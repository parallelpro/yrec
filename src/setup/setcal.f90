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
! that chkcal.f90's Newton correction (applied on the even NK runs
! following the odd-NK rescaling run) has fixed X/alpha evaluation
! points to iterate against.
subroutine setcal(age_scale_factor)

      use const_lib
      implicit none

      double precision, intent(in) :: age_scale_factor

! MHP 6/13 ADD OPTION TO CALIBRATE SOLAR Z/X, SOLAR Z/X, SOLAR AGE
! common/cals2/: only target_solar_age is used here. Naming matches
! chkcal.f90.
      double precision :: luminosity_tolerance, radius_tolerance, &
           zx_tolerance
      logical :: calibrate_solar_model, calibrate_solar_zx
      double precision :: target_solar_zx, target_solar_age
      common/cals2/ luminosity_tolerance, radius_tolerance, zx_tolerance, &
           calibrate_solar_model, calibrate_solar_zx, target_solar_zx, &
           target_solar_age



! common/newxym/: initial_x_array/initial_z_array/mixing_length_array/
! has_senv0_array/senv0_array, all used here. Naming matches
! chkcal.f90.
      double precision :: initial_x_array(50), initial_z_array(50), &
           mixing_length_array(50)
      logical :: has_senv0_array(50)
      double precision :: senv0_array(50)
      common /newxym/ initial_x_array, initial_z_array, &
           mixing_length_array, has_senv0_array, senv0_array

! common/calsun/: only solar_calibration_active is used here. Naming
! matches chkcal.f90.
      double precision :: dlum_dx, drad_dx, dlum_dalpha, drad_dalpha, &
           log_l_prev, log_r_prev, delta_x, delta_alpha
      logical :: solar_calibration_active
      common/calsun/ dlum_dx, drad_dx, dlum_dalpha, drad_dalpha, &
           log_l_prev, log_r_prev, delta_x, delta_alpha, &
           solar_calibration_active

! MHP 8/25 Removed character file names from common block
! common/zramp/: rsclzc/rsclzm1/rsclzm2 are used here; the remaining
! members are unused placeholders. Naming matches gtlaol2.f90.
      double precision :: rsclzc(50), rsclzm1(50), rsclzm2(50)
      integer :: iolaol2, ioopal2, nk
      logical :: use_z_ramp
      common/zramp/ rsclzc, rsclzm1, rsclzm2, iolaol2, ioopal2, nk, &
           use_z_ramp

      save

! --- locals ---
      integer :: i, j

! SET UP RUN TO CALIBRATE A SOLAR MODEL.
! THIS CONSISTS OF SETTING THE NUMBER OF RUNS TO THE MAXIMUM (50),
! AND COPYING THE RELEVANT PARAMETERS FROM THE FIRST TWO RUNS TO
! THE NEXT SERIES OF 24 CALIBRATING RUNS.
! mhp 5/96 changed to do solar models in 3 runs rather than 2.
      num_runs = 48
      do i = 2,48
         initial_x_array(i) = initial_x_array(1)
         initial_z_array(i) = initial_z_array(1)
         mixing_length_array(i) = mixing_length_array(1)
         has_senv0_array(i) = has_senv0_array(1)
         senv0_array(i) = senv0_array(1)
      end do
      do i = 4,46,3
         rescale_kind(i) = rescale_kind(1)
         first_call_flag(i) = .true.
         num_models(i) = num_models(1)
         rsclzc(i) = rsclzc(1)
         rsclzm1(i) = rsclzm1(1)
         rsclzm2(i) = rsclzm2(1)
         do j = 1,4
            rescale_params(j,i) = rescale_params(j,1)
         end do
      end do
! MHP 06/13 HARDWIRE RUN 2 TO 1D8 YEARS AND RUN3 TO CALSOLAGE YEARS
      target_end_age(2) = 1.0d8*age_scale_factor
      timestep_override(2) = timestep_override(2)*age_scale_factor
      target_end_age(3) = target_solar_age*age_scale_factor
      timestep_override(3) = timestep_override(3)*age_scale_factor
      do i = 5,47,3
         rescale_kind(i) = 1
         first_call_flag(i) = .false.
         num_models(i) = num_models(2)
         target_end_age(i) = target_end_age(2)
         end_age_stop_active(i) = end_age_stop_active(2)
         timestep_override(i) = timestep_override(2)
         timestep_override_active(i) = timestep_override_active(2)
      end do
      do i = 6,48,3
         rescale_kind(i) = 1
         first_call_flag(i) = .false.
         num_models(i) = num_models(3)
         target_end_age(i) = target_end_age(3)
         end_age_stop_active(i) = end_age_stop_active(3)
         timestep_override(i) = timestep_override(3)
         timestep_override_active(i) = timestep_override_active(3)
      end do
      solar_calibration_active = .false.
      return
end subroutine setcal

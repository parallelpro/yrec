!----------------------------------------------------------------------
! calibration
!----------------------------------------------------------------------
! The calibration loops' two halves: the setup before a calibration run
! (setup_solar_calibration, setup_star_calibration) and the check after
! each run that converges or sets up the next trial (check_solar_calibration,
! check_star_calibration).
module calibration
      implicit none
      private
      public :: check_solar_calibration, check_star_calibration, setup_solar_calibration, setup_star_calibration

contains

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! check_solar_calibration (formerly chkcal)
!----------------------------------------------------------------------
! Checks whether the current model is a converged solar calibration
! (log L and log R, and optionally log Z/X, within tolerance of the
! solar values) and, if not, applies a fixed 2x2 Newton correction to
! the envelope hydrogen fraction (X) and mixing-length alpha using
! pre-determined empirical partial derivatives of log L and log R with
! respect to X and alpha (star%dlum_dx/star%drad_dx/star%dlum_dalpha/star%drad_dalpha).
! Called after every completed calibration TRIPLE (mod(nk,3) == 0;
! the protocol comment in run_yrec's run loop is the reference):
! card 3k+1 rescales, 3k+2 settles to 1e8 yr, 3k+3 evolves to the
! solar age. On a miss, the corrected (X[, Z], alpha) are written
! into the NEXT triple's cards (run_index+1 .. +3), reading the
! previous guess from run_index-2 (this triple's rescale card).
subroutine check_solar_calibration(log_l_lsun, log_r_rsun, run_index, current_zx)

      use star_info_lib, only: star
      use io_units
      use math_lib
      implicit none

      double precision, intent(in) :: log_l_lsun, log_r_rsun
      integer, intent(in) :: run_index
      double precision, intent(in) :: current_zx

! --- locals ---
      double precision :: log_zx_mismatch, delta_z

! CHECK TO SEE IF THE MODEL IS A CALIBRATED SOLAR MODEL.
! IF NOT, ESTIMATE CORRECTIONS TO X AND ALPHA.
      if(abs(log_l_lsun).lt.star%ctrl%luminosity_tolerance .and. &
           abs(log_r_rsun).lt.star%ctrl%radius_tolerance)then
! ADD CHECK FOR Z/X
         if(star%ctrl%calibrate_solar_zx .and.star%ctrl%target_solar_zx.gt.0.0d0)then
            log_zx_mismatch = log10(current_zx)-log10(star%ctrl%target_solar_zx)
            if(abs(log_zx_mismatch).lt.star%ctrl%zx_tolerance)then
               star%solar_calibration_active = .true.
               return
            endif
         else
! CALIBRATED SOLAR MODEL.  SET UP OUTPUT FLAGS AND EXIT
            star%solar_calibration_active = .true.
            return
         endif
      endif

! Added code to use pre-determined partial derivatives    LLP  6/28/09
      star%dlum_dx = -3.78d0            ! empirical result:  -3.783    RMS error .070
      star%drad_dx = -0.89d0             ! empirical result:  -0.890    RMS error .048
      star%dlum_dalpha = 0.0139d0            ! empirical result:  +0.139    RMS error .0022
      star%drad_dalpha = -0.050d0            ! empirical result:  -0.0504   RMS error .0059
! (2026 goto campaign: an unconditional `goto 1234` here bypassed the
! original run_index-based partial-derivative code, which was therefore
! dead; it was removed. See git history for the bypassed block.)
!     USE DERIVATIVES OF L AND R WITH RESPECT TO X AND ALPHA TO
!     GET IMPROVED GUESSES FOR ALPHA AND X.
         star%delta_alpha = ((log_l_lsun*star%drad_dx/star%dlum_dx-log_r_rsun)/ &
              (star%drad_dalpha-star%dlum_dalpha*star%drad_dx/star%dlum_dx))
         star%delta_x = -(log_l_lsun + star%dlum_dalpha*star%delta_alpha)/star%dlum_dx
         star%job%mixing_length_array(run_index+1) = &
              star%job%mixing_length_array(run_index-2)+star%delta_alpha
         star%job%mixing_length_array(run_index+2) = star%job%mixing_length_array(run_index+1)
         star%job%mixing_length_array(run_index+3) = star%job%mixing_length_array(run_index+1)
         if(star%ctrl%calibrate_solar_zx)then
            star%job%rescale_params(3,run_index+1) = &
                 star%job%rescale_params(3,run_index-2)*star%ctrl%target_solar_zx/current_zx
            delta_z = star%job%rescale_params(3,run_index+1) - &
                 star%job%rescale_params(3,run_index-2)
            star%job%initial_z_array(run_index+1) = star%job%rescale_params(3,run_index+1)
            star%job%initial_z_array(run_index+2) = star%job%rescale_params(3,run_index+1)
            star%job%initial_z_array(run_index+3) = star%job%rescale_params(3,run_index+1)
         endif
         star%job%rescale_params(2,run_index+1) = star%job%rescale_params(2,run_index-2)+star%delta_x
         star%job%initial_x_array(run_index+1) = star%job%rescale_params(2,run_index+1)
         star%job%initial_x_array(run_index+2) = star%job%rescale_params(2,run_index+1)
         star%job%initial_x_array(run_index+3) = star%job%rescale_params(2,run_index+1)
         write(terminal_unit,*) "New BL, Old BL, Delta BL: ", &
             log_l_lsun, star%log_l_prev, log_l_lsun-star%log_l_prev
         write(terminal_unit,*) "New RL, Old RL, Delta RL: ", &
             log_r_rsun, star%log_r_prev, log_r_rsun-star%log_r_prev
         write(terminal_unit,*) "New X, Old X, DX: ", &
             star%job%rescale_params(2,run_index+1), star%job%rescale_params(2,run_index-2), &
             star%delta_x
         write(terminal_unit,*) "New A, Old A, DA: ", &
             star%job%mixing_length_array(run_index+1), &
             star%job%mixing_length_array(run_index-2), star%delta_alpha
         if(star%ctrl%calibrate_solar_zx)then
            write(terminal_unit,*) "New Z, Old Z, DZ: ", &
             star%job%rescale_params(3,run_index+1), star%job%rescale_params(3,run_index-2), &
             delta_z
         endif
         star%log_l_prev = log_l_lsun
         star%log_r_prev = log_r_rsun
end subroutine check_solar_calibration

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! check_star_calibration (formerly chkscal)
!----------------------------------------------------------------------
! Checks whether the evolving track has passed a target stellar radius
! (R*), interpolates the luminosity and age at that radius, and either
! signals a converged calibration (log L/Lsun at R* matches the target
! luminosity to within tolerance) or sets up the envelope hydrogen
! fraction X for the next trial run via linear/secant interpolation in
! X vs. log L at R*.
! ONLY CALLED FOR EVEN NK, ASSUMES RESCALING ON ODD NK AND EVOLVING
! ON EVEN NK
subroutine check_star_calibration(log_l_lsun, log_teff, current_age, run_index)

      use star_info_lib, only: star
      use io_units
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: log_l_lsun, log_teff, current_age
      integer, intent(in) :: run_index

! locals
      double precision :: teff_current, r_rsun_current
      double precision :: dlogl_dlogr, dage_dlogr
      double precision :: new_x, prev_x, dx_dlogl
! the hydrogen-fraction step used to seed the second calibration run
! (the "Y=Y+0.01" step in the messages below)
      double precision, parameter :: x_calibration_step = 0.01d0

!     star%star_found_flag                 T - have got a star at R* and L*
!     star%just_passed_target_radius_flag  T - on this run have just passed R*
!     ctrl%target_luminosity_lsun          luminosity (L/Lsun) wanted by adjusting X
!     ctrl%target_star_luminosity_tolerance  tolerance on that luminosity
!     job%target_radius_rsun               R* (R/Rsun, LINEAR, despite the
!                                          log_r_* names of the model values)
!     r_rsun_current                   R/Rsun of current model (linear)
!     star%log_r_prev_model                R/Rsun of previous model (linear)
!     current_age / star%age_prev_model    age of current / previous model (Gyr)
!     star%age_at_target_radius            age of model at R*
!     log_l_lsun / star%log_l_prev_model   log(L/Lsun) of current / previous model
!     star%log_l_at_target_radius          log(L/Lsun) at R*
!     star%log_l_at_target_radius_prev_run same, previous run
!     prev_x                               X of previous run = rescale_params(2, NK-1)
!
!     Check if star has passed R*.
!     If not store L and age and return.
      star%just_passed_target_radius_flag=.false.
      teff_current = exp10(log_teff)
      r_rsun_current = sqrt((exp10(log_l_lsun))*star%solar_luminosity_cgs/ &
           (c4pi*csig))/(teff_current*teff_current*star%solar_radius_cgs)
      if(r_rsun_current.gt.star%log_r_prev_model) then
         if(.not.(r_rsun_current.gt.star%job%target_radius_rsun.and. &
              star%log_r_prev_model.le.star%job%target_radius_rsun)) then
              star%log_l_prev_model = log_l_lsun
              star%age_prev_model = current_age
              star%log_r_prev_model = r_rsun_current
              return
          end if
      else
          if (.not.(r_rsun_current.lt.star%job%target_radius_rsun.and. &
               star%log_r_prev_model.ge.star%job%target_radius_rsun)) then
              star%log_l_prev_model = log_l_lsun
              star%age_prev_model = current_age
              star%log_r_prev_model = r_rsun_current
              return
           end if
      endif
!
      write(*,*) ' Just passed R*'
      write(run_log_unit,*) '#Just passed R*'
      star%just_passed_target_radius_flag = .true.
!
!     Have previous L,Age and current L,Age (one before R* and
!     one after R*).  Interpolate to get L,Age at R*
      dlogl_dlogr = (log_l_lsun-star%log_l_prev_model)/ &
           (r_rsun_current-star%log_r_prev_model)
      star%log_l_at_target_radius = log_l_lsun + &
           dlogl_dlogr*(star%job%target_radius_rsun-r_rsun_current)
      dage_dlogr = (current_age-star%age_prev_model)/ &
           (r_rsun_current-star%log_r_prev_model)
      star%age_at_target_radius = current_age + &
           dage_dlogr*(star%job%target_radius_rsun-r_rsun_current)
      write(*,*) ' X, LogL/Lsun at R* =', star%job%rescale_params(2,run_index-1), &
           star%log_l_at_target_radius
      write(run_log_unit,*) '#X, LogL/Lsun at R* =', &
           star%job%rescale_params(2,run_index-1), star%log_l_at_target_radius
      if (abs(exp10(star%log_l_at_target_radius)-star%ctrl%target_luminosity_lsun) &
           .le. star%ctrl%target_star_luminosity_tolerance) then
!        Get here then have track that passes through specified
!        L and R. Use age at R for final
!        run to stop at that age. Do one more run
!        stopping at interpolated age.
         star%star_found_flag=.true.
       star%job%end_age_stop_active(run_index+1) = .true.
       star%job%target_end_age(run_index+1) = star%age_at_target_radius*1.0d9
       star%job%end_age_stop_active(run_index+2) = .true.
       star%job%target_end_age(run_index+2) = star%age_at_target_radius*1.0d9
       star%job%rescale_params(2,run_index+1) = star%job%rescale_params(2,run_index-1)
       star%job%initial_x_array(run_index+1) = star%job%rescale_params(2,run_index+1)
       star%job%initial_x_array(run_index+2) = star%job%initial_x_array(run_index+1)
!
      write(*,*) ' Have hit R* & L*, prepare final run to age:', &
        star%age_at_target_radius
      write(run_log_unit,*)'#Have hit R* & L*, prepare final run to age:', &
        star%age_at_target_radius
         return
      else
         if (run_index .eq. 2) then
!           First time through. Save L and X at R*.
!           Add 0.01 to Y. Start next run.
            star%log_l_at_target_radius_prev_run = star%log_l_at_target_radius
          prev_x = star%job%rescale_params(2,run_index-1)
          new_x = prev_x - x_calibration_step
          star%job%rescale_params(2,run_index+1) = new_x
          star%job%initial_x_array(run_index+1) = new_x
          star%job%initial_x_array(run_index+2) = new_x
!
      write(*,*) ' NK=2, Y=Y+0.01, Setup next run, X=', new_x
      write(run_log_unit,*) '#NK=2, Y=Y+0.01, Setup next run, X=', new_x
            return
       else
!           If NK=4,6,8,... (second and more times through) then
!           Use current and previous values of L at R and X to calculate
!           dX/dlogL. Save L.  Start next run.
            new_x = star%job%rescale_params(2, run_index-1)
          prev_x = star%job%rescale_params(2, run_index-3)
            dx_dlogl = (new_x-prev_x)/ &
                 (star%log_l_at_target_radius-star%log_l_at_target_radius_prev_run)
          new_x = dx_dlogl*(log10(star%ctrl%target_luminosity_lsun)- &
               star%log_l_at_target_radius)+new_x
!
      write(*,*) ' Setup next run, NK, X =', run_index, new_x
      write(run_log_unit, *) ' Setup next run, NK, X =', run_index, new_x
          star%log_l_at_target_radius_prev_run = star%log_l_at_target_radius
          star%job%rescale_params(2,run_index+1) = new_x
          star%job%initial_x_array(run_index+1) = new_x
          star%job%initial_x_array(run_index+2) = new_x
       end if
      endif
end subroutine check_star_calibration

! ---------------------------------------------------------------
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
      use controls_buffer, only: max_runs
      implicit none

      double precision, intent(in) :: age_scale_factor
! --- locals ---
      integer :: i, j
! 16 three-run calibration cycles (the run arrays in controls_buffer are
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

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! setscal
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setscal.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Sets up the run-list to evolve a model to a target point on the
! HR diagram (target luminosity plus either target Teff or target
! R/Rsun, per star%ctrl%specify_teff_flag): expands the two
! seed runs into 50 runs by copying their rescaling/mixing-length
! parameters forward (odd runs from 3 rerun the seed rescaling; even
! runs from 4 reuse run 2's age/timestep stop criteria), mirroring
! setup_solar_calibration.f90's solar-calibration run-list expansion.
subroutine setup_star_calibration

      use star_info_lib, only: star
      use controls_buffer, only: max_runs
      use io_units
      use physical_constants
      use math_lib
      implicit none
! --- locals ---
      integer :: i, j

      star%star_found_flag = .false.
      star%just_passed_target_radius_flag = .false.
      if (star%ctrl%specify_teff_flag) then
         star%job%target_radius_rsun = sqrt(star%ctrl%target_luminosity_lsun* &
              star%solar_luminosity_cgs/(c4pi*csig))/(star%job%target_teff*star%job%target_teff* &
              star%solar_radius_cgs)
      else
         star%job%target_teff = pow((star%ctrl%target_luminosity_lsun*star%solar_luminosity_cgs)/ &
              (c4pi*csig*star%job%target_radius_rsun*star%job%target_radius_rsun* &
              star%solar_radius_cgs*star%solar_radius_cgs), 0.25d0)
      end if
      star%log_r_prev_model = 0
!     SET UP RUN TO EVOLVE TO L, Teff IN HR-DIAGRM.
!     THIS CONSISTS OF SETTING THE NUMBER OF RUNS TO THE MAXIMUM
!     (max_runs = 50), AND COPYING THE RELEVANT PARAMETERS FROM THE
!     FIRST TWO RUNS TO THE 48 CALIBRATING RUNS THAT FOLLOW.
      star%job%num_runs = max_runs
      do i = 2,max_runs
         star%job%initial_x_array(i) = star%job%initial_x_array(1)
         star%job%initial_z_array(i) = star%job%initial_z_array(1)
         star%job%mixing_length_array(i) = star%job%mixing_length_array(1)
         star%job%has_senv0_array(i) = star%job%has_senv0_array(1)
         star%job%senv0_array(i) = star%job%senv0_array(1)
      end do
      do i = 3,max_runs-1,2
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
      do i = 4,max_runs,2
         star%job%rescale_kind(i) = 1
         star%job%first_call_flag(i) = .false.
         star%job%num_models(i) = star%job%num_models(2)
         star%job%target_end_age(i) = star%job%target_end_age(2)
         star%job%end_age_stop_active(i) = star%job%end_age_stop_active(2)
         star%job%timestep_override(i) = star%job%timestep_override(2)
         star%job%timestep_override_active(i) = star%job%timestep_override_active(2)
      end do
      write(*,*) ' Evolve to R*, L* = ', star%job%target_radius_rsun, &
           star%ctrl%target_luminosity_lsun
      write(run_log_unit,*) '#Evolve to R*, L* = ', star%job%target_radius_rsun, &
           star%ctrl%target_luminosity_lsun
      return
end subroutine setup_star_calibration

end module calibration

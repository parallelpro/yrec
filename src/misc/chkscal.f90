!----------------------------------------------------------------------
! chkscal
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original chkscal.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Checks whether the evolving track has passed a target stellar radius
! (R*), interpolates the luminosity and age at that radius, and either
! signals a converged calibration (log L/Lsun at R* matches the target
! luminosity to within tolerance) or sets up the envelope hydrogen
! fraction X for the next trial run via linear/secant interpolation in
! X vs. log L at R*.
! ONLY CALLED FOR EVEN NK, ASSUMES RESCALING ON ODD NK AND EVOLVING
! ON EVEN NK
subroutine chkscal(log_l_lsun, log_teff, current_age, run_index)

      use luout_lib
      use const_lib
      implicit none

      double precision, intent(in) :: log_l_lsun, log_teff, current_age
      integer, intent(in) :: run_index




! common/ckind/: rescale_params is used here (RESCAL); num_models/
! rescale_kind/first_call_flag/num_runs are unused placeholders. Naming
! matches wrthead.f90/wrtmonte.f90/chkcal.f90.
      double precision :: rescale_params(4,50)
      integer :: num_models(50), rescale_kind(50)
      logical :: first_call_flag(50)
      integer :: num_runs
      common /ckind/ rescale_params, num_models, rescale_kind, &
           first_call_flag, num_runs

!      COMMON/SETT/ENDAGE(50),SETDT(50),LENDAG(50),LSETDT(50)

! common/newxym/: only initial_x_array is used here. Naming matches
! wrtmonte.f90/chkcal.f90.
      double precision :: initial_x_array(50), initial_z_array(50), &
           mixing_length_array(50)
      logical :: has_senv0_array(50)
      double precision :: senv0_array(50)
      common /newxym/ initial_x_array, initial_z_array, &
           mixing_length_array, has_senv0_array, senv0_array

! common/calstar/: target-star calibration state; all members used
! here. Not referenced in any already-converted file.
      double precision :: target_luminosity_lsun, luminosity_tolerance, &
           target_teff, target_radius_rsun, log_l_prev_model, &
           log_r_prev_model, age_at_target_radius, &
           log_l_at_target_radius, log_l_at_target_radius_prev_run, &
           age_prev_model
      logical :: star_found_flag, specify_teff_flag, &
           just_passed_target_radius_flag, calibrate_star_flag
      common/calstar/ target_luminosity_lsun, luminosity_tolerance, &
           target_teff, target_radius_rsun, &
           log_l_prev_model, log_r_prev_model, age_at_target_radius, &
           log_l_at_target_radius, log_l_at_target_radius_prev_run, &
           age_prev_model, star_found_flag, specify_teff_flag, &
           just_passed_target_radius_flag, calibrate_star_flag


      save

! locals
      double precision :: teff_current, log_r_rsun_current
      double precision :: dlogl_dlogr, dage_dlogr
      double precision :: new_x, prev_x, dx_dlogl

!     LSTAR     T - have got a star at Teff and L
!     LPASSR    T - on run have just passed Teff
!     XLS       Luminosity (L/Lsun) wanted by adjusting Y
!     XLSTOL    tolerance wanted for luminosity
!     LTEFF     T - specify Teff for star
!               F - specify R/Rsun for star
!     STEFF     Effective temperature of star (K) or...
!     SR        Radius of star (R/Rsun)
!     TEFF      Teff of current model
!     ALR       log(R/Rsun) of current model
!     ALRI      log(R/Rsun) of previous model
!     DAGE      age of current model (Gyr)
!     AGEI      age of previous model (Gyr)
!     AGER      age of model at R*
!     BL        luminosity of current model log(L/Lsun)
!     BLI       luminosity of previous model
!     BLR       luminosity of model at R
!     BLRP      luminosity of model at R* of previous run
!     XP        X of previous run = RESCAL(2, NK-1)
!
!     Check if star has passed R*.
!     If not store L and age and return.
      just_passed_target_radius_flag=.false.
      teff_current = 10.0d0**log_teff
      log_r_rsun_current = sqrt((10.0d0**log_l_lsun)*solar_luminosity_cgs/ &
           (c4pi*csig))/(teff_current*teff_current*solar_radius_cgs)
      if(log_r_rsun_current.gt.log_r_prev_model) then
         if(.not.(log_r_rsun_current.gt.target_radius_rsun.and. &
              log_r_prev_model.le.target_radius_rsun)) then
              log_l_prev_model = log_l_lsun
              age_prev_model = current_age
              log_r_prev_model = log_r_rsun_current
              return
          end if
      else
          if (.not.(log_r_rsun_current.lt.target_radius_rsun.and. &
               log_r_prev_model.ge.target_radius_rsun)) then
              log_l_prev_model = log_l_lsun
              age_prev_model = current_age
              log_r_prev_model = log_r_rsun_current
              return
           end if
      endif
!
!     Check if track has passed through Teff the right number of
!     times. If not store L and age and return.
! ZZZ
      write(*,*) ' Just passed R*'
      write(itrack,*) '#Just passed R*'
      just_passed_target_radius_flag = .true.
!
!     Have previous L,Age and current L,Age (one before R* and
!     one after R*).  Interpolate to get L,Age at R*
      dlogl_dlogr = (log_l_lsun-log_l_prev_model)/ &
           (log_r_rsun_current-log_r_prev_model)
      log_l_at_target_radius = log_l_lsun + &
           dlogl_dlogr*(target_radius_rsun-log_r_rsun_current)
      dage_dlogr = (current_age-age_prev_model)/ &
           (log_r_rsun_current-log_r_prev_model)
      age_at_target_radius = current_age + &
           dage_dlogr*(target_radius_rsun-log_r_rsun_current)
      write(*,*) ' X, LogL/Lsun at R* =', rescale_params(2,run_index-1), &
           log_l_at_target_radius
      write(itrack,*) '#X, LogL/Lsun at R* =', &
           rescale_params(2,run_index-1), log_l_at_target_radius
      if (abs(10.0d0**log_l_at_target_radius-target_luminosity_lsun) &
           .le. luminosity_tolerance) then
!        Get here then have track that passes through specified
!        L and R. Use age at R for final
!        run to stop at that age. Do one more run
!        stopping at interpolated age.
         star_found_flag=.true.
       end_age_stop_active(run_index+1) = .true.
       target_end_age(run_index+1) = age_at_target_radius*1.0d9
       end_age_stop_active(run_index+2) = .true.
       target_end_age(run_index+2) = age_at_target_radius*1.0d9
       rescale_params(2,run_index+1) = rescale_params(2,run_index-1)
       initial_x_array(run_index+1) = rescale_params(2,run_index+1)
       initial_x_array(run_index+2) = initial_x_array(run_index+1)
!
      write(*,*) ' Have hit R* & L*, prepare final run to age:', &
        age_at_target_radius
      write(itrack,*)'#Have hit R* & L*, prepare final run to age:', &
        age_at_target_radius
         return
      else
         if (run_index .eq. 2) then
!           First time through. Save L and X at R*.
!           Add 0.01 to Y. Start next run.
            log_l_at_target_radius_prev_run = log_l_at_target_radius
          prev_x = rescale_params(2,run_index-1)
          new_x = prev_x - 0.01d0
          rescale_params(2,run_index+1) = new_x
          initial_x_array(run_index+1) = new_x
          initial_x_array(run_index+2) = new_x
!
      write(*,*) ' NK=2, Y=Y+0.01, Setup next run, X=', new_x
      write(itrack,*) '#NK=2, Y=Y+0.01, Setup next run, X=', new_x
            return
       else
!           If NK=4,6,8,... (second and more times through) then
!           Use current and previous values of L at R and X to calculate
!           dX/dlogL. Save L.  Start next run.
            new_x = rescale_params(2, run_index-1)
          prev_x = rescale_params(2, run_index-3)
            dx_dlogl = (new_x-prev_x)/ &
                 (log_l_at_target_radius-log_l_at_target_radius_prev_run)
          new_x = dx_dlogl*(log10(target_luminosity_lsun)- &
               log_l_at_target_radius)+new_x
!
      write(*,*) ' Setup next run, NK, X =', run_index, new_x
      write(itrack, *) ' Setup next run, NK, X =', run_index, new_x
          log_l_at_target_radius_prev_run = log_l_at_target_radius
          rescale_params(2,run_index+1) = new_x
          initial_x_array(run_index+1) = new_x
          initial_x_array(run_index+2) = new_x
       end if
      endif
      return
end subroutine chkscal

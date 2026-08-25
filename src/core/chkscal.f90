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

      use star_info_lib, only: star
      use luout_lib
      use const_lib
      implicit none

      double precision, intent(in) :: log_l_lsun, log_teff, current_age
      integer, intent(in) :: run_index





!      COMMON/SETT/ENDAGE(50),SETDT(50),LENDAG(50),LSETDT(50)
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
      log_r_rsun_current = sqrt((10.0d0**log_l_lsun)*star%solar_luminosity_cgs/ &
           (c4pi*csig))/(teff_current*teff_current*star%solar_radius_cgs)
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
      write(*,*) ' X, LogL/Lsun at R* =', star%ctrl%rescale_params(2,run_index-1), &
           log_l_at_target_radius
      write(itrack,*) '#X, LogL/Lsun at R* =', &
           star%ctrl%rescale_params(2,run_index-1), log_l_at_target_radius
      if (abs(10.0d0**log_l_at_target_radius-star%ctrl%target_luminosity_lsun) &
           .le. star%ctrl%target_star_luminosity_tolerance) then
!        Get here then have track that passes through specified
!        L and R. Use age at R for final
!        run to stop at that age. Do one more run
!        stopping at interpolated age.
         star_found_flag=.true.
       star%ctrl%end_age_stop_active(run_index+1) = .true.
       star%ctrl%target_end_age(run_index+1) = age_at_target_radius*1.0d9
       star%ctrl%end_age_stop_active(run_index+2) = .true.
       star%ctrl%target_end_age(run_index+2) = age_at_target_radius*1.0d9
       star%ctrl%rescale_params(2,run_index+1) = star%ctrl%rescale_params(2,run_index-1)
       star%ctrl%initial_x_array(run_index+1) = star%ctrl%rescale_params(2,run_index+1)
       star%ctrl%initial_x_array(run_index+2) = star%ctrl%initial_x_array(run_index+1)
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
          prev_x = star%ctrl%rescale_params(2,run_index-1)
          new_x = prev_x - 0.01d0
          star%ctrl%rescale_params(2,run_index+1) = new_x
          star%ctrl%initial_x_array(run_index+1) = new_x
          star%ctrl%initial_x_array(run_index+2) = new_x
!
      write(*,*) ' NK=2, Y=Y+0.01, Setup next run, X=', new_x
      write(itrack,*) '#NK=2, Y=Y+0.01, Setup next run, X=', new_x
            return
       else
!           If NK=4,6,8,... (second and more times through) then
!           Use current and previous values of L at R and X to calculate
!           dX/dlogL. Save L.  Start next run.
            new_x = star%ctrl%rescale_params(2, run_index-1)
          prev_x = star%ctrl%rescale_params(2, run_index-3)
            dx_dlogl = (new_x-prev_x)/ &
                 (log_l_at_target_radius-log_l_at_target_radius_prev_run)
          new_x = dx_dlogl*(log10(star%ctrl%target_luminosity_lsun)- &
               log_l_at_target_radius)+new_x
!
      write(*,*) ' Setup next run, NK, X =', run_index, new_x
      write(itrack, *) ' Setup next run, NK, X =', run_index, new_x
          log_l_at_target_radius_prev_run = log_l_at_target_radius
          star%ctrl%rescale_params(2,run_index+1) = new_x
          star%ctrl%initial_x_array(run_index+1) = new_x
          star%ctrl%initial_x_array(run_index+2) = new_x
       end if
      endif
      return
end subroutine chkscal

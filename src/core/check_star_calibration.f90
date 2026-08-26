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
subroutine check_star_calibration(log_l_lsun, log_teff, current_age, run_index)

      use star_info_lib, only: star
      use luout_lib
      use phys_const_lib
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
      star%just_passed_target_radius_flag=.false.
      teff_current = 10.0d0**log_teff
      log_r_rsun_current = sqrt((10.0d0**log_l_lsun)*star%solar_luminosity_cgs/ &
           (c4pi*csig))/(teff_current*teff_current*star%solar_radius_cgs)
      if(log_r_rsun_current.gt.star%log_r_prev_model) then
         if(.not.(log_r_rsun_current.gt.star%job%target_radius_rsun.and. &
              star%log_r_prev_model.le.star%job%target_radius_rsun)) then
              star%log_l_prev_model = log_l_lsun
              star%age_prev_model = current_age
              star%log_r_prev_model = log_r_rsun_current
              return
          end if
      else
          if (.not.(log_r_rsun_current.lt.star%job%target_radius_rsun.and. &
               star%log_r_prev_model.ge.star%job%target_radius_rsun)) then
              star%log_l_prev_model = log_l_lsun
              star%age_prev_model = current_age
              star%log_r_prev_model = log_r_rsun_current
              return
           end if
      endif
!
!     Check if track has passed through Teff the right number of
!     times. If not store L and age and return.
! ZZZ
      write(*,*) ' Just passed R*'
      write(short_file_unit,*) '#Just passed R*'
      star%just_passed_target_radius_flag = .true.
!
!     Have previous L,Age and current L,Age (one before R* and
!     one after R*).  Interpolate to get L,Age at R*
      dlogl_dlogr = (log_l_lsun-star%log_l_prev_model)/ &
           (log_r_rsun_current-star%log_r_prev_model)
      star%log_l_at_target_radius = log_l_lsun + &
           dlogl_dlogr*(star%job%target_radius_rsun-log_r_rsun_current)
      dage_dlogr = (current_age-star%age_prev_model)/ &
           (log_r_rsun_current-star%log_r_prev_model)
      star%age_at_target_radius = current_age + &
           dage_dlogr*(star%job%target_radius_rsun-log_r_rsun_current)
      write(*,*) ' X, LogL/Lsun at R* =', star%job%rescale_params(2,run_index-1), &
           star%log_l_at_target_radius
      write(short_file_unit,*) '#X, LogL/Lsun at R* =', &
           star%job%rescale_params(2,run_index-1), star%log_l_at_target_radius
      if (abs(10.0d0**star%log_l_at_target_radius-star%ctrl%target_luminosity_lsun) &
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
      write(short_file_unit,*)'#Have hit R* & L*, prepare final run to age:', &
        star%age_at_target_radius
         return
      else
         if (run_index .eq. 2) then
!           First time through. Save L and X at R*.
!           Add 0.01 to Y. Start next run.
            star%log_l_at_target_radius_prev_run = star%log_l_at_target_radius
          prev_x = star%job%rescale_params(2,run_index-1)
          new_x = prev_x - 0.01d0
          star%job%rescale_params(2,run_index+1) = new_x
          star%job%initial_x_array(run_index+1) = new_x
          star%job%initial_x_array(run_index+2) = new_x
!
      write(*,*) ' NK=2, Y=Y+0.01, Setup next run, X=', new_x
      write(short_file_unit,*) '#NK=2, Y=Y+0.01, Setup next run, X=', new_x
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
      write(short_file_unit, *) ' Setup next run, NK, X =', run_index, new_x
          star%log_l_at_target_radius_prev_run = star%log_l_at_target_radius
          star%job%rescale_params(2,run_index+1) = new_x
          star%job%initial_x_array(run_index+1) = new_x
          star%job%initial_x_array(run_index+2) = new_x
       end if
      endif
      return
end subroutine check_star_calibration

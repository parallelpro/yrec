!----------------------------------------------------------------------
! pdist
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original pdist.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Calculates the distance the star has travelled in the HR diagram
! (in a weighted log-L/log-Teff/age metric) since the last pulsation
! model was written; if far enough, opens a fresh set of pulsation
! model/envelope/atmosphere output files (named from pulse_mod_path/
! pulse_env_path/pulse_atm_path plus a run/model-number suffix) and
! sets pulsation_output_active so the caller writes to them; otherwise
! closes the previous set of files.
! MHP 8/25 passed file names directly and removed common block with
! file names.
subroutine pdist(prev_log_l, prev_log_teff, prev_age, path_length_sq, &
     log_luminosity, log_teff, model_number, pulse_atm_path, &
     pulse_env_path, pulse_mod_path)
      use star_info_lib, only: star
      use const_lib
      use luout_lib
      implicit none

      double precision, intent(inout) :: prev_log_l, prev_log_teff, &
           prev_age, path_length_sq
      double precision, intent(in) :: log_luminosity, log_teff
      integer, intent(in) :: model_number
      character(len=256), intent(in) :: pulse_atm_path, pulse_env_path, &
           pulse_mod_path

      character(len=256) :: output_path, temp_string
      double precision :: delta_log_l, delta_log_teff, delta_age
      double precision :: weighted_l_sq, weighted_teff_sq, weighted_age_sq
      integer :: trim_col

      delta_log_l = log_luminosity-prev_log_l
      delta_log_teff = log_teff-prev_log_teff
      delta_age = star%run%dage-prev_age
      weighted_l_sq = po_weight_l*delta_log_l*po_weight_l*delta_log_l
      weighted_teff_sq = po_weight_teff*delta_log_teff*po_weight_teff* &
           delta_log_teff
      weighted_age_sq = po_weight_age*delta_age*po_weight_age*delta_age
      path_length_sq = path_length_sq+weighted_l_sq+weighted_teff_sq+ &
           weighted_age_sq
      prev_log_l = log_luminosity
      prev_log_teff = log_teff
      prev_age = star%run%dage
      if (path_length_sq .ge. po_max_len_sq) then
          pulsation_output_active = .true.
          path_length_sq = 0.0d0
          trim_col = index(pulse_mod_path,' ')-1
          if (model_number.lt.10000) then
             write(temp_string,'(I2.2,''_'',I4.4)') nk, model_number
             output_path = pulse_mod_path(1:trim_col) // temp_string(1:7)
          else
             write(temp_string,'(I2.2,''_'',I5.5)') nk, model_number
             output_path = pulse_mod_path(1:trim_col) // temp_string(1:8)
          end if
          open(opal_model_unit,file=output_path,status='NEW',form='FORMATTED')
          trim_col = index(pulse_env_path,' ')-1
          if (model_number.lt.10000) then
             write(temp_string,'(I2.2,''_'',I4.4)') nk, model_number
             output_path = pulse_env_path(1:trim_col) // temp_string(1:7)
          else
             write(temp_string,'(I2.2,''_'',I5.5)') nk, model_number
             output_path = pulse_env_path(1:trim_col) // temp_string(1:8)
          end if
          open(opal_envelope_unit,file=output_path,status='NEW',form='FORMATTED')
          trim_col = index(pulse_atm_path,' ')-1
          if (model_number.lt.10000) then
             write(temp_string,'(I2.2,''_'',I4.4)') nk, model_number
             output_path = pulse_atm_path(1:trim_col) // temp_string(1:7)
          else
             write(temp_string,'(I2.2,''_'',I5.5)') nk, model_number
             output_path = pulse_atm_path(1:trim_col) // temp_string(1:8)
          end if
          open(opal_atm_unit,file=output_path,status='NEW',form='FORMATTED')
       else
          close(opal_model_unit)
          close(opal_envelope_unit)
          close(opal_atm_unit)
          pulsation_output_active = .false.
       end if
       write(iowr,276) path_length_sq,delta_log_l,delta_log_teff,delta_age
       write(short_file_unit,276) path_length_sq,weighted_l_sq, &
            weighted_teff_sq,weighted_age_sq
 276   format(5X,' LEN SQ=',1PE10.3, '  D L=',1PE10.3, &
            '  D Te=',1PE10.3,'  D AGE=',1PE10.3)

       return
end subroutine pdist

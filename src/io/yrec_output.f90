!----------------------------------------------------------------------
! yrec_output
!----------------------------------------------------------------------
! Added 2026 (MESA-style output centralization). The output-mode
! branch (controls_lib's use_legacy_output) has exactly three
! decision points, all routed here or clearly wrapped at the source:
!
!   parmin's legacy open-block   -- wrapped in use_legacy_output; its
!                                   MESA else-branch calls
!                                   output_init_mesa and forces off
!                                   the flags behind physics-time
!                                   legacy streams (envint's profile
!                                   blocks, legacy OPAL pulse files).
!   output_run_header            -- per run: legacy wrthead banner or
!                                   nothing.
!   output_write_model           -- per converged model: legacy wrtout
!                                   or MESA-mode gettau +
!                                   write_history + GYRE trigger.
!
! run_yrec's final-model store (wrtlst/putstore) and .snu calibration
! writes are single-purpose legacy blocks wrapped in place with the
! same flag. MESA-mode output today: CASE.history (MESA history.data
! layout), CASE.log (all short_file_unit diagnostics -- ~40 files
! across the code write to that unit, so it is retargeted rather than
! left unopened), GYRE pulse files when pulse_gyre_interval > 0.
module yrec_output
      use const_lib
      implicit none
      private
      public :: output_init_mesa, output_run_header, output_write_model

contains

! MESA-mode job initialization: no legacy output files exist; the
! shared diagnostics unit becomes CASE.log (derived from the .short
! path, which every deck defines).
subroutine output_init_mesa(fshort)
      use luout_lib
      character(len=*), intent(in) :: fshort
      character(len=256) :: log_path
      integer :: n

      n = index(fshort, '.short')
      if (n > 0) then
         log_path = fshort(1:n) // 'log'
      else
         log_path = trim(fshort) // '.log'
      end if
      open(short_file_unit, file=log_path, form='FORMATTED', &
           status='REPLACE')
end subroutine output_init_mesa

subroutine output_run_header(star_mass_msun)
      double precision, intent(in) :: star_mass_msun

      if (use_legacy_output) then
         call wrthead(star_mass_msun)
      end if
end subroutine output_run_header

subroutine output_write_model(timestep_yr, log_gravity, has_h_shell, &
     h_shell_begin_index, h_shell_end_index, h_shell_mid_index, &
     trial_sign_flag, punch_pending_flag, total_angular_momentum, &
     total_rotational_kinetic_energy)
      use star_info_lib, only: star
      use luout_lib
      double precision, intent(in) :: timestep_yr
      double precision, intent(out) :: log_gravity
      logical, intent(in) :: has_h_shell
      integer, intent(in) :: h_shell_begin_index, h_shell_end_index, &
           h_shell_mid_index
      double precision, intent(in) :: trial_sign_flag
      logical, intent(inout) :: punch_pending_flag
      double precision, intent(in) :: total_angular_momentum, &
           total_rotational_kinetic_energy

      if (use_legacy_output) then
         call wrtout(timestep_yr, log_gravity, has_h_shell, &
              h_shell_begin_index, h_shell_end_index, h_shell_mid_index, &
              trial_sign_flag, punch_pending_flag, total_angular_momentum, &
              total_rotational_kinetic_energy)
      else
! Every history quantity was computed by update_output_diagnostics
! (which also ran gettau in MESA mode) and stored in star_info;
! the writers below are pure readers. wrtout computes log_gravity as
! an output on the legacy path; hand back the stored value here.
         log_gravity = star%run%log_g_surface
         call write_history()
! The GYRE writer is the MESA-mode pulsation mechanism (wrtout
! triggers it on the legacy path).
         if (pulse_gyre_interval > 0) then
            if (mod(star%model_number, pulse_gyre_interval) == 0) then
               call write_gyre_pulse(star%nz, star%model_number, star%m, &
                    star%logRho, star%luminosity_lsun, star%logP, &
                    star%logR, star%logT, star%omega)
            end if
         end if
! Keep the log live during the run, like the history file.
         flush(short_file_unit)
      end if
end subroutine output_write_model

end module yrec_output

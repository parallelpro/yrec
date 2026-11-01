!----------------------------------------------------------------------
! run_log
!----------------------------------------------------------------------
! The compact MESA-style progress line and the gate for the solver
! diagnostics. The run log holds a run header, warnings and events,
! and one progress line every terminal_interval models (the final
! model of a card always prints); the physical summary is in the
! history file, and the solver diagnostics (Henyey iteration trace,
! envelope/atmosphere integration statistics, triangle/rezoning
! bookkeeping) print only with report_solver_diagnostics.
module run_log
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: solver_diagnostics, set_solver_diagnostics, log_model_line, &
           log_final_model_line
      public :: log_run_summary
      public :: progress_lines

! Reprint the column header every so many printed lines (MESA's
! write_header_frequency).
      integer, parameter :: header_every = 10

! The progress lines' state, per job (part of yrec_output's
! output_state, which run_yrec owns): the lines printed since the
! header, and the last model printed, so the same model never prints
! twice (the final-model call may coincide with an interval line).
      type :: progress_lines
         integer :: lines_since_header = header_every
         integer :: last_printed_model = -1
      end type progress_lines

! the control report_solver_diagnostics (set_solver_diagnostics), which
! star-free code asks through solver_diagnostics()
      logical :: solver_diagnostics_on = .false.

contains

! The control report_solver_diagnostics, as the controls read installs
! it (set_solver_diagnostics), so that star-free code can ask.
subroutine set_solver_diagnostics(on)
      logical, intent(in) :: on
      solver_diagnostics_on = on
end subroutine set_solver_diagnostics

! The gate for verbose solver forensics in the run log.
logical function solver_diagnostics()
      solver_diagnostics = solver_diagnostics_on
end function solver_diagnostics

! The compact per-model progress line, to the terminal and the run
! log, throttled by terminal_interval (0 = off; final always prints).
!   star reads: model_number, physics%terminal_interval
!   star writes: none
subroutine log_model_line(s, progress)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(progress_lines), intent(inout) :: progress
      if (s%physics%terminal_interval <= 0) return
      if (mod(s%model_number, s%physics%terminal_interval) /= 0) return
      call write_model_line(s, progress)
end subroutine log_model_line

! The end-of-card call: the last converged model always prints.
!   star reads: none
!   star writes: none
subroutine log_final_model_line(s, progress)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(progress_lines), intent(inout) :: progress
      call write_model_line(s, progress)
end subroutine log_final_model_line

!   star reads: age_gyr, log_L_lsun, log_R_surface_rsun, log_Teff,
!     model_number, newton_iterations, nz, timestep_yr, xa
!   star writes: none
subroutine write_model_line(s, progress)
      use star_info_lib, only: star_info, i_h1
      use io_units, only: run_log_unit, terminal_unit
      type(star_info), intent(in) :: s
      type(progress_lines), intent(inout) :: progress
      character(len=*), parameter :: header = &
           '   model  zones         age_Gyr        dt_yr   log_Teff' // &
           '      log_L      log_R center_h1  iters'
      character(len=*), parameter :: line_format = &
           '(1x,i7,1x,i6,1x,es15.8,1x,es12.5,3(1x,f10.6),2x,f8.6,1x,i6)'
      if (s%model_number == progress%last_printed_model) return
      if (progress%lines_since_header >= header_every) then
         write(terminal_unit,'(a)') header
         write(run_log_unit,'(a)') header
         progress%lines_since_header = 0
      end if
      write(terminal_unit,line_format) s%model_number, s%nz, s%age_gyr, &
           s%timestep_yr, s%log_Teff, s%log_L_lsun, &
           s%log_R_surface_rsun, s%xa(i_h1,1), s%newton_iterations
      write(run_log_unit,line_format) s%model_number, s%nz, s%age_gyr, &
           s%timestep_yr, s%log_Teff, s%log_L_lsun, &
           s%log_R_surface_rsun, s%xa(i_h1,1), s%newton_iterations
      progress%lines_since_header = progress%lines_since_header + 1
      progress%last_printed_model = s%model_number
end subroutine write_model_line

! The end-of-run summary: why the run ended (star%termination_reason,
! set by stop_conditions and the calibration verdict) and the final
! model. The wall-clock line differs from run to run, so it goes to the
! terminal only, which keeps the run log byte-pinnable.
!   star reads: age_gyr, model_number, termination_reason
!   star writes: none
subroutine log_run_summary(s, wall_seconds)
      use star_info_lib, only: star_info
      use io_units, only: run_log_unit, terminal_unit
      type(star_info), intent(in) :: s
      real(dp), intent(in) :: wall_seconds
      character(len=*), parameter :: reason_format = "(1x,'run finished: ',a)", &
           model_format = "(3x,'final model',i7,'   age ',es15.8,' Gyr')", &
           clock_format = "(3x,'wall-clock ',f9.1,' s')"
      write(terminal_unit,'(a)') ''
      write(run_log_unit,'(a)') ''
      write(terminal_unit,reason_format) trim(s%termination_reason)
      write(run_log_unit,reason_format) trim(s%termination_reason)
      write(terminal_unit,model_format) s%model_number, s%age_gyr
      write(run_log_unit,model_format) s%model_number, s%age_gyr
      write(terminal_unit,clock_format) wall_seconds
end subroutine log_run_summary

end module run_log

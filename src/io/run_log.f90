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
      implicit none
      private
      public :: solver_diagnostics, log_model_line, log_final_model_line, log_reset
      public :: log_run_summary

! Reprint the column header every so many printed lines (MESA's
! write_header_frequency), and never print the same model twice
! (the final-model call may coincide with an interval line).
      integer, parameter :: header_every = 10
      integer :: lines_since_header = header_every
      integer :: last_printed_model = -1

contains

! Reset the per-run formatting state (in-process re-entry: called
! from output_init_mesa at every job start).
subroutine log_reset()
      lines_since_header = header_every
      last_printed_model = -1
end subroutine log_reset

! The gate for verbose solver forensics in the run log.
logical function solver_diagnostics()
      use star_info_lib, only: star
      solver_diagnostics = star%ctrl%report_solver_diagnostics
end function solver_diagnostics

! The compact per-model progress line, to the terminal and the run
! log, throttled by terminal_interval (0 = off; final always prints).
subroutine log_model_line()
      use star_info_lib, only: star
      if (star%ctrl%terminal_interval <= 0) return
      if (mod(star%model_number, star%ctrl%terminal_interval) /= 0) return
      call write_model_line()
end subroutine log_model_line

! The end-of-card call: the last converged model always prints.
subroutine log_final_model_line()
      call write_model_line()
end subroutine log_final_model_line

subroutine write_model_line()
      use star_info_lib, only: star, i_h1
      use io_units
      character(len=*), parameter :: header = &
           '   model  zones         age_Gyr        dt_yr   log_Teff' // &
           '      log_L      log_R center_h1  iters'
      if (star%model_number == last_printed_model) return
      if (lines_since_header >= header_every) then
         write(terminal_unit,'(a)') header
         write(run_log_unit,'(a)') header
         lines_since_header = 0
      end if
      write(terminal_unit,10) star%model_number, star%nz, star%dage, &
           star%timestep_yr, star%log_Teff, star%log_L, &
           star%log_R_surface, star%xa(i_h1,1), star%newton_iterations
      write(run_log_unit,10) star%model_number, star%nz, star%dage, &
           star%timestep_yr, star%log_Teff, star%log_L, &
           star%log_R_surface, star%xa(i_h1,1), star%newton_iterations
   10 format(1x,i7,1x,i6,1x,es15.8,1x,es12.5,3(1x,f10.6),2x,f8.6,1x,i6)
      lines_since_header = lines_since_header + 1
      last_printed_model = star%model_number
end subroutine write_model_line

! End-of-run summary (2026, run-log verbosity item): why the run
! ended (star%termination_reason, set by stop_conditions and the
! calibration verdict) plus the final model. The wall-clock line is
! nondeterministic and goes to the TERMINAL ONLY, keeping the run
! log byte-pinnable.
subroutine log_run_summary(wall_seconds)
      use star_info_lib, only: star
      use io_units
      double precision, intent(in) :: wall_seconds
      write(terminal_unit,'(a)') ''
      write(run_log_unit,'(a)') ''
      write(terminal_unit,10) trim(star%termination_reason)
      write(run_log_unit,10) trim(star%termination_reason)
      write(terminal_unit,20) star%model_number, star%dage
      write(run_log_unit,20) star%model_number, star%dage
      write(terminal_unit,30) wall_seconds
   10 format(1x,'run finished: ',a)
   20 format(3x,'final model',i7,'   age ',es15.8,' Gyr')
   30 format(3x,'wall-clock ',f9.1,' s')
end subroutine log_run_summary

end module run_log

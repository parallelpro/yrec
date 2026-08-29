!----------------------------------------------------------------------
! run_log_lib
!----------------------------------------------------------------------
! Added 2026 (run-log redesign): the compact MESA-style progress
! line and the verbosity gate for the solver forensics.
!
! The default run log is small: a run header, warnings/events, and
! one progress line per terminal_interval models (the final model of
! a card always prints). Everything the old .short printed per model
! beyond that either lives in history.data (the physical summary --
! deleted outright) or is solver forensics (Henyey iteration trace,
! envelope/atmosphere integration statistics, triangle/rezoning
! bookkeeping) now gated behind report_solver_diagnostics.
module run_log_lib
      implicit none
      private
      public :: solver_diagnostics, log_model_line, log_final_model_line

! Reprint the column header every so many printed lines (MESA's
! write_header_frequency), and never print the same model twice
! (the final-model call may coincide with an interval line).
      integer, parameter :: header_every = 10
      integer :: lines_since_header = header_every
      integer :: last_printed_model = -1

contains

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
      use luout_lib
      character(len=*), parameter :: header = &
           '   model  zones         age_Gyr        dt_yr   log_Teff' // &
           '      log_L      log_R center_h1  iters'
      if (star%model_number == last_printed_model) return
      if (lines_since_header >= header_every) then
         write(iowr,'(a)') header
         write(short_file_unit,'(a)') header
         lines_since_header = 0
      end if
      write(iowr,10) star%model_number, star%nz, star%dage, &
           star%timestep_yr, star%log_Teff, star%log_L, &
           star%log_R_surface, star%xa(i_h1,1), star%newton_iterations
      write(short_file_unit,10) star%model_number, star%nz, star%dage, &
           star%timestep_yr, star%log_Teff, star%log_L, &
           star%log_R_surface, star%xa(i_h1,1), star%newton_iterations
   10 format(1x,i7,1x,i6,1x,es15.8,1x,es12.5,3(1x,f10.6),2x,f8.6,1x,i6)
      lines_since_header = lines_since_header + 1
      last_printed_model = star%model_number
end subroutine write_model_line

end module run_log_lib

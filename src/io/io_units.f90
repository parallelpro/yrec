!----------------------------------------------------------------------
! io_units
!----------------------------------------------------------------------
! The units of the run's files, and check_io, the report for a failed
! open, read or write, and report_error, the report of any other error.
! The three fixed units written from everywhere, assigned at each read
! (run_files' assign_unit_numbers):
!   run_log_unit    -- the run log (log_output_file, default run.log)
!   terminal_unit   -- stdout (always 6; kept as a named unit only
!                      because ~20 files write through it)
!   last_model_unit -- the final .mod model (last_model_file)
! The files the controls read opens (run_files' open_run_files, with
! newunit; -1 until then):
!   starting_model_unit -- the starting model (first_model_file)
!   dynamics_unit       -- the Monte Carlo samples (dynamics_file)
!   mc_summary_unit, mc_model_unit -- the Monte Carlo summary and model
!                          files (monte_carlo_file1_path, _file2_path)
!   isochrone_unit      -- the isochrone file (isochrone_file)
! The table readers open their tables themselves, with newunit.
module io_units
      implicit none
      private
      public :: last_model_unit, run_log_unit, terminal_unit, check_io, &
           report_error
      public :: starting_model_unit, dynamics_unit, mc_summary_unit, &
           mc_model_unit, isochrone_unit

      integer :: last_model_unit, run_log_unit, terminal_unit
      integer :: starting_model_unit = -1, dynamics_unit = -1, &
           mc_summary_unit = -1, mc_model_unit = -1, isochrone_unit = -1

contains

! ---------------------------------------------------------------
! After an open, read or write with iostat=ios and iomsg=message: on
! failure, report "cannot <action> <path>: <message>" on the terminal,
! and in the run log once it is open; close close_unit if given; and
! set ierr = 1.
subroutine check_io(ios, message, action, path, ierr, close_unit)
      integer, intent(in) :: ios
      character(len=*), intent(in) :: message, action, path
      integer, intent(inout) :: ierr
      integer, intent(in), optional :: close_unit

      if (ios == 0) return
      if (present(close_unit)) close(close_unit)
      call report_error('cannot ' // action // ' ' // trim(path) // ': ' // trim(message))
      ierr = 1
end subroutine check_io

! ---------------------------------------------------------------
! An error: "ERROR: <text>" on the terminal, and in the run log once it
! is open.
subroutine report_error(text)
      character(len=*), intent(in) :: text
      logical :: log_open

      write(terminal_unit,'(1x,2a)') 'ERROR: ', text
      inquire(unit=run_log_unit, opened=log_open)
      if (log_open .and. run_log_unit /= terminal_unit) &
           write(run_log_unit,'(1x,2a)') 'ERROR: ', text
end subroutine report_error

end module io_units

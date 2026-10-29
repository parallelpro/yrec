!----------------------------------------------------------------------
! io_units
!----------------------------------------------------------------------
! The three fixed units written from everywhere, assigned at each read
! (run_files' assign_unit_numbers), and check_io, the report for a
! failed open, read or write, and report_error, the report of any other
! error. Newer readers and writers use
! open(newunit=...) instead of fixed units.
!   run_log_unit    -- the run log (log_output_file, default run.log)
!   terminal_unit   -- stdout (always 6; kept as a named unit only
!                      because ~20 files write through it)
!   last_model_unit -- the final .mod model (last_model_file)
module io_units
      implicit none
      private
      public :: last_model_unit, run_log_unit, terminal_unit, check_io, &
           report_error

      integer :: last_model_unit, run_log_unit, terminal_unit

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

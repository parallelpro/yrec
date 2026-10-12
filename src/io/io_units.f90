!----------------------------------------------------------------------
! io_units
!----------------------------------------------------------------------
! The three fixed units written from everywhere, assigned at each read
! (run_files' assign_unit_numbers), and check_io, the report for a
! failed open, read or write. Newer readers and writers use
! open(newunit=...) instead of fixed units.
!   run_log_unit    -- the run log (log_output_file, default run.log)
!   terminal_unit   -- stdout (always 6; kept as a named unit only
!                      because ~20 files write through it)
!   last_model_unit -- the final .mod model (last_model_file)
module io_units
      implicit none

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
      logical :: log_open

      if (ios == 0) return
      if (present(close_unit)) close(close_unit)
      write(terminal_unit,'(1x,6a)') 'ERROR: cannot ', action, ' ', trim(path), &
           ': ', trim(message)
      inquire(unit=run_log_unit, opened=log_open)
      if (log_open .and. run_log_unit /= terminal_unit) &
           write(run_log_unit,'(1x,6a)') 'ERROR: cannot ', action, ' ', trim(path), &
           ': ', trim(message)
      ierr = 1
end subroutine check_io

end module io_units

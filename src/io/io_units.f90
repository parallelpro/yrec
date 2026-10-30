!----------------------------------------------------------------------
! io_units
!----------------------------------------------------------------------
! The units of the run's files, and io's messages, each on the terminal
! and in the run log; a message from before the run log opens (the
! controls read's) waits for it, and open_run_log writes the waiting
! lines first (flush_pending_log), so the run log holds them in order:
!   check_io        after an open, read or write: "ERROR: cannot <action>
!                   <path>: <reason>" on failure
!   report_error    "ERROR: <text>", an error the run stops on
!   report_warning  "WARNING: <text>", a probable mistake the run goes on
!                   with (an inert control, a combination that cancels)
!   report_note     "<text>", what the read did (a mixture applied)
!   log_line        a line for the run log only (its header, the run plan)
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
           report_error, report_warning, report_note, log_line, &
           flush_pending_log, discard_pending_log
      public :: starting_model_unit, dynamics_unit, mc_summary_unit, &
           mc_model_unit, isochrone_unit

      integer :: last_model_unit, run_log_unit, terminal_unit
      integer :: starting_model_unit = -1, dynamics_unit = -1, &
           mc_summary_unit = -1, mc_model_unit = -1, isochrone_unit = -1

! The run log's lines written before it opened, in order (the most a
! read writes is its warnings, the header and a line or three per run
! card: well under max_pending_lines).
      integer, parameter :: max_pending_lines = 1000, max_line_length = 1024
      character(len=max_line_length) :: pending_line(max_pending_lines)
      integer :: pending_length(max_pending_lines)
      integer :: num_pending = 0

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
      call report('ERROR: ', text)
end subroutine report_error

! ---------------------------------------------------------------
! A warning: "WARNING: <text>", on the terminal and in the run log once
! it is open.
subroutine report_warning(text)
      character(len=*), intent(in) :: text
      call report('WARNING: ', text)
end subroutine report_warning

! ---------------------------------------------------------------
! A note: "<text>", on the terminal and in the run log once it is open.
subroutine report_note(text)
      character(len=*), intent(in) :: text
      call report('', text)
end subroutine report_note

! prefix // text, after one blank, on the terminal and in the run log
subroutine report(prefix, text)
      character(len=*), intent(in) :: prefix, text

      write(terminal_unit,'(1x,2a)') prefix, trim(text)
      call log_line(' ' // prefix // trim(text))
end subroutine report

! ---------------------------------------------------------------
! One line of the run log, as given: written if the log is open, else
! kept for flush_pending_log.
subroutine log_line(text)
      character(len=*), intent(in) :: text

      if (run_log_is_open()) then
         write(run_log_unit,'(a)') text
      else if (num_pending < max_pending_lines) then
         num_pending = num_pending + 1
         pending_line(num_pending) = text
         pending_length(num_pending) = min(len(text), max_line_length)
      end if
end subroutine log_line

! ---------------------------------------------------------------
! The lines kept while the run log was closed, written to it (open_run_log
! calls this as soon as it has opened the log).
subroutine flush_pending_log
      integer :: k

      if (.not. run_log_is_open()) return
      do k = 1, num_pending
         write(run_log_unit,'(a)') pending_line(k)(1:pending_length(k))
      end do
      num_pending = 0
end subroutine flush_pending_log

! ---------------------------------------------------------------
! The kept lines dropped (each controls read starts without any).
subroutine discard_pending_log
      num_pending = 0
end subroutine discard_pending_log

logical function run_log_is_open()
      inquire(unit=run_log_unit, opened=run_log_is_open)
      if (run_log_unit == terminal_unit) run_log_is_open = .false.
end function run_log_is_open

end module io_units

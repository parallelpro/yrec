!----------------------------------------------------------------------
! core_messages
!----------------------------------------------------------------------
! core's messages, in one place. Each route writes the lines it is given
! to its destinations:
!   core_report        to the terminal and to the run log: an error the
!                      run stops on, or a warning the run goes on with;
!   core_terminal      to the terminal only;
!   core_log           to the run log only (the per-model notes);
!   core_stdout        to standard output only: what core's write(*)
!                      meant, which is standard output even before io
!                      assigns its units (the terminal is standard error
!                      until then; run_messages);
!   core_stdout_and_log
!                      to standard output and to the run log.
! The terminal and the run log are run_messages' terminal_unit and
! run_log_unit. Before io assigns them they are standard error, as for
! every domain's messages.
!
! The caller formats the message into an array of lines with an
! internal write (one record per element), after setting every element
! to unwritten_line. The lines the write did not reach are not printed,
! so the output has the records the format makes, as a write to a unit
! would. Each line is written without its trailing blanks.
!
! A domain's computational kernels print nothing: they return a status,
! and the driver that called them reports through these routes
! (GUIDELINES "Code style" 7; tools/check_boundaries.py's no-I/O check).
module core_messages
      use, intrinsic :: iso_fortran_env, only: output_unit
      use run_messages, only: run_log_unit, terminal_unit
      implicit none
      private
      public :: core_report, core_terminal, core_log, core_stdout, &
           core_stdout_and_log

! The length of one message line. A line may hold a path of up to 256
! characters and its context.
      integer, parameter, public :: message_length = 400
      character(len=*), parameter, public :: unwritten_line = achar(0)

contains

subroutine core_report(lines)
      character(len=*), intent(in) :: lines(:)
      call write_core_lines(terminal_unit, lines)
      call write_core_lines(run_log_unit, lines)
end subroutine core_report

subroutine core_terminal(lines)
      character(len=*), intent(in) :: lines(:)
      call write_core_lines(terminal_unit, lines)
end subroutine core_terminal

subroutine core_log(lines)
      character(len=*), intent(in) :: lines(:)
      call write_core_lines(run_log_unit, lines)
end subroutine core_log

subroutine core_stdout(lines)
      character(len=*), intent(in) :: lines(:)
      call write_core_lines(output_unit, lines)
end subroutine core_stdout

subroutine core_stdout_and_log(lines)
      character(len=*), intent(in) :: lines(:)
      call write_core_lines(output_unit, lines)
      call write_core_lines(run_log_unit, lines)
end subroutine core_stdout_and_log

! ---------------------------------------------------------------
! The messages of envelope/, structure/ and step/ that have a
! procedure of their own.

! ---------------------------------------------------------------
! The messages of model/, run/ and diagnostics/ that have a
! procedure of their own.

! ---------------------------------------------------------------
subroutine write_core_lines(unit, lines)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: lines(:)
      integer :: k
      do k = 1, size(lines)
         if (lines(k)(1:1) == unwritten_line) cycle
         write(unit,'(a)') trim(lines(k))
      end do
end subroutine write_core_lines

end module core_messages

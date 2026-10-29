!----------------------------------------------------------------------
! rotation_messages
!----------------------------------------------------------------------
! Rotation's messages, in one place:
!   rotation_report  an error the run stops on, or a cut of the
!                    sub-step: to the terminal and to the run log;
!   rotation_note    a report of each sub-step (the dynamical-shear
!                    interfaces, the transport's convergence): to the
!                    terminal, and to the run log when solver diagnostics
!                    are on (run_log's solver_diagnostics).
! The caller formats the message into an array of lines with an
! internal write (one record per element) after setting every element
! to unwritten_line; the lines the write did not reach are not printed,
! so the output has the records the format makes, as a write to a unit
! would. The terminal is standard output (what rotation's write(6) and
! write(*) meant), whether or not the run has assigned io's units.
module rotation_messages
      use, intrinsic :: iso_fortran_env, only: output_unit
      use io_units, only: run_log_unit
      use run_log, only: solver_diagnostics
      implicit none
      private
      public :: rotation_report, rotation_note
      integer, parameter, public :: message_length = 160
      character(len=*), parameter, public :: unwritten_line = achar(0)

contains

subroutine rotation_report(lines)
      character(len=*), intent(in) :: lines(:)
      call write_lines(output_unit, lines)
      call write_lines(run_log_unit, lines)
end subroutine rotation_report

subroutine rotation_note(lines)
      character(len=*), intent(in) :: lines(:)
      call write_lines(output_unit, lines)
      if (solver_diagnostics()) call write_lines(run_log_unit, lines)
end subroutine rotation_note

subroutine write_lines(unit, lines)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: lines(:)
      integer :: k
      do k = 1, size(lines)
         if (lines(k)(1:1) == unwritten_line) cycle
         write(unit,'(a)') trim(lines(k))
      end do
end subroutine write_lines

end module rotation_messages

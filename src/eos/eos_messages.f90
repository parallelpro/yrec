!----------------------------------------------------------------------
! eos_messages
!----------------------------------------------------------------------
! The equation of state's messages, in one place (rotation_messages is
! the model):
!   eos_log     a warning or a diagnostic the run continues after (a
!               Saha or density iteration that did not converge, an
!               OPAL lookup outside its table): to the run log;
!   eos_report  an error eos returns ierr for: to the terminal and to
!               the run log;
!   eos_log_once  a warning the run log gets once after each eos_init
!               (reset_eos_warnings), by its id (once_*); eos_warned
!               tells whether it has been written.
! The caller formats the message into an array of lines with an
! internal write (one record per element) after setting every element
! to unwritten_line; the lines the write did not reach are not printed,
! so the output has the records the format makes, as a write to a unit
! would. The terminal is standard output.
module eos_messages
      use, intrinsic :: iso_fortran_env, only: output_unit
      use run_messages, only: run_log_unit
      implicit none
      private
      public :: eos_log, eos_report, eos_log_once, reset_eos_warnings, eos_warned
      integer, parameter, public :: message_length = 200
      character(len=*), parameter, public :: unwritten_line = achar(0)
! the warnings written once per eos_init: the fully ionized gas
! extrapolating the Fermi-Dirac table
      integer, parameter, public :: once_fermi_extrapolation = 1
      integer, parameter :: n_once = 1
      logical :: warned(n_once) = .false.

contains

subroutine eos_log(lines)
      character(len=*), intent(in) :: lines(:)
      call write_message_lines(run_log_unit, lines)
end subroutine eos_log

subroutine eos_report(lines)
      character(len=*), intent(in) :: lines(:)
      call write_message_lines(output_unit, lines)
      call write_message_lines(run_log_unit, lines)
end subroutine eos_report

subroutine eos_log_once(id, lines)
      integer, intent(in) :: id
      character(len=*), intent(in) :: lines(:)
      if (warned(id)) return
      warned(id) = .true.
      call eos_log(lines)
end subroutine eos_log_once

subroutine reset_eos_warnings()
      warned = .false.
end subroutine reset_eos_warnings

logical function eos_warned(id)
      integer, intent(in) :: id
      eos_warned = warned(id)
end function eos_warned

subroutine write_message_lines(unit, lines)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: lines(:)
      integer :: k
      do k = 1, size(lines)
         if (lines(k)(1:1) == unwritten_line) cycle
         write(unit,'(a)') trim(lines(k))
      end do
end subroutine write_message_lines

end module eos_messages

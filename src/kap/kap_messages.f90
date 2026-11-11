!----------------------------------------------------------------------
! kap_messages
!----------------------------------------------------------------------
! The opacities' messages, in one place (rotation_messages and
! eos_messages are the models). kap writes no other way:
!   kap_report    an error kap returns ierr for, found at kap_init (a
!                 setting that cannot work, a table that cannot be
!                 opened or read): to the terminal and to the run log;
!   kap_terminal  a table reader's or a lookup's diagnostic that has
!                 always gone to the terminal only;
!   kap_log       a table reader's or a lookup's diagnostic that has
!                 always gone to the run log only.
! kap_terminal and kap_log keep each older message where it went, word
! for word; making every error a kap_report, in sentence case, is the
! messages' rewording (a later commit that re-pins the logs).
! The caller formats the message into an array of lines with an
! internal write (one record per element) after setting every element
! to unwritten_line; the lines the write did not reach are not printed,
! so the output has the records the format makes, as a write to a unit
! would. A list-directed record that ends in a number goes through
! kap_terminal_values or kap_log_values instead, which write it to the
! unit as it always was: the internal write would lose the blanks the
! number's field ends with. The terminal is standard output (what kap's
! write(*) and print* meant), whether or not the run has assigned io's
! units.
module kap_messages
      use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
      use run_messages, only: run_log_unit
      implicit none
      private
      public :: kap_report, kap_terminal, kap_log, kap_terminal_values, &
           kap_log_values
! long enough for the longest record kap formats (a table path of 256
! characters and its message), which an external unit would not split
      integer, parameter, public :: message_length = 512
      character(len=*), parameter, public :: unwritten_line = achar(0)

contains

subroutine kap_report(lines)
      character(len=*), intent(in) :: lines(:)
      call write_kap_lines(output_unit, lines)
      call write_kap_lines(run_log_unit, lines)
end subroutine kap_report

subroutine kap_terminal(lines)
      character(len=*), intent(in) :: lines(:)
      call write_kap_lines(output_unit, lines)
end subroutine kap_terminal

subroutine kap_log(lines)
      character(len=*), intent(in) :: lines(:)
      call write_kap_lines(run_log_unit, lines)
end subroutine kap_log

! The list-directed record text, values (or values alone): to the
! terminal, or to the run log.
subroutine kap_terminal_values(values, text)
      real(dp), intent(in) :: values(:)
      character(len=*), intent(in), optional :: text
      call write_kap_values(output_unit, values, text)
end subroutine kap_terminal_values

subroutine kap_log_values(values, text)
      real(dp), intent(in) :: values(:)
      character(len=*), intent(in), optional :: text
      call write_kap_values(run_log_unit, values, text)
end subroutine kap_log_values

subroutine write_kap_values(unit, values, text)
      integer, intent(in) :: unit
      real(dp), intent(in) :: values(:)
      character(len=*), intent(in), optional :: text
      if (present(text)) then
         write(unit,*) text, values
      else
         write(unit,*) values
      end if
end subroutine write_kap_values

subroutine write_kap_lines(unit, lines)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: lines(:)
      integer :: k
      do k = 1, size(lines)
         if (lines(k)(1:1) == unwritten_line) cycle
         write(unit,'(a)') trim(lines(k))
      end do
end subroutine write_kap_lines

end module kap_messages

!----------------------------------------------------------------------
! kap_messages
!----------------------------------------------------------------------
! The opacities' messages, in one place (rotation_messages and
! eos_messages are the models). kap writes only errors, and only
! through kap_report, which writes to the terminal and to the run log.
! The run log's lines go through run_messages' log_line, which knows
! whether the log is open: before it opens (an embedding host, a test
! program) a line waits for it, and nothing is written to an unconnected
! unit (bug_sweep M112).
! The table readers report their own errors at kap_init; the lookups and
! the tables' builds write nothing and return a one-line reason with
! ierr, which kap_lib reports with the entry point that met it
! (kap_init, kap_get or kap_update_surface_tables).
! The caller formats the message into an array of lines with an
! internal write (one record per element) after setting every element
! to unwritten_line; the lines the write did not reach are not printed,
! so the output has the records the format makes, as a write to a unit
! would. The terminal is standard output, whether or not the run has
! assigned io's units.
module kap_messages
      use, intrinsic :: iso_fortran_env, only: output_unit
      use run_messages, only: log_line
      implicit none
      private
      public :: kap_report
! long enough for the longest record kap formats (a table path of 256
! characters and its message), which an external unit would not split
      integer, parameter, public :: message_length = 512
      character(len=*), parameter, public :: unwritten_line = achar(0)

contains

! The written lines of a message, on the terminal (standard output) and
! in the run log.
subroutine kap_report(lines)
      character(len=*), intent(in) :: lines(:)
      integer :: k
      do k = 1, size(lines)
         if (lines(k)(1:1) == unwritten_line) cycle
         write(output_unit,'(a)') trim(lines(k))
      end do
      do k = 1, size(lines)
         if (lines(k)(1:1) == unwritten_line) cycle
         call log_line(trim(lines(k)))
      end do
end subroutine kap_report

end module kap_messages

!----------------------------------------------------------------------
! kap_messages
!----------------------------------------------------------------------
! The opacities' messages, in one place (rotation_messages and
! eos_messages are the models). kap writes only errors, and only
! through kap_report, which writes to the terminal and to the run log.
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
      use run_messages, only: run_log_unit
      implicit none
      private
      public :: kap_report
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

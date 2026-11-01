!----------------------------------------------------------------------
! net_messages
!----------------------------------------------------------------------
! net/'s messages, in one place. Only the configuration prints (the
! checks of set_net_settings and the legacy cross-section map); the
! kernels print nothing and return ierr.
!   net_report  an error the run stops on, or a warning about the
!               configuration: to the terminal, and to the run log once
!               it is open.
! The caller formats the message into an array of lines (one record per
! element) with an internal write.
module net_messages
      use, intrinsic :: iso_fortran_env, only: output_unit
      use run_messages, only: run_log_unit
      implicit none
      private
      public :: net_report
      integer, parameter, public :: message_length = 160

contains

subroutine net_report(lines)
      character(len=*), intent(in) :: lines(:)
      logical :: log_open
      integer :: k

      do k = 1, size(lines)
         write(output_unit,'(a)') trim(lines(k))
      end do
      inquire(unit=run_log_unit, opened=log_open)
      if (.not. log_open .or. run_log_unit == output_unit) return
      do k = 1, size(lines)
         write(run_log_unit,'(a)') trim(lines(k))
      end do
end subroutine net_report

end module net_messages

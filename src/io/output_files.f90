!----------------------------------------------------------------------
! output_files
!----------------------------------------------------------------------
! The files io writes, opened and closed in one way, with a failure
! reported through io_units' check_io ("ERROR: cannot write <path>:
! <reason>", ierr = 1, which stops the run with exit 1):
!   open_output_file   a new file at path (any old one replaced), for
!                      writing
!   close_output_file  after the writes: the first failed write
!                      reported, the file closed, and its size on disk
!                      checked
! A writer writes with iostat=ios and iomsg=msg and skips every write
! after the first that fails (`if (ios == 0) write ...`); the file is
! closed and the failure reported once, by close_output_file. A file
! kept open across models (the history, the .mod, the isochrone file)
! has each write checked the same way, by its writer.
!
! The size check is what catches a full disk under gfortran: its runtime
! drops the bytes a full disk refuses without reporting it on the
! write, the flush or the close (measured 2026-10-31 on a 2 MB disk
! image, and with RLIMIT_FSIZE: every iostat 0, the file cut short). The
! runtime still counts the bytes it was given, which inquire by unit
! returns; after the close, inquire by file gives the size on disk, and
! a file shorter on disk than that was cut.
module output_files
      use, intrinsic :: iso_fortran_env, only: int64
      use io_units, only: check_io, report_error
      implicit none
      private
      public :: open_output_file, close_output_file

contains

! ---------------------------------------------------------------
! A new formatted file at path on unit u (newunit), for writing.
subroutine open_output_file(path, u, ierr)
      character(len=*), intent(in) :: path
      integer, intent(out) :: u
      integer, intent(inout) :: ierr
      character(len=256) :: msg
      integer :: ios

      open(newunit=u, file=path, status='REPLACE', action='WRITE', form='FORMATTED', &
           iostat=ios, iomsg=msg)
      call check_io(ios, msg, 'write', path, ierr)
end subroutine open_output_file

! ---------------------------------------------------------------
! The file at path on unit u, written: ios and msg are its writes'
! status (the first failed write's, or 0). The unit is closed in every
! case; a failed write, a failed close, or a file shorter on disk than
! what was written to it (the module header) is reported.
subroutine close_output_file(u, path, ios, msg, ierr)
      integer, intent(in) :: u, ios
      character(len=*), intent(in) :: path, msg
      integer, intent(inout) :: ierr
      character(len=256) :: flush_msg, close_msg
      integer :: flush_ios, close_ios
      integer(int64) :: written, on_disk

      flush(u, iostat=flush_ios, iomsg=flush_msg)
      inquire(unit=u, size=written)
      close(u, iostat=close_ios, iomsg=close_msg)
! the first failure: a write, the flush or the close; else the size
      if (ios /= 0) then
         call check_io(ios, msg, 'write', path, ierr)
      else if (flush_ios /= 0) then
         call check_io(flush_ios, flush_msg, 'write', path, ierr)
      else if (close_ios /= 0) then
         call check_io(close_ios, close_msg, 'write', path, ierr)
      else
         inquire(file=path, size=on_disk)
         if (on_disk < written) then
            call report_error('cannot write ' // trim(path) // ': the file on disk is ' // &
                 'shorter than what was written to it (is the disk full?)')
            ierr = 1
         end if
      end if
end subroutine close_output_file

end module output_files

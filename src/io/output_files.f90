!----------------------------------------------------------------------
! output_files
!----------------------------------------------------------------------
! The files io writes, opened, flushed and closed in one way, with a
! failure reported through run_messages' check_io ("ERROR: cannot write
! <path>: <reason>", ierr = 1, which stops the run with exit 1):
!   open_output_file   a new file at path (any old one replaced), for
!                      writing
!   flush_output_file  the writes so far handed to the system: the first
!                      failed write reported, and the file's size on
!                      disk checked
!   close_output_file  after the writes: flush_output_file, then the
!                      file closed
! A writer writes with iostat=ios and iomsg=msg and skips every write
! after the first that fails (`if (ios == 0) write ...`); the failure is
! reported once, by flush_output_file. A file written in one go ends with
! close_output_file. A file kept open across models has each write
! checked the same way, by its writer; the history file is flushed with
! flush_output_file after each row (history_output).
!
! The size check is what catches a full disk under gfortran: its runtime
! drops the bytes a full disk refuses without reporting it on the
! write, the flush or the close (measured 2026-10-31 on a 2 MB disk
! image, and with RLIMIT_FSIZE: every iostat 0, the file cut short). The
! runtime still counts the bytes it was given, which inquire by unit
! returns: its own count, not the file's size (after another process
! cut the file to 190 bytes, it still says 130000). After the flush, the
! file the unit has open is measured through its descriptor, and a file
! shorter than what the unit wrote was cut.
!
! Another run may write the same path (several decks in one directory
! all write output/history.data; M174). So the file measured is the
! unit's own, not the one the path names now, which the other run may
! have put there. And it is measured right after the flush that hands
! the bytes over, not later: the other run's open with status='replace'
! truncates the same file in place (gfortran opens an existing file with
! O_TRUNC, keeping its inode), so a check long after this run's last
! write would see that run's truncation as this run's write failing.
! That is why the history, open for the run, is checked at each row and
! not again at its close.
module output_files
      use, intrinsic :: iso_fortran_env, only: int64
      use, intrinsic :: iso_c_binding, only: c_int, c_int64_t
      use run_messages, only: check_io, report_error
      implicit none
      private
      public :: open_output_file, flush_output_file, close_output_file

! lseek's whence: from the start, from the current offset, from the end
! (the same values on macOS and Linux)
      integer(c_int), parameter :: seek_set = 0, seek_cur = 1, seek_end = 2

      interface
! off_t lseek(int fd, off_t offset, int whence); off_t is 64 bits on
! macOS and on 64-bit Linux
         integer(c_int64_t) function c_lseek(fd, offset, whence) bind(c, name='lseek')
            import :: c_int, c_int64_t
            integer(c_int), value :: fd, whence
            integer(c_int64_t), value :: offset
         end function c_lseek
      end interface

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
! The writes so far on unit u (the file at path) handed to the system:
! ios and msg are their status (the first failed write's, or 0). A
! failed write, a failed flush, or a file shorter on disk than what was
! written to it (the module header) is reported.
subroutine flush_output_file(u, path, ios, msg, ierr)
      integer, intent(in) :: u, ios
      character(len=*), intent(in) :: path, msg
      integer, intent(inout) :: ierr
      character(len=256) :: flush_msg
      integer :: flush_ios
      integer(int64) :: written, on_disk

      if (ios /= 0) then
         call check_io(ios, msg, 'write', path, ierr)
         return
      end if
      flush(u, iostat=flush_ios, iomsg=flush_msg)
      if (flush_ios /= 0) then
         call check_io(flush_ios, flush_msg, 'write', path, ierr)
         return
      end if
      inquire(unit=u, size=written)
      on_disk = size_on_disk(u)
      if (on_disk >= 0 .and. on_disk < written) then
         call report_error('cannot write ' // trim(path) // ': the file on disk is ' // &
              'shorter than what was written to it (is the disk full?)')
         ierr = 1
      end if
end subroutine flush_output_file

! ---------------------------------------------------------------
! The file at path on unit u, written: ios and msg are its writes'
! status (the first failed write's, or 0). The writes are flushed and
! checked (flush_output_file) and the unit closed in every case; a
! failed close is reported when nothing failed before it.
subroutine close_output_file(u, path, ios, msg, ierr)
      integer, intent(in) :: u, ios
      character(len=*), intent(in) :: path, msg
      integer, intent(inout) :: ierr
      character(len=256) :: close_msg
      integer :: close_ios, flush_ierr

      flush_ierr = 0
      call flush_output_file(u, path, ios, msg, flush_ierr)
      close(u, iostat=close_ios, iomsg=close_msg)
      if (flush_ierr /= 0) then
         ierr = flush_ierr
      else
         call check_io(close_ios, close_msg, 'write', path, ierr)
      end if
end subroutine close_output_file

! ---------------------------------------------------------------
! The size in bytes of the file open on unit u, from its descriptor
! (gfortran's fnum; lseek to the end, then back to the offset the
! runtime left, which it relies on), or -1 where it cannot be measured:
! a file that cannot seek, or a compiler other than gfortran, whose
! runtime the check is for. (gfortran's fstat would give the size too,
! but in a 32-bit integer: wrong past 2 GiB.)
function size_on_disk(u) result(bytes)
      integer, intent(in) :: u
      integer(int64) :: bytes
#ifdef __GFORTRAN__
      integer(c_int) :: fd
      integer(c_int64_t) :: offset

      bytes = -1
      fd = int(fnum(u), c_int)
      offset = c_lseek(fd, 0_c_int64_t, seek_cur)
      if (offset < 0) return
      bytes = c_lseek(fd, 0_c_int64_t, seek_end)
      if (c_lseek(fd, offset, seek_set) /= offset) bytes = -1
#else
      bytes = -1
#endif
end function size_on_disk

end module output_files

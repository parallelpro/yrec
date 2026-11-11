!----------------------------------------------------------------------
! opacity_table_files
!----------------------------------------------------------------------
! Opening and reading kap's table files, for every table reader. A table
! that cannot be opened, or a record that cannot be read (a short,
! truncated or malformed file), is reported with the table's name and
! path (kap_report: terminal and run log) and returns ierr = 1, so that
! kap_init returns the error and the run stops with a message instead of
! a Fortran runtime error (L5's rule, which io's files follow).
!     open_table    open a table for reading: status 'old', so a missing
!                   path is an error and never created
!     read_failed   after a read with iostat and iomsg: .true. (with the
!                   report made, the unit closed and ierr = 1) when the
!                   read failed
module opacity_table_files
      use kap_messages, only: kap_report, message_length
      implicit none
      private
      public :: open_table, read_failed

contains

subroutine open_table(path, label, unit, ierr)
      character(len=*), intent(in) :: path, label
      integer, intent(out) :: unit, ierr
      integer :: ios
      character(len=256) :: iomsg
      character(len=message_length) :: line

      ierr = 0
      open(newunit=unit, file=path, status='old', action='read', iostat=ios, iomsg=iomsg)
      if (ios /= 0) then
         line = 'kap_init: cannot open the ' // label // ' table ' // trim(path) // ': ' // trim(iomsg)
         call kap_report([line])
         ierr = 1
      end if
end subroutine open_table

logical function read_failed(ios, iomsg, label, path, unit, ierr)
      integer, intent(in) :: ios, unit
      character(len=*), intent(in) :: iomsg, label, path
      integer, intent(inout) :: ierr
      character(len=message_length) :: line

      read_failed = ios /= 0
      if (.not. read_failed) return
      close(unit)
      if (ios < 0) then
         line = 'kap_init: the ' // label // ' table ' // trim(path) // ' ends early'
      else
         line = 'kap_init: cannot read the ' // label // ' table ' // trim(path) // ': ' // trim(iomsg)
      end if
      call kap_report([line])
      ierr = 1
end function read_failed

end module opacity_table_files

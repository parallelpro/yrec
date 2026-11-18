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
!     table_z_differs   after reading the Z a table file states: .true.
!                   (with the report made, the unit closed and ierr = 1)
!                   when it is not the Z of the control the table is
!                   for
module opacity_table_files
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use kap_messages, only: kap_report, message_length, unwritten_line
      implicit none
      private
      public :: open_table, read_failed, table_z_differs

! How far the Z a table file states may be from its control's: the files
! state Z to four or five decimals, rounded (Z01895.OPALN states 0.0190,
! Z0190AG91.DBGLAOL 0.01903), so a control may differ from its table's
! header in the fourth decimal; a larger difference is another table.
      real(dp), parameter :: table_z_tolerance = 1.0d-4

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
      character(len=message_length) :: lines(1)

      read_failed = ios /= 0
      if (.not. read_failed) return
      if (ios < 0) then
         lines(1) = 'kap_init: the ' // label // ' table ' // trim(path) // ' ends early'
      else
         lines(1) = 'kap_init: cannot read the ' // label // ' table ' // trim(path) // ': ' // trim(iomsg)
      end if
      call refuse_table(unit, lines, ierr)
end function read_failed

! The Z a table file states (file_z) against the Z of the control the
! table is for (control_z, of the control named control_name), which
! blends and helium cores take as the table's Z: beyond
! table_z_tolerance, an error naming both.
logical function table_z_differs(file_z, control_z, control_name, label, path, unit, ierr)
      real(dp), intent(in) :: file_z, control_z
      character(len=*), intent(in) :: control_name, label, path
      integer, intent(in) :: unit
      integer, intent(inout) :: ierr
      character(len=message_length) :: lines(1)

      table_z_differs = abs(file_z - control_z) > table_z_tolerance
      if (.not. table_z_differs) return
      lines = unwritten_line
      write(lines,'(a,f8.5,a,f8.5)') 'kap_init: the ' // label // ' table ' // trim(path) // &
           ' is at Z =', file_z, ', not at ' // control_name // ' =', control_z
      call kap_report(lines)
      ierr = 1
      call refuse_table(unit, lines, ierr)
end function table_z_differs

! A table refused while it is being read: the unit closed, the reason
! reported, ierr = 1.
subroutine refuse_table(unit, lines, ierr)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: lines(:)
      integer, intent(inout) :: ierr

      close(unit)
end subroutine refuse_table

end module opacity_table_files

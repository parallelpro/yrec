!----------------------------------------------------------------------
! kurucz_reader
!----------------------------------------------------------------------
! Reading Kurucz 1990 table files into kurucz_table_sets (the family's
! only I/O; kurucz_opacity splines and looks them up).
!     read_kurucz_tables   the table at kurucz_table_z1 and, with
!                          two_z_tables, the one at kurucz_table_z2
module kurucz_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: kurucz_max_num_densities, kurucz_max_num_temps, &
           kurucz_table_set
      use kap_messages, only: kap_report, message_length, unwritten_line
      use opacity_table_files, only: open_table, read_failed
      use math_lib
      implicit none
      private
      public :: read_kurucz_tables

contains

!----------------------------------------------------------------------
! read_kurucz_tables
!----------------------------------------------------------------------
! Reads the Kurucz table (and, with two_z_tables, the second Z's) into
! tables(1) (and (2)) with the internal read_kurucz_table. A table that
! cannot be opened or read, or that does not fit kurucz_table_set,
! returns ierr = 1.
subroutine read_kurucz_tables(tables, kurucz_table_path, kurucz_table2_path, &
     two_z_tables, ierr)

! (1) at kurucz_table_z1 and, with two_z_tables, (2) at kurucz_table_z2
      type(kurucz_table_set), intent(inout) :: tables(:)
      character(len=256), intent(in) :: kurucz_table_path, kurucz_table2_path
      logical, intent(in) :: two_z_tables
      integer, intent(out) :: ierr
      integer :: num_read, density_index
      integer :: read_status
      real(dp) :: prev_grid_temp, grid_temp, pressure, &
           log10_opacity0, log10_opacity1, log10_opacity2, &
           log10_opacity4, log10_opacity8, electron_density, &
           atom_density, density, unused_col

      ierr = 0
! The range check of compute_kurucz_opacity (log rho <= -3, log T <= 4.1)
! is the first table's only, as in the F77 code.
      tables(2)%check_range = .false.

      call read_kurucz_table(kurucz_table_path, &
           tables(1), ierr)
      if (ierr /= 0) return

      if (two_z_tables) then
         call read_kurucz_table(kurucz_table2_path, &
              tables(2), ierr)
         if (ierr /= 0) return
      end if

      return

contains

! read_kurucz_table: read one Kurucz table file into tbl: after its two
! header lines (three in some files), one row per (T, rho), each T's
! rows together; a row's log10 kappa (at no microturbulence) is stored
! as kappa.
subroutine read_kurucz_table(path, tbl, ierr)
      character(len=256), intent(in) :: path
      type(kurucz_table_set), intent(inout) :: tbl
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success
      integer :: unit
      character(len=message_length) :: lines(1)
      integer :: ios
      character(len=256) :: iomsg

    1 format(2f5.2,5f7.3,3f9.5,f8.3)
      call open_table(path, 'Kurucz 1990', unit, ierr)
      if (ierr /= 0) return
      read(unit,'(/)',iostat=ios,iomsg=iomsg)
      if (read_failed(ios, iomsg, 'Kurucz 1990', path, unit, ierr)) return

      num_read = 0
      prev_grid_temp = 0.0d0
      table_read: do
      read(unit,1,iostat=read_status) grid_temp, pressure, &
           log10_opacity0, log10_opacity1, log10_opacity2, &
           log10_opacity4, log10_opacity8, electron_density, &
           atom_density, density, unused_col
      if (read_status < 0) exit table_read    ! was end=120
! a record that does not parse is a header line before the first row
! (some tables have three), and an error after it
      if (read_status > 0 .and. num_read == 0) cycle table_read
      if (read_status > 0) then
         lines = unwritten_line
         write(lines,'(a,i0)') 'kap_init: the Kurucz 1990 table ' // trim(path) // &
              ' has an unreadable row after temperature row ', num_read
         call kap_report(lines)
         close(unit)
         ierr = 1
         return
      end if
      if (prev_grid_temp /= grid_temp) then
         num_read = num_read+1
         if (num_read > kurucz_max_num_temps) then
            lines = unwritten_line
            write(lines,'(a,i0,a)') 'kap_init: the Kurucz 1990 table ' // trim(path) // &
                 ' has more than ', kurucz_max_num_temps, ' temperature rows'
            call kap_report(lines)
            close(unit)
            ierr = 1
            return
         end if
         tbl%grid_logt(num_read) = grid_temp
         prev_grid_temp = grid_temp
         density_index = 1
      end if
      if (num_read < 1) then
         lines = unwritten_line
         write(lines,'(a)') 'kap_init: the Kurucz 1990 table ' // trim(path) // &
              ' has no temperature row before its data'
         call kap_report(lines)
         close(unit)
         ierr = 1
         return
      end if
      if (density_index > kurucz_max_num_densities) then
         lines = unwritten_line
         write(lines,'(a,i0,a)') 'kap_init: the Kurucz 1990 table ' // trim(path) // &
              ' has more than ', kurucz_max_num_densities, ' densities at one temperature'
         call kap_report(lines)
         close(unit)
         ierr = 1
         return
      end if
      tbl%log10_rho(num_read, density_index) = density
      tbl%opacity(num_read, density_index) = exp10(log10_opacity0)
      density_index = density_index+1
      end do table_read

      tbl%num_temps = num_read
!     Close the table we have read
      close(unit,err=99)
      return
   99 continue
      lines = unwritten_line
      write(lines,'(a)') 'kap_init: cannot close the Kurucz 1990 table ' // trim(path)
      call kap_report(lines)
      ierr = 1
end subroutine read_kurucz_table

end subroutine read_kurucz_tables

end module kurucz_reader

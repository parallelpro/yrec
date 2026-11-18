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
! tables(1) (and (2)) with read_kurucz_table. A table that cannot be
! opened or read, or that does not fit kurucz_table_set, returns ierr
! (1).
subroutine read_kurucz_tables(tables, kurucz_table_path, kurucz_table2_path, &
     two_z_tables, ierr)
! (1) at kurucz_table_z1 and, with two_z_tables, (2) at kurucz_table_z2
      type(kurucz_table_set), intent(inout) :: tables(:)
      character(len=256), intent(in) :: kurucz_table_path, kurucz_table2_path
      logical, intent(in) :: two_z_tables
      integer, intent(out) :: ierr

      ierr = 0
! The range check of compute_kurucz_opacity (log rho <= -3, log T <= 4.1)
! is the first table's only, as in the F77 code.
      tables(2)%check_range = .false.
      call read_kurucz_table(kurucz_table_path, tables(1), ierr)
      if (ierr /= 0) return
      if (two_z_tables) then
         call read_kurucz_table(kurucz_table2_path, tables(2), ierr)
      end if
end subroutine read_kurucz_tables

! One Kurucz table file, at path, into tbl: after its two header lines
! (three in some files), one row per (T, rho), each T's rows together.
! A row gives log T, log P, log10 kappa at five microturbulent
! velocities (0, 1, 2, 4 and 8 km/s), log Ne, log Na, log rho and a last
! column; log T, log rho and kappa at no microturbulence are kept.
subroutine read_kurucz_table(path, tbl, ierr)
      character(len=256), intent(in) :: path
      type(kurucz_table_set), intent(inout) :: tbl
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success
      integer :: unit, ios, num_temps, density_index
! a row as read (only log T, log10 kappa at 0 km/s and log rho are used)
      real(dp) :: logt, previous_logt, log_pressure, log10_opacity(5), log_electron_density, &
           log_atom_density, log_rho, last_column
      character(len=message_length) :: lines(1)
      character(len=256) :: iomsg

      call open_table(path, 'Kurucz 1990', unit, ierr)
      if (ierr /= 0) return
      read(unit,'(/)',iostat=ios,iomsg=iomsg)
      if (read_failed(ios, iomsg, 'Kurucz 1990', path, unit, ierr)) return

      num_temps = 0
      previous_logt = 0.0d0
      lines = unwritten_line
      do
         read(unit,'(2f5.2,5f7.3,3f9.5,f8.3)',iostat=ios) logt, log_pressure, log10_opacity, &
              log_electron_density, log_atom_density, log_rho, last_column
         if (ios < 0) exit
! a record that does not parse is a header line before the first row
! (some tables have three), and an error after it
         if (ios > 0 .and. num_temps == 0) cycle
         if (ios > 0) then
            write(lines,'(a,i0)') 'kap_init: the Kurucz 1990 table ' // trim(path) // &
                 ' has an unreadable row after temperature row ', num_temps
            exit
         end if
         if (previous_logt /= logt) then
            num_temps = num_temps+1
            if (num_temps > kurucz_max_num_temps) then
               write(lines,'(a,i0,a)') 'kap_init: the Kurucz 1990 table ' // trim(path) // &
                    ' has more than ', kurucz_max_num_temps, ' temperature rows'
               exit
            end if
            tbl%grid_logt(num_temps) = logt
            previous_logt = logt
            density_index = 1
         end if
         if (num_temps < 1) then
            write(lines,'(a)') 'kap_init: the Kurucz 1990 table ' // trim(path) // &
                 ' has no temperature row before its data'
            exit
         end if
         if (density_index > kurucz_max_num_densities) then
            write(lines,'(a,i0,a)') 'kap_init: the Kurucz 1990 table ' // trim(path) // &
                 ' has more than ', kurucz_max_num_densities, ' densities at one temperature'
            exit
         end if
         tbl%log10_rho(num_temps, density_index) = log_rho
         tbl%opacity(num_temps, density_index) = exp10(log10_opacity(1))
         density_index = density_index+1
      end do
      if (lines(1)(1:1) == unwritten_line .and. num_temps == 0) then
         write(lines,'(a)') 'kap_init: the Kurucz 1990 table ' // trim(path) // &
              ' has no table row'
      end if
      if (lines(1)(1:1) /= unwritten_line) then
         call kap_report(lines)
         close(unit)
         ierr = 1
         return
      end if
      tbl%num_temps = num_temps
      close(unit,iostat=ios)
      if (ios /= 0) then
         write(lines,'(a)') 'kap_init: cannot close the Kurucz 1990 table ' // trim(path)
         call kap_report(lines)
         ierr = 1
      end if
end subroutine read_kurucz_table

end module kurucz_reader

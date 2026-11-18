!----------------------------------------------------------------------
! laol_reader
!----------------------------------------------------------------------
! Reading the LAOL89 table files into laol_table_sets and the pure-Z
! table into a pure_z_table_set (the family's only I/O; laol_opacity
! splines and looks them up).
!     read_laol_tables    the table at laol_table_z1 and, with
!                         two_z_tables, the one at laol_table_z2
!     read_pure_z_table   the pure-Z table (Z = 1, C and O)
module laol_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: laol_table_set, n_laol_rho, n_laol_t, n_laol_x, &
           pure_z_table_set
      use kap_messages, only: kap_report, message_length, unwritten_line
      use opacity_table_files, only: open_table, read_failed, table_z_differs
      use math_lib
      implicit none
      private
      public :: read_laol_tables, read_pure_z_table

contains

!----------------------------------------------------------------------
! read_laol_tables
!----------------------------------------------------------------------
! Reads the LAOL89 table (and, with two_z_tables, the second Z's) into
! tables(1) (and (2)) with read_laol_table, and records the
! extrapolation distance. A table that cannot be opened or read, does
! not fit laol_table_set, or states another Z than its control's
! (table_z) returns ierr = 1.
subroutine read_laol_tables(tables, laol_table_path, laol_table2_path, two_z_tables, &
     table_z, extrapolation_tolerance, ierr)
! (1) at laol_table_z1 and, with two_z_tables, (2) at laol_table_z2
      type(laol_table_set), intent(inout) :: tables(:)
      character(len=256), intent(in) :: laol_table_path, laol_table2_path
      logical, intent(in) :: two_z_tables
! laol_table_z1 and laol_table_z2
      real(dp), intent(in) :: table_z(:)
! laol_extrapolation_tolerance (its ln is the distance in log rho)
      real(dp), intent(in) :: extrapolation_tolerance
      integer, intent(out) :: ierr

      ierr = 0
      tables(2)%table_number = 2
      tables(1)%log_extrapolation_tolerance = log(extrapolation_tolerance)
      tables(2)%log_extrapolation_tolerance = log(extrapolation_tolerance)
      call read_laol_table(laol_table_path, tables(1), table_z(1), 'laol_table_z1', ierr)
      if (ierr /= 0) return
      if (two_z_tables) then
         call read_laol_table(laol_table2_path, tables(2), table_z(2), 'laol_table_z2', ierr)
      end if
end subroutine read_laol_tables

! One LAOL89 table file, at path, into tbl: its sizes, its grids in X,
! density and T, and the opacities. The header's two metal mixtures are
! read and not kept: the Cox & Stewart mixture of its zlot part (C, N,
! O, Ne, Na, Mg, Al, Si, Ar, Fe by weight; the star takes its surface
! mixture from the controls) and the 18 metals of its zhit part, for a
! Debye-Hueckel correction (eos takes those from its own envelope
! mixture, M37). The Z the header states for the two parts (zlot, zhit)
! must be the Z of the table's control (control_z, of control_name).
subroutine read_laol_table(path, tbl, control_z, control_name, ierr)
      character(len=256), intent(in) :: path
      type(laol_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: control_z
      character(len=*), intent(in) :: control_name
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success
      integer :: laol_unit, ix, ir, ios
      real(dp) :: zlot, zhit, cox_stewart_mixture(10), debye_huckel_mixture(18)
      character(len=message_length) :: lines(1)
      character(len=256) :: iomsg

      call open_table(path, 'LAOL 1989', laol_unit, ierr)
      if (ierr /= 0) return
!     Read in array sizes
      read(laol_unit,'(/,18x,i2,9x,i3,14x,i3)',iostat=ios,iomsg=iomsg) tbl%num_x,tbl%num_rho,tbl%num_t
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
      if (tbl%num_x > n_laol_x-1.or.tbl%num_rho > n_laol_rho.or.tbl%num_t > n_laol_t) then
         lines = unwritten_line
         write(lines,'(a,3(i0,a),3(i0,a))') 'kap_init: the LAOL 1989 table ' // trim(path) // &
              ' has ', tbl%num_x, ' X, ', tbl%num_rho, ' rho and ', tbl%num_t, &
              ' T values, more than ', n_laol_x-1, ', ', n_laol_rho, ' and ', n_laol_t
         call kap_report(lines)
         close(laol_unit)
         ierr = 1
         return
      end if
!     zlot, the metal abundance of the table's low-temperature part, and
!     its Cox & Stewart mixture, scaled to sum to zlot
      read(laol_unit,'(47x,f8.5,//,1p6e12.5,/,1p4e12.5)',iostat=ios,iomsg=iomsg) zlot, &
           cox_stewart_mixture
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
      if (table_z_differs(zlot, control_z, control_name, 'LAOL 1989', path, laol_unit, ierr)) return
!     zhit, the metal abundance of the table's high-temperature part, and
!     its 18 metals, scaled to sum to zhit
      read(laol_unit,'(54x,f8.5,/////,1p6e12.5,/,1p6e12.5,/1p6e12.5)',iostat=ios,iomsg=iomsg) zhit, &
           debye_huckel_mixture
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
      if (table_z_differs(zhit, control_z, control_name, 'LAOL 1989', path, laol_unit, ierr)) return
!     Read in h mass fractions of table
      read(laol_unit,'(/,(1p6e12.5))',iostat=ios,iomsg=iomsg) tbl%grid_x(1:tbl%num_x)
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
!     Read in density grid of table
      read(laol_unit,'(/,(1p6e12.5))',iostat=ios,iomsg=iomsg) tbl%grid_rho(1:tbl%num_rho)
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
!     Read in temperature grid of table
      read(laol_unit,'(/,(1p6e12.5))',iostat=ios,iomsg=iomsg) tbl%grid_t(1:tbl%num_t)
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
!     Read in opacities: after a line, for each X and density a line and
!     then the opacities at every T
      read(laol_unit,'(1x)',iostat=ios,iomsg=iomsg)
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
      do ix=1,tbl%num_x
         do ir=1,tbl%num_rho
            read(laol_unit,'(1x)',iostat=ios,iomsg=iomsg)
            if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
            read(laol_unit,'(1p6e12.5)',iostat=ios,iomsg=iomsg) tbl%opacity(ix,ir,1:tbl%num_t)
            if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
         end do
      end do
      close(laol_unit)
end subroutine read_laol_table

!----------------------------------------------------------------------
! read_pure_z_table
!----------------------------------------------------------------------
! Reads the pure-Z table (Z = 1, one X; the LAOL89 file layout) into tbl
! and records the extrapolation distance. Its header's Z, mixtures and X
! are read and not kept. A table that cannot be opened or read, or does
! not fit pure_z_table_set, returns ierr = 1.
subroutine read_pure_z_table(tbl, pure_z_table_path, extrapolation_tolerance, ierr)
      type(pure_z_table_set), intent(inout) :: tbl
      character(len=256), intent(in) :: pure_z_table_path
! laol_extrapolation_tolerance (its ln is the distance in log rho)
      real(dp), intent(in) :: extrapolation_tolerance
      integer, intent(out) :: ierr
      integer :: pure_z_unit, num_x, ir, ios
      real(dp) :: header_z, cox_stewart_mixture(10), header_x
      character(len=message_length) :: lines(1)
      character(len=256) :: iomsg

      tbl%log_extrapolation_tolerance = log(extrapolation_tolerance)
      call open_table(pure_z_table_path, 'pure-Z', pure_z_unit, ierr)
      if (ierr /= 0) return
!     Read in array sizes
      read(pure_z_unit,'(/,18x,i2,9x,i3,14x,i3)',iostat=ios,iomsg=iomsg) num_x,tbl%num_rho,tbl%num_t
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
      if (num_x /= 1.or.tbl%num_rho > n_laol_rho.or.tbl%num_t > n_laol_t) then
         lines = unwritten_line
         write(lines,'(a,3(i0,a),2(i0,a))') 'kap_init: the pure-Z table ' // &
              trim(pure_z_table_path) // ' has ', num_x, ' X, ', tbl%num_rho, ' rho and ', &
              tbl%num_t, ' T values, where it may have 1, up to ', n_laol_rho, ' and ', n_laol_t
         call kap_report(lines)
         close(pure_z_unit)
         ierr = 1
         return
      end if
!     the header's Z and mixture of its low-temperature part, the Z of
!     its high-temperature part (with its mixture's lines skipped)
      read(pure_z_unit,'(47x,f8.5,//,1p6e12.5,/,1p4e12.5)',iostat=ios,iomsg=iomsg) header_z, &
           cox_stewart_mixture
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
      read(pure_z_unit,'(54x,f8.5,///////)',iostat=ios,iomsg=iomsg) header_z
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
!     Read in h mass fractions of table
      read(pure_z_unit,'(/,(1p6e12.5))',iostat=ios,iomsg=iomsg) header_x
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
!     Read in density grid of table
      read(pure_z_unit,'(/,(1p6e12.5))',iostat=ios,iomsg=iomsg) tbl%rho_grid(1:tbl%num_rho)
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
!     Read in temperature grid of table
      read(pure_z_unit,'(/,(1p6e12.5))',iostat=ios,iomsg=iomsg) tbl%logt_grid(1:tbl%num_t)
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
!     Read in pure Z opacities: after a line, for each density a line
!     and then the opacities at every T
      read(pure_z_unit,'(1x)',iostat=ios,iomsg=iomsg)
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
      do ir=1,tbl%num_rho
         read(pure_z_unit,'(1x)',iostat=ios,iomsg=iomsg)
         if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
         read(pure_z_unit,'(1p6e12.5)',iostat=ios,iomsg=iomsg) tbl%opacity(ir,1:tbl%num_t)
         if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
      end do
      close(pure_z_unit)
end subroutine read_pure_z_table

end module laol_reader

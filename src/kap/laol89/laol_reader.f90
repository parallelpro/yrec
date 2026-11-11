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
      use opacity_table_files, only: open_table, read_failed
      use math_lib
      implicit none
      private
      public :: read_laol_tables, read_pure_z_table

contains

!----------------------------------------------------------------------
! read_laol_tables (formerly rdlaol)
!----------------------------------------------------------------------
! Reads the LAOL89 table (and, with two_z_tables, the second Z's) into
! tables(1) (and (2)) with the internal read_laol_table, and records the
! extrapolation distance. A table that cannot be opened or read, or does
! not fit laol_table_set, returns ierr = 1.
subroutine read_laol_tables(tables, laol_work_array, laol_table_path, &
     laol_table2_path, two_z_tables, extrapolation_tolerance, ierr)
      use composition_layout, only: ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, &
           ix_o, ix_n, ix_ar, ix_ne

! (1) at laol_table_z1 and, with two_z_tables, (2) at laol_table_z2
      type(laol_table_set), intent(inout) :: tables(:)
      logical, intent(in) :: two_z_tables
! laol_extrapolation_tolerance (its ln is the distance in log rho)
      real(dp), intent(in) :: extrapolation_tolerance
      integer, intent(out) :: ierr
      real(dp), intent(inout) :: laol_work_array(:)
      character(len=256), intent(in) :: laol_table_path, laol_table2_path

! work_array2: the second table's Cox & Stewart mixture, read and not
! used
      real(dp) :: work_array2(12)
      integer :: ii, ix, ir, it

      ierr = 0
      tables(2)%table_number = 2
      tables(1)%log_extrapolation_tolerance = log(extrapolation_tolerance)
      tables(2)%log_extrapolation_tolerance = log(extrapolation_tolerance)
      call read_laol_table(laol_table_path, tables(1), &
           laol_work_array, ierr)
      if (ierr /= 0) return

      if (two_z_tables) then
         call read_laol_table(laol_table2_path, tables(2), &
              work_array2, ierr)
         if (ierr /= 0) return
      end if
      return

contains

! read_laol_table: read one LAOL89 table file at `path` into tbl. The
! header's Cox & Stewart metal mixture goes to work_array (the
! caller's laol_work_array for table 1, a discarded local for table
! 2); its 18-element mixture for the EOS's Debye-Hueckel correction is
! read and not kept (eos takes those metals from its own envelope
! mixture, M37).
subroutine read_laol_table(path, tbl, work_array, ierr)
      character(len=256), intent(in) :: path
      type(laol_table_set), intent(inout) :: tbl
      real(dp), intent(inout) :: work_array(:)
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success
      integer :: laol_unit
!     Zlot/zhit and the Debye-Hueckel mixture are read-list targets
!     that are never used afterwards.
      real(dp) :: zlot, zhit, debye_huckel_z(18)
      character(len=message_length) :: lines(1)
      integer :: ios
      character(len=256) :: iomsg

      call open_table(path, 'LAOL 1989', laol_unit, ierr)
      if (ierr /= 0) return
!     Read in array sizes
      read(laol_unit,100,iostat=ios,iomsg=iomsg) tbl%num_x,tbl%num_rho,tbl%num_t
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
  100 format(/,18x,i2,9x,i3,14x,i3)
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
!     the relative abundances by weight of the metals C, N, O, Ne, Na,
!     Mg, Si, Ar, Fe (the Cox & Stewart mixture) of the table's zlot
!     part, scaled to sum to zlot
      read(laol_unit,120,iostat=ios,iomsg=iomsg) zlot, work_array(ix_c), work_array(ix_n), &
           work_array(ix_o), work_array(ix_ne), work_array(ix_na), &
           work_array(ix_mg), work_array(ix_al), work_array(ix_si), &
           work_array(ix_ar), work_array(ix_fe)
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
  120 format(47x,f8.5,//,1p6e12.5,/,1p4e12.5)
!     zhit, the metal abundance of the table's hit part, and the 18
!     metals C, N, O, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, Ca, Ti, Cr, Mn,
!     Fe, Ni scaled to sum to zhit (for a Debye-Hueckel correction)
      read(laol_unit,130,iostat=ios,iomsg=iomsg) zhit, debye_huckel_z
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
  130 format(54x,f8.5,/////,1p6e12.5,/,1p6e12.5,/1p6e12.5)
!     Read in h mass fractions of table
      read(laol_unit,140,iostat=ios,iomsg=iomsg) (tbl%grid_x(ii),ii=1,tbl%num_x)
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
  140 format(/,(1p6e12.5))
!     Read in density grid of table
      read(laol_unit,150,iostat=ios,iomsg=iomsg) (tbl%grid_rho(ii),ii=1,tbl%num_rho)
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
  150 format(/,(1p6e12.5))
!     Read in temperature grid of table
      read(laol_unit,160,iostat=ios,iomsg=iomsg) (tbl%grid_t(ii),ii=1,tbl%num_t)
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
  160 format(/,(1p6e12.5))
!     Read in opacities
      read(laol_unit,170,iostat=ios,iomsg=iomsg)
      if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
  170 format(1x)
      do ix=1,tbl%num_x
         do ir=1,tbl%num_rho
            read(laol_unit,200,iostat=ios,iomsg=iomsg)
            if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
  200       format(1x)
            read(laol_unit,210,iostat=ios,iomsg=iomsg) (tbl%opacity(ix,ir,it),it=1,tbl%num_t)
            if (read_failed(ios, iomsg, 'LAOL 1989', path, laol_unit, ierr)) return
  210       format(1p6e12.5)
         end do
      end do
      close(laol_unit)
end subroutine read_laol_table

end subroutine read_laol_tables

!----------------------------------------------------------------------
! read_pure_z_table (formerly rdzlaol)
!----------------------------------------------------------------------
! Reads the pure-Z table (Z = 1, one X) into tbl and records the
! extrapolation distance. A table that cannot be opened or read, or does
! not fit pure_z_table_set, returns ierr = 1.
subroutine read_pure_z_table(tbl, pure_z_table_path, extrapolation_tolerance, ierr)

      type(pure_z_table_set), intent(inout) :: tbl
! laol_extrapolation_tolerance (its ln is the distance in log rho)
      real(dp), intent(in) :: extrapolation_tolerance
      integer, intent(out) :: ierr
      character(len=256), intent(in) :: pure_z_table_path
      integer :: n, i, ii, ir, it
      integer :: pure_z_unit
      real(dp) :: dummy(104)
      character(len=message_length) :: lines(1)
      integer :: ios
      character(len=256) :: iomsg

      tbl%log_extrapolation_tolerance = log(extrapolation_tolerance)
      call open_table(pure_z_table_path, 'pure-Z', pure_z_unit, ierr)
      if (ierr /= 0) return
!     Read in array sizes
      read(pure_z_unit,100,iostat=ios,iomsg=iomsg) n,tbl%num_rho,tbl%num_t
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
  100 format(/,18x,i2,9x,i3,14x,i3)
      if (n /= 1.or.tbl%num_rho > n_laol_rho.or.tbl%num_t > n_laol_t) then
         lines = unwritten_line
         write(lines,'(a,3(i0,a),2(i0,a))') 'kap_init: the pure-Z table ' // &
              trim(pure_z_table_path) // ' has ', n, ' X, ', tbl%num_rho, ' rho and ', &
              tbl%num_t, ' T values, where it may have 1, up to ', n_laol_rho, ' and ', n_laol_t
         call kap_report(lines)
         close(pure_z_unit)
         ierr = 1
         return
      end if
      read(pure_z_unit,120,iostat=ios,iomsg=iomsg) (dummy(i),i=1,11)
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
  120 format(47x,f8.5,//,1p6e12.5,/,1p4e12.5)
      read(pure_z_unit,131,iostat=ios,iomsg=iomsg) dummy(1)
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
  131 format(54x,f8.5,///////)
!     Read in h mass fractions of table
      read(pure_z_unit,140,iostat=ios,iomsg=iomsg) dummy(1)
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
  140 format(/,(1p6e12.5))
!     Read in density grid of table
      read(pure_z_unit,150,iostat=ios,iomsg=iomsg) (tbl%rho_grid(ii),ii=1,tbl%num_rho)
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
  150 format(/,(1p6e12.5))
!     Read in temperature grid of table
      read(pure_z_unit,160,iostat=ios,iomsg=iomsg) (tbl%logt_grid(ii),ii=1,tbl%num_t)
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
  160 format(/,(1p6e12.5))
!     Read in pure Z opacities
      read(pure_z_unit,170,iostat=ios,iomsg=iomsg)
      if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
  170 format(1x)
      do ir=1,tbl%num_rho
         read(pure_z_unit,200,iostat=ios,iomsg=iomsg)
         if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
  200    format(1x)
         read(pure_z_unit,210,iostat=ios,iomsg=iomsg) (tbl%opacity(ir,it),it=1,tbl%num_t)
         if (read_failed(ios, iomsg, 'pure-Z', pure_z_table_path, pure_z_unit, ierr)) return
  210    format(1p6e12.5)
         end do
      close(pure_z_unit)

end subroutine read_pure_z_table

end module laol_reader

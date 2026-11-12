!----------------------------------------------------------------------
! scv_table_reader
!----------------------------------------------------------------------
! The SCV tables' reader, read_scv_tables (eos_init calls it, through
! SCV's table-source entry): the files' I/O, kept apart from the SCV EOS
! that interpolates them (scv_eos prints nothing). The tables live in
! scv_eos_tables.
module scv_table_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use scv_eos_tables, only: iscv_log10_p, nptsx, scv_nt, tablex, tabley, tablez, tlogx, &
           tablenv_hydrogen_fraction, tablenv_metal_fraction, scv_tables_loaded
      use eos_messages, only: eos_report, message_length, unwritten_line
      implicit none
      private
      public :: read_scv_tables

contains

!----------------------------------------------------------------------
! read_scv_tables
!----------------------------------------------------------------------
! Reads the SCV tables (columns 1-11 of the pure-H and pure-He tables,
! 1-13 of the pure-metal one, computed with the Prather EOS) into
! tablex, tabley and tablez. The three share one grid, the H table's: T
! rows tlogx, point counts nptsx and each point's log P_gas (column 1).
! The envelope table, built from the tables read before, no longer
! stands (its composition is reset; the next eos_set_mixture builds it
! again). A table that cannot be opened, a read that fails (a truncated
! or corrupt table, or a row longer than the arrays), and a He or Z
! table on another grid are ierr = 1 with a message; every table opened
! is closed.
subroutine read_scv_tables(scv_h_table_path, scv_he_table_path, scv_z_table_path, &
     ierr)
      character(len=*), intent(in) :: scv_h_table_path, scv_he_table_path, &
           scv_z_table_path
      integer, intent(out) :: ierr
! the units of the H, He and Z tables (0: not open)
      integer :: units(3), k

      tablenv_hydrogen_fraction = -1.0d0
      tablenv_metal_fraction = -1.0d0
      scv_tables_loaded = .false.
      units = 0
      call open_scv_table(scv_h_table_path, units(1), ierr)
      if (ierr == 0) call open_scv_table(scv_he_table_path, units(2), ierr)
      if (ierr == 0) call open_scv_table(scv_z_table_path, units(3), ierr)
      if (ierr == 0) call read_scv_rows(units, scv_h_table_path, scv_he_table_path, &
           scv_z_table_path, ierr)
      scv_tables_loaded = ierr == 0
      do k = 1, 3
         if (units(k) /= 0) close(units(k))
      end do
end subroutine read_scv_tables

! The rows of the three tables, on units(1:3) (H, He, Z; their paths
! for the messages): for each T row its log T and point count, then its
! points; the metal table lists them from the highest pressure down. The
! He and Z tables' log T, point counts and log P_gas must be the H
! table's (off_grid names the first that is not).
subroutine read_scv_rows(units, h_path, he_path, z_path, ierr)
      integer, intent(in) :: units(:)
      character(len=*), intent(in) :: h_path, he_path, z_path
      integer, intent(out) :: ierr
      character(len=message_length) :: eos_lines(4)
      integer :: t_idx, p_idx, col_idx, ios, off_grid
! the table being read (1 H, 2 He, 3 Z): the one named when a read fails
      integer :: reading
! the He and Z tables' log T and point count of the row
      real(dp) :: copy_log10_t(2:3)
      integer :: copy_npts(2:3)

      ierr = 0
      off_grid = 0
      reading = 1
      do t_idx = 1, scv_nt
         reading = 1
         read(units(1),1,iostat=ios) tlogx(t_idx),nptsx(t_idx)
         if (ios == 0 .and. (nptsx(t_idx) < 4 .or. nptsx(t_idx) > size(tablex, 2))) ios = -1
         if (ios == 0) reading = 2
         if (ios == 0) read(units(2),1,iostat=ios) copy_log10_t(2),copy_npts(2)
         if (ios == 0) reading = 3
         if (ios == 0) read(units(3),1,iostat=ios) copy_log10_t(3),copy_npts(3)
    1    format(f5.2,i4)
         if (ios /= 0) exit
         off_grid = first_off_grid(copy_log10_t /= tlogx(t_idx) .or. copy_npts /= nptsx(t_idx))
         if (off_grid /= 0) exit
         do p_idx = 1, nptsx(t_idx)
            reading = 1
            read(units(1),2,iostat=ios) (tablex(t_idx,p_idx,col_idx),col_idx=1,11)
            if (ios == 0) reading = 2
            if (ios == 0) read(units(2),2,iostat=ios) (tabley(t_idx,p_idx,col_idx),col_idx=1,11)
    2       format(f6.2,1p2e13.5,0p,8f9.4)
            if (ios /= 0) exit
         end do
         if (ios /= 0) exit
         reading = 3
         do p_idx = nptsx(t_idx),1,-1
            read(units(3),3,iostat=ios)(tablez(t_idx,p_idx,col_idx),col_idx=1,13)
    3       format(f6.2,12f9.4)
            if (ios /= 0) exit
         end do
         if (ios /= 0) exit
         off_grid = first_off_grid([any(tabley(t_idx,1:nptsx(t_idx),iscv_log10_p) /= &
              tablex(t_idx,1:nptsx(t_idx),iscv_log10_p)), any(tablez(t_idx,1:nptsx(t_idx),iscv_log10_p) /= &
              tablex(t_idx,1:nptsx(t_idx),iscv_log10_p))])
         if (off_grid /= 0) exit
      end do
      if (ios /= 0) then
         eos_lines = unwritten_line
         write(eos_lines,'(a,a,a,i0,a,a)') 'eos_init: the SCV ', trim(scv_table_name(reading)), &
              ' table ends early or is corrupt, at T row ', t_idx, ': ', trim(scv_table_path(reading))
         call eos_report(eos_lines)
         ierr = 1
      else if (off_grid /= 0) then
         eos_lines = unwritten_line
         if (off_grid == 2) then
            write(eos_lines,'(a,i0,a,a)') 'eos_init: an SCV table''s grid (log T, point count or ' // &
                 'log P_gas) is not the H table''s, at T row ', t_idx, ': ', trim(he_path)
         else
            write(eos_lines,'(a,i0,a,a)') 'eos_init: an SCV table''s grid (log T, point count or ' // &
                 'log P_gas) is not the H table''s, at T row ', t_idx, ': ', trim(z_path)
         end if
         call eos_report(eos_lines)
         ierr = 1
      end if

contains

! the name and the path of table k (1 H, 2 He, 3 Z)
pure function scv_table_name(k) result(name)
      integer, intent(in) :: k
      character(len=2) :: name
      name = 'H '
      if (k == 2) name = 'He'
      if (k == 3) name = 'Z '
end function scv_table_name

function scv_table_path(k) result(path)
      integer, intent(in) :: k
      character(len=max(len(h_path), len(he_path), len(z_path))) :: path
      path = h_path
      if (k == 2) path = he_path
      if (k == 3) path = z_path
end function scv_table_path

! 2 or 3 for the He or the Z table (differs(1) or (2)), the first that
! differs; 0 for neither
pure integer function first_off_grid(differs)
      logical, intent(in) :: differs(:)
      first_off_grid = 0
      if (differs(2)) first_off_grid = 3
      if (differs(1)) first_off_grid = 2
end function first_off_grid

end subroutine read_scv_rows

subroutine open_scv_table(path, table_unit, ierr)
      character(len=*), intent(in) :: path
      integer, intent(out) :: table_unit, ierr
      character(len=message_length) :: eos_lines(4)
      open(newunit=table_unit, file=path, status='OLD', iostat=ierr)
      if (ierr /= 0) then
         eos_lines = unwritten_line
         write(eos_lines,*) 'eos_init: cannot open the SCV table ', trim(path)
         call eos_report(eos_lines)
         table_unit = 0
         ierr = 1
      end if
end subroutine open_scv_table

end module scv_table_reader

!----------------------------------------------------------------------
! opal_table_reader
!----------------------------------------------------------------------
! The OPAL EOS tables' reader: the files' I/O, kept apart from the
! interpolation (opal_eos prints nothing of its own; it calls readcoeos
! on a vintage's first use).
!     open_opal_eos_table  open the selected vintage's table (eos_init,
!                          through OPAL's table-source entry)
!     readcoeos            read it into the vintage's instance, with its
!                          messages (report_*)
module opal_table_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opal_eos_tables, only: n_eos_mv, n_eos_mx, opal_layout, &
           opal_shared, opal_eos_vintage, opal_flag_set, opal_empty_entry
      use eos_messages, only: eos_report, message_length, unwritten_line
      implicit none
      private
      public :: open_opal_eos_table, readcoeos

contains

!----------------------------------------------------------------------
! open_opal_eos_table
!----------------------------------------------------------------------
! Opens the OPAL EOS table of the vintage the settings select on a new
! unit (opal_shared%table_unit), closing one left open by an earlier
! eos_init, and marks the vintage's table unread, so that its reader
! reads this file on first use (refilling the table first) whatever an
! earlier eos_init loaded.
subroutine open_opal_eos_table(vintage, opal_eos_table_path, ierr)
      integer, intent(in) :: vintage
      character(len=*), intent(in) :: opal_eos_table_path
      integer, intent(out) :: ierr
      character(len=message_length) :: eos_lines(4)
      logical :: is_open

      if (opal_shared%table_unit /= 0) then
         inquire(unit=opal_shared%table_unit, opened=is_open)
         if (is_open) close(opal_shared%table_unit)
      end if
      opal_shared%table_loaded_flag(vintage) = 0
      opal_shared%readco_init_flag(vintage) = 0
      open(newunit=opal_shared%table_unit, file=opal_eos_table_path, status='OLD', &
           iostat=ierr)
      if (ierr /= 0) then
         eos_lines = unwritten_line
         write(eos_lines,*) 'eos_init: cannot open the OPAL EOS table ', trim(opal_eos_table_path)
         call eos_report(eos_lines)
         ierr = 1
      end if
end subroutine open_opal_eos_table

!----------------------------------------------------------------------
! readcoeos (formerly readco, readcoeos01 and readcoeos06)
!----------------------------------------------------------------------
! Reads an OPAL EOS table (one block per tabulated X) into v, then
! derives the grid spacings interpolate_opal_eos uses and the table's
! edge in density at each T6 row (find_opal_ramp, 1995 and 2001). Called
! once per vintage from interpolate_opal_eos, on the unit eos_init's
! open_opal_eos_table opened. The vintage's row of opal_layout gives the
! table's size, its variable count and its formats; the header's moles
! per gram, mean molecular weight and species fractions, and each row's
! log10 R (1995, 2001) or amu and log10 n_e (2006), are read and not
! kept. Every X block must be on the first block's grid: its density at
! each column and the T6 of each row it reads (the grids are kept from
! the first block). A read that fails, a block out of its layout or off
! the first block's grid (M291) and a zero T6 are ierr = 1, with a
! message, the unit closed.
subroutine readcoeos(v, vintage, ierr)
      integer, intent(in) :: vintage

      type(opal_eos_vintage), intent(inout) :: v
      integer, intent(out) :: ierr

      integer, parameter :: mx = n_eos_mx, mv = n_eos_mv
! the vintage's table size
      integer :: nr, nt

      character(len=1) :: blank_line

! --- locals ---
      integer :: x_loop_index
      integer :: x_idx, var_idx, t6_idx, r_idx, t6_scan_idx
      integer :: density_row, t6_row, record_number
      integer :: column_idx, table_var_idx
      real(dp) :: unused_field, unused_header(2), unused_species(6), unused_columns(2)
! a T6 row's temperature as read, before it is compared with the first
! block's
      real(dp) :: t6_read
      integer :: ios
      logical :: layout_error

      nr = opal_layout(vintage)%nr
      nt = opal_layout(vintage)%nt
      ierr = 0
      if (opal_shared%table_unit == 0) then
         call report_no_table(vintage)
         ierr = 1
         return
      end if

      blank_line = ' '

      if (opal_shared%readco_init_flag(vintage) /= opal_flag_set) then
         do x_idx = 1, mx
            do var_idx = 1, mv
               do t6_idx = 1, nt
                  do r_idx = 1, nr
                     v%eos_table(x_idx, var_idx, t6_idx, r_idx) = opal_empty_entry
                  end do
               end do
            end do
         end do
         do var_idx = 1, mv
            v%eos_output(var_idx) = 1.0d0
         end do
         opal_shared%readco_init_flag(vintage) = opal_flag_set
      end if

! .....read  tables (on the unit open_opal_eos_table opened); a read
! that fails or a block out of its layout ends the read with ierr = 1,
! the unit closed
      call read_table_blocks(ios, layout_error)
      if (ios /= 0 .or. layout_error) then
         if (ios /= 0) call report_bad_table(vintage, x_loop_index, density_row)
         close (opal_shared%table_unit)
         opal_shared%table_unit = 0
         ierr = 1
         return
      end if

      do t6_scan_idx = 1, nt
         if (v%t6_list(1,t6_scan_idx) == 0.0d0) then
            call report_zero_t6(vintage, t6_scan_idx)
            close (opal_shared%table_unit)
            opal_shared%table_unit = 0
            ierr = 1
            return
         end if
         v%t6_grid(t6_scan_idx) = v%t6_list(1,t6_scan_idx)
      end do
      do t6_idx = 2, nt
         v%t6_grid_spacing_inv(t6_idx) = 1.0d0/(v%t6_grid(t6_idx) - &
              v%t6_grid(t6_idx-1))
      end do
      v%density_grid(1) = v%density_grid_table(1,1)
      do r_idx = 2, nr
         v%density_grid(r_idx) = v%density_grid_table(1,r_idx)
         v%density_grid_spacing_inv(r_idx) = 1.0d0/(v%density_grid(r_idx) - &
              v%density_grid(r_idx-1))
      end do
      do x_idx = 2, mx
         v%x_grid_spacing_inv(x_idx) = 1.0d0/(v%x_grid_copy(x_idx) - v%x_grid_copy(x_idx-1))
      end do
      close (opal_shared%table_unit)
      opal_shared%table_unit = 0
! MHP 7/2003 added ramp between OPAL and other EOS: the edge of the
! table at high rho, fixed T (find_opal_ramp reads it for 2001)
      v%t_row_index = 1
      do t6_idx = 1, nt
         v%density_edge_at_t(t6_idx) = v%density_grid(opal_layout(vintage)%density_index_edge(t6_idx))
      end do

contains

! Each X block of the table: its header, then for each density column
! its record (which must be the column's number, with the first block's
! density and at least the T6 rows the table needs) and its T6 rows in
! the vintage's row format, each read row's T6 the first block's (a
! later block's is compared with the one before it, which the first
! block's run of checks makes the first block's).
subroutine read_table_blocks(ios, layout_error)
      integer, intent(out) :: ios
      logical, intent(out) :: layout_error
      integer :: u
      u = opal_shared%table_unit
      layout_error = .false.
      density_row = 0
      do x_loop_index = 1, mx
         read (u,opal_layout(vintage)%header_format,iostat=ios) &
              v%hydrogen_fraction_header(x_loop_index), &
              v%z_table(x_loop_index), unused_header(1), unused_header(2)
         if (ios /= 0) return
         read (u,'(21X,E14.7,4X,E14.7,3X,E11.4,3X,E11.4,3X,E11.4, &
              &4X,E11.4)',iostat=ios) unused_species
         if (ios /= 0) return
         read (u,'(A)',iostat=ios) blank_line
         if (ios /= 0) return
         do density_row = 1, nr
            read (u,'(2I5,2F12.7,17X,E15.7)',iostat=ios) record_number, &
                 v%temperature_count_used(x_loop_index,density_row), &
                 unused_field, unused_field, &
                 v%density_grid_table(x_loop_index,density_row)
            if (ios /= 0) return
            if (record_number /= density_row) then
               call report_record_number(vintage, record_number, density_row)
               layout_error = .true.
               return
            end if
            if (v%density_grid_table(x_loop_index,density_row) /= &
                 v%density_grid_table(1,density_row)) then
               call report_opal_grid_mismatch(vintage, x_loop_index, density_row, 0, &
                    v%density_grid_table(x_loop_index,density_row), &
                    v%density_grid_table(1,density_row))
               layout_error = .true.
               return
            end if
            read(u,'(A)',iostat=ios) blank_line
            if (ios == 0) read(u,'(A)',iostat=ios) blank_line
            if (ios /= 0) return
            if (v%temperature_count_used(x_loop_index,density_row) < &
                 opal_layout(vintage)%t6_index_lo(density_row)) then
               call report_short_column(vintage, v%hydrogen_fraction_header(x_loop_index), &
                    v%density_grid_table(x_loop_index,density_row))
               layout_error = .true.
               return
            end if
            do t6_row = 1, v%temperature_count_used(x_loop_index,density_row)
               if (t6_row > opal_layout(vintage)%t6_index_lo(density_row)) then
                  read (u,'(A)',iostat=ios) blank_line
                  if (ios /= 0) return
                  cycle
               end if
               read (u, opal_layout(vintage)%row_format,iostat=ios) t6_read, &
                    (unused_columns(column_idx), column_idx=1,opal_layout(vintage)%n_skipped_columns), &
                    (v%eos_table(x_loop_index,table_var_idx,t6_row,density_row), &
                    table_var_idx=1,opal_layout(vintage)%n_vars)
               if (ios /= 0) return
               if (x_loop_index > 1 .and. t6_read /= v%t6_list(density_row,t6_row)) then
                  call report_opal_grid_mismatch(vintage, x_loop_index, density_row, t6_row, &
                       t6_read, v%t6_list(density_row,t6_row))
                  layout_error = .true.
                  return
               end if
               v%t6_list(density_row,t6_row) = t6_read
            end do
            read(u,'(A)',iostat=ios) blank_line
            if (ios == 0) read(u,'(A)',iostat=ios) blank_line
            if (ios == 0) read(u,'(A)',iostat=ios) blank_line
            if (ios /= 0) return
         end do
         read(u,'(A)',iostat=ios) blank_line
         if (ios /= 0) return
      end do
end subroutine read_table_blocks

end subroutine readcoeos

!----------------------------------------------------------------------
! The reader's messages, one text for every vintage after its name
! (opal_layout's eos_name), to the terminal and the run log (eos_report).
!----------------------------------------------------------------------
! the reader: a density block out of order
subroutine report_record_number(vintage, record_number, density_row)
      integer, intent(in) :: vintage, record_number, density_row
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,i0,a,i0,a)') trim(opal_layout(vintage)%eos_name), &
           ': the table is out of its layout: record ', record_number, ' where density column ', &
           density_row, ' was expected'
      call eos_report(eos_lines)
end subroutine report_record_number

! the reader: a density column with fewer T6 rows than the table needs
subroutine report_short_column(vintage, hydrogen_fraction, density)
      integer, intent(in) :: vintage
      real(dp), intent(in) :: hydrogen_fraction, density
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,f6.4,a,es14.4)') trim(opal_layout(vintage)%eos_name), &
           ': a density column of the table has fewer T6 rows than it needs: X ', &
           hydrogen_fraction, ', density', density
      call eos_report(eos_lines)
end subroutine report_short_column

! the reader: an X block off the first block's grid, at a density column
! (t6_row 0: its density) or at a T6 row of it (that row's T6), with the
! value read and the first block's
subroutine report_opal_grid_mismatch(vintage, x_block, density_column, t6_row, value, first_value)
      integer, intent(in) :: vintage, x_block, density_column, t6_row
      real(dp), intent(in) :: value, first_value
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      if (t6_row == 0) then
         write(eos_lines,'(a,a,i0,a,i0,a,es15.7,a,es15.7)') trim(opal_layout(vintage)%eos_name), &
              ': the table''s X block ', x_block, ' is not on the first block''s grid: ' // &
              'the density of column ', density_column, ' is', value, ' against', first_value
      else
         write(eos_lines,'(a,a,i0,a,i0,a,i0,a,f11.6,a,f11.6)') trim(opal_layout(vintage)%eos_name), &
              ': the table''s X block ', x_block, ' is not on the first block''s grid: ' // &
              'the T6 of density column ', density_column, ', row ', t6_row, ' is', value, &
              ' against', first_value
      end if
      call eos_report(eos_lines)
end subroutine report_opal_grid_mismatch

! the reader: no table open (eos_init opens the selected vintage's; a
! vintage installed without it, or a table read before, has none)
subroutine report_no_table(vintage)
      integer, intent(in) :: vintage
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a)') trim(opal_layout(vintage)%eos_name), &
           ': no table is open to read (eos_init opens it)'
      call eos_report(eos_lines)
end subroutine report_no_table

! the reader: a read that failed (the table ends early or is corrupt)
subroutine report_bad_table(vintage, x_block, density_row)
      integer, intent(in) :: vintage, x_block, density_row
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,i0,a,i0)') trim(opal_layout(vintage)%eos_name), &
           ': the table ends early or is corrupt: X block ', x_block, ', density column ', &
           density_row
      call eos_report(eos_lines)
end subroutine report_bad_table

! the reader: a zero T6
subroutine report_zero_t6(vintage, t6_idx)
      integer, intent(in) :: vintage, t6_idx
      character(len=message_length) :: eos_lines(4)
      eos_lines = unwritten_line
      write(eos_lines,'(a,a,i0,a)') trim(opal_layout(vintage)%eos_name), ': the table''s T6 row ', &
           t6_idx, ' is zero'
      call eos_report(eos_lines)
end subroutine report_zero_t6

end module opal_table_reader

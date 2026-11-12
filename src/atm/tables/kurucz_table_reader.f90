!----------------------------------------------------------------------
! kurucz_table_reader
!----------------------------------------------------------------------
! Reads a surface-pressure table, Kurucz (1993) or Castelli & Kurucz
! (2003), at atm_init:
!     read_surface_pressure_table  the file into a surface_pressure_table,
!                                  and each row's valid log g columns
! The layout (fixed formats, values in 1pe16.8): a header line with the
! log metal abundance; "LOG TEFF:" and the num_teff log Teff rows, four
! to a line; "LOG G:" and the 11 log g columns; then, for each row, its
! 11 log P on three lines. -999 (any log P <= 0) marks a point the model
! grid does not have, at the low- or high-gravity end of a row.
module kurucz_table_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use atm_tables, only: surface_pressure_table, surface_pressure_num_g, &
           kurucz_num_teff
      use atm_messages, only: report_table_open_failure, report_table_read_failure
      implicit none
      private
      public :: read_surface_pressure_table

      integer, parameter :: ng = surface_pressure_num_g
! the formats: the header's abundance; the Kurucz table's 57 and the
! Castelli/Kurucz table's 76 log Teff; the log g columns; a row's log P
      character(len=*), parameter :: abundance_format = '(1x,10x,1pe16.8,/)'
      character(len=*), parameter :: kurucz_teff_format = '(1p4e16.8,13(/1p4e16.8),/1pe16.8,/)'
      character(len=*), parameter :: castelli_teff_format = '(1p4e16.8,18(/1p4e16.8),/)'
      character(len=*), parameter :: log_g_format = '(1p4e16.8,/1p4e16.8,/1p3e16.8,/)'
      character(len=*), parameter :: row_format = '(1p4e16.8,/1p4e16.8,/1p3e16.8)'

contains

!----------------------------------------------------------------------
! read_surface_pressure_table
!----------------------------------------------------------------------
! Reads the table at table_path, of num_teff rows (kurucz_num_teff or
! castelli_num_teff), into tbl, then finds each row's first and last
! valid log g column. A file that cannot be opened or read returns
! ierr = 1 after its message.
subroutine read_surface_pressure_table(table_path, num_teff, tbl, ierr)
      character(len=*), intent(in) :: table_path
      integer, intent(in) :: num_teff
      type(surface_pressure_table), intent(out) :: tbl
      integer, intent(out) :: ierr

      integer :: table_unit, row, ios
      real(dp) :: log_metal_abundance
      character(len=256) :: message

      ierr = 0
      allocate(tbl%log_teff(num_teff), tbl%log_g(ng), tbl%log10_pressure(num_teff,ng), &
           tbl%first_column(num_teff), tbl%last_column(num_teff))
      tbl%log_teff = 0.0d0
      tbl%log_g = 0.0d0
      tbl%log10_pressure = 0.0d0
      tbl%first_column = 0
      tbl%last_column = 0
      tbl%num_teff = num_teff

      open(newunit=table_unit, file=table_path, status='OLD', iostat=ios)
      if (ios /= 0) then
         call report_table_open_failure(table_path)
         ierr = 1
         return
      end if
      read(table_unit, abundance_format, iostat=ios, iomsg=message) log_metal_abundance
      if (ios == 0) then
         if (num_teff == kurucz_num_teff) then
            read(table_unit, kurucz_teff_format, iostat=ios, iomsg=message) tbl%log_teff
         else
            read(table_unit, castelli_teff_format, iostat=ios, iomsg=message) tbl%log_teff
         end if
      end if
      if (ios == 0) read(table_unit, log_g_format, iostat=ios, iomsg=message) tbl%log_g
      do row = 1, num_teff
         if (ios /= 0) exit
         read(table_unit, row_format, iostat=ios, iomsg=message) tbl%log10_pressure(row,:)
      end do
      close(table_unit)
      if (ios /= 0) then
         call report_table_read_failure(table_path, message)
         ierr = 1
         return
      end if

      call find_valid_columns(tbl)
      tbl%loaded = .true.
end subroutine read_surface_pressure_table

! Each row's last valid log g column (the highest with log P > 0) and
! its first (the lowest of the valid run that ends there; 1 when the run
! reaches the first column). A row without a valid entry gets first
! column ng and last column 0, which the lookup refuses.
pure subroutine find_valid_columns(tbl)
      type(surface_pressure_table), intent(inout) :: tbl
      integer :: row, column
      logical :: found_valid_pressure

      do row = 1, tbl%num_teff
         found_valid_pressure = .false.
         do column = ng, 1, -1
            if (tbl%log10_pressure(row,column) <= 0.0d0) then
               if (found_valid_pressure) then
                  tbl%first_column(row) = column + 1
                  exit
               end if
            else
               if (.not. found_valid_pressure) tbl%last_column(row) = column
               found_valid_pressure = .true.
            end if
         end do
         if (column < 1) tbl%first_column(row) = 1
         if (.not. found_valid_pressure) tbl%first_column(row) = ng
      end do
end subroutine find_valid_columns

!======================================================================
! Sources: the F77's SETUPS read (SURFP's table) and its Castelli/Kurucz
! twin (JNT 6/14); G. Somers (5/15) added the per-row valid
! columns.

end module kurucz_table_reader

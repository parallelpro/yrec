!----------------------------------------------------------------------
! output_columns
!----------------------------------------------------------------------
! What the history and profile writers share: the column record (a
! name and how it is written), and the selection of columns from a
! columns file. Format-blind: each writer hands parse_columns its own
! column table and default selection.
!
! A columns file has one column name per line, '!' comments, blank
! lines ignored; unknown names are fatal (config error), with the
! valid names listed in the run log. A blank/absent control selects
! the compiled-in default: default_columns.inc is GENERATED from the
! authoritative defaults/{history,profile}_columns.list (uncommented
! names, file order) by tools/gen_default_columns.py -- edit a .list
! and rebuild (make reruns the generator) to change the defaults.
module output_columns
      implicit none
      private
      public :: max_cols, column_name_len, parse_columns
      public :: output_column, real_column, integer_column
      public :: column_integer_format, column_real_format, column_numbers_format, &
           column_names_format
      public :: hist_default_names, n_hist_default
      public :: prof_default_names, n_prof_default

      integer, parameter :: max_cols = 128
! the longest column name (gen_default_columns.py's WIDTH)
      integer, parameter :: column_name_len = 24
! How a column is written: real_column as es40.16e3, integer_column as
! i40 (the id columns, the zone number, the mixing type).
      integer, parameter :: real_column = 1, integer_column = 2
! MESA's layout of a history or profile file: every value in a field of
! 40 after a blank, an integer value, a real value, and a line of
! column numbers or names (the global block's and the columns').
      character(len=*), parameter :: column_integer_format = '(1x,i40)', &
           column_real_format = '(1x,es40.16e3)', &
           column_numbers_format = '(999(1x,i40))', &
           column_names_format = '(999(1x,a40))'

! One column of a history or profile file: its name in the file and how
! it is written. Each writer's table holds one per column, at the
! column's named index (ih_* in history_output, ip_* in stitched_model).
      type :: output_column
         character(len=column_name_len) :: name = ' '
         integer :: kind = real_column
      end type output_column

      include 'default_columns.inc'

contains

! ---------------------------------------------------------------
subroutine parse_columns(fname, columns, default_names, sel, nsel, &
           label, ierr)
      use run_messages, only: check_io, report_error
      character(len=*), intent(in) :: fname, label
      type(output_column), intent(in) :: columns(:)
      character(len=column_name_len), intent(in) :: default_names(:)
      integer, intent(out) :: sel(max_cols), nsel
      integer, intent(out) :: ierr
      character(len=256) :: line, msg
      integer :: u, ios, i

      ierr = 0
      nsel = 0
      if (len_trim(fname) == 0) then
! blank columns file: the compiled-in default selection (generated
! from defaults/<label>_columns.list -- see default_columns.inc)
         do i = 1, size(default_names)
            call append_column(default_names(i), columns, sel, nsel, &
                 label, ierr)
            if (ierr /= 0) return
         end do
         return
      end if
      open(newunit=u, file=fname, status='OLD', action='READ', iostat=ios, &
           iomsg=msg)
      call check_io(ios, msg, 'open the ' // label // '_columns_file', fname, ierr)
      if (ierr /= 0) return
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         i = index(line, '!')
         if (i > 0) line = line(1:i-1)
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         call append_column(line, columns, sel, nsel, label, ierr)
         if (ierr /= 0) then
            close(u)
            return
         end if
      end do
      close(u)
      if (nsel == 0) then
         call report_error('the ' // label // '_columns_file ' // trim(fname) // &
              ' selects no columns')
         ierr = 1
         return
      end if
end subroutine parse_columns

! ---------------------------------------------------------------
! Look a column name up in the writer's table and append its index
! to the selection; unknown names are fatal, with the valid list
! written to the run log.
subroutine append_column(name, columns, sel, nsel, label, ierr)
      use run_messages, only: run_log_unit, report_error
      character(len=*), intent(in) :: name, label
      type(output_column), intent(in) :: columns(:)
      integer, intent(inout) :: sel(max_cols), nsel
      integer, intent(out) :: ierr
      character(len=80) :: text
      integer :: j

      ierr = 0
      if (nsel >= max_cols) then
         write(text, '(a,i0,3a)') 'more than ', max_cols, ' ', label, ' columns selected'
         call report_error(trim(text))
         ierr = 1
         return
      end if
      do j = 1, size(columns)
         if (trim(name) == trim(columns(j)%name)) then
            nsel = nsel + 1
            sel(nsel) = j
            return
         end if
      end do
      call report_error('unknown ' // label // ' column: ' // trim(name) // &
           ' (the run log lists the valid ones)')
      write(run_log_unit,'(3a)') ' valid ', label, ' columns:'
      do j = 1, size(columns)
         write(run_log_unit,'(2x,a)') trim(columns(j)%name)
      end do
      ierr = 1
end subroutine append_column

end module output_columns

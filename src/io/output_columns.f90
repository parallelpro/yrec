!----------------------------------------------------------------------
! output_columns_lib
!----------------------------------------------------------------------
! Column-selection machinery shared by the history and profile
! writers (split out of yrec_output, 2026). Format-blind: callers
! hand parse_columns their own name table and default selection.
!
! A columns file has one column name per line, '!' comments, blank
! lines ignored; unknown names are fatal (config error), with the
! valid names listed in the run log. A blank/absent control selects
! the compiled-in default: default_columns.inc is GENERATED from the
! authoritative defaults/{history,profile}_columns.list (uncommented
! names, file order) by tools/gen_default_columns.py -- edit a .list
! and rebuild (make reruns the generator) to change the defaults.
module output_columns_lib
      implicit none
      private
      public :: max_cols, parse_columns
      public :: hist_default_names, n_hist_default
      public :: prof_default_names, n_prof_default

      integer, parameter :: max_cols = 128
      include 'default_columns.inc'

contains

! ---------------------------------------------------------------
subroutine parse_columns(fname, names, ncol, default_names, ndefault, &
           sel, nsel, label, ierr)
      use luout_lib
      character(len=*), intent(in) :: fname, label
      integer, intent(out) :: ierr
      integer, intent(in) :: ncol, ndefault
      character(len=24), intent(in) :: names(ncol), default_names(ndefault)
      integer, intent(out) :: sel(max_cols), nsel
      character(len=256) :: line
      integer :: u, ios, i

      ierr = 0
      nsel = 0
      if (len_trim(fname) == 0) then
! blank columns file: the compiled-in default selection (generated
! from defaults/<label>_columns.list -- see default_columns.inc)
         do i = 1, ndefault
            call append_column(default_names(i), names, ncol, sel, nsel, &
                 label, ierr)
            if (ierr /= 0) return
         end do
         return
      end if
      open(newunit=u, file=fname, status='OLD', action='READ', iostat=ios)
      if (ios /= 0) then
         write(*,*) 'cannot open ', trim(label), '_columns_file: ', &
              trim(fname)
         write(run_log_unit,*) 'cannot open ', trim(label), &
              '_columns_file: ', trim(fname)
! 2026 io-writer stops -> ierr (stage-3 pattern): config error
! returned to output_init_mesa -> read_input -> read_controls.
         ierr = 1
         return
      end if
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         i = index(line, '!')
         if (i > 0) line = line(1:i-1)
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         call append_column(line, names, ncol, sel, nsel, label, ierr)
         if (ierr /= 0) return
      end do
      close(u)
      if (nsel == 0) then
         write(*,*) trim(label), '_columns_file selected no columns: ', &
              trim(fname)
         ierr = 1
         return
      end if
end subroutine parse_columns

! ---------------------------------------------------------------
! Look a column name up in the writer's table and append its index
! to the selection; unknown names are fatal, with the valid list
! written to the run log.
subroutine append_column(name, names, ncol, sel, nsel, label, ierr)
      use luout_lib
      character(len=*), intent(in) :: name, label
      integer, intent(in) :: ncol
      character(len=24), intent(in) :: names(ncol)
      integer, intent(inout) :: sel(max_cols), nsel
      integer, intent(out) :: ierr
      integer :: j

      ierr = 0
      do j = 1, ncol
         if (trim(name) == trim(names(j))) then
            nsel = nsel + 1
            sel(nsel) = j
            return
         end if
      end do
      write(*,*) 'unknown ', trim(label), ' column: ', trim(name)
      write(run_log_unit,*) 'unknown ', trim(label), &
           ' column: ', trim(name)
      write(run_log_unit,*) 'valid ', trim(label), ' columns:'
      do j = 1, ncol
         write(run_log_unit,'(2x,a)') trim(names(j))
      end do
      ierr = 1
end subroutine append_column

end module output_columns_lib

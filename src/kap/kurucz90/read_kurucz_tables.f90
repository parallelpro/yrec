!----------------------------------------------------------------------
! read_kurucz_tables
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setkrz.f; only variable names, source form, and comment
! style were updated.
!
! Reads the Kurucz opacity table (and, if a second Z table is
! requested, a second Kurucz table) into opacity_table%kurucz(1)/(2),
! then calls build_kurucz_splines on each for the spline
! interpolation coefficients used later by kurucz.f90.
!
! 2026 wave 3 (R5): the two verbatim read loops (first table on
! star%ctrl%kurucz_table_unit, second on kurucz_table2_unit) became
! the internal read_kurucz_table below, called once per set; the
! order read 1, read 2, splines 1, splines 2 is the original's. The
! loop's working scalars stay host-associated so that whatever value
! they carried from the first table into the second is carried
! exactly as before.
subroutine read_kurucz_tables(kurucz_table_path, kurucz_table2_path, ierr)
      use star_info_lib, only: star

      use opacity_table_lib
      use math_lib
      implicit none

      character(len=256), intent(in) :: kurucz_table_path, kurucz_table2_path
      integer, intent(out) :: ierr
      integer :: num_read, density_index
      integer :: read_status
      double precision :: prev_grid_temp, grid_temp, pressure, &
           log10_opacity0, log10_opacity1, log10_opacity2, &
           log10_opacity4, log10_opacity8, electron_density, &
           atom_density, density, unused_col

      ierr = 0
! The early-out range check in kurucz() belongs to the first table
! only (the former kurucz2.f90 never had it).
      opacity_table%kurucz(2)%check_range = .false.

      call read_kurucz_table(star%ctrl%kurucz_table_unit, kurucz_table_path, &
           opacity_table%kurucz(1), ierr)
      if (ierr /= 0) return

! DBG 12/95 read in second Z table if requested
      if (star%use_two_z_tables) then
         call read_kurucz_table(star%ctrl%kurucz_table2_unit, kurucz_table2_path, &
              opacity_table%kurucz(2), ierr)
         if (ierr /= 0) return
      end if

      call build_kurucz_splines(opacity_table%kurucz(1), ierr)
      if (ierr /= 0) return
! DBG 12/95 second Z table
      if (star%use_two_z_tables) then
         call build_kurucz_splines(opacity_table%kurucz(2), ierr)
      end if

      return

contains

! read_kurucz_table: read one Kurucz table file on `unit` into tbl.
subroutine read_kurucz_table(unit, path, tbl, ierr)
      integer, intent(in) :: unit
      character(len=256), intent(in) :: path
      type(kurucz_table_set), intent(inout) :: tbl
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success

    1 format(2f5.2,5f7.3,3f9.5,f8.3)
!     OPEN TABLE
      open(unit, file=path, status='OLD')
      read(unit,'(/)')

      num_read = 0
      prev_grid_temp = 0.0d0
      table_read: do
      read(unit,1,iostat=read_status) grid_temp, pressure, &
           log10_opacity0, log10_opacity1, log10_opacity2, &
           log10_opacity4, log10_opacity8, electron_density, &
           atom_density, density, unused_col
      if (read_status .gt. 0) cycle table_read   ! was err=110: reread next record
      if (read_status .lt. 0) exit table_read    ! was end=120
      if (prev_grid_temp.ne.grid_temp) then
         num_read = num_read+1
         if (num_read.gt.kurucz_max_num_temps) then
            write(*,*) 'read_kurucz_tables: too many temperature rows in ', &
                 trim(path)
            ierr = 1
            return
         end if
         tbl%grid_logt(num_read) = grid_temp
         prev_grid_temp = grid_temp
         density_index = 1
      endif
      if (num_read.lt.1) then
         write(*,*) 'read_kurucz_tables: no header row read from ', &
              trim(path)
         ierr = 1
         return
      end if
      tbl%log10_rho(num_read, density_index) = density
      tbl%log10_opacity(num_read, density_index) = exp10(log10_opacity0)
      density_index = density_index+1
      end do table_read

      tbl%num_temps = num_read
!     CLOSE THE TABLE WE HAVE READ
      close(unit,err=99)
      return
   99 continue
      write(*,*) 'read_kurucz_tables: error closing table file'
      ierr = 1
end subroutine read_kurucz_table

end subroutine read_kurucz_tables

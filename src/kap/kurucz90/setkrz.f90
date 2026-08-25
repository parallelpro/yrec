!----------------------------------------------------------------------
! setkrz
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setkrz.f; only variable names, source form, and comment
! style were updated.
!
! Reads the Kurucz opacity table (and, if a second Z table is
! requested, a second Kurucz table) and builds the spline
! interpolation coefficients used later by kurucz.f90/kurucz2.f90.
subroutine setkrz(kurucz_table_path, kurucz_table2_path)
      use star_info_lib, only: star

      use opacity_table_lib
      implicit none
      integer, parameter :: max_num_temps = 60
      integer, parameter :: max_num_densities = 50
      integer, parameter :: num_x_tables = 1
      integer, parameter :: num_x_temp_entries = max_num_temps*num_x_tables

! MHP 8/25 Remove unused variables
      character(len=256), intent(in) :: kurucz_table_path, kurucz_table2_path
      integer :: x_table_count, num_read, density_index
      integer :: read_status
      double precision :: prev_grid_temp, grid_temp, pressure, &
           log10_opacity0, log10_opacity1, log10_opacity2, &
           log10_opacity4, log10_opacity8, electron_density, &
           atom_density, density, unused_col

      x_table_count = num_x_tables
    1 format(2f5.2,5f7.3,3f9.5,f8.3)
!     OPEN TABLE
      open(star%ctrl%kurucz_table_unit, file=kurucz_table_path, status='OLD')
      read(star%ctrl%kurucz_table_unit,'(/)')

      num_read = 0
      prev_grid_temp = 0.0d0
      table_read: do
      read(star%ctrl%kurucz_table_unit,1,iostat=read_status) grid_temp, pressure, &
           log10_opacity0, log10_opacity1, log10_opacity2, &
           log10_opacity4, log10_opacity8, electron_density, &
           atom_density, density, unused_col
      if (read_status .gt. 0) cycle table_read   ! was err=110: reread next record
      if (read_status .lt. 0) exit table_read    ! was end=120
      if (prev_grid_temp.ne.grid_temp) then
         num_read = num_read+1
         if (num_read.gt.max_num_temps) stop ' KURUCZ INPUT ERROR'
         opacity_table%kurucz_grid_logt(num_read) = grid_temp
         prev_grid_temp = grid_temp
         density_index = 1
      endif
      if (num_read.lt.1) stop ' KURUCZ INPUT ERROR'
      opacity_table%kurucz_log10_rho(num_read, density_index) = density
      opacity_table%kurucz_log10_opacity(num_read, density_index) = 10.0d0**log10_opacity0
      density_index = density_index+1
      end do table_read

      opacity_table%kurucz_num_temps = num_read
!     CLOSE THE TABLE WE HAVE READ
      close(star%ctrl%kurucz_table_unit,err=99)

! DBG 12/95 read in second Z table if requested
      if (star%use_two_z_tables) then
!        OPEN TABLE
         open(star%ctrl%ikur2, file=kurucz_table2_path, status='OLD')
         read(star%ctrl%ikur2,'(/)')

         num_read = 0
         prev_grid_temp = 0.0d0
         table2_read: do
         read(star%ctrl%ikur2,1,iostat=read_status) grid_temp, pressure, &
              log10_opacity0, log10_opacity1, log10_opacity2, &
              log10_opacity4, log10_opacity8, electron_density, &
              atom_density, density, unused_col
         if (read_status .gt. 0) cycle table2_read   ! was err=210
         if (read_status .lt. 0) exit table2_read    ! was end=220
         if (prev_grid_temp.ne.grid_temp) then
            num_read = num_read+1
            if (num_read.gt.max_num_temps) stop ' KURUCZ INPUT ERROR'
            opacity_table%kurucz2_grid_logt(num_read) = grid_temp
            prev_grid_temp = grid_temp
            density_index = 1
         endif
         if (num_read.lt.1) stop ' KURUCZ INPUT ERROR'
         opacity_table%kurucz2_log10_rho(num_read, density_index) = density
         opacity_table%kurucz2_log10_opacity(num_read, density_index) = 10.0d0**log10_opacity0
         density_index = density_index+1
         end do table2_read

         opacity_table%kurucz2_num_temps = num_read
!        CLOSE THE TABLE WE HAVE READ
         close(star%ctrl%ikur2,err=99)
      end if

      call ykoeff

      return
   99 stop 'ERROR IN FILE CLOSING'
end subroutine setkrz

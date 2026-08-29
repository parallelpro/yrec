!----------------------------------------------------------------------
! setllo
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setllo.f; only variable names, source form, and comment
! style were updated. FORMAT strings and the READ/OPEN statements
! that must byte-for-byte match the on-disk table file are preserved
! verbatim.
!
! DBG 5/94 Modified to read in second opacity table at different Z.
! Reads the OPAL92 (Lawrence Livermore) opacity table(s) and builds
! the spline interpolation coefficients used later by opal92_surface_table.f90 and
! the OPAL92 lookup routines (yllo3d/yllo3d2, not part of this batch).
subroutine read_opal92_tables(opal92_table_path, opal92_table2_path)
      use star_info_lib, only: star

      use opacity_table_lib
      implicit none
! runtime-allocated unit for the second (Z-ramp) OPAL92 table
! (formerly luout_lib's fixed ioopal2 = 64)
      integer :: opal92b_unit
      integer, parameter :: num_t = 50
      integer, parameter :: num_d = 17
      integer, parameter :: num_x = 3
      integer, parameter :: num_xt = num_t*num_x
! CONST is declared but not referenced anywhere in the original
! setllo.f body; preserved here as an unused placeholder.
      double precision, parameter :: const_unused = 11604.5d0

! MHP 8/25 removed variables not used in subroutine
      character(len=256), intent(in) :: opal92_table_path, opal92_table2_path
      double precision :: local_grid_y(num_x), local_grid_z(num_x)
      integer :: i, k, density_index, num_temps_read
      double precision :: grid_temp_k

!     OPEN TABLE
      open(unit=star%ctrl%laol_table_unit,file=opal92_table_path)
      do i=1,num_x
!        READ GRID POINT FOR ABUNDANCE
!        READ NUMBER OF GRIDS FOR DENSITY, AND TEMPERATURE
        read(star%ctrl%laol_table_unit,190,end=97) opacity_table%opal92_grid_x(i), local_grid_z(i)
        local_grid_y(i)=1.0d0-opacity_table%opal92_grid_x(i)-local_grid_z(i)
  190   format(33x,f7.4,2x,f7.4)
         read(star%ctrl%laol_table_unit,'()')
!        READ  LOG(DENSITY/TEMPERATURE**3)
            read(star%ctrl%laol_table_unit, 200) (opacity_table%opal92_grid_logr(density_index), density_index=1, num_d)
  200   format (6x, 17f7.1)
!        READ GRID VALUES FOR TEMPERATURE, AND OPACITY TABLE
         do k=1, num_t
         read(star%ctrl%laol_table_unit,196,end=93) grid_temp_k, &
              (opacity_table%opal92_log10_opacity(k+(i-1)*num_t,density_index),density_index=1,num_d)
         opacity_table%opal92_grid_logt(k)=dlog10(grid_temp_k)
         end do
   93    num_temps_read=k-1
  196    format(18f7.3)
!
      end do
!     CLOSE THE TABLE WE HAVE READ
   97 close(star%ctrl%laol_table_unit,err=99)
      opacity_table%opal92_num_temps = num_temps_read
      opacity_table%opal92_num_x=i-1

! DBG 5/94 Second Opacity Table read here
      if (star%use_two_z_tables) then
         open(newunit=opal92b_unit,file=opal92_table2_path)
         do i=1,num_x
            read(opal92b_unit,190,end=597) opacity_table%opal92_grid_x_z2(i), local_grid_z(i)
            local_grid_y(i)=1.0d0-opacity_table%opal92_grid_x_z2(i)-local_grid_z(i)
            read(opal92b_unit,'()')
            read(opal92b_unit, 200) (opacity_table%opal92_grid_logr_z2(density_index), density_index=1, num_d)
            do k=1, num_t
               read(opal92b_unit,196,end=593) grid_temp_k, &
                    (opacity_table%opal92_log10_opacity_z2(k+(i-1)*num_t,density_index),density_index=1,num_d)
               opacity_table%opal92_grid_logt_z2(k)=dlog10(grid_temp_k)
            end do
  593       num_temps_read=k-1
         end do
  597    close(opal92b_unit,err=99)
         opacity_table%opal92_num_temps_z2 = num_temps_read
         opacity_table%opal92_num_x_z2=i-1
      end if
!
      call opal92_table_prep

      return
   99 stop 'ERROR IN FILE CLOSING'
end subroutine read_opal92_tables

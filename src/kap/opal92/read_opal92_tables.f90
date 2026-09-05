!----------------------------------------------------------------------
! read_opal92_tables
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setllo.f; only variable names, source form, and comment
! style were updated. FORMAT strings and the READ/OPEN statements
! that must byte-for-byte match the on-disk table file are preserved
! verbatim.
!
! DBG 5/94 Modified to read in second opacity table at different Z.
! Reads the OPAL92 (Lawrence Livermore) opacity table(s) and, via
! opal92_table_prep, builds the spline interpolation coefficients used
! later by opal92_surface_table.f90 and the OPAL92 lookup routine
! opal92_interp3d.
!
! 2026 wave 3 (R5): the verbatim second-Z read loop (into
! opacity_table%opal92(2)) became the internal read_opal92_table
! below, called once per set in the original order. The two OPEN
! statements stay here as they were (the first table on the fixed
! star%ctrl%opal92_table_unit, the second on a NEWUNIT); the READ,
! FORMAT and CLOSE statements moved into the helper verbatim, and the
! loop's working scalars stay host-associated.
subroutine read_opal92_tables(opal92_table_path, opal92_table2_path, ierr)
      use star_info_lib, only: star

      use opacity_table_lib
      use math_lib
      implicit none
! runtime-allocated unit for the second (Z-ramp) OPAL92 table
! (formerly luout_lib's fixed ioopal2 = 64)
      integer :: opal92b_unit

      character(len=256), intent(in) :: opal92_table_path, opal92_table2_path
      integer, intent(out) :: ierr
!     local_grid_z IS A READ-LIST TARGET THAT IS NEVER USED AFTERWARDS.
      double precision :: local_grid_z(n_opal92_x)
      integer :: i, k, density_index, num_temps_read
      double precision :: grid_temp_k

!     OPEN TABLE
      ierr = 0
      open(unit=star%ctrl%opal92_table_unit,file=opal92_table_path)
      call read_opal92_table(star%ctrl%opal92_table_unit, opacity_table%opal92(1), ierr)
      if (ierr /= 0) return

! DBG 5/94 Second Opacity Table read here
      if (star%use_two_z_tables) then
         open(newunit=opal92b_unit,file=opal92_table2_path)
         call read_opal92_table(opal92b_unit, opacity_table%opal92(2), ierr)
         if (ierr /= 0) return
      end if
!
      call opal92_table_prep(ierr)

      return

contains

! read_opal92_table: read the already-opened table on `unit` into tbl
! and close the unit.
subroutine read_opal92_table(unit, tbl, ierr)
      integer, intent(in) :: unit
      type(opal92_table_set), intent(inout) :: tbl
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success

      do i=1,n_opal92_x
!        READ GRID POINT FOR ABUNDANCE
!        READ NUMBER OF GRIDS FOR DENSITY, AND TEMPERATURE
        read(unit,190,end=97) tbl%grid_x(i), local_grid_z(i)
  190   format(33x,f7.4,2x,f7.4)
         read(unit,'()')
!        READ  LOG(DENSITY/TEMPERATURE**3)
            read(unit, 200) (tbl%grid_logr(density_index), density_index=1, n_opal92_d)
  200   format (6x, 17f7.1)
!        READ GRID VALUES FOR TEMPERATURE, AND OPACITY TABLE
         do k=1, n_opal92_t
         read(unit,196,end=93) grid_temp_k, &
              (tbl%log10_opacity(k+(i-1)*n_opal92_t,density_index),density_index=1,n_opal92_d)
         tbl%grid_logt(k)=log10(grid_temp_k)
         end do
   93    num_temps_read=k-1
  196    format(18f7.3)
!
      end do
!     CLOSE THE TABLE WE HAVE READ
   97 close(unit,err=99)
      tbl%num_temps = num_temps_read
      tbl%num_x=i-1
      return
   99 continue
      write(*,*) 'read_opal92_tables: error closing table file'
      ierr = 1
end subroutine read_opal92_table

end subroutine read_opal92_tables

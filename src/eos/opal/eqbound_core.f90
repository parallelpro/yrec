!----------------------------------------------------------------------
! eqbound_core
!----------------------------------------------------------------------
! Readability W3 (2026): the table-edge search and ramp-factor
! calculation shared by eqbound.f90 (OPAL 1995) and eqbound01.f90 (OPAL
! 2001). Those two routines were identical statement for statement
! from the T6 row search onward, differing only in the table they
! address and its row count nt; each wrapper keeps its own density
! pre-check (5.0 for 1995, 7.0 for 2001) and then calls this. Not used
! by eqbound06.f90, whose edge search is a different algorithm.
!
! Arguments: the vintage's row/column counts nt/nr, t6_grid (nt
! descending T6 rows), the caller's cached row index t_row_index
! (updated here), density_edge_at_t (table-edge density per row),
! density_index_edge (column of that edge per row) and density_grid
! (the nr density columns), then the point (t6, log10_density) and the
! eqbound outputs.
subroutine eqbound_core(nt, nr, t6_grid, t_row_index, density_edge_at_t, &
     density_index_edge, density_grid, t6, log10_density, ramp_factor, &
     in_opal_table, needs_ramp, ierr)

      use math_lib
      implicit none

      integer, intent(in) :: nt, nr
      double precision, intent(in) :: t6_grid(nt), density_edge_at_t(nt)
      integer, intent(in) :: density_index_edge(nt)
      double precision, intent(in) :: density_grid(nr)
      integer, intent(inout) :: t_row_index

      double precision, intent(in) :: t6, log10_density
      double precision, intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp
! --- locals ---
      double precision :: table_edge_density, ramp_start_density
      double precision :: t6_top_of_table
      integer :: t6_scan_idx

      integer, intent(out) :: ierr

      ierr = 0

!     find nearest table element in t.
      if (t6.lt.t6_grid(t_row_index)) then
         do t6_scan_idx = t_row_index+1, nt
            if (t6.ge.t6_grid(t6_scan_idx)) then
               t_row_index = t6_scan_idx - 1
               exit
            end if
         end do
         if (t6_scan_idx > nt) then
         t6_top_of_table = t6_grid(nt)
!        caller should have stopped outside table bounds; error exit
         write(*,5) t6, t6_top_of_table, t_row_index
    5    format(' ERROR IN OPAL EOS: OUTSIDE TABLE IN T6',2F10.6,I5)
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
         end if
      else
         do t6_scan_idx = t_row_index, 1, -1
            if (t6.le.t6_grid(t6_scan_idx)) then
               t_row_index = t6_scan_idx
               exit
            end if
         end do
         if (t6_scan_idx < (1)) then
         t6_top_of_table = t6_grid(1)
!        caller should have stopped outside table bounds; error exit
         write(*,5) t6, t6_top_of_table, t_row_index
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
         end if
      end if
!     define table edge in rho by linear interpolation.
      table_edge_density = density_edge_at_t(t_row_index+1)
      table_edge_density = log10(table_edge_density)
!     define beginning of ramp in the same fashion -
!     ramp is defined as one table element wide.
      ramp_start_density = density_grid(density_index_edge(t_row_index+1)-1)
      ramp_start_density = log10(ramp_start_density)
!     check if within table bounds in rho
      if (log10_density.gt.table_edge_density) then
         in_opal_table = .false.
         needs_ramp = .true.
         ramp_factor = 0d0
         return
      end if

!     If we get here, the point is in the table.
      in_opal_table = .true.

!     Now we check if ramping is needed.
!     First we check if ramping in temperature is needed.
      if (t6.le.t6_grid(nt-1)) then
         needs_ramp = .true.
         ramp_factor = (t6-t6_grid(nt))/(t6_grid(nt-1)-t6_grid(nt))
      else if (t6.ge.t6_grid(2)) then
         needs_ramp = .true.
         ramp_factor = (t6_grid(1)-t6)/(t6_grid(1)-t6_grid(2))
      else if (log10_density.ge.ramp_start_density) then
!        If we get here, ramping in density is needed.
         needs_ramp = .true.
         ramp_factor = (table_edge_density-log10_density)/ &
              (table_edge_density-ramp_start_density)
      else
         needs_ramp = .false.
      end if

      return
end subroutine eqbound_core

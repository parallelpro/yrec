!----------------------------------------------------------------------
! eqbound
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eqbound.f; only variable names, source form, and comment
! style were updated.
!
! MHP 8/98: checks whether a (T, log10 rho) point falls within the
! OPAL 1995 EOS table, and if it is near the table's edge, computes
! the ramp factor used to blend the OPAL result with the Yale/SCV
! result (see eqstat2.f90's use_opal95_eos branch, which calls this
! immediately after oeqos).
subroutine eqbound(temperature, log10_density, ramp_factor, &
     in_opal_table, needs_ramp)

      implicit none

      integer, parameter :: mx = 5, mv = 10, nr = 77, nt = 56

      double precision, intent(in) :: temperature, log10_density
      double precision, intent(out) :: ramp_factor
      logical, intent(out) :: in_opal_table, needs_ramp

! common/a/: the main OPAL 1995 EOS table and its interpolation
! grids/scratch arrays. Only density_grid, t6_grid,
! density_index_edge_at_t (via common/rmpopeos/, see below) are used
! here; the rest are placeholders. Naming matches readco.f90/esac.f90.
      double precision :: eos_table(mx,mv,nt,nr), t6_list(nr,nt), &
           density_grid(nr), t6_grid(nt), x_interp_result(nt,nr), &
           x_interp_result_alt(nt,nr), x_grid_spacing_inv(mx), &
           t6_grid_spacing_inv(nt), density_grid_spacing_inv(nr), &
           x_grid(mx)
      integer :: x_loop_index, x_index_lo
      common/a/ eos_table, t6_list, density_grid, t6_grid, &
           x_interp_result, x_interp_result_alt, x_grid_spacing_inv, &
           t6_grid_spacing_inv, density_grid_spacing_inv, x_grid, &
           x_loop_index, x_index_lo

! common/rmpopeos/: edge-of-table ramp data. Naming matches readco.f90.
      double precision :: density_edge_at_t(nt)
      integer :: density_index_edge_at_t(nt), t_row_index
      common/rmpopeos/ density_edge_at_t, density_index_edge_at_t, &
           t_row_index

      save

! --- locals ---
      double precision :: t6, table_edge_density, ramp_start_density
      double precision :: t6_top_of_table, t_fraction
      integer :: t6_scan_idx

      t6 = temperature*1.0d-6
!     exit if outside table in rho
      if (log10_density.lt.-14d0 .or. log10_density.gt.5.0d0) then
         goto 9999        ! Out of Table in density. Go to Error exit
      end if
!     find nearest table element in t.
      if (t6.lt.t6_grid(t_row_index)) then
         do t6_scan_idx = t_row_index+1, nt
            if (t6.ge.t6_grid(t6_scan_idx)) then
               t_row_index = t6_scan_idx - 1
               goto 10
            end if
         end do
         t6_top_of_table = t6_grid(nt)
!        sr call should have been stopped outside table bounds; stop code
         write(*,5) t6, t6_top_of_table, t_row_index
    5    format(' ERROR IN OPAL EOS: OUTSIDE TABLE IN T6',2F10.6,I5)
         stop
   10    continue
      else
         do t6_scan_idx = t_row_index, 1, -1
            if (t6.le.t6_grid(t6_scan_idx)) then
               t_row_index = t6_scan_idx
               goto 20
            end if
         end do
         t6_top_of_table = t6_grid(1)
!        sr call should have been stopped outside table bounds; stop code
         write(*,5) t6, t6_top_of_table, t_row_index
         stop
   20    continue
      end if
      t_fraction = (t6 - t6_grid(t_row_index+1))/ &
           (t6_grid(t_row_index)-t6_grid(t_row_index+1))
!     define table edge in rho by linear interpolation.
      table_edge_density = density_edge_at_t(t_row_index+1)
      table_edge_density = log10(table_edge_density)
!     define beginning of ramp in the same fashion -
!     ramp is defined as one table element wide.
      ramp_start_density = density_grid(density_index_edge_at_t(t_row_index+1)-1)
      ramp_start_density = log10(ramp_start_density)
!     check if within table bounds in rho
      if (log10_density.gt.table_edge_density) then
         goto 9999     ! Out of table in density. Go to error exit
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
!        If we get here, we are in the middle of the table and
!        no ramping is needed.
         needs_ramp = .false.
      end if

      return

!     Error Exit.
 9999 continue
      in_opal_table = .false.       ! Not in table
      needs_ramp = .true.           ! Turn on ramping
      ramp_factor = 0d0             ! Set ramping factor to zero
!                            This way, out of table results are ramped to zero.
      return
end subroutine eqbound

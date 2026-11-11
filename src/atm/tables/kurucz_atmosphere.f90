!----------------------------------------------------------------------
! kurucz_atmosphere
!----------------------------------------------------------------------
! The surface boundary condition from a surface-pressure table, Kurucz
! (1993) or Castelli & Kurucz (2003): log P at T = Teff, at (log Teff,
! log g).
!     interpolate_surface_pressure  the lookup: the table in, log P and a
!                                   status out (a star-blind kernel; it
!                                   prints nothing)
! kurucz_table_reader reads the tables, atm_lib keeps them and reports a
! failed lookup (atm_messages).
module kurucz_atmosphere
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use atm_tables, only: surface_pressure_table, surface_pressure_num_g
      use numerics_lib, only: inter3, natural_spline, splint
      implicit none
      private
      public :: interpolate_surface_pressure
      public :: surface_pressure_outside

! The ierr of a point the table cannot give (a spline's failure returns
! its numerics_status instead).
      integer, parameter :: surface_pressure_outside = 1

! The tables' fixed lower edges, below which a point is refused (about
! one row and one column below the first), and the log Teff above which
! the tables hold only their top three log g columns (3-point Lagrange in
! log g) or top two (linear).
      real(dp), parameter :: table_log_teff_min = 3.5d0, table_log_g_min = -0.5d0
      real(dp), parameter :: log_teff_three_columns = 4.5d0, log_teff_two_columns = 4.55d0

! The table's last column, and the spline stencils' size.
      integer, parameter :: ng = surface_pressure_num_g
      integer, parameter :: num_nodes = 4

contains

!----------------------------------------------------------------------
! interpolate_surface_pressure
!----------------------------------------------------------------------
! log P at (log10_teff, log10_gravity): in each of the four log Teff rows
! about the point, a natural cubic spline in log g through the row's
! four valid columns about it (above log Teff 4.5, 3-point Lagrange or
! linear in the top columns); then a natural cubic spline in log Teff
! through the four rows' values. Outside the table, ierr =
! surface_pressure_outside: below log Teff 3.5 or log g -0.5, more than
! one table step above the last row or column, or a stencil that reaches
! a row without the columns it needs (fewer than four valid ones, or not
! the top ones above log Teff 4.5). A spline's failure returns its
! status in ierr and spline_status. log10_pressure is set only on
! success.
pure subroutine interpolate_surface_pressure(tbl, log10_teff, log10_gravity, &
     log10_pressure, ierr, spline_status)
      type(surface_pressure_table), intent(in) :: tbl
      real(dp), intent(in) :: log10_teff, log10_gravity
      real(dp), intent(inout) :: log10_pressure
      integer, intent(out) :: ierr, spline_status

      real(dp) :: teff_nodes(num_nodes), pressure_at_nodes(num_nodes), &
           teff_spline_deriv(num_nodes), interpolated_value
      integer :: row, first_row, node

      ierr = 0
      spline_status = 0
      if (outside_table(tbl, log10_teff, log10_gravity)) then
         ierr = surface_pressure_outside
         return
      end if

! the four rows about log Teff
      do row = 1, tbl%num_teff
         if (log10_teff <= tbl%log_teff(row)) exit
      end do
      if (row > tbl%num_teff) row = tbl%num_teff
      first_row = max(1, row-2)
      first_row = min(tbl%num_teff-3, first_row)
      teff_nodes = tbl%log_teff(first_row:first_row+3)

! log P at log g in each
      do node = 1, num_nodes
         call pressure_in_row(tbl, first_row+node-1, log10_gravity, &
              pressure_at_nodes(node), ierr, spline_status)
         if (ierr /= 0) return
      end do

! then in log Teff
      call natural_spline(teff_nodes, pressure_at_nodes, num_nodes, teff_spline_deriv)
      call splint(teff_nodes, pressure_at_nodes, num_nodes, teff_spline_deriv, &
           log10_teff, interpolated_value, spline_status)
      if (spline_status /= 0) then
         ierr = spline_status
         return
      end if
      log10_pressure = interpolated_value
end subroutine interpolate_surface_pressure

! The point is outside what the table can give: below its fixed lower
! edges, or more than one table step above its last row or column.
pure logical function outside_table(tbl, log10_teff, log10_gravity)
      type(surface_pressure_table), intent(in) :: tbl
      real(dp), intent(in) :: log10_teff, log10_gravity
      real(dp) :: log_teff_max, log_g_max
      integer :: n

      n = tbl%num_teff
      log_teff_max = tbl%log_teff(n) + (tbl%log_teff(n) - tbl%log_teff(n-1))
      log_g_max = tbl%log_g(ng) + (tbl%log_g(ng) - tbl%log_g(ng-1))
      outside_table = log10_teff < table_log_teff_min .or. log10_gravity < table_log_g_min .or. &
           log10_teff > log_teff_max .or. log10_gravity > log_g_max
end function outside_table

! log P at log10_gravity in one row. Above log Teff 4.55, linear in the
! top two columns; above 4.5, 3-point Lagrange in the top three; below,
! a natural cubic spline through the row's four valid columns about the
! point: its top four from the second-to-last valid column up, its first
! four below its second valid column. A row without the columns its
! branch needs refuses the point (the Castelli/Kurucz tables' last row
! holds only the top column; an all -999 row has first_column = ng and
! last_column = 0).
pure subroutine pressure_in_row(tbl, row, log10_gravity, value, ierr, spline_status)
      type(surface_pressure_table), intent(in) :: tbl
      integer, intent(in) :: row
      real(dp), intent(in) :: log10_gravity
      real(dp), intent(out) :: value
      integer, intent(inout) :: ierr, spline_status

      real(dp) :: weights(3), unused_dweights(3), fraction
      integer :: first_valid, last_valid, first, k

      value = 0.0d0
      first_valid = tbl%first_column(row)
      last_valid = tbl%last_column(row)
      if (tbl%log_teff(row) > log_teff_three_columns) then
         if (tbl%log_teff(row) > log_teff_two_columns) then
            if (first_valid > ng-1 .or. last_valid < ng) then
               ierr = surface_pressure_outside
               return
            end if
            fraction = (log10_gravity - tbl%log_g(ng-1))/(tbl%log_g(ng) - tbl%log_g(ng-1))
            value = tbl%log10_pressure(row,ng-1) + &
                 fraction*(tbl%log10_pressure(row,ng) - tbl%log10_pressure(row,ng-1))
         else
            if (first_valid > ng-2 .or. last_valid < ng) then
               ierr = surface_pressure_outside
               return
            end if
            call inter3(tbl%log_g(ng-2:ng), weights, unused_dweights, log10_gravity)
            value = tbl%log10_pressure(row,ng-2)*weights(1) + &
                 tbl%log10_pressure(row,ng-1)*weights(2) + &
                 tbl%log10_pressure(row,ng)*weights(3)
         end if
         return
      end if

      if (last_valid - first_valid + 1 < num_nodes) then
         ierr = surface_pressure_outside
         return
      end if
      if (log10_gravity >= tbl%log_g(last_valid-1)) then
         first = last_valid - 3
      else
! the interval [log_g(k+1), log_g(k+2)) that holds the point, and the
! stencil k..k+3 about it; below the second valid column, the first four
         first = first_valid
         do k = last_valid-3, first_valid, -1
            if (log10_gravity < tbl%log_g(k+2) .and. log10_gravity >= tbl%log_g(k+1)) then
               first = k
               exit
            end if
         end do
      end if
      call spline_in_log_g(tbl, row, first, log10_gravity, value, spline_status)
      if (spline_status /= 0) ierr = spline_status
end subroutine pressure_in_row

! The natural cubic spline through columns first..first+3 of a row, at
! log10_gravity.
pure subroutine spline_in_log_g(tbl, row, first, log10_gravity, value, spline_status)
      type(surface_pressure_table), intent(in) :: tbl
      integer, intent(in) :: row, first
      real(dp), intent(in) :: log10_gravity
      real(dp), intent(out) :: value
      integer, intent(out) :: spline_status

      real(dp) :: nodes(num_nodes), values(num_nodes), deriv(num_nodes)

      nodes = tbl%log_g(first:first+3)
      values = tbl%log10_pressure(row,first:first+3)
      call natural_spline(nodes, values, num_nodes, deriv)
      call splint(nodes, values, num_nodes, deriv, log10_gravity, value, spline_status)
end subroutine spline_in_log_g

!======================================================================
! Sources: the F77's SURFP and KCSURFP (one routine since 2026). G.
! Somers (5/15) bounded the gravity search by each row's valid columns,
! so that rows ending in -999 entries are not interpolated through them.
! The checks of the top columns, the upper range and the row's own top
! four are docs/bug_sweep.md's M231, M232 and M234 (2026-11).

end module kurucz_atmosphere

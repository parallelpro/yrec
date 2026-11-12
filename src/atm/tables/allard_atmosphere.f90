!----------------------------------------------------------------------
! allard_atmosphere
!----------------------------------------------------------------------
! The surface boundary condition from the Allard model atmospheres
! (NextGen, AMES-Cond, AMES-Dusty) at one [Fe/H] and [alpha/Fe]: log P
! at T = Teff, or log P and log T at tau = 100, at (log Teff, log g).
!     interpolate_allard_surface  the lookup: the tables in, an
!                                 allard_lookup_result out (a star-blind
!                                 kernel; it prints nothing)
!     allard_lookup_result        the result, or why the point is outside
! allard_table_reader reads the tables, atm_lib keeps them and reports a
! lookup outside them (atm_messages).
module allard_atmosphere
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use atm_tables, only: allard_table
      use numerics_lib, only: locate, polint
      implicit none
      private
      public :: interpolate_allard_surface, allard_lookup_result

! A lookup's outcome: the result, or what put the point outside the
! tables. teff_above and gravity_outside (either or both) and
! row_outside are a point the caller takes elsewhere (the star layer
! falls back to the gray atmosphere); teff_below and a polint failure
! are errors.
      type :: allard_lookup_result
! log P (at T = Teff, or at tau = 100) and log T (log Teff, or at
! tau = 100); set only when the point is inside
         real(dp) :: log10_pressure = 0.0d0, log10_temperature = 0.0d0
! log Teff at or above the tables' log_teff_max; log g at or beyond
! their log_g_min or log_g_max
         logical :: teff_above = .false., gravity_outside = .false.
! log Teff at or below log_teff_min
         logical :: teff_below = .false.
! the log Teff row of the stencil whose log g range does not hold log g
! (0: none)
         integer :: row_outside = 0
! a polint's failure (numerics_status; 0: none)
         integer :: polint_status = 0
      end type allard_lookup_result

! the Lagrange stencils' size, in log g and in log Teff
      integer, parameter :: num_nodes = 4

contains

!----------------------------------------------------------------------
! interpolate_allard_surface
!----------------------------------------------------------------------
! The tables at (log_teff, log_g), in four steps:
!   1. the tables' range (find_table_ranges): log Teff below
!      log_teff_max and above log_teff_min, one row beyond the first and
!      the last row; log g strictly inside (log_g_min, log_g_max), the
!      widest of the rows' ranges, each of which reaches four column
!      widths below its row's first valid log g and one above its last
!      (M233);
!   2. the stencil: the four log Teff rows about log Teff (clamped to the
!      table) and, in each, the four valid log g columns about log g
!      (clamped to the row's valid ones); log g must lie in each row's
!      own range (row_log_g_min, row_log_g_max);
!   3. 4-point Lagrange interpolation (polint) in log g along each of the
!      four rows, for log P, and log P and log T at tau = 100;
!   4. 4-point Lagrange interpolation in log Teff of the rows' values:
!      log P at T = Teff, with log T = log Teff, or with use_tau100 log P
!      and log T at tau = 100.
pure subroutine interpolate_allard_surface(tbl, use_tau100, log_teff, log_g, lookup)
      type(allard_table), intent(in) :: tbl
      logical, intent(in) :: use_tau100
      real(dp), intent(in) :: log_teff, log_g
      type(allard_lookup_result), intent(out) :: lookup

      integer :: first_row, first_column(num_nodes)
      real(dp) :: pressure_col(num_nodes), pressure_tau100_col(num_nodes), &
           temperature_tau100_col(num_nodes)

! 1. the tables' range
      lookup%teff_above = log_teff >= tbl%log_teff_max
      lookup%gravity_outside = log_g >= tbl%log_g_max .or. log_g <= tbl%log_g_min
      if (lookup%teff_above .or. lookup%gravity_outside) return
      if (log_teff <= tbl%log_teff_min) then
         lookup%teff_below = .true.
         return
      end if
! 2. the stencil
      call find_stencil(tbl, log_teff, log_g, first_row, first_column, lookup%row_outside)
      if (lookup%row_outside > 0) return
! 3. along each row, in log g
      call interpolate_rows(tbl, first_row, first_column, log_g, pressure_col, &
           pressure_tau100_col, temperature_tau100_col, lookup%polint_status)
      if (lookup%polint_status /= 0) return
! 4. across the rows, in log Teff
      call interpolate_across_rows(tbl, first_row, use_tau100, log_teff, pressure_col, &
           pressure_tau100_col, temperature_tau100_col, lookup)
end subroutine interpolate_allard_surface

! The first of the four log Teff rows about log_teff (locate gives the
! row below the point: one more below it, clamped to the table), and in
! each row the first of the four valid log g columns about log_g
! (likewise, clamped to the row's valid columns). row_outside: the first
! row whose log g range (row_log_g_min, row_log_g_max) does not hold
! log_g, else 0.
pure subroutine find_stencil(tbl, log_teff, log_g, first_row, first_column, row_outside)
      type(allard_table), intent(in) :: tbl
      real(dp), intent(in) :: log_teff, log_g
      integer, intent(out) :: first_row, first_column(:), row_outside

      integer :: i, row, first_valid, last_valid, below

      row_outside = 0
      first_column = 0
      call locate(tbl%log_teff(1:tbl%num_teff), tbl%num_teff, log_teff, below)
      first_row = below - 1
      if (first_row < 1) first_row = 1
      if (first_row > tbl%num_teff - 3) first_row = tbl%num_teff - 3

      do i = 1, num_nodes
         row = first_row + i - 1
         first_valid = tbl%first_column(row)
         last_valid = tbl%last_column(row)
         call locate(tbl%log_g(first_valid:last_valid), last_valid - first_valid + 1, &
              log_g, below)
         first_column(i) = first_valid + below - 2
         if (first_column(i) < first_valid) first_column(i) = first_valid
         if (first_column(i) > last_valid - 3) first_column(i) = last_valid - 3
         if (log_g < tbl%row_log_g_min(row) .or. log_g > tbl%row_log_g_max(row)) then
            row_outside = row
            return
         end if
      end do
end subroutine find_stencil

! Along each of the four rows, 4-point Lagrange interpolation in log g
! of log P and of log P and log T at tau = 100. status: the first
! polint's failure, else 0.
pure subroutine interpolate_rows(tbl, first_row, first_column, log_g, pressure_col, &
     pressure_tau100_col, temperature_tau100_col, status)
      type(allard_table), intent(in) :: tbl
      integer, intent(in) :: first_row, first_column(:)
      real(dp), intent(in) :: log_g
      real(dp), intent(out) :: pressure_col(:), pressure_tau100_col(:), &
           temperature_tau100_col(:)
      integer, intent(out) :: status

      integer :: i, row, first
      real(dp) :: log_g_nodes(num_nodes), unused_dy

      pressure_col = 0.0d0
      pressure_tau100_col = 0.0d0
      temperature_tau100_col = 0.0d0
      do i = 1, num_nodes
         row = first_row + i - 1
         first = first_column(i)
         log_g_nodes = tbl%log_g(first:first+3)
         call polint(log_g_nodes, tbl%log10_pressure(row,first:first+3), num_nodes, log_g, &
              pressure_col(i), unused_dy, status)
         if (status /= 0) return
         call polint(log_g_nodes, tbl%log10_pressure_tau100(row,first:first+3), num_nodes, log_g, &
              pressure_tau100_col(i), unused_dy, status)
         if (status /= 0) return
         call polint(log_g_nodes, tbl%log10_temperature_tau100(row,first:first+3), num_nodes, &
              log_g, temperature_tau100_col(i), unused_dy, status)
         if (status /= 0) return
      end do
end subroutine interpolate_rows

! 4-point Lagrange interpolation in log Teff of the rows' log P (and log
! T = log Teff), or with use_tau100 of their log P and log T at
! tau = 100, into lookup.
pure subroutine interpolate_across_rows(tbl, first_row, use_tau100, log_teff, pressure_col, &
     pressure_tau100_col, temperature_tau100_col, lookup)
      type(allard_table), intent(in) :: tbl
      integer, intent(in) :: first_row
      logical, intent(in) :: use_tau100
      real(dp), intent(in) :: log_teff
      real(dp), intent(in) :: pressure_col(:), pressure_tau100_col(:), &
           temperature_tau100_col(:)
      type(allard_lookup_result), intent(inout) :: lookup

      real(dp) :: log_teff_nodes(num_nodes), unused_dy

      log_teff_nodes = tbl%log_teff(first_row:first_row+3)
      if (.not. use_tau100) then
         call polint(log_teff_nodes, pressure_col, num_nodes, log_teff, lookup%log10_pressure, &
              unused_dy, lookup%polint_status)
         if (lookup%polint_status /= 0) return
         lookup%log10_temperature = log_teff
      else
         call polint(log_teff_nodes, pressure_tau100_col, num_nodes, log_teff, &
              lookup%log10_pressure, unused_dy, lookup%polint_status)
         if (lookup%polint_status /= 0) return
         call polint(log_teff_nodes, temperature_tau100_col, num_nodes, log_teff, &
              lookup%log10_temperature, unused_dy, lookup%polint_status)
      end if
end subroutine interpolate_across_rows

!======================================================================
! Sources: the F77's ALSURFP (the Allard atmospheres, MHP 6/97; the
! 2007 NextGen layout, 2008), with its four steps; the model
! atmospheres of F. Allard, P. H. Hauschildt and collaborators.

end module allard_atmosphere

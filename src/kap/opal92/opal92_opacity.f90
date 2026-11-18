!----------------------------------------------------------------------
! opal92_opacity
!----------------------------------------------------------------------
! The atomic opacities of the 1992 OPAL (Livermore) tables, one file per
! Z with tables in X, log T and log R; a second file at another Z when
! the run blends two Z tables.
!     prepare_opal92_tables       each (X, T) row's spline in log R, from
!                                 the table(s) opal92_reader read
!     build_opal92_surface_table  the fixed-X table at the envelope
!                                 composition, per Z table
!     compute_opal92_opacity      the opacity and its log derivatives by
!                                 interpolation in X, T and R
! The tables, their splines, the surface tables and the lookup's cursors
! are opal92_table_sets (opacity_tables), which kap_lib holds and passes
! in. The sources are at the end of the file.
module opal92_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: n_opal92_d, n_opal92_t, n_opal92_x, &
           opal92_missing_opacity, opal92_table_set, opal92_x_match_tol, &
           kap_settings_t, kap_table_paths_t, opacity_table_state
      use opacity_source, only: kap_point, kap_result, kap_source
      use opacity_blending, only: blend_results_in_z
      use numerics_lib, only: findex, natural_spline, splint, &
           compute_spline_coefficients, numerics_status_text
      use opal92_reader, only: read_opal92_tables
      use math_lib
      implicit none
      private
      public :: opal92_source
      public :: prepare_opal92_tables, build_opal92_surface_table, compute_opal92_opacity

! The X table a lookup reads is one of the file's n_opal92_x tables or,
! by this index, the surface table built at the envelope X
      integer, parameter :: opal92_surface_slot = n_opal92_x + 1
! Below the second X table a point is interpolated with the first one
! only up to this log R; above it, it is extrapolated in X from the
! second and third tables (the F77 code's rule)
      real(dp), parameter :: first_x_table_logr_max = -1.0d0
! the spline coefficient work array's interval count, more than any
! row's log R columns
      integer, parameter :: max_spline_intervals = 100

contains


!----------------------------------------------------------------------
! prepare_opal92_tables
!----------------------------------------------------------------------
! Builds the natural cubic splines in log R of every (X, T) row of the
! table(s) (spline_opal92_rows), as coefficients per interval: the
! second table's with two_z_tables. A row that does not start at the
! first log R returns ierr with the reason in message.
pure subroutine prepare_opal92_tables(tables, two_z_tables, message, ierr)
      type(opal92_table_set), intent(inout) :: tables(:)
      logical, intent(in) :: two_z_tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      call spline_opal92_rows(tables(1), message, ierr)
      if (ierr /= 0) return
      if (two_z_tables) call spline_opal92_rows(tables(2), message, ierr)
end subroutine prepare_opal92_tables

! One table set's splines in log R: for each (X, T) row, the cells with
! data (above opal92_missing_opacity), which must start at the first
! log R; their count and the row's natural cubic spline, as four
! coefficients per interval (tbl%spline_coeffs). A row is addressed as
! T row + (X index - 1) x the number of T rows: by the rows read for
! the spline's arrays and by n_opal92_t in log10_opacity, where the
! reader stored it (the two are equal: a short table is refused).
pure subroutine spline_opal92_rows(tbl, message, ierr)
      type(opal92_table_set), intent(inout) :: tbl
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr
      real(dp) :: spline_work(4,max_spline_intervals), logr_nodes(n_opal92_d), &
           row_opacity(n_opal92_d), row_second_derivs(n_opal92_d)
      integer :: ix, t_row, row, column, num_filled, interval, i
      real(dp) :: log10_opacity

      ierr = 0
      message = ''
      do ix = 1,tbl%num_x
         do t_row = 1,tbl%num_temps
            row = t_row + (ix-1)*tbl%num_temps
            num_filled = 0
            do column = 1,n_opal92_d
               log10_opacity = tbl%log10_opacity(t_row+n_opal92_t*(ix-1),column)
! (an empty cell)
               if (log10_opacity <= opal92_missing_opacity) cycle
               if (num_filled <= 0) then
                  tbl%density_start_index(row) = column
                  if (column /= 1) then
                     write(message,'(a,f6.3,a)') 'the OPAL 1992 table''s row at log T6', &
                          tbl%grid_logt(t_row), ' does not start at its first log R'
                     ierr = 1
                     return
                  end if
               end if
               num_filled = num_filled + 1
               logr_nodes(num_filled) = tbl%grid_logr(column)
               row_opacity(num_filled) = log10_opacity
            end do
            tbl%density_count(row) = num_filled
            if (num_filled <= 1) cycle
! the row's natural cubic spline in log R, as coefficients per interval
            call natural_spline(logr_nodes, row_opacity, num_filled, row_second_derivs)
            call compute_spline_coefficients(logr_nodes, row_opacity, &
                 row_second_derivs, num_filled, spline_work)
            do interval = 1,num_filled
               do i = 1,4
                  tbl%spline_coeffs(row,i + (interval-1)*4) = spline_work(i,interval)
               end do
            end do
         end do
      end do
end subroutine spline_opal92_rows

!----------------------------------------------------------------------
! build_opal92_surface_table
!----------------------------------------------------------------------
! Builds each table set's surface table (its "4th X table", at the
! envelope X, so that lookups there need no X interpolation) with
! build_opal92_surface_rows: the second set's with two_z_tables. The
! density cursor carries from the first set into the second. A row that
! does not start at the first log R returns ierr with the reason in
! message.
pure subroutine build_opal92_surface_table(tables, hydrogen_fraction, table_z, &
     two_z_tables, message, ierr)
! (1) at table_z(1) and, with two_z_tables, (2) at table_z(2)
! (opal_table_z1 and opal_table_z2)
      type(opal92_table_set), intent(inout) :: tables(:)
      real(dp), intent(in) :: hydrogen_fraction, table_z(:)
      logical, intent(in) :: two_z_tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr
! findex's starting guesses: none (0)
      integer :: x_index, logr_cursor

      ierr = 0
      message = ''
      x_index = 0
      logr_cursor = 0
      call build_opal92_surface_rows(tables(1), hydrogen_fraction, table_z(1), x_index, &
           logr_cursor, message, ierr)
      if (ierr /= 0) return
      if (two_z_tables) then
         call build_opal92_surface_rows(tables(2), hydrogen_fraction, table_z(2), x_index, &
              logr_cursor, message, ierr)
      end if
end subroutine build_opal92_surface_table

! One table set's surface table at hydrogen_fraction (and Z = table_z,
! recorded): log kappa at every (T, log R) node interpolated linearly in
! X between the two X tables around the surface X (the last two above
! the grid, the first two below it), then splined in log R (only
! log kappa is kept; interpolate_opal92_2d's derivatives are not used
! here). x_index and logr_cursor are findex's cursors, the caller's.
pure subroutine build_opal92_surface_rows(tbl, hydrogen_fraction, table_z, x_index, logr_cursor, &
     message, ierr)
      type(opal92_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: hydrogen_fraction, table_z
      integer, intent(inout) :: x_index, logr_cursor
      character(len=*), intent(out) :: message
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success
      real(dp) :: spline_work(4,max_spline_intervals), row_opacity(n_opal92_d), &
           row_second_derivs(n_opal92_d)
      integer :: t_row, row, column, first_column, last_column, interval, i
      real(dp) :: log_t6, log_r, x_fraction
! the lookups in the two X tables (kappa and the derivatives are not used)
      real(dp) :: opacity, log10_opacity_lo, log10_opacity_hi, dlnkap_dlnrho, dlnkap_dlnt

      message = ''
!     the surface table's composition
      tbl%surface_x=hydrogen_fraction
      tbl%surface_z=table_z
      call findex(tbl%grid_x, n_opal92_x, hydrogen_fraction, x_index)
      if (x_index < 0) x_index=-x_index
      if (x_index >= 3) x_index=2
      if (x_index <= 0) then
         write(message,'(a,f8.5,a)') 'the surface X =', hydrogen_fraction, &
              ' is outside the OPAL 1992 tables'' X grid'
         ierr = 1
         return
      end if
      tbl%surface_x_index=x_index
      do t_row=1,tbl%num_temps
         x_index=tbl%surface_x_index
         log_t6=tbl%grid_logt(t_row)
         row=t_row+(x_index-1)*tbl%num_temps
         first_column=tbl%density_start_index(row)
         if (first_column /= 1) then
            write(message,'(a,f6.3,a)') 'the OPAL 1992 table''s row at log T6', log_t6, &
                 ' does not start at its first log R'
            ierr = 1
            return
         end if
         last_column=first_column+tbl%density_count(row)-1
         do column=first_column,last_column
            log_r=tbl%grid_logr(column)
            call interpolate_opal92_2d(tbl,log_t6,log_r,x_index,t_row,logr_cursor, &
                 opacity,log10_opacity_lo,dlnkap_dlnrho,dlnkap_dlnt, message, ierr)
            if (ierr /= 0) return
            call interpolate_opal92_2d(tbl,log_t6,log_r,x_index+1,t_row,logr_cursor, &
                 opacity,log10_opacity_hi,dlnkap_dlnrho,dlnkap_dlnt, message, ierr)
            if (ierr /= 0) return
            x_fraction=(hydrogen_fraction-tbl%grid_x(x_index))/ &
                 (tbl%grid_x(x_index+1)-tbl%grid_x(x_index))
            row_opacity(column)=(log10_opacity_hi-log10_opacity_lo)*x_fraction+log10_opacity_lo
         end do
         call natural_spline(tbl%grid_logr, row_opacity, last_column, row_second_derivs)
         call compute_spline_coefficients(tbl%grid_logr, row_opacity, &
              row_second_derivs, last_column, spline_work)
         do interval=1,last_column
            do i=1,4
               tbl%surface_spline_coeffs(t_row,i+(interval-1)*4)=spline_work(i,interval)
            end do
         end do
      end do
end subroutine build_opal92_surface_rows

!----------------------------------------------------------------------
! compute_opal92_opacity
!----------------------------------------------------------------------
! The opacity at (log rho, log T, X) and its log derivatives from table
! set tbl (the first or the second Z's): in the table's coordinates
! (log T6 and log R = log rho - 3 log T6), by cubic splines in log R and
! log T (interpolate_opal92_2d) in one X table where X matches it (the
! surface table, or a table of the grid, within opal92_x_match_tol),
! else in the two X tables around X, interpolated (or extrapolated)
! linearly in X. A point off the tables returns ierr = 1 with the reason
! in message. tbl%index_x, index_t and index_rho are the warm-start
! cursors.
pure subroutine compute_opal92_opacity(tbl, log10_density, log10_temperature, hydrogen_fraction, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, message, ierr)
      type(opal92_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      logical :: single_x_table
      real(dp) :: log_t6, log_r, x_fraction
      integer :: ix
! the lookups in the X table (and the next one): kappa, log10 kappa and
! its derivatives at constant log R
      real(dp) :: opacity_lo, log10_opacity_lo, dlnkap_dlnr_lo, dlnkap_dlnt_lo
      real(dp) :: opacity_hi, log10_opacity_hi, dlnkap_dlnr_hi, dlnkap_dlnt_hi
      real(dp) :: dlnkap_dlnr, dlnkap_dlnt_at_r

      ierr = 0
      message = ''
      single_x_table = .false.
! the table's coordinates: log T6 = log10(T/1e6 K) and log R =
! log10(rho/T6**3)
      log_t6 = log10_temperature - 6.0d0
      log_r = log10_density - 3.0d0*log_t6

      if (abs(tbl%surface_x-hydrogen_fraction) <= opal92_x_match_tol) then
         tbl%index_x = opal92_surface_slot
         single_x_table = .true.
      end if
      if (.not. single_x_table) then
         do ix = 1,n_opal92_x
            if (abs(tbl%grid_x(ix)-hydrogen_fraction) <= opal92_x_match_tol) then
               tbl%index_x = ix
               single_x_table = .true.
               exit
            end if
         end do
      end if
      if (.not. single_x_table) then
         call findex(tbl%grid_x, n_opal92_x, hydrogen_fraction, tbl%index_x)
         if (tbl%index_x < 0) tbl%index_x = -tbl%index_x
         if (tbl%index_x <= 1.and.log_r > first_x_table_logr_max) tbl%index_x = 2
         if (tbl%index_x >= 3) tbl%index_x = 2
         if (tbl%index_x <= 0) then
            write(message,'(a,f8.5,a)') 'X =', hydrogen_fraction, &
                 ' is outside the OPAL 1992 tables'' X grid'
            ierr = 1
            return
         end if
      end if
      call findex(tbl%grid_logt, tbl%num_temps, log_t6, tbl%index_t)
      if (tbl%index_t < 0.and.tbl%grid_logt(tbl%num_temps) == log_t6) tbl%index_t = -tbl%index_t
      if (tbl%index_t < 0) then
         write(message,'(a,f9.6,a)') 'log T =', log10_temperature, &
              ' is outside the OPAL 1992 tables'
         ierr = 1
         return
      end if
      call interpolate_opal92_2d(tbl, log_t6, log_r, tbl%index_x, tbl%index_t, tbl%index_rho, &
           opacity_lo, log10_opacity_lo, dlnkap_dlnr_lo, dlnkap_dlnt_lo, message, ierr)
      if (ierr /= 0) return
      if (single_x_table) then
! Use only one X table
         log10_opacity = log10_opacity_lo
         opacity = opacity_lo
         dlnkap_dlnr = dlnkap_dlnr_lo
         dlnkap_dlnt_at_r = dlnkap_dlnt_lo
      else
! Linear interpolation (or extrapolation) in X
         call interpolate_opal92_2d(tbl, log_t6, log_r, tbl%index_x+1, tbl%index_t, tbl%index_rho, &
              opacity_hi, log10_opacity_hi, dlnkap_dlnr_hi, dlnkap_dlnt_hi, message, ierr)
         if (ierr /= 0) return
         x_fraction = (hydrogen_fraction-tbl%grid_x(tbl%index_x))/ &
              (tbl%grid_x(tbl%index_x+1)-tbl%grid_x(tbl%index_x))
         log10_opacity = (log10_opacity_hi-log10_opacity_lo)*x_fraction + log10_opacity_lo
         dlnkap_dlnr = (dlnkap_dlnr_hi-dlnkap_dlnr_lo)*x_fraction + dlnkap_dlnr_lo
         dlnkap_dlnt_at_r = (dlnkap_dlnt_hi-dlnkap_dlnt_lo)*x_fraction + dlnkap_dlnt_lo
         opacity = exp10(log10_opacity)
      end if
! from the derivative at constant log R to that at constant rho
      dlnkap_dlnrho = dlnkap_dlnr
      dlnkap_dlnt = dlnkap_dlnt_at_r - 3.0d0*dlnkap_dlnrho

end subroutine compute_opal92_opacity

!----------------------------------------------------------------------
! interpolate_opal92_2d
!----------------------------------------------------------------------
! The opacity in one X table of tbl (x_table; opal92_surface_slot is the
! surface table) at (log T6, log R): on up to six T rows around t_index,
! log kappa and its log R derivative from the row's spline in log R
! (extrapolated linearly outside the row), then natural splines in
! log T6 through those, for log kappa, dlnkap_dlnt (at constant log R)
! and dlnkap_dlnrho. logr_cursor is the density cursor. A row that does
! not start at the first density, a log T6 outside the rows or a failed
! spline returns ierr with the reason in message.
pure subroutine interpolate_opal92_2d(tbl, log_t6, log_r, x_table, t_index, &
     logr_cursor, opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, message, ierr)
      type(opal92_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: log_t6, log_r
      integer, intent(in) :: x_table, t_index
      integer, intent(inout) :: logr_cursor
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

! the rows used: their log T6, and log kappa and its log R derivative
! at log R on each; the T splines' second derivatives
      real(dp) :: row_log_t6(n_opal92_t), row_log10_opacity(n_opal92_t)
      real(dp) :: row_dlnkap_dlnr(n_opal92_t), second_derivs(n_opal92_t)
      integer :: num_rows, t_row, first_row, last_row
      integer :: x_index, row, first_column, last_column, coeff_base
      real(dp) :: delta_logr, c1, c2, c3, c4

      ierr = 0
      message = ''
! For six grid points of temperature
      first_row = t_index - 2
      if (first_row <= 0) first_row = 1
      last_row = t_index + 3
      if (last_row > tbl%num_temps) last_row = tbl%num_temps
      if (x_table < opal92_surface_slot) then
         x_index = x_table
      else
         x_index = tbl%surface_x_index
      end if
      num_rows = 0
      do t_row = first_row,last_row
         row = t_row + (x_index-1)*tbl%num_temps
         first_column = tbl%density_start_index(row)
         if (first_column /= 1) then
            write(message,'(a,f6.3,a)') 'the OPAL 1992 table''s row at log T6', &
                 tbl%grid_logt(t_row), ' does not start at its first log R'
            ierr = 1
            return
         end if
         last_column = first_column + tbl%density_count(row) - 1
         call findex(tbl%grid_logr, last_column, log_r, logr_cursor)
         num_rows = num_rows + 1
         row_log_t6(num_rows) = tbl%grid_logt(t_row)
         if (logr_cursor < 0) then
! outside the row: linear extrapolation from its end
            logr_cursor = -logr_cursor
            coeff_base = 4*(logr_cursor - first_column)
            delta_logr = log_r - tbl%grid_logr(logr_cursor)
            if (x_table < opal92_surface_slot) then
               c1 = tbl%spline_coeffs(row,coeff_base+1)
               c2 = tbl%spline_coeffs(row,coeff_base+2)
            else
               c1 = tbl%surface_spline_coeffs(t_row,coeff_base+1)
               c2 = tbl%surface_spline_coeffs(t_row,coeff_base+2)
            end if
            row_log10_opacity(num_rows) = c2*delta_logr + c1
            row_dlnkap_dlnr(num_rows) = c2
         else
! inside: the spline's cubic on the interval that holds log R, and its
! derivative
            coeff_base = 4*(logr_cursor - first_column)
            delta_logr = log_r - tbl%grid_logr(logr_cursor)
            if (x_table < opal92_surface_slot) then
               c1 = tbl%spline_coeffs(row,coeff_base+1)
               c2 = tbl%spline_coeffs(row,coeff_base+2)
               c3 = tbl%spline_coeffs(row,coeff_base+3)
               c4 = tbl%spline_coeffs(row,coeff_base+4)
            else
               c1 = tbl%surface_spline_coeffs(t_row,coeff_base+1)
               c2 = tbl%surface_spline_coeffs(t_row,coeff_base+2)
               c3 = tbl%surface_spline_coeffs(t_row,coeff_base+3)
               c4 = tbl%surface_spline_coeffs(t_row,coeff_base+4)
            end if
            row_log10_opacity(num_rows) = ((c4*delta_logr+c3)*delta_logr+c2)*delta_logr + c1
            row_dlnkap_dlnr(num_rows) = (3.0d0*c4*delta_logr+2.0d0*c3)*delta_logr + c2
         end if
      end do
      if (row_log_t6(1) > log_t6.or.row_log_t6(num_rows) < log_t6) then
         write(message,'(a,f10.6,a)') 'log T6 =', log_t6, &
              ' is outside the OPAL 1992 table''s rows around it'
         ierr = 1
         return
      end if
! log kappa and its T derivative from the spline in log T6
      call natural_spline(row_log_t6, row_log10_opacity, num_rows, second_derivs)
      call splint(row_log_t6, row_log10_opacity, num_rows, second_derivs, log_t6, log10_opacity, &
           ierr, dy=dlnkap_dlnt)
      if (ierr /= 0) then
         message = 'the OPAL 1992 table''s spline in log T failed: ' // &
              trim(numerics_status_text(ierr))
         return
      end if
      opacity = exp10(log10_opacity)
! and the log R derivative, from its own spline in log T6
      call natural_spline(row_log_t6, row_dlnkap_dlnr, num_rows, second_derivs)
      call splint(row_log_t6, row_dlnkap_dlnr, num_rows, second_derivs, log_t6, dlnkap_dlnrho, ierr)
      if (ierr /= 0) then
         message = 'the OPAL 1992 table''s spline in log T failed: ' // &
              trim(numerics_status_text(ierr))
         return
      end if

end subroutine interpolate_opal92_2d

!----------------------------------------------------------------------
! The family's descriptor (opacity_source) and its procedures
!----------------------------------------------------------------------
! An interior family of one table per Z: with use_two_z_tables, a second
! table at opal_table_z2, blended linearly in Z with the first; in a
! helium-burning core, the first table alone, at opal_table_z1.
pure function opal92_source(settings) result(source)
      type(kap_settings_t), intent(in) :: settings
      type(kap_source) :: source

      source%name = 'OPAL 1992'
      source%blends_two_z = .true.
      source%table_z = [settings%opal_table_z1, settings%opal_table_z2]
      source%table_z_controls = [character(len=16) :: 'opal_table_z1', 'opal_table_z2']
      source%selected = settings%use_opal92_tables
      source%load => load_opal92
      source%refresh => refresh_opal92
      source%lookup => lookup_opal92
      source%core_lookup => core_lookup_opal92
end function opal92_source

subroutine load_opal92(settings, paths, envelope_hydrogen_fraction, tables, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(kap_table_paths_t), intent(in) :: paths
      real(dp), intent(in) :: envelope_hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

! (the reader reports its own errors)
      message = ''
      call read_opal92_tables(tables%opal92, paths%opal92, paths%opal92_z2, &
           settings%use_two_z_tables, [settings%opal_table_z1, settings%opal_table_z2], ierr)
      if (ierr /= 0) return
      call prepare_opal92_tables(tables%opal92, settings%use_two_z_tables, message, ierr)
      if (ierr /= 0) return
      call build_opal92_surface_table(tables%opal92, envelope_hydrogen_fraction, &
           [settings%opal_table_z1, settings%opal_table_z2], &
           settings%use_two_z_tables, message, ierr)
end subroutine load_opal92

pure subroutine refresh_opal92(settings, hydrogen_fraction, tables, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      real(dp), intent(in) :: hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      call build_opal92_surface_table(tables%opal92, hydrogen_fraction, &
           [settings%opal_table_z1, settings%opal_table_z2], &
           settings%use_two_z_tables, message, ierr)
end subroutine refresh_opal92

pure subroutine lookup_opal92(settings, tables, point, res, in_table, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: in_table
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr
      type(kap_result) :: res_z2

      in_table = .true.
      call compute_opal92_opacity(tables%opal92(1), point%log10_density, &
           point%log10_temperature, point%hydrogen_fraction, res%opacity, &
           res%log10_opacity, res%dlnkap_dlnrho, res%dlnkap_dlnt, message, ierr)
      if (ierr /= 0) return
      if (settings%use_two_z_tables) then
         call compute_opal92_opacity(tables%opal92(2), point%log10_density, &
              point%log10_temperature, point%hydrogen_fraction, res_z2%opacity, &
              res_z2%log10_opacity, res_z2%dlnkap_dlnrho, res_z2%dlnkap_dlnt, message, ierr)
         if (ierr /= 0) return
         call blend_results_in_z(res, res_z2, point%metal_fraction, &
              settings%opal_table_z1, settings%opal_table_z2)
      end if
end subroutine lookup_opal92

pure subroutine core_lookup_opal92(settings, tables, point, res, table_z, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      real(dp), intent(out) :: table_z
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      call compute_opal92_opacity(tables%opal92(1), point%log10_density, &
           point%log10_temperature, point%hydrogen_fraction, res%opacity, &
           res%log10_opacity, res%dlnkap_dlnrho, res%dlnkap_dlnt, message, ierr)
      table_z = settings%opal_table_z1
end subroutine core_lookup_opal92

!----------------------------------------------------------------------
! Sources
!----------------------------------------------------------------------
! The tables: the 1992 OPAL opacities (Rogers & Iglesias 1992, ApJS 79,
! 507), one file per Z in input/opacity/atomic/legacy/opal92/. The code:
! the F77 SETLLO, LL4TH, YLLO2D and their companions, with DBG's second
! Z table (5/94) and MHP's array dimensions (10/02), converted to Fortran
! 90 in 2026: the two tables' copies one routine on a table set passed
! in, the opens with status 'old' and a short table an error (M220), the
! findex guesses defined (M104).

end module opal92_opacity

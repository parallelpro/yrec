!----------------------------------------------------------------------
! laol_opacity
!----------------------------------------------------------------------
! The Los Alamos (LAOL89) atomic opacity tables, one per Z (a second at
! another Z when the run blends two), and the LAOL pure-Z table (Z = 1,
! a pure mixture of equal numbers of C and O, PURECO.DBGLAOL) that
! helium-burning cores blend toward.
!     read_laol_tables, build_laol_splines    load and spline the
!                                             LAOL89 tables
!     compute_laol_opacity                    the opacity and its log
!                                             derivatives from them
!     read_pure_z_table, build_pure_z_splines,
!     compute_pure_z_opacity                  the same for the pure-Z
!                                             table
! Both interpolate a column (laol_column_opacity) with cubic splines in
! log rho and in T, extrapolating up to ln(laol_extrapolation_tolerance)
! in log rho past a row's ends (a natural log on a log10 axis, as the F77
! code had it: 2.3 dex for the default 10; bug_sweep M227), and the
! LAOL89 tables linearly in X. The tables are laol_table_sets and a
! pure_z_table_set (opacity_tables), which kap_lib holds and passes in.
! The sources are at the end of the file.
module laol_opacity
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: laol_log10_opacity_cap, laol_opacity_cap, laol_table_set, n_laol_rho, &
           n_laol_t, n_laol_x, pure_z_table_set, kap_settings_t, kap_table_paths_t, &
           opacity_table_state
      use opacity_source, only: kap_point, kap_result, kap_source
      use opacity_blending, only: blend_results_in_z
      use numerics_lib, only: natural_spline, locate, splint, numerics_status_text
      use laol_reader, only: read_laol_tables, read_pure_z_table
      use math_lib
      implicit none
      private
      public :: laol89_source, pure_z_source
! laol_column_opacity's outcomes
      integer, parameter :: laol_column_ok = 0, laol_column_outside_density = 1, &
           laol_column_outside_temperature = 2, laol_column_splint_failed = 3
      public :: build_laol_splines, compute_laol_opacity, &
           build_pure_z_splines, compute_pure_z_opacity

contains


!----------------------------------------------------------------------
! build_laol_splines (formerly sulaol)
!----------------------------------------------------------------------
! Converts each table's T grid to log10 T and builds the natural cubic
! spline in log rho of every (X, T) row (the internal spline_laol_table),
! for compute_laol_opacity; the second table with two_z_tables.
subroutine build_laol_splines(tables, two_z_tables)
      type(laol_table_set), intent(inout) :: tables(:)
      logical, intent(in) :: two_z_tables

      real(dp) :: row_log10_opacity(n_laol_rho), row_log_rho(n_laol_rho), &
           row_d2opacity(n_laol_rho)
      integer :: it, ix, ir, num_valid_rho

      do it=1, tables(1)%num_t
         tables(1)%grid_t(it) = log10(tables(1)%grid_t(it))
      end do
      call spline_laol_table(tables(1))
      if (two_z_tables) then
       do it=1, tables(2)%num_t
          tables(2)%grid_t(it) = log10(tables(2)%grid_t(it))
       end do
         call spline_laol_table(tables(2))
      end if
      return

contains

! spline_laol_table: the natural cubic spline in log rho of every
! (X, T) row of tbl with at least 4 points (fewer: none, num_points 0).
subroutine spline_laol_table(tbl)
      type(laol_table_set), intent(inout) :: tbl

      do ix=1, tbl%num_x
         do it=1, tbl%num_t
            num_valid_rho=0
            do ir=1, tbl%num_rho
                tbl%slaol_opacity(ix,ir,it) = 0.0d0
                tbl%slaol_log_rho(ix,ir,it) = 0.0d0
                tbl%slaol_d2opacity(ix,ir,it) = 0.0d0
                if (tbl%opacity(ix,ir,it) /= 0.0d0) then
                   num_valid_rho = num_valid_rho+1
                   row_log10_opacity(num_valid_rho) = log10(tbl%opacity(ix,ir,it))
                   row_log_rho(num_valid_rho) = log10(tbl%grid_rho(ir))
                end if
            end do
            if (num_valid_rho >= 4) then
               tbl%slaol_num_points(ix,it)=num_valid_rho
               call natural_spline(row_log_rho, row_log10_opacity, num_valid_rho, &
                    row_d2opacity)
               do ir=1,num_valid_rho
                   tbl%slaol_opacity(ix,ir,it) = row_log10_opacity(ir)
                   tbl%slaol_log_rho(ix,ir,it) = row_log_rho(ir)
                   tbl%slaol_d2opacity(ix,ir,it) = row_d2opacity(ir)
               end do
            else
               tbl%slaol_num_points(ix,it) = 0
            end if
         end do
      end do
end subroutine spline_laol_table

end subroutine build_laol_splines

!----------------------------------------------------------------------
! compute_laol_opacity (formerly gtlaol)
!----------------------------------------------------------------------
! The LAOL89 opacity at (log rho, log T, X) from table tbl (laol(1) or
! laol(2)): the two X columns bracketing X each give log kappa and its
! derivatives (laol_column_opacity), interpolated linearly in X (or,
! outside the X grid, extrapolated). log10 kappa is capped at
! laol_log10_opacity_cap. A point no column covers returns ierr = 1 with
! the reason in message, which names the table (the second Z's too).
subroutine compute_laol_opacity(tbl, log10_density, log10_temperature, hydrogen_fraction, &
     opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, message, ierr)
      type(laol_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr
      real(dp) :: opacity_by_x(2), dlnkap_dlnrho_by_x(2), dlnkap_dlnt_by_x(2), x_values(2)
      integer :: x_grid_index, num_valid_x, x_loop_index, status
      real(dp) :: log10_opacity_value, slope
      character(len=16) :: table_name

      ierr = 0
      message = ''
      call locate(tbl%grid_x, tbl%num_x, hydrogen_fraction, x_grid_index)
      if (x_grid_index == tbl%num_x) then
          x_grid_index = tbl%num_x-1
      end if
      if (x_grid_index == 0) then
          x_grid_index = 1
      end if
      num_valid_x = 0
      do x_loop_index = x_grid_index, x_grid_index + 1
         num_valid_x = num_valid_x + 1
         call laol_column_opacity(tbl%slaol_opacity(x_loop_index,:,:), &
              tbl%slaol_log_rho(x_loop_index,:,:), tbl%slaol_d2opacity(x_loop_index,:,:), &
              tbl%slaol_num_points(x_loop_index,:), tbl%grid_t, tbl%num_t, &
              tbl%log_extrapolation_tolerance, log10_density, log10_temperature, &
              opacity_by_x(num_valid_x), dlnkap_dlnrho_by_x(num_valid_x), &
              dlnkap_dlnt_by_x(num_valid_x), status, ierr)
         if (status /= laol_column_ok) then
            table_name = 'LAOL 1989'
            if (tbl%table_number == 2) table_name = 'second LAOL 1989'
            call laol_column_message(status, ierr, table_name, log10_density, &
                 log10_temperature, message)
            if (status /= laol_column_splint_failed) ierr = 1
            return
         end if
         x_values(num_valid_x) = tbl%grid_x(x_loop_index)
      end do
! linear in X between the two columns
      slope = (opacity_by_x(2)-opacity_by_x(1))/(x_values(2)-x_values(1))
      log10_opacity_value = opacity_by_x(1)+slope*(hydrogen_fraction-x_values(1))
      slope = (dlnkap_dlnt_by_x(2)-dlnkap_dlnt_by_x(1))/(x_values(2)-x_values(1))
      dlnkap_dlnt = dlnkap_dlnt_by_x(1) + slope*(hydrogen_fraction-x_values(1))
      slope = (dlnkap_dlnrho_by_x(2)-dlnkap_dlnrho_by_x(1))/(x_values(2)-x_values(1))
      dlnkap_dlnrho = dlnkap_dlnrho_by_x(1)+slope*(hydrogen_fraction-x_values(1))
      call cap_laol_opacity(log10_opacity_value, opacity, log10_opacity)
end subroutine compute_laol_opacity


!----------------------------------------------------------------------
! build_pure_z_splines (formerly zsulaol)
!----------------------------------------------------------------------
! Converts the pure-Z table's T grid to log10 T and builds the natural
! cubic spline in log rho of each T row, for compute_pure_z_opacity.
subroutine build_pure_z_splines(tbl)
      type(pure_z_table_set), intent(inout) :: tbl

      real(dp) :: row_log10_opacity(n_laol_rho), row_log_rho(n_laol_rho), &
           row_d2opacity(n_laol_rho)
      integer :: it, ir, num_valid_rho

      do it=1, tbl%num_t
         tbl%logt_grid(it) = log10(tbl%logt_grid(it))
      end do
        do it=1, tbl%num_t
            num_valid_rho=0
            do ir=1, tbl%num_rho
                  tbl%spline_opacity(ir,it) = 0.0d0
                  tbl%spline_log_rho(ir,it) = 0.0d0
                  tbl%spline_d2opacity(ir,it) = 0.0d0
                  if (tbl%opacity(ir,it) /= 0.0d0) then
                     num_valid_rho = num_valid_rho+1
                     row_log10_opacity(num_valid_rho) = log10(tbl%opacity(ir,it))
                     row_log_rho(num_valid_rho) = log10(tbl%rho_grid(ir))
                  end if
            end do
            if (num_valid_rho >= 4) then
               tbl%spline_num_points(it)=num_valid_rho
               call natural_spline(row_log_rho, row_log10_opacity, num_valid_rho, &
                    row_d2opacity)
               do ir=1,num_valid_rho
                     tbl%spline_opacity(ir,it) = row_log10_opacity(ir)
                     tbl%spline_log_rho(ir,it) = row_log_rho(ir)
                     tbl%spline_d2opacity(ir,it) = row_d2opacity(ir)
               end do
            else
               tbl%spline_num_points(it) = 0
            end if
        end do
end subroutine build_pure_z_splines

!----------------------------------------------------------------------
! compute_pure_z_opacity (formerly gtpurz)
!----------------------------------------------------------------------
! The pure-Z table's opacity at (log rho, log T): its one column
! (laol_column_opacity), log10 kappa capped as compute_laol_opacity's.
subroutine compute_pure_z_opacity(tbl, log10_density, log10_temperature, opacity, &
     log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, message, ierr)
      type(pure_z_table_set), intent(in) :: tbl
      real(dp), intent(in) :: log10_density, log10_temperature
      real(dp), intent(out) :: opacity, log10_opacity, &
           dlnkap_dlnrho, dlnkap_dlnt
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr
      real(dp) :: log10_opacity_value
      integer :: status

      ierr = 0
      message = ''
      call laol_column_opacity(tbl%spline_opacity, tbl%spline_log_rho, tbl%spline_d2opacity, &
           tbl%spline_num_points, tbl%logt_grid, tbl%num_t, tbl%log_extrapolation_tolerance, &
           log10_density, log10_temperature, log10_opacity_value, dlnkap_dlnrho, &
           dlnkap_dlnt, status, ierr)
      if (status /= laol_column_ok) then
         call laol_column_message(status, ierr, 'pure-Z', log10_density, &
              log10_temperature, message)
         if (status /= laol_column_splint_failed) ierr = 1
         return
      end if
      call cap_laol_opacity(log10_opacity_value, opacity, log10_opacity)
end subroutine compute_pure_z_opacity

!----------------------------------------------------------------------
! laol_column_opacity
!----------------------------------------------------------------------
! One column of an LAOL table (one X; the pure-Z table's only one), from
! its splines in log rho, one per T row (spline_opacity etc. are
! (density, T)): on the four T rows around log T (select_four_point_window)
! log kappa at log rho from the row's spline, or within
! log_extrapolation_tolerance (in log rho) of a row's ends linearly from
! its two end knots, with d log kappa/d log rho the secant of the knots
! that bracket log rho; then a natural spline in log T through the rows
! that gave a value, for log kappa, and the secants of its bracketing
! rows for the derivatives (dln kappa/dln rho interpolated linearly in T
! between theirs). status: laol_column_ok; laol_column_outside_density
! (a row with fewer than 4 points), laol_column_outside_temperature
! (fewer than 4 rows reach log rho), or laol_column_splint_failed with
! splint's code in jerr.
subroutine laol_column_opacity(spline_opacity, spline_log_rho, spline_d2opacity, &
     spline_num_points, logt_grid, num_t, log_extrapolation_tolerance, log10_density, &
     log10_temperature, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt, status, jerr)
      real(dp), intent(in) :: spline_opacity(:,:), spline_log_rho(:,:), spline_d2opacity(:,:)
      integer, intent(in) :: spline_num_points(:), num_t
      real(dp), intent(in) :: logt_grid(:), log_extrapolation_tolerance, log10_density, &
           log10_temperature
      real(dp), intent(out) :: log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
      integer, intent(out) :: status, jerr
      real(dp) :: row_log10_opacity(n_laol_rho), row_log_rho(n_laol_rho), &
           row_d2opacity(n_laol_rho), dlnkap_dlnrho_by_t(n_laol_t)
      real(dp) :: logt_interp_opacity(n_laol_t), logt_values(n_laol_t), &
           logt_d2opacity(n_laol_t)
      integer :: t_locate_guess, t_index, t_range_lo, t_range_hi, num_valid_t, &
           num_valid_rho, rho_loop_index, spline_index_lo, spline_index_hi
      real(dp) :: log10_opacity_value, slope

      status = laol_column_ok
      jerr = 0
      call locate(logt_grid, num_t, log10_temperature, t_locate_guess)
      num_valid_t = 0
      call select_four_point_window(t_locate_guess, num_t, t_range_lo, t_range_hi)
      do t_index = t_range_lo, t_range_hi
         num_valid_rho = spline_num_points(t_index)
         if (num_valid_rho < 4) then
            status = laol_column_outside_density
            return
         end if
         do rho_loop_index = 1, num_valid_rho
            row_log10_opacity(rho_loop_index) = spline_opacity(rho_loop_index,t_index)
            row_log_rho(rho_loop_index) = spline_log_rho(rho_loop_index,t_index)
            row_d2opacity(rho_loop_index) = spline_d2opacity(rho_loop_index,t_index)
         end do
         if (log10_density > row_log_rho(1) .and. &
              log10_density < row_log_rho(num_valid_rho)) then
            call splint(row_log_rho, row_log10_opacity, num_valid_rho, &
                 row_d2opacity, log10_density, log10_opacity_value, &
                 jerr, spline_index_lo, spline_index_hi)
            if (jerr /= 0) then
               status = laol_column_splint_failed
               return
            end if
            num_valid_t = num_valid_t+1
            logt_interp_opacity(num_valid_t) = log10_opacity_value
            logt_values(num_valid_t) = logt_grid(t_index)
            dlnkap_dlnrho_by_t(num_valid_t) = &
                 (row_log10_opacity(spline_index_hi)-row_log10_opacity(spline_index_lo))/ &
                 (row_log_rho(spline_index_hi)-row_log_rho(spline_index_lo))
         else if (log10_density > row_log_rho(1)-log_extrapolation_tolerance .and. &
              log10_density <= row_log_rho(1)) then
            slope = (row_log10_opacity(2)-row_log10_opacity(1))/ &
                 (row_log_rho(2)-row_log_rho(1))
            log10_opacity_value = row_log10_opacity(1)+slope*(log10_density-row_log_rho(1))
            num_valid_t = num_valid_t+1
            logt_interp_opacity(num_valid_t) = log10_opacity_value
            logt_values(num_valid_t) = logt_grid(t_index)
            dlnkap_dlnrho_by_t(num_valid_t) = slope
         else if (log10_density >= row_log_rho(num_valid_rho) .and. &
              log10_density < row_log_rho(num_valid_rho)+log_extrapolation_tolerance) then
            slope = (row_log10_opacity(num_valid_rho-1)-row_log10_opacity(num_valid_rho))/ &
                 (row_log_rho(num_valid_rho-1)-row_log_rho(num_valid_rho))
            log10_opacity_value = row_log10_opacity(num_valid_rho)+ &
                 slope*(log10_density-row_log_rho(num_valid_rho))
            num_valid_t = num_valid_t+1
            logt_interp_opacity(num_valid_t) = log10_opacity_value
            logt_values(num_valid_t) = logt_grid(t_index)
            dlnkap_dlnrho_by_t(num_valid_t) = slope
         end if
      end do
      if (num_valid_t < 4) then
         status = laol_column_outside_temperature
         return
      end if
      call natural_spline(logt_values, logt_interp_opacity, num_valid_t, &
           logt_d2opacity)
      call splint(logt_values, logt_interp_opacity, num_valid_t, &
           logt_d2opacity, log10_temperature, log10_opacity_value, &
           jerr, spline_index_lo, spline_index_hi)
      if (jerr /= 0) then
         status = laol_column_splint_failed
         return
      end if
      log10_opacity = log10_opacity_value
      slope=(dlnkap_dlnrho_by_t(spline_index_hi)-dlnkap_dlnrho_by_t(spline_index_lo))/ &
           (logt_values(spline_index_hi)-logt_values(spline_index_lo))
      dlnkap_dlnrho = dlnkap_dlnrho_by_t(spline_index_lo)+ &
           slope*(log10_temperature-logt_values(spline_index_lo))
      dlnkap_dlnt = (logt_interp_opacity(spline_index_hi)-logt_interp_opacity(spline_index_lo))/ &
           (logt_values(spline_index_hi)-logt_values(spline_index_lo))
end subroutine laol_column_opacity

! kappa and log10 kappa from log10_opacity_value, capped at
! laol_opacity_cap (log10 kappa laol_log10_opacity_cap) where it overflows
! The reason a column lookup failed (status, laol_column_opacity's), for
! the table named: a spline's numerics status (jerr), or the point off
! the table's rows in log rho or its columns in log T, beyond the
! extrapolation tolerance.
subroutine laol_column_message(status, jerr, table_name, log10_density, &
     log10_temperature, message)
      integer, intent(in) :: status, jerr
      character(len=*), intent(in) :: table_name
      real(dp), intent(in) :: log10_density, log10_temperature
      character(len=*), intent(out) :: message

      if (status == laol_column_splint_failed) then
         message = 'the ' // table_name // ' table''s spline in log T failed: ' // &
              trim(numerics_status_text(jerr))
      else if (status == laol_column_outside_density) then
         write(message,'(a,f8.4,a,f8.4,a)') 'log rho =', log10_density, ', log T =', &
              log10_temperature, ' is outside the ' // table_name // ' table in log rho'
      else
         write(message,'(a,f8.4,a,f8.4,a)') 'log rho =', log10_density, ', log T =', &
              log10_temperature, ' is outside the ' // table_name // ' table in log T'
      end if
end subroutine laol_column_message

pure subroutine cap_laol_opacity(log10_opacity_value, opacity, log10_opacity)
      real(dp), intent(in) :: log10_opacity_value
      real(dp), intent(out) :: opacity, log10_opacity

      if (log10_opacity_value > laol_log10_opacity_cap) then
         opacity = laol_opacity_cap
         log10_opacity = laol_log10_opacity_cap
      else
         opacity = exp10(log10_opacity_value)
         log10_opacity = log10_opacity_value
      end if
end subroutine cap_laol_opacity

!----------------------------------------------------------------------
! select_four_point_window (formerly util/xrng4)
!----------------------------------------------------------------------
! The four table points window_start..window_end around grid point
! grid_index of a grid of grid_size points, for a 4-point interpolation:
! the first or last four points near the ends, else grid_index-1 ..
! grid_index+2.
subroutine select_four_point_window(grid_index, grid_size, window_start, window_end)
      integer, intent(in) :: grid_index, grid_size
      integer, intent(out) :: window_start, window_end
      if (grid_index <= 2) then
         window_start = 1
         window_end = 4
      else if (grid_index >= grid_size-2) then
         window_start = grid_size-3
         window_end = grid_size
      else
         window_start = grid_index - 1
         window_end = grid_index + 2
      end if
end subroutine select_four_point_window

!----------------------------------------------------------------------
! The family's descriptor (opacity_source) and its procedures
!----------------------------------------------------------------------
! An interior family of one table per Z, as OPAL 1992 (laol_table_z1,
! laol_table_z2).
function laol89_source(settings) result(source)
      type(kap_settings_t), intent(in) :: settings
      type(kap_source) :: source

      source%name = 'LAOL 1989'
      source%blends_two_z = .true.
      source%table_z = [settings%laol_table_z1, settings%laol_table_z2]
      source%table_z_controls = [character(len=16) :: 'laol_table_z1', 'laol_table_z2']
      source%selected = settings%use_laol89_tables
      source%load => load_laol89
      source%lookup => lookup_laol89
      source%core_lookup => core_lookup_laol89
end function laol89_source

subroutine load_laol89(settings, paths, envelope_hydrogen_fraction, tables, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(kap_table_paths_t), intent(in) :: paths
      real(dp), intent(in) :: envelope_hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr
! the table header's 12-element metal mixture (read, not used: the
! star takes its surface mixture from the controls)
      real(dp) :: laol_header_mixture(12)

! (the reader reports its own errors)
      message = ''
      laol_header_mixture = 0.0d0
      call read_laol_tables(tables%laol, laol_header_mixture, paths%laol, &
           paths%laol_z2, settings%use_two_z_tables, &
           settings%laol_extrapolation_tolerance, ierr)
      if (ierr /= 0) return
      call build_laol_splines(tables%laol, settings%use_two_z_tables)
end subroutine load_laol89

subroutine lookup_laol89(settings, tables, point, res, in_table, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: in_table
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr
      type(kap_result) :: res_z2

      in_table = .true.
      call compute_laol_opacity(tables%laol(1), point%log10_density, &
           point%log10_temperature, point%hydrogen_fraction, res%opacity, &
           res%log10_opacity, res%dlnkap_dlnrho, res%dlnkap_dlnt, message, ierr)
      if (ierr /= 0) return
      if (settings%use_two_z_tables) then
         call compute_laol_opacity(tables%laol(2), point%log10_density, &
              point%log10_temperature, point%hydrogen_fraction, res_z2%opacity, &
              res_z2%log10_opacity, res_z2%dlnkap_dlnrho, res_z2%dlnkap_dlnt, message, ierr)
         if (ierr /= 0) return
         call blend_results_in_z(res, res_z2, point%metal_fraction, &
              settings%laol_table_z1, settings%laol_table_z2)
      end if
end subroutine lookup_laol89

subroutine core_lookup_laol89(settings, tables, point, res, table_z, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      real(dp), intent(out) :: table_z
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      call compute_laol_opacity(tables%laol(1), point%log10_density, &
           point%log10_temperature, point%hydrogen_fraction, res%opacity, &
           res%log10_opacity, res%dlnkap_dlnrho, res%dlnkap_dlnt, message, ierr)
      table_z = settings%laol_table_z1
end subroutine core_lookup_laol89

! The pure-Z table (Z = 1, pure C and O), the helium-burning core's
! table: kap_eval blends an interior family's opacity toward it.
function pure_z_source(settings) result(source)
      type(kap_settings_t), intent(in) :: settings
      type(kap_source) :: source

      source%name = 'pure Z (LAOL)'
      source%selected = settings%use_pure_z_table
      source%load => load_pure_z
      source%lookup => lookup_pure_z
end function pure_z_source

subroutine load_pure_z(settings, paths, envelope_hydrogen_fraction, tables, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(kap_table_paths_t), intent(in) :: paths
      real(dp), intent(in) :: envelope_hydrogen_fraction
      type(opacity_table_state), intent(inout) :: tables
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

! (the reader reports its own errors)
      message = ''
      call read_pure_z_table(tables%pure_z, paths%pure_z, &
           settings%laol_extrapolation_tolerance, ierr)
      if (ierr /= 0) return
      call build_pure_z_splines(tables%pure_z)
end subroutine load_pure_z

subroutine lookup_pure_z(settings, tables, point, res, in_table, message, ierr)
      type(kap_settings_t), intent(in) :: settings
      type(opacity_table_state), intent(inout) :: tables
      type(kap_point), intent(in) :: point
      type(kap_result), intent(inout) :: res
      logical, intent(out) :: in_table
      character(len=*), intent(out) :: message
      integer, intent(out) :: ierr

      in_table = .true.
      call compute_pure_z_opacity(tables%pure_z, point%log10_density, &
           point%log10_temperature, res%opacity, res%log10_opacity, &
           res%dlnkap_dlnrho, res%dlnkap_dlnt, message, ierr)
end subroutine lookup_pure_z

!----------------------------------------------------------------------
! Sources
!----------------------------------------------------------------------
! The tables: the Los Alamos opacities (Huebner et al. 1977; the LAOL89
! files in input/opacity/atomic/legacy/opal92/) and the pure-Z table
! PURECO.DBGLAOL (input/opacity/atomic/lanl/: "a pure mixture of equal
! numbers of C and O", its first line). The code: the F77 RDLAOL, SULAOL,
! GTLAOL, RDZLAOL, ZSULAOL, GTPURZ and XRNG4, with DBG's second Z table
! (4/94) and the Debye-Hueckel mixture's read (7/92), converted to
! Fortran 90 in 2026: the two tables' copies one routine on a table set
! passed in, the second table's T grid converted over its own rows
! (M228), the lookups' common column one kernel, the reads checked
! (M220).

end module laol_opacity
